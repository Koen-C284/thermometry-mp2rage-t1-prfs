% compute_T1.m - Look-Locker T1 Map Estimation Function
%

%
% Syntax:
%   [t1_map, t1_map_image] = compute_T1(img, info)
%
% Inputs:
%   img  - 6D or 7D complex/magnitude MRI dataset [x, y, z, coils, dynamics, inversion_times, (real/imag)]
%   info - Struct containing image header definitions (must contain info.imgdef.trigger_time.uniq)
%
% Outputs:
%   t1_map       - 3D Look-Locker corrected T1 map in seconds [x, y, z]
%   t1_map_image - 3D signal intensity reference image [x, y, z]

function [t1_map, t1_map_image] = compute_T1(img, info)


%% === Define vectors ===
if isfield(info,'imgdef') && isfield(info.imgdef,'trigger_time') && isfield(info.imgdef.trigger_time,'uniq')
    TI_ms = double(info.imgdef.trigger_time.uniq(:));
else
    error('Could not find trigger_time in info.imgdef.trigger_time.uniq');
end

tVec = TI_ms * 1e-3;

%% === Define Dimensions ===
[nx, ny, nz, ~, nDyn, nTI, ~] = size(img);

if length(size(img)) == 6
    data = img;
elseif length(size(img)) == 7
    data = img(:,:,:,:,:,:,1) .* exp(1i * img(:,:,:,:,:,:,2));
else
    error("Unknown data size for T1 map")
end
data = reshape(data, [nx, ny, nz, nDyn, nTI]);
data = mean(data, 4);
t1_map_image = mean(abs(data), 5);
%% === Build LUT for T1* ===
% Choose a practical grid in seconds (inversion recovery typically 0.1–4 s)
T1min = 0.05; T1max = 4.0; nGrid = 4000;
T1grid = linspace(T1min, T1max, nGrid).'; % column
nlsS = buildNlsStruct(tVec, T1grid);

%% === Make a simple mask (skip background) ===
% Use max across TI as signal proxy
sigMax = max(abs(data), [], ndims(data));
% sigMax has shape [x y z dyn]; build a consistent 4D logical mask
mask = sigMax > (0.1 * max(sigMax(:))); % 10% of global max

%% === Pre-allocate outputs ===
T1star = nan(nx,ny,nz, 'double');
Amap   = nan(nx,ny,nz, 'double');
Bmap   = nan(nx,ny,nz, 'double');
Res    = nan(nx,ny,nz, 'double');

%% === Fit per dynamic (keeps memory reasonable) ===
fprintf('Fitting Look-Locker T1*... (nx=%d, ny=%d, nz=%d, dyn=%d, TI=%d)\n', nx,ny,nz,nDyn,nTI);

% Extract [x y z TI] for this dynamic
Y = squeeze(data(:,:,:,:)); % -> [nx ny nz nTI]

Y = reshape(Y, [], nTI);      % -> [Nvox nTI]
m = reshape(mask(:,:,:), [], 1);

% Pre-alloc local
T1s_local = nan(size(m), 'double');
A_local   = nan(size(m), 'double');
B_local   = nan(size(m), 'double');
R_local   = nan(size(m), 'double');

idx = find(m);
for ii = 1:numel(idx)
    v = idx(ii);
    y = double(squeeze(Y(v, :))).'; % column [nTI x 1]

    % Skip zeros/NaNs
    if ~any(isfinite(y)) || max(abs(y)) < eps
        continue;
    end

    % Decide which fitting routine to use based on complex-valued input
    if ~isreal(y)
        % Use rdNls which handles complex inputs (returns T1, b, a, res)
        [T1s, b, a, res] = rdNls(y, nlsS);
    else
        % Use rdNlsPr (polarity-restoration variant) for magnitude/real data
        [T1s, b, a, res] = rdNlsPr(y, nlsS);
    end

    T1s_local(v) = double(T1s);
    A_local(v)   = double(a);
    B_local(v)   = double(b);
    R_local(v)   = double(res);
end

% Reshape back
T1star(:,:,:) = reshape(T1s_local, nx,ny,nz);
Amap(:,:,:)   = reshape(A_local,   nx,ny,nz);
Bmap(:,:,:)   = reshape(B_local,   nx,ny,nz);
Res(:,:,:)    = reshape(R_local,   nx,ny,nz);


%% === Look-Locker correction ===
% Model fitted above: S(t) = a + b*exp(-t/T1*). Look-Locker T1 correction for the segmented (TFE) readout.
ratio = -(Bmap ./ Amap - 1);
t1_map = mask.*abs(T1star .* ratio);

% Clean unreasonable values
t1_map(~isfinite(t1_map)) = 0;
t1_map(t1_map < T1min)        = 0;
t1_map(t1_map > T1max)        = 0; % clip to [T1min, T1max]

fprintf("Done fitting LL T1 maps ... \n")
end


%% === Helper functions ===
% Reduced-dimension nonlinear least-squares T1 fitting (buildNlsStruct, rdNls, rdNlsPr)
% adapted from Barral JK et al., Magn Reson Med 2010;64(4):1057-1067.

    function nlsS = buildNlsStruct(tVec, T1grid)
        % Precompute exponentials and projection helpers for LUT
        tVec = tVec(:);                  % [nTI x 1]
        alpha = 1 ./ T1grid(:).';        % [1 x nGrid]
        theExp = exp(-tVec * alpha);     % [nTI x nGrid]

        nlsS.tVec       = tVec;
        nlsS.N          = numel(tVec);
        nlsS.T1Vec      = T1grid(:);
        nlsS.theExp     = theExp;        % [N x G]
        nlsS.sumExp     = sum(theExp, 1).';  % [G x 1]
        nlsS.rhoNormVec = sum(theExp.^2, 1).' - (1/nlsS.N) * (nlsS.sumExp).^2; % [G x 1]
    end

    function [T1Est, bEst, aEst, res] = rdNls(data, nlsS)

        y = data(:);
        ySum = sum(y);

        rhoTyVec = (y.'*nlsS.theExp).' - ySum/nlsS.N*sum(nlsS.theExp,1)';
        theExp = nlsS.theExp;

        [~,ind] = max( abs(rhoTyVec).^2./nlsS.rhoNormVec);

        T1Vec = nlsS.T1Vec;

        % The estimated parameters
        T1Est = T1Vec(ind);
        bEst = rhoTyVec(ind) / nlsS.rhoNormVec(ind);
        aEst = 1/nlsS.N*(ySum - bEst*sum(theExp(:,ind)));

        % Compute the residual
        model = aEst + bEst*exp(-nlsS.tVec/T1Est);
        res = 1/sqrt(nlsS.N) * norm(1 - model./data);
    end

    function [T1Est, bEst, aEst, res] = rdNlsPr(data, nlsS)

        aEstTmp = zeros(1,2);
        bEstTmp = zeros(1,2);
        T1EstTmp = zeros(1,2);
        resTmp = zeros(1,2);

        y = data(:);

        [~, minInd] = min(y);

        for ii = 1:2

            if ii == 1
                dataTmp = y.*[-ones(minInd,1); ones(nlsS.N - minInd,1)];
            elseif ii == 2
                dataTmp = y.*[-ones(minInd-1,1); ones(nlsS.N - (minInd-1),1)];
            end

            ySum = sum(dataTmp);
            rhoTyVec = (dataTmp.'*nlsS.theExp).' - 1./nlsS.N*sum(nlsS.theExp,1)'*ySum;
            rhoNormVec = nlsS.rhoNormVec;

            [~,ind] = max(abs(rhoTyVec).^2./rhoNormVec);
            T1Vec = nlsS.T1Vec;

            T1EstTmp(ii) = T1Vec(ind);
            bEstTmp(ii) = rhoTyVec(ind)/ rhoNormVec(ind);
            aEstTmp(ii) = 1/nlsS.N*(ySum - bEstTmp(ii)*sum(nlsS.theExp(:,ind)));
            modelValue = aEstTmp(ii) + bEstTmp(ii)*exp(-nlsS.tVec/T1EstTmp(ii));
            resTmp(ii) = 1/sqrt(nlsS.N) * norm(1 - modelValue./dataTmp);
        end

        [res,ind] = min(resTmp);
        aEst = aEstTmp(ind);
        bEst = bEstTmp(ind);
        T1Est = T1EstTmp(ind);

    end
