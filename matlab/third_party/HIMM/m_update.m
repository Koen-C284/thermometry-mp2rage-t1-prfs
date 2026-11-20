function m = m_update(imgsw,Zw,m,tes,algp,mask,body_mask,muscle_mask,R,drift_fit,onlypositive)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
% HIMM: Harmonic Initialized Model-based Multi-echo
% m_update: function to fit PRF shift due to heating

% --- EXTRACT FAT MASK FROM ALGP ---
% We extract the fat mask here to avoid changing the function signature
if isfield(algp, 'fat_mask')
    fat_mask = algp.fat_mask;
else
    fat_mask = []; % Fallback if not provided
end

% reshape images and model into vectors
innprod = Zw.*conj(imgsw);

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;

hesssum = 0; gradsum = 0;
for ii = 1:length(tes)
    hesssum = hesssum + t4(:,:,ii)./t2(:,:,ii)*tes(ii)^2;
    gradsum = gradsum + t4(:,:,ii)*tes(ii);
end

if exist('R','var') % add roughness penalty
    if ~isempty(R)
        hesssum = hesssum + R.denom(R,m);
        gradsum = gradsum + R.cgrad(R,m);
    end
end

if exist('mask','var') % apply mask
    if ~isempty(mask)
        gradsum = gradsum.*mask;
    end
end

if exist('signal_mask','var') % apply secondary signal mask
    if ~isempty(signal_mask)
        gradsum = gradsum.*abs(1-signal_mask);
    end
end

% --- DEFENSE 1: ZERO OUT THE GRADIENT IN FAT VOXELS ---
if ~isempty(fat_mask)
    % Multiply by the logical NOT of the fat mask. 
    % (1 - 1 = 0), so the gradient in fat becomes exactly 0.
    gradsum = gradsum .* (1 - double(fat_mask));
end

% update m
dm = -gradsum./hesssum;
dm(isnan(dm)) = 0;
m = m + dm;

% Enforce limits
if drift_fit == 0
    muscle_mask = imerode(muscle_mask,strel('diamond', 0)); % 3
    m(muscle_mask == 0) = 0;
end

% --- DEFENSE 2: HARD ZERO THE TEMPERATURE IN FAT VOXELS ---
if ~isempty(fat_mask)
    m(fat_mask == 1) = 0;
end

% Enforce negativity (temp increase -> freq decrease)
if onlypositive
    m(m > 0) = 0;
end
end