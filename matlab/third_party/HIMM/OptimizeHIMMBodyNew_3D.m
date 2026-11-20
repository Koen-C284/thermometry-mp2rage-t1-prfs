function [m,wts,Ac,c,A,phi,Ac_original] = OptimizeHIMMBodyNew(imgs, algp, scanp, lib, Acswitch,mswitch,phiswitch,minit,phiinit)
%% Multi-echo MR Thermometry around the human knee at 3T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
%
% Function to solve for heating, field shifts, and select a baseline in the
% presence of fat. 
%
% Algorithm code for the iterative solver was developed by Megan Poorman 
% and William Grissom. The code presented in this version does not use the 
% regularization term anymore as it is not focused on applications using
% localized heating. Alternarively, the drift field (Ac) was initialized 
% using near-harmonic 2D reconstruction (proposed by Salomir et al.)
%
% Creators: Megan Poorman, William Grissom
% Updated by: Mathijs Kikken
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x slices x Necho)
%   ----- algorithm params
%     algp.order = order;
%     algp.nmiter = number of m fit iters;
%     algp.nciter = number of c fit iters;
%     algp.stopthresh = fit tolerance;
%     algp.suppress_info = whether updates will be visualized on the command window
%
%   ----- scan parameters
%     scanp.dim = [rows cols scans];    %% LET OP: bevat nu ook als 3e dimensie 'scans'
%     scanp.b0 = b0;
%     scanp.tes = tes;
%     scanp.prc = prc;
% 
%   ----- baseline library
%     lib.Wlib = Wlib;
%     lib.Flib = Flib;
%     lib.dw0lib = dw0lib;
%     lib.R2starlib = R2starlib;
%     lib.fatmodel = fatmodel;        % minit: initial m (frequency shift map) [3D]
%
% Output: 
%   m:      frequency shift due to heating
%   wts:    baseline weight
%   c:      polynomial coefficients fit
%   A:      polynomial basis functions
%   Ac:     estimated drift field
%   phi:    Tx/Rx gain estimated
%   Ac_original:   original drift field at that time, used to calculate residuals


%% --- Initialization --- %%
%--- construct 2nd order finite difference roughness penalty object
% keyboard
if isfield('algp','beta')
  R = Robject(ones(scanp.dim(1),scanp.dim(2)),'order',2,'beta',algp.beta,'type_denom','matlab');
else
  R = [];
end

%--- Initialize vars for backward compatibility
if ~isfield(algp,'usequadprog')
   algp.usequadprog = 1;
end
if ~isfield(algp,'nlib')
    algp.nlib = size(lib.Wlib,3);
end

%--- Create polynomial matrix of given order (and degree)
if strcmp(algp.method,'poly')
    A = get_polynomial_function(algp.order,scanp.dim,3);     %% LET OP: dimensies van 2 -> 3 gezet
elseif strcmp(algp.method,'spher')
    A = get_spherical_function(algp.order, scanp.dim);
else
    message = {'The function with coefficients that are to be '
        'optimized was specified incorrection, '
        'please set algo.method to the correct string. '};
    errorStruct.message = sprintf('%s\n',message{:});
    errorStruct.identifier = 'MyFunction:parameterIncorrecctlyDefined';
    error(errorStruct)
end

%--- Initialize Ac, Ac_salomir, m, and phi
if isfield(lib,'cinit')
    c = lib.cinit;
else
    c = zeros(size(A,2),1);
end

Ac = reshape(A*c, scanp.dim);    % assumes equal number of pixels in the x & y directions

if ~exist('minit','var')
    m = zeros(scanp.dim(1),scanp.dim(2),scanp.dim(3));    %LET OP scanpdim(3) toegevoegd voor 3D maken!
else
    m = minit;
end

if ~exist('phiinit','var')
    phi = zeros(scanp.dim(1),scanp.dim(2),scanp.dim(3));    %LET OP scanpdim(3) toegevoegd voor 3D maken!
else
    phi = phiinit;
end 

% Define Ac initialization: to be used in salomir approach
if isfield(lib,'Acinit')
    Ac_salomir = lib.Acinit;
else
    Ac_salomir = zeros(scanp.dim);
end

if length(scanp.dim) == 1 %assume square
    scanp.dim = [scanp.dim scanp.dim];
end


%% --- Start HIMM iterations --- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Summary of the current algorithm code

% The code consists of 3 separate sections
% 1. Near Harmonic 2D reconstruction (initialization phase)
%    phase accumulation between the current dynamic and the reference 
%    dynamic(s) is determined. This phase accumulation is used for fat-based
%    near-harmomic 2D reconstruction. 
%
% 2. Optimization (with positivity constrain)
%    A jointly optimization of drift fields and temperature fields is 
%    performed. Here, temperature changes are constrained to be positive 
%    such that drift fields converge in the right direction. 
%
% 3. Optimization (without positivity constrain)
%    A jointly optimization of drift fields and temperature fields is 
%    performed. Here, temperature changes are not constrained to be positive anymore.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization factor (median of heating dynamics * num echoes)
algp.norm = abs(median(imgs,'all') * length(scanp.tes));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. initialize and perform near-harmonic 2D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run f_update with current m, Ac, phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First iteration
if algp.nlib < size(lib.Wlib,3)
    % compress the library down to the nlib entries that best 
    % correlate with the dynamic image
    [wts,algp.libinds] = f_update(abs(imgs),abs(Z),algp.usequadprog,algp.nlib);
    lib.Wlib = lib.Wlib(:,:,:,algp.libinds);
    lib.Flib = lib.Flib(:,:,:,algp.libinds);
    lib.dw0lib = lib.dw0lib(:,:,:,algp.libinds);
    lib.R2starlib = lib.R2starlib(:,:,:,algp.libinds);
else
    % For the first few iterations we will initialize with complex values from only the final echo
    wts = f_update(imgs(:,:,:,1:end),lib.drift_lib.imgs_ref(:,:,1:end,:),algp.usequadprog);           % : toegevoegd in imgs, omdat z er nu ook in zit
end

% Update normalization factor: median of heating dynamics * num echoes
[body_mask, fat_idx, muscle_mask] = get_masks(algp.body_mask,lib.Wlib,lib.Flib,wts,scanp.dim,algp.th_fat,algp.thickness_fat,algp.body_erode);
% size(imgs)NH2Dr_fit_mask
% size(body_mask)
body_masks = repmat(body_mask,1,1,1,length(scanp.tes));
% size(body_masks)
imgsdyn_masked = imgs; imgsdyn_masked(body_masks ~= 1) = nan; 
imgsdyn_masked = imgsdyn_masked(:);
imgsdyn_masked = imgsdyn_masked(~isnan(imgsdyn_masked));
algp.norm = abs(median(imgsdyn_masked) * length(scanp.tes));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run c_update with current f, m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if Acswitch == 1

    % Near-harmonic 2D reconstruction is performed in the fatty region only
    % QR decomposition is performed in both fat and muscle regions
    NH2DR_fit_mask = fat_idx;
    QR_fit_mask = fat_idx + muscle_mask; QR_fit_mask(QR_fit_mask > 1) = 1;

    for kk = 1:algp.nciter
            
        % Calculate and apply weights to baseline images
        Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);       %scanp.dim(3) toegevoegd

        % We want to get a first estimate of Ac by performing near-harmonic 2D
        % reconstruction on the fat. For that, we first need to calculate the
        % phase accumulation in fat
        for jj = 1:size(lib.Wlib,4)
            % Compute model values
            Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
             scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
             lib.R2starlib(:,:,:,jj),zeros(size(m)),lib.fatmodel);
            Ztot = Ztot + Z*wts(jj);
        end
        Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]);        %% LET OP: size(Z,3) veranderd naar 1, omdat phi zelf al 3D is 
    
        % Update m
        Ac = m_update(imgs(:,:,:,scanp.echo_init_ix:end),Ztot(:,:,:,scanp.echo_init_ix:end),Ac,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,muscle_mask,R,Acswitch,0);
        Ac_original = Ac;
        Ac = Ac.*NH2DR_fit_mask;
    end


    % Apply the 2D smoothing filter, by looping over the slices
    for zz = 1:size(Ac,3)
        Ac_slice = Ac(:,:,zz);
        NH2DR_fit_mask_slice = NH2DR_fit_mask(:, :, zz);

        Ac_slice(NH2DR_fit_mask_slice == 0) = NaN;
        Ac_smooth = smoothdata2(Ac_slice, "movmedian", algp.median_filter_size, "omitnan");
        Ac_smooth(NH2DR_fit_mask_slice == 0) = 0;

        Ac(:,:,zz) = Ac_smooth;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run Salomir with recently acquired Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('        Performing Near-Harmonic 2D reconstruction\n');
    %figure; imagesc(Ac.*NH2DR_fit_mask); title('Drift field prior to near-harmonic 2D reconstruction'); colorbar; axis image; axis off; drawnow;
    [Ac,c,Ac_salomir] = NH2DReconExecution(A,Ac,Ac_salomir,fat_idx,body_mask,QR_fit_mask,algp,scanp);  % AC erin is alleen rondje vet, AC eruit is volledige mask van driftveld
    % figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4]);
    % for slice = 1:5
    %     subplot(1,5,slice); imagesc(Ac_salomir(:,:,slice).*body_mask(:,:,slice), [15 25]); axis image; axis off; colorbar;
    % end 
    % drawnow;
    % 
    % figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4]);
    % for slice = 1:5
    %     subplot(1,5,slice); imagesc(Ac(:,:,slice).*body_mask(:,:,slice), [15 25]); axis image; axis off; colorbar;
    % end    
    % drawnow;    
        % figure; imagesc(Ac_salomir.*body_mask); title('Drift field after near-harmonic 2D reconstruction'); colorbar; axis image; axis off; drawnow;
        % figure; imagesc(Ac.*body_mask); title('Drift field after near-harmonic 2D reconstruction'); colorbar; axis image; axis off; drawnow;
end


if Acswitch == 0
    Ac = zeros(size(Ac));
    c = zeros(size(c));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run m_update with current f, Ac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mswitch
    for kk = 1:algp.nmiter

        % Build Z matrix using true values of Ac and m
        Zw = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
        Zf = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
        
        for jj = 1:size(lib.Wlib,4)
            % Compute model values of water and fat separately
            [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,:,jj).*exp(1i*phi),lib.Flib(:,:,:,jj).*exp(1i*phi), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ... 
                lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
            Zw = Zw + foow*wts(jj);
            Zf = Zf + foof*wts(jj);
        end

        % Update m
        m = m_update(imgs(:,:,:,scanp.echo_init_ix:end)-Zf(:,:,:,scanp.echo_init_ix:end),Zw(:,:,:,scanp.echo_init_ix:end),m,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,muscle_mask,R,Acswitch,0);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. start the iterative process of fitting with spherical harmonics
%    Executed with a positivity constrain in temperature
%    Ac is initialized properly through near-harmonic 2D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if algp.positivity_constrain == 1

    % calculate a terrible cost that is guaranteed to get iterations going
    cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp); 
    cost = 2*cost;
    ii = 0; % iteration counter
    
    % The optimization algorithm is run for at least 2 iterations
    converged = false;
    while converged == false || ii <= 1
    
        iter_stop = 9; % maximum of 10 iterations
        if ii > iter_stop
            break
        end   
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run c_update with current f, m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Evaluate in voxels where the body is present, such that fitting is insensitive to background noise
        fit_mask = fat_idx + muscle_mask; fit_mask(fit_mask > 1) = 1;
        
        if Acswitch && ii > 0 
            for kk = 1:algp.nciter
            
                % Calculate and apply weights to baseline images
                Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
                
                for jj = 1:size(lib.Wlib,4)
                    % Compute model values
                    Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
                        lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                    Ztot = Ztot + Z*wts(jj);
                end
                Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]);
                
                % Evaluate on pixels in exclusively body
                % Voxels in the background will not be evaluated
                fit_masks_TEs = repmat(fit_mask,1,1,1,length(scanp.tes));
                fit_masks_TEs(fit_masks_TEs>1) = 1;
                imgs_temp = imgs;
                tes_drift = scanp.tes;
        
                Ztot(fit_masks_TEs ~= 1) = 0;
                imgs_temp(fit_masks_TEs ~= 1) = 0;
                
                % Update c
                c = c_update(imgs_temp(:,:,:,scanp.echo_init_ix:end),Ztot(:,:,:,scanp.echo_init_ix:end),A,c,tes_drift(scanp.echo_init_ix:end));
                
                % Recalculate frequency shift map (assumes equal number of pixels in the x & y directions)
                Ac = reshape(A*c,scanp.dim);
            end
        end

        if Acswitch == 0
            Ac = zeros(size(Ac));
            c = zeros(size(c));
        end
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run phi_update with current f, m, Ac
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This accounts for any universal DC offsets due to reciever gain etc.
        if phiswitch
            for kk = 1:algp.nciter
                
                % Calculate and apply weights to baseline images
                Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
     
                for jj = 1:size(lib.Wlib,4)
                    % Compute model values
                    Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
                        lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                    Ztot = Ztot + Z*wts(jj);
                end
                Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]); 
    
                % Update phi
                phi = phi_update(imgs,Ztot,phi(:),scanp.tes);
    
                % Recalculate frequency shift map
                phi = reshape(phi,[scanp.dim(1) scanp.dim(2) scanp.dim(3)]);
    
            end
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run m_update with current f, Ac
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if mswitch
            for kk = 1:algp.nmiter
    
                % Build Z matrix using true values of Ac and m
                Zw = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
                Zf = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
                
                for jj = 1:size(lib.Wlib,4)
                    % Compute model values of water and fat separately
                    [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,:,jj).*exp(1i*phi),lib.Flib(:,:,:,jj).*exp(1i*phi), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ... 
                        lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                    Zw = Zw + foow*wts(jj);
                    Zf = Zf + foof*wts(jj);
                end
    
                % Update m
                m = m_update(imgs(:,:,:,scanp.echo_init_ix:end)-Zf(:,:,:,scanp.echo_init_ix:end),Zw(:,:,:,scanp.echo_init_ix:end),m,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,muscle_mask,R,Acswitch,1);
            end
        end
        
        % report cost + iteration
        costOld = cost; ii = ii + 1;
        cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp);
        converged = costOld - cost < algp.stopthresh*costOld;
    
        if algp.suppress_info ~= 1
            fprintf('        Optimization iteration %d (positivity constrain in dT): Cost = %0.2d | difference = %0.2d (%d)\n',ii,cost,cost-costOld,converged);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. start the iterative process without positivity constrain in dT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine temperature without positivity constrain
if mswitch
    for kk = 1:algp.nmiter

        % Build Z matrix using true values of Ac and m
        Zw = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
        Zf = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
        
        for jj = 1:size(lib.Wlib,4)
            % Compute model values of water and fat separately
            [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,:,jj).*exp(1i*phi),lib.Flib(:,:,:,jj).*exp(1i*phi), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ... 
                lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
            Zw = Zw + foow*wts(jj);
            Zf = Zf + foof*wts(jj);
        end

        % Update m
        m = m_update(imgs(:,:,:,scanp.echo_init_ix:end)-Zf(:,:,:,scanp.echo_init_ix:end),Zw(:,:,:,scanp.echo_init_ix:end),m,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,muscle_mask,R,Acswitch,0);
    end
end

% calculate a terrible cost that is guaranteed to get iterations going
cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp); 
cost = 2*cost;
ii = 0; % iteration counter

% The optimization algorithm is run for at least 2 iterations
converged = false;
while converged == false || ii <= 1

    iter_stop = 9; % maximum of 10 iterations
    if ii > iter_stop
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run c_update with current f, m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Evaluate in voxels where the body is present, such that fitting is insensitive to background noise
    fit_mask = fat_idx + muscle_mask; fit_mask(fit_mask > 1) = 1;
    
    if Acswitch && ii > 0
        for kk = 1:algp.nciter
        
            % Calculate and apply weights to baseline images
            Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
            
            for jj = 1:size(lib.Wlib,4)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
                    lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]);
            
            % Evaluate on pixels in exclusively body
            % Voxels in the background will not be evaluated
            fit_masks_TEs = repmat(fit_mask,1,1,1,length(scanp.tes));
            fit_masks_TEs(fit_masks_TEs>1) = 1;
            imgs_temp = imgs;
            tes_drift = scanp.tes;
    
            Ztot(fit_masks_TEs ~= 1) = 0;
            imgs_temp(fit_masks_TEs ~= 1) = 0;
            
            % Update c
            c = c_update(imgs_temp(:,:,:,scanp.echo_init_ix:end),Ztot(:,:,:,scanp.echo_init_ix:end),A,c,tes_drift(scanp.echo_init_ix:end));
            
            % Recalculate frequency shift map (assumes equal number of pixels in the x & y directions)
            Ac = reshape(A*c, scanp.dim);   
        end
    end

    if Acswitch == 0
        Ac = zeros(size(Ac));
        c = zeros(size(c));
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run phi_update with current f, m, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This accounts for any universal DC offsets due to reciever gain etc.
    if phiswitch
        for kk = 1:algp.nciter
            
            % Calculate and apply weights to baseline images
            Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
 
            for jj = 1:size(lib.Wlib,4)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
                    lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]); 

            % Update phi
            phi = phi_update(imgs,Ztot,phi(:),scanp.tes);

            % Recalculate frequency shift map
            phi = reshape(phi,[scanp.dim(1) scanp.dim(2) scanp.dim(3)]);

        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run m_update with current f, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mswitch
        for kk = 1:algp.nmiter

            % Build Z matrix using true values of Ac and m
            Zw = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
            Zf = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
            
            for jj = 1:size(lib.Wlib,4)
                % Compute model values of water and fat separately
                [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,:,jj).*exp(1i*phi),lib.Flib(:,:,:,jj).*exp(1i*phi), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ... 
                    lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
                Zw = Zw + foow*wts(jj);
                Zf = Zf + foof*wts(jj);
            end

            % Update m
            m = m_update(imgs(:,:,:,scanp.echo_init_ix:end)-Zf(:,:,:,scanp.echo_init_ix:end),Zw(:,:,:,scanp.echo_init_ix:end),m,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,muscle_mask,R,Acswitch,0);
        end
    end
    
    % report cost + iteration
    costOld = cost;ii = ii + 1;
    cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp);
    converged = costOld - cost < algp.stopthresh*costOld;

    if algp.suppress_info ~= 1
        fprintf('        Optimization iteration %d (no positivity constrain in dT): Cost = %0.2d | difference = %0.2d (%d)\n',ii,cost,cost-costOld,converged);
    end
end

% Perform temperature determination again, but now without the masks
m = m_update(imgs(:,:,:,scanp.echo_init_ix:end)-Zf(:,:,:,scanp.echo_init_ix:end),Zw(:,:,:,scanp.echo_init_ix:end),m,scanp.tes(scanp.echo_init_ix:end),algp,[],body_mask,ones(size(muscle_mask)),R,Acswitch,0);

%% Cost function
% Function to evaluate the overall cost function
function cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp)
    Ztot = zeros([scanp.dim(1) scanp.dim(2) scanp.dim(3) length(scanp.tes)]);
    for jj = 1:size(lib.Wlib,4)
        Z = calcmeimgs(lib.Wlib(:,:,:,jj),lib.Flib(:,:,:,jj), ...
            scanp.tes,scanp.b0,lib.dw0lib(:,:,:,jj)+Ac, ...
            lib.R2starlib(:,:,:,jj),m,lib.fatmodel);
        Ztot = Ztot + Z*wts(jj);
    end

    Ztot = Ztot.*repmat(exp(1i*phi),[1 1 1 size(Z,4)]);
    cost = norm((imgs(:)-Ztot(:))/algp.norm)^2;

