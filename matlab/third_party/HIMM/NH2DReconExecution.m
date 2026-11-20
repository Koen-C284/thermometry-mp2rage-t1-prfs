function [Ac,c,Ac_salomir] = NH2DReconExecution(A,Ac,Ac_salomir,Sal_mask,body_mask,brain_mask,algp,scanp)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% NH2DReconExecution: perform near harmonic 2D reconstruction
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs:               Dynamic multi-echo images (Nx x Ny x Necho)
%   m:                  frequency shift due to heating
%   A:                  basis used to model the frequency shift due to drift
%   Ac:                 frequency shift due to drift determined in the previous iteration
%   Ac_salomir:         frequency shift due to drift, estimated through
%                       near-harmonic 2D reconstruction in the previous 
%                       iteration
%   wts:                baseline weights
%   phase_motion:       frequency shift due to motion
%   phi:                Tx/Rx gain
%   Sal_mask:           fat mask
%   signal_mask:        regions with low signal mask
%   body_mask:          mask of the regions that have tissues
%   brain_mask:         brain mask
%   R:                  2nd order finite difference roughness penalty object
%   algp:               algorithm params structure
%   scanp:              scan params structure
%   lib:                baseline library structure
%
% Output: 
%   Ac:                 frequency shift due to drift, estimated using
%                       near-harmonic 2D reconstruction and QR decomposition 
%                       in the current iteration
%   c:                  spherical harmonics coefficients fit
%   Ac_salomir:         frequency shift due to drift, estimated through
%                       near-harmonic 2D reconstruction in the current 
%                       iteration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the fat mask and perform near-harmonic 2D recon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Mask based on the body mask, fat mask, and regions that have been corrected for motion
    mask_B0 = body_mask;
    Sal_mask = Sal_mask.*mask_B0;

    % Describe the background field as a harmonic function
    Ac_pre_fit = Ac;
    Ac = Salomir_fit(Ac,Ac_salomir,Sal_mask,algp.Nstd_ol,algp.neighbours,algp.neighbours);
    %Ac = Salomir_fit_Dirichelet_article(Ac,Ac_salomir,Sal_mask,algp.Nstd_ol,algp.neighbours,algp.neighbours);

    % Make sure the next near-harmonic 2D reconstruction is initialized with the previous
    Ac_salomir = Ac;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtering using a median filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Perform median filter, where values outside motion-corrected regions are excluded
%     size_kernel = algp.median_filter_size_fit;
%     Ac_premedfilter = Ac.*body_mask;%Ac_premedfilter = Ac.*brain_mask;
%     Ac_premedfilter(body_mask == 0) = nan;%Ac_premedfilter(brain_mask == 0) = nan;
%     Ac = nanmedfilt2(Ac_premedfilter,size_kernel);
%     mask_Acfit = ~isnan(Ac);
%     mask_Acfit = mask_Acfit.*mask_B0;
%     Ac(isnan(Ac)) = 0;
% 
%     figure; 
%     subplot(1,3,1); imagesc(Ac_pre_fit); axis image off; drawnow; title('Ac prior to Salomir fit');
%     %subplot(1,3,2); imagesc(Ac_premedfilter.*brain_mask); axis image off; drawnow; title('Ac corrected');
%     subplot(1,3,2); imagesc(Ac_premedfilter.*body_mask); axis image off; drawnow; title('Ac corrected');
%     subplot(1,3,3); imagesc(Ac.*mask_Acfit); axis image off; title('Ac corrected + filtered'); drawnow; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run c_fit to update c to approximate the Salomir map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parameters
    % indepvar: the independent variable is the physical mask with ones in
    %           regions where off-resonance is assumed to be estimated correctly
    % depvar:   the estimated off-resonace in the masked body region (dependent variable)
    %mask_Acfit = ~isnan(Ac);
    %mask_Acfit = mask_Acfit.*mask_B0;
    %mask_fit = mask_Acfit;
    mask_fit = brain_mask;
    indepvar = imerode(mask_fit,strel('diamond', 0));
    depvar = Ac(indepvar == 1); 

    % Perform the fit with the desired order (and degree) as provided in
    % the polynomial / spherical harmomics array (A)
    % This fit makes use of QR decomposition to decompose into the
    % basic functions of polynomials / spherical harmonics
    c = c_fit(indepvar,depvar,scanp.dim,A);
    if scanp.dim(1)==scanp.dim(2)
        Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
    else
        Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
    end
end