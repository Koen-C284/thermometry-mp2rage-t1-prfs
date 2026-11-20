function phase = phaseacc_update(imgs,Z,phase,tes,prc)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% phaseacc_update: function to fit phase accumulation
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs:  Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:     Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   phase: initial guess
%   tes:   echo times
%   prc:   direction of precession
%
% Output: 
%   phi: dim x dim fit


phase = phase(:);
innprod = Z.*conj(imgs);
innprod = permute(innprod,[3 1 2]); innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
    hesssum = hesssum + t4(:,ii)./t2(:,ii)*(tes(ii))^2;
    gradsum = gradsum - (1-2*prc)* t4(:,ii)*(tes(ii));
end

% update phase
hesssum(isnan(hesssum)) = 0;

dphase = (gradsum)./hesssum; 
dphase(isnan(dphase) | isinf(dphase)) = 0;

phase = phase + dphase;