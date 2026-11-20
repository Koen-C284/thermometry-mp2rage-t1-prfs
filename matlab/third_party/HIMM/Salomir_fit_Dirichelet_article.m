function [image] = Salomir_fit_Dirichelet_article(image,prev,mask,Nstd_ol,neighbours,min_nb)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% Salomir_fit_Dirichelet_article: function that performs the near-harmonic 2D reconstruction
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Method for performing a background fit
% Adapted from Reference-Free PRFS.. Salomir et al Feb 2012
% 1) Define the starting point of iteration: 
%    - At the fat border, experimental values are used (at each iteration)
%    - At other voxels, the result of the fitting at the previous time step is used
% 2) Iterate to find background using eq. 13 from paper Salomir 
%
% Inputs:
%   time_step:              index of the image
%   image:                  2D unwrapped (subtracted) phase image
%   prev:                   result of the fit at the previous time step
%   mask:                   fat mask
%   Nstd_ol and neighbours: parameters used to filter out local extreme values in the fat border
%   min_nb:                 parameter used to eliminate voxels at the fat-tissue 
%                           interface from the fat mask) is one at voxels where the body is
%
% Outputs:
%   image:                  result of the background fit 



%% Get the border pixels
bordermask = mask - imerode(mask,strel('diamond', 1)'); % border of muscle
bordermask(1,:)=1; bordermask(end,:)=1;                 % border of image
bordermask(:,1)=1; bordermask(:,end)=1;                 % border of image


%% Start the background fitting

% Put everything within the ROI to zero
roi = mask+bordermask; roi(roi>1)=1;     % sets all values inside the roi to 1
image = image.*roi;                      % ensures all values outside the roi are 0

image(1,:)=0; image(end,:)=0;            % ensures that the edges of the image are explicitly set to zero
image(:,1)=0; image(:,end)=0;

image=image.*roi;

% Use the previous scan as initial input
image = image + prev.*(1-roi);   %prev is a function input and often squeeze(background(k-1,:,:))

% Find the size of the mask to only have for-loops over that dimensions
[xm,ym]=find(roi==0);
kmin=min(xm);kmax=max(xm);
lmin=min(ym);lmax=max(ym);


%% Save the starting image as a gif
% h=figure;
% filename = 'Dirichelet.gif';
% 
% imagesc(image)
% axis image off
% set(h,'Color',[1,1,1]);
% drawnow
% frame = getframe(1);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif', 'Loopcount',1,'DelayTime',1);

for iterations=1:200000

    oldimage=image; %image changes every iteration.  

    for k=kmin:kmax
        for l=lmin:lmax
            if roi(k,l)==0
                image(k,l)=(image(k-1,l)+image(k+1,l)+image(k,l+1)+image(k,l-1))/4 -eps/4;
                % eps=Floating-point relative accuracy=2.2204e-16  
                % This is the first order Laplacian, which should be zero because we fit spherical harmonics.
                % We fit them starting at the boudary and moving inward.
            end
        end
    end

%% Plot during fitting
%     imagesc(image)
%     axis image off
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%    %Save the images as animated gif
%     if mod(iterations,10)==0
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
%     end

    %pause(0.005)

    oldimage=oldimage-image;
    if sum(oldimage(:).^2)<1*10^-9
        disp(['Number of iterations: ',num2str(iterations)])
        break
    end
end

end