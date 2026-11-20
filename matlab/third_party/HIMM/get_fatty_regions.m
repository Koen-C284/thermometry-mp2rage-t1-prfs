function [body_mask,mask_F,mask_no_signal,muscle_mask] = get_fatty_regions(manual_body_mask, Wlib, Flib, wts, dims, th_f_w, thickness_fat, th_no_signal, body_erode, magn_img, factor_noise)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% get_fatty_regions: function to create the masks used in the HIMM algorithm
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   Wlib:           water image
%   Flib:           fat image
%   wts:            coefficients to combine the library imaging into the image 
%                   best suited to the current dynamic
%   dims:           dimensions of the current data
%   th_f_w:         threshold applied to the relative water/fat image used to 
%                   construct the fat mask
%   th_no_signal:   list with 2 values used to construct the low signal mask 
%                   (through thresholding)
%   SE: structure   element used to erode make the fat mask
%   body_factor:    factor used on top of the Otsu's threshold to define the body mask
%   body_erode:     factor used to construct the structure element to erode the body mask
%
% Output: 
%   mask_F:         fat mask
%   body_mask:      mask of the regions that have tissues
%   mask_no_signal: regions with low signal have value 1

    % Body segmentation of pixels with high Fat / Water content
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify region that can be seen as being human body (separate background)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create 'in phase' (IP) image by summing water and fat images over the image library
    IP = zeros([dims(1) dims(2)]);
    for jj = 1:size(Wlib,3)
        IP = IP + Wlib(:,:,jj)*wts(jj) + Flib(:,:,jj)*wts(jj);
    end

    % Otsu's method
    IP = mat2gray(abs(IP));             % Rescale to [0 1] range
    [th_body] = multithresh(IP,2);      % Use Otsu's method to select a number of thresholds

    body_mask_threshold = zeros(size(IP));
    body_mask_threshold(IP > 0.5*mean(IP,'all')) = 1;
    body_mask_threshold = imclose(body_mask_threshold,strel("disk",5));
    body_mask = body_mask_threshold.*manual_body_mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get fat & water content and apply processing to get fatty region mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    water_content = zeros([dims(1) dims(2)]);
    fat_content = zeros([dims(1) dims(2)]);
    for jj = 1:size(Wlib,3)
        water_content = water_content + abs(Wlib(:,:,jj))*wts(jj);
        fat_content = fat_content + abs(Flib(:,:,jj))*wts(jj);
    end

    % Get fat/water content
    rel_fw = abs(fat_content) ./ (abs(water_content) + abs(fat_content)); 

    % Post processing
    body_mask = imerode(body_mask,strel('diamond', body_erode));
    rel_fw(isnan(rel_fw)) = 0;          % Set NaN values to background
    rel_fw(body_mask ~= 1) = 0;         % Set all non-human-body values to background
    mask_F = zeros([dims(1) dims(2)]);  % Create an empty mask that is to be filled with fatty regions
    mask_F(rel_fw > th_f_w) = 1;        % All pixels with fat/water content above threshold are set to 1 
    %mask_F = imdilate(mask_F,strel('disk',2));

%     % Remove masked pixels in thick regions and keep in thin regions
%     mask_F = bwareaopen(mask_F, ceil(dims(1)*dims(2)*0.005)); % Remove objects containing fewer than dim(1)*dim(2)*0.005 pixels
%     mask_F_erode = imerode(mask_F, SE);         % Removal of superficial layers
%     mask_F_dilate = imdilate(mask_F_erode, SE);	% Grow regions back
%     mask_F_thin = mask_F - mask_F_dilate;       % Regions that did not grow back are thin
%     mask_F = mask_F_erode + mask_F_thin;        % Keep the eroded and thin regions 
    
    % Get a mask of regions that do not provide (a lot of) signal
    water_scale = mat2gray(abs(water_content)); % Rescale to [0 1] range
    fat_scale = mat2gray(abs(fat_content));     % Rescale to [0 1] range
    mask_no_signal = zeros([dims(1) dims(2)]);
    mask_no_signal((water_scale < th_no_signal(1)) & (fat_scale < th_no_signal(2))) = 1;
    mask_no_signal(body_mask == 0) = 1;

    % Make muscle mask (using Otsu's method) and change fat mask accordingly
    fat_temp_mask = zeros([dims(1) dims(2)]);
    muscle_mask = zeros([dims(1) dims(2)]);
    F = mat2gray(abs(fat_content));       % Rescale to [0 1] range
    W = mat2gray(abs(water_content));       % Rescale to [0 1] range
    [th_fat] = multithresh(F,2);         % Use Otsu's method to select a number of thresholds
    [th_muscle] = multithresh(W,2);         % Use Otsu's method to select a number of thresholds
    fat_temp_mask(F > 2.0*th_fat(1)) = 1;
    muscle_mask(W > 1.5*th_muscle(1)) = 1;
    muscle_mask(mask_no_signal == 1) = 0;
    mask_F(mask_no_signal == 1) = 0;
    mask_F(muscle_mask == 1) = 0;
    muscle_mask(fat_temp_mask == 1) = 0;

    % Only select the outer rim of fat 
    % If there is a vessel inside the fat, it might cause drift errors due to cardiac pulsation
    body_dilated = imdilate(body_mask,strel('disk',2));
    body_dilated = abs(body_mask - body_dilated);
    edge_body = imdilate(body_dilated,strel('disk',2 + thickness_fat));
    edge_body(body_mask == 0) = 0;

    muscle_mask_dilated = imdilate(muscle_mask,strel('disk',6));
    body_mask_eroded = imerode(body_mask,strel('disk',20));
    bone_mask = mask_F;
    bone_mask(muscle_mask_dilated == 1) = 0;
    bone_mask(body_mask_eroded == 0) = 0;

    %mask_F(edge_body == 0) = 0;
    %mask_F = mask_F + bone_mask;
    mask_F(mask_F > 1) = 1;

end