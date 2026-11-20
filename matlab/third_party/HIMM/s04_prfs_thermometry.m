%% Multi-echo MR Thermometry at 3T using near-harmonic 2D reconstruction for initialization
% 
%  HIMM: Harmonic Initialized Model-based Multi-echo  
%  Main script to extract temperature and drift fields from measured data
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% The basis of this algorithm was created by Megan Poorman & William Grissom
% (https://github.com/poormanme/waterFatSeparated_MRThermometry)
% The code provided here does not use the l1 norm for regularization, but 
% rather initializes the drift field (spatial B0 distortion) using the
% algorithm of near-harmonic 2D reconstruction (proposed by Salomir et al.)
%
% Updated from: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Do not reproduce, distribute, or modify without proper citation according
% to license files


if exist('scan_data_dir', 'var')
    saved_scan_data_dir = scan_data_dir;
    clearvars -except saved_scan_data_dir;
    scan_data_dir = saved_scan_data_dir;
    clear saved_scan_data_dir;
else
    clear all;
end
clc; close all;

%%
addpath(pwd);

%% Define scanner and sequence parameters

b0 = 3;                     % Tesla
gyro = 42.5778;             % MHz/T
frequency_system = b0*gyro; % Frequency measured by used MR system [Hz]
prc = +1;                   % Direction of precession (scanner-dependent)
                            % +1 = heating induces negative apparent freq shift
                            % -1 = heating induces positive apparent freq shift
phi0shift = 0;              % Tx/Rx gain
alpha = -0.01;              % ppm/deg C
rad2degC = prc/(2*pi*b0*alpha*gyro);  % Factor to convert radians to temperature


%% Load dynamic (baseline/heating) data
% Scan definition
sequence     = 'heating';        % 'baseline' or 'heating'
num_stacks   = 1;                % Number of transverse images acquired in a single dynamic
discard_dyns = 5;                % Number of dynamics that are discarded  ( 5 )
num_avgs     = 5;                % Number of dynamics that will be averaged ( 5 )
slices       = 1;                % ADJUST to number of slices


%% Load scan and define which (averaged) dynamics should be used as reference 
if ~exist('scan_data_dir', 'var') || isempty(scan_data_dir)
    project_root = fullfile('..', '..');
    reg_dir = fullfile(project_root, 'data', '02_interim', '03_registered');
    out_dir = fullfile(project_root, 'data', '03_processed', 'prfs_thermometry');
else
    reg_dir = fullfile(scan_data_dir, '02_interim', '03_registered');
    out_dir = fullfile(scan_data_dir, '03_processed', 'prfs_thermometry');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Specify which echoes are used for the reconstruction, and the appendix of
mp2rage_path = fullfile(reg_dir, 'mp2rage_image.nii.gz');
mp2rage_path_info = fullfile(reg_dir, 'mp2rage_image.json');

scan_raw = niftiread(mp2rage_path);
json_str = fileread(mp2rage_path_info);
GRE_info = jsondecode(json_str);
tes = GRE_info.imgdef.echo_time.uniq.*1e-3; 



scan = double(scan_raw (:,:,:,:,:,2,1) .* exp(1i* scan_raw (:,:,:,:,:,2,2)));

%% Plot images for all echoes
figure('Visible', 'off');
num_echoes = size(scan, 4);
num_col = ceil(num_echoes / 3); % Adjust number of rows based on 3 columns
max_intensity = max(squeeze(sum(abs(scan(:,:,1,:,1)), 5)),[], 'all');
for echo_ix = 1:num_echoes
    subplot(3, num_col, echo_ix);
    intensity_image = squeeze(sum(abs(scan(:,:,1,echo_ix,1)), 5));
    scaled_image = intensity_image / max_intensity; % Scale intensity to [0, 1]
    imagesc(scaled_image, [0, 1]); % Display scaled image
    axis image off;
    title(['Echo ' num2str(echo_ix)]);
    colormap gray;
end

%% Addition: remove motion
% 1. Define the directory and full file paths
disp_path = fullfile(reg_dir, 'displacement.nii.gz');

% 2. Check if the files exist
if isfile(disp_path)
    % Load the NIfTI files and convert to double for processing
    disp_data = double(niftiread(disp_path));
    X_disp_nii = disp_data(:,:,:,1);
    Y_disp_nii = disp_data(:,:,:,2);
    
    % Initialize the DisplacementField array to match your previous format
    % Format: [X, Y, 2 (X&Y), Dynamics, Slices]
    DisplacementField = zeros([size(scan,1), size(scan,2), 2, size(scan,5), max(slices)]);
    
    % 3. Loop through slices to populate the array
    for slice = slices
        % Assuming standard 4D NIfTI dimensions: [X, Y, Slices, Dynamics]
        if ndims(X_disp_nii) >= 4
            DisplacementField(:,:,1,:,slice) = reshape(X_disp_nii(:,:,slice,:), [size(scan,1), size(scan,2), 1, size(scan,5)]);
            DisplacementField(:,:,2,:,slice) = reshape(Y_disp_nii(:,:,slice,:), [size(scan,1), size(scan,2), 1, size(scan,5)]);
        else
            % Fallback just in case the NIfTI is only 3D [X, Y, Slices] with no dynamics
            DisplacementField(:,:,1,:,slice) = X_disp_nii;
            DisplacementField(:,:,2,:,slice) = Y_disp_nii;
        end
    end
    
    % 4. Compute the mask as before
    DF_mask = squeeze(sum(sum(DisplacementField,4),3)) ~= 0;
    
else
    warning('Displacement NIfTI files not found in the specified directory. Initializing zero displacement field.');
    DisplacementField = zeros([size(scan,1), size(scan,2), 2, size(scan,5), max(slices)]);
    DF_mask = true([size(scan,1), size(scan,2)]);
end

%%
% Compute temporal SD for the magnitude images te evaluate residual motion - Standard deviation: for each voxel (or pixel), you're measuring how much the signal fluctuates over time
figure('Visible', 'off'); 
body_masks = zeros([size(scan,1) size(scan,2) size(scan,3)]);
for slice_ix = 1:size(scan,3)
    abs_data_echo_all = squeeze(sum(abs(scan(:,:,slice_ix,:,:)),4));
    abs_data_echo_all_norm = abs_data_echo_all./ max(abs_data_echo_all,[],'all');
    sd_scan = std(squeeze(abs_data_echo_all_norm),0,3);
    sum_scan = sum(squeeze(abs_data_echo_all_norm),3);

    mask_temp = zeros([size(scan,1) size(scan,2)]);
    mask_temp(sum_scan > 0.50*mean(sum_scan,"all")) = 1;             % Lowered threshold to ensure full body/phantom coverage for drift field estimation
    body_masks(:,:,slice_ix) = imclose(mask_temp,strel("disk",10));

    subplot(2,size(scan,3),slice_ix); imagesc(sd_scan, [0 0.10]); axis image off; colormap gray;% colorbar; 
    subplot(2,size(scan,3),slice_ix+size(scan,3)); imagesc(body_masks(:,:,slice_ix)); axis image off; colormap gray;% colorbar; 
end


%% Identify Displacement over time

ylim_plot_DF = 3;
ylim_plot_DC = 1;

% Select motion threshold: if motion becomes larger than this value, a new reference is selected
median_smooth_size_DF = 1;    % smoothing kernel size applied on max. deformation field | 1
dx = 3;                       % span of dynamics over which the slope is calculate      | 9
motion_threshhold = 0.045;    % slope of the amount of motion in pixels                 | 0.075
slope_opening_kernel = 5;     % opening kernel size applied on thesholded slope         | 5
slope_dilation_kernel = 5;    % dilation kernel size applied on thesholded slope        | 5
slope_closing_kernel = 5;     % closing kernel size applied on thesholded slope         | 5

% Determine maximum displacement per slice in every dynamic
XY_displacement = squeeze(sqrt(DisplacementField(:,:,1,:,:).^2+DisplacementField(:,:,2,:,:).^2));


max_displacement_over_t = zeros([size(XY_displacement,3) size(XY_displacement,4)]);
for slice_ix = 1:size(XY_displacement,4)
    for dyn_ix = 1:size(XY_displacement,3)
        XY_displacement_masked = XY_displacement(:,:,dyn_ix,slice_ix);
        max_displacement_over_t(dyn_ix,slice_ix) = max(abs(XY_displacement_masked),[],'all');
    end
end

% Determine displacement per voxel in every dynamic
all_displacement_over_t = zeros([size(XY_displacement,3) sum(DF_mask,'all')]);
for dyn_ix = 1:size(XY_displacement,3)
    XY_displacement_masked = squeeze(XY_displacement(:,:,dyn_ix,:));
    XY_displacement_masked = XY_displacement_masked(DF_mask == 1);
    all_displacement_over_t(dyn_ix,:) = XY_displacement_masked;
end


% Smooth over time and determine slope voxel slice
max_displacement_over_t_smooth = zeros([size(XY_displacement,3) sum(DF_mask,'all')]);
slopes = zeros([size(XY_displacement,3) sum(DF_mask,'all')]);
for voxel_ix = 1:size(all_displacement_over_t,2)
    max_displacement_over_t_smooth(:,voxel_ix) = smoothdata(all_displacement_over_t(1:end,voxel_ix),"movmean",median_smooth_size_DF);
    x = linspace(1,size(max_displacement_over_t_smooth,1),size(max_displacement_over_t_smooth,1))';
    for dyn_ix = 1:size(max_displacement_over_t_smooth,1)
        i_start = dyn_ix - round(dx/2);
        i_end = dyn_ix + round(dx/2);
        if i_start < 1 
            i_start = 1;
        end
        if i_end > size(max_displacement_over_t_smooth,1)
            i_end = size(max_displacement_over_t_smooth,1);
        end
        slopes(dyn_ix,voxel_ix) = (max_displacement_over_t_smooth(i_end,voxel_ix)-max_displacement_over_t_smooth(i_start,voxel_ix))/dx;
    end
end

% Check which values of the slope exceed the motion threshold value and apply median filter
clear WFSBeginCondition;
exceed_motion = abs(max(slopes,[],2))>motion_threshhold;
exceed_motion_positive_slope = abs(max(slopes,[],2))>motion_threshhold;
exceed_motion_negative_slope = abs(min(slopes,[],2))>motion_threshhold;
exceed_motion_positive_slope(1:discard_dyns+1) = 1;
exceed_motion_negative_slope(1:discard_dyns+1) = 1;
exceed_motion_filtered_positive_slope = imopen(exceed_motion_positive_slope, true(slope_opening_kernel, 1));
exceed_motion_filtered_positive_slope = imdilate(exceed_motion_filtered_positive_slope, true(slope_dilation_kernel, 1));
exceed_motion_filtered_negative_slope = imopen(exceed_motion_negative_slope, true(slope_opening_kernel, 1));
exceed_motion_filtered_negative_slope = imdilate(exceed_motion_filtered_negative_slope, true(slope_dilation_kernel, 1));
exceed_motion_filtered = exceed_motion_filtered_positive_slope + exceed_motion_filtered_negative_slope;
exceed_motion_filtered(exceed_motion_filtered > 1) = 1;
exceed_motion_filtered = imclose(exceed_motion_filtered, true(slope_closing_kernel, 1));
WFSBeginCondition = [discard_dyns+1];
for dyn_ix = discard_dyns+2:size(exceed_motion,1)

    if exceed_motion_filtered(dyn_ix) == 1
        WFSBeginCondition(end+1) = dyn_ix;
    end
    %end
end

% Provide visualization of displacement over time as well as selection of references
figure('units','normalized','outerposition',[0 0 1 1], 'Visible', 'off');
ix_plots = linspace(1,size(XY_displacement,3),5);
plot_ix = 1;
for slice = slices
    for ix = 1:5
        subplot(length(slices),5,plot_ix);
        img_plot = XY_displacement(:,:,round(ix_plots(ix)),slice);
        img_plot(body_masks(:,:,slice) == 0) = nan;
        img_obj = imagesc(img_plot, [0 ylim_plot_DF]); axis image off; 
        set(img_obj, 'AlphaData', ~isnan(img_plot));
        if slice == 1
            title(['Dynamic ' num2str(round(ix_plots(ix)))],'color','w','FontSize',16);
        end
        plot_ix = plot_ix + 1;
    end
end
ax = gca; ax.FontSize = 16; ax.XColor = 'w'; set(gcf,'Color','k');
if strcmp(sequence,'heating')
    saveas(gcf,['\MotionOverTime_distributions_',sequence, ...
        '_smooth',num2str(median_smooth_size_DF), ...
        '_dx',num2str(dx), ...
        '_threshold',num2str(motion_threshhold), ...
        '_opening',num2str(slope_opening_kernel), ...
        '_dilation',num2str(slope_dilation_kernel), ...
        '_closing',num2str(slope_closing_kernel),'.svg']);
else
    saveas(gcf,['\MotionOverTime_distributions_',sequence, ...
        '_smooth',num2str(median_smooth_size_DF), ...
        '_dx',num2str(dx), ...
        '_threshold',num2str(motion_threshhold), ...
        '_opening',num2str(slope_opening_kernel), ...
        '_dilation',num2str(slope_dilation_kernel), ...
        '_closing',num2str(slope_closing_kernel),'.svg']);
end

figure('units','normalized','outerposition',[0 0 1 1], 'Visible', 'off');
yyaxis left
colors = sky(5);
for kk = 1:size(max_displacement_over_t,2)
    plot(max_displacement_over_t(:,kk),'LineWidth',4,'Color',[colors(kk,:) 0.8],'LineStyle','-','Marker', 'none'); hold on;
end
xlim([0 size(max_displacement_over_t,1)]); ylim([0 ylim_plot_DF]);
ax = gca; ax.FontSize = 24; ax.XColor = 'w'; set(gcf,'Color','k');

yyaxis right
plot(max(slopes,[],2),'LineWidth',4,'Color',"#D95319",'LineStyle','-'); hold on;
plot(min(slopes,[],2),'LineWidth',4,'Color',"#D95319",'LineStyle','-'); hold on;
xlim([0 size(max_displacement_over_t,1)]); ylim([-ylim_plot_DC-0.01 ylim_plot_DC+0.01]);
ax = gca; ax.FontSize = 24; ax.XColor = 'w'; set(gcf,'Color','k');% ax.YColor = 'w'; 

yyaxis left
for ix = 0:discard_dyns
    plot([ix,ix],[-0.5,ylim_plot_DF+0.5], ...
        'LineWidth',4, 'Color',[0 0 0 0.5],'LineStyle','-'); hold on; 
end
for ix = 1:length(WFSBeginCondition)
    plot([WFSBeginCondition(ix),WFSBeginCondition(ix)],[-0.5,ylim_plot_DF+0.5], ...
        'LineWidth',4, 'Color',[0 0 0 0.5],'LineStyle','-'); hold on; 
end
if strcmp(sequence,'heating')
    saveas(gcf,['\MotionOverTime_plots_',sequence, ...
        '_smooth',num2str(median_smooth_size_DF), ...
        '_dx',num2str(dx), ...
        '_threshold',num2str(motion_threshhold), ...
        '_opening',num2str(slope_opening_kernel), ...
        '_dilation',num2str(slope_dilation_kernel), ...
        '_closing',num2str(slope_closing_kernel),'.svg']);
else
    saveas(gcf,['\MotionOverTime_plots_',sequence, ...
        '_smooth',num2str(median_smooth_size_DF), ...
        '_dx',num2str(dx), ...
        '_threshold',num2str(motion_threshhold), ...
        '_opening',num2str(slope_opening_kernel), ...
        '_dilation',num2str(slope_dilation_kernel), ...
        '_closing',num2str(slope_closing_kernel),'.svg']);
end

% Loop over begin conditions and only select the ones which do not have a new reference in the very next dynamic
if length(WFSBeginCondition) > 1
    WFSBeginCondition_first = [WFSBeginCondition(1)];
    WFSBeginCondition_last = [];
    for ix = 2:length(WFSBeginCondition)
        if WFSBeginCondition(ix)-1 ~= WFSBeginCondition(ix-1)
            WFSBeginCondition_first(end+1) = WFSBeginCondition(ix);
        end
    end
    for ix = 1:length(WFSBeginCondition)-1
        if WFSBeginCondition(ix)+1 ~= WFSBeginCondition(ix+1)
            WFSBeginCondition_last(end+1) = WFSBeginCondition(ix);
        end
    end
    WFSBeginCondition_last(end+1) = WFSBeginCondition(end);
else
    WFSBeginCondition_first = WFSBeginCondition(1);
    WFSBeginCondition_last = WFSBeginCondition(1);
end

% Determine how Water/Fat separation partitions are made based on averaging
clear WFlib_partitions;
for condition = 1:length(WFSBeginCondition)
    initIx = 1+num_stacks*(0+WFSBeginCondition(condition)-1);
    WFlib_partitions(:,condition) = reshape(initIx:initIx+num_stacks*num_avgs-1,num_stacks*num_avgs,[]);
end
WFlib_partitions(WFlib_partitions>size(scan,5)) = size(scan,5);

%% Express the first dynamic(s) into water and fat components (W, F, R2*, and dw0)

% Set water/fat recon parameters
waterfatparams.water_shift_ppm = -0.22;               % Shift of the water peak with respect to 4.7 ppm (check Ren et al, JLR 2008) | -0.22 seems to work
waterfatparams.R2star_range = [10 301 301];            % [minimum | maximum | number_for_quantization]
waterfatparams.off_resonance_range = [-600 600 601];  % [minimum | maximum | number_for_quantization]
waterfatparams.itr = 40;                              % Number of graph cut iterations
waterfatparams.lambda = 0.01;                         % Regularization parameter (smoothening of static off-resonance)
waterfatparams.threshold = 0.01;
waterfatparams.num_fat_peaks = 7;                 	  % Choose an x-peak model, options are:
                                                      % [ 1 | 3 | 4 | 5 | 6 | 7 | 9 ]

%% Apply water/fat separation to the first dynamic of every stack
clear Wlib Flib dw0lib R2starlib

for slice = 1:size(scan,3)
    disp(['Performing W/F separation for slice ' num2str(slice) ' out of ' num2str(size(scan,3))]);
    for REPPConditionIx = 1:size(WFlib_partitions,2)
        disp(['   reference ' num2str(REPPConditionIx) ' out of ' num2str(size(WFlib_partitions,2))]);
        [Wlib(:,:,slice,REPPConditionIx),Flib(:,:,slice,REPPConditionIx),dw0lib(:,:,slice,REPPConditionIx),R2starlib(:,:,slice,REPPConditionIx),fatmodel] = perform_water_fat_separation( ...
            mean(squeeze(scan(:,:,slice,:,WFlib_partitions(:,REPPConditionIx)')),4),GRE_info,frequency_system,b0,prc,waterfatparams);
    end
end


%% Plot the fat figures and percentages
figure('Visible', 'off');

% Calculate fat percentage
fat_percentage = mean(abs(Flib), 3) ./ (mean(abs(Wlib), 3) + mean(abs(Flib), 3)) * 100;

% Create a 1 x 3 subplot layout
subplot(1, 3, 1);
imagesc(mean(abs(Wlib(:,:,1)), 3)); % Average over slices for visualization using absolute values
title('Water Component (Wlib)');
xlabel('Readout');
ylabel('Phase Encoding');
colorbar;
axis equal tight; % Ensure equal scaling for the axes and make axis tight

% Plot Flib as a subplot
subplot(1, 3, 2);
imagesc(mean(abs(Flib(:,:,1)), 3)); % Average over slices for visualization using absolute values
title('Fat Component (Flib)');
xlabel('Readout');
ylabel('Phase Encoding');
colorbar;
axis equal tight; % Ensure equal scaling for the axes and make axis tight

% Plot fat percentage as a subplot
subplot(1, 3, 3);
imagesc(fat_percentage(:,:,1).*body_masks); % Display fat percentage
title('Fat Percentage');
xlabel('Readout');
ylabel('Phase Encoding');
colorbar;
axis equal tight; % Ensure equal scaling for the axes and make axis tight
%% Masking - use manually created body mask to generate fat, signal, and muscle masks 
% The drift field will be initialized using near harmonic 2D reconstruction, 
% this algorithm requires a mask of fatty regions from which the harmonic
% reconstruction can be initialized.
% In addition to the fat mask, we also create masks for the body and for
% regions with negligable signal
% Three parameters are to be optimized for the fat, body and muscle masks:
th_fat = 0.60;                   % Higher value results in smaller 'fat' mask (0.80)
thickness_fat = 0;               % Dilation value that determines the fay layer thickness (from edge to inside)  (if 0, subcutaneous fat is also used)
body_erode = 0;                  % Integer value that determined how much the body mask is eroded (0)

% Compute the fat, body and no signal mask with current parameters
figure('Visible', 'off'); 
for slice = 1:size(scan,3)
    [body_mask(:,:,slice),fat_mask(:,:,slice),muscle_mask(:,:,slice)] = get_masks(body_masks(:,:,slice), Wlib(:,:,slice,1), Flib(:,:,slice,1), 1, [size(scan,1) size(scan,2)], th_fat, thickness_fat, body_erode);

    % Provide a visualization
    subplot(size(scan,3),3,(slice-1)*3+1); imagesc(fat_mask(:,:,slice)); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); %title('fat mask');
    subplot(size(scan,3),3,(slice-1)*3+2); imagesc(body_mask(:,:,slice)); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); %title('body mask');
    subplot(size(scan,3),3,(slice-1)*3+3); imagesc(muscle_mask(:,:,slice)); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); %title('muscle mask');
end

%% Currently, no averaging is applied

scan_T_avg = scan;

%% Algorithm parameters

% Only algorithm parameters may be altered, which does not hold for scanner 
% and sequence since those parameters are constant

% Settings for drift, temperature, and motion field fit
algp.Acswitch = 1;                % whether drift will be fitted for or not
algp.order = 3;                   % polynomial / spherical background phase shift order 
algp.method = 'spher';             % 'poly' for polynomial and 'spher' for spherical harmonics
algp.positivity_constrain = 0;    % have some iterations with a positivity constrain in temperature

% Parameter structure: algorithm parameters
algp.num_salomir_fits = 1;        %1
algp.nmiter = 5;                  % # m iterations
algp.nciter = 5;                  % # c iterations
algp.stopthresh = 1e-5;           % algorithm stopping threshold (lowered from 1e-4 for precise drift field prediction)
algp.rad2degC = rad2degC;         % Conversion between degree C and radians
algp.suppress_info = 0;           % Whether command window information is shown per iteration or not
initWithPreviousEstimate = 0;     % Parameter to speed up the algorithm

% Parameter structure: algorithm masking parameters
algp.th_fat = th_fat;               % threshold for fat definition for when l1 norm is only calculated in fats
algp.thickness_fat = thickness_fat; % maximum thickness of the fat layer in voxels
algp.body_erode = body_erode;       % Value used to make the body mask slightly smaller
algp.Nstd_ol = 2;                   % Filter out local extreme values in the fat border
algp.neighbours = 2;                % Define local neighbourhood
algp.min_nb = 0;                    % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)

% Parameter structure: scan parameters
scanp.dim = [size(scan,1) size(scan,2)];
scanp.b0 = b0;
scanp.tes = tes;
scanp.prc = prc;
% Determine slice thickness field name and assign to scanp.slice_th
if isfield(GRE_info.imgdef,'slice_thickness') && isfield(GRE_info.imgdef.slice_thickness,'uniq')
    scanp.slice_th = GRE_info.imgdef.slice_thickness.uniq;
elseif isfield(GRE_info.imgdef,'slice_thickness_in_mm') && isfield(GRE_info.imgdef.slice_thickness_in_mm,'uniq')
    scanp.slice_th = GRE_info.imgdef.slice_thickness_in_mm.uniq;
else
    error('slice thickness field not found in GRE_info.imgdef');
end


scanp.voxel_size = GRE_info.imgdef.pixel_spacing_x_y.uniq;

%% Algorithm parameters

% Only algorithm parameters may be altered, which does not hold for scanner 
% and sequence since those parameters are constant

% % Settings for drift, temperature, and motion field fit
% algp.Acswitch = 1;                % whether drift will be fitted for or not
% algp.order = 5;                   % polynomial / spherical background phase shift order 
% algp.method = 'poly';             % 'poly' for polynomial and 'spher' for spherical harmonics
% algp.positivity_constrain = 0;    % have some iterations with a positivity constrain in temperature
% 
% % Parameter structure: algorithm parameters
% algp.num_salomir_fits = 1;   %1
% algp.nmiter = 5;                  % # m iterations
% algp.nciter = 5;                  % # c iterations
% algp.stopthresh = 1e-5;           % algorithm stopping threshold
% algp.rad2degC = rad2degC;         % Conversion between degree C and radians
% algp.suppress_info = 0;           % Whether command window information is shown per iteration or not
% initWithPreviousEstimate = 0;     % Parameter to speed up the algorithm
% 
% % Parameter structure: algorithm masking parameters
% algp.th_fat = th_fat;               % threshold for fat definition for when l1 norm is only calculated in fats
% algp.thickness_fat = thickness_fat; % maximum thickness of the fat layer in voxels
% algp.body_erode = body_erode;       % Value used to make the body mask slightly smaller
% algp.Nstd_ol = 2;                   % Filter out local extreme values in the fat border
% algp.neighbours = 2;                % Define local neighbourhood
% algp.min_nb = 0;                    % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)
% 
% % Parameter structure: scan parameters
% scanp.dim = [size(scan,1) size(scan,2)];
% scanp.b0 = b0;
% scanp.tes = tes;
% scanp.prc = prc;
% % Determine slice thickness field name and assign to scanp.slice_th
% if isfield(GRE_info.imgdef,'slice_thickness') && isfield(GRE_info.imgdef.slice_thickness,'uniq')
%     scanp.slice_th = GRE_info.imgdef.slice_thickness.uniq;
% elseif isfield(GRE_info.imgdef,'slice_thickness_in_mm') && isfield(GRE_info.imgdef.slice_thickness_in_mm,'uniq')
%     scanp.slice_th = GRE_info.imgdef.slice_thickness_in_mm.uniq;
% else
%     error('slice thickness field not found in GRE_info.imgdef');
% end
% 
% 
% scanp.voxel_size = GRE_info.imgdef.pixel_spacing_x_y.uniq;
%% Perform iterative optimization to separate the temperature and drift fields
clear lib cfunc mfunc Acfunc wts phi PhaseAcc mfuncComb AcfuncLib AcfuncLib phase_motionLib;
disp('Optimization of multi-echo fat-suppressed MR thermometry has started..')

% Provide which dynamics are to be implemented in the optimization
show_ix = size(scan,5)-discard_dyns;   % number of dynamics that will be evaluated | size(scan,5)-discard_dyns
init_ix = discard_dyns + 1;            % first dynamic to be evaluated | discard_dyns + 1
step = 1;                              % step size in dynamics

% Initialize temperature and drift to zeros in the discarded dynamics
mfunc = zeros([scanp.dim size(slices,2) discard_dyns]);
Acfunc = zeros([scanp.dim size(slices,2) discard_dyns]);
if strcmp(algp.method,'poly')
    A = get_polynomial_function(scanp.dim,algp.order);
elseif strcmp(algp.method,'spher')
    A = get_spherical_function(scanp.dim,algp.order);
end
cfunc = zeros([size(A,2) size(slices,2) discard_dyns]);

% Update median filter range (prior to near harmonic 2D reconstruction)
algp.median_filter_size = 5;

% Update echo times
scanp.echo_init_ix = 1; % 1 use all echoes; 2 use 2nd echo onwards; 3 use 3rd echo onwards; etc. ( 1 )
scanp.tes = tes;

tic 
for ii = init_ix:step:init_ix+show_ix-1

    % Check to which REPP condition the current dynamic corresponds
    REPPCondition_dyns = WFlib_partitions(1,:); % WFlib_partitions(1,:) | REPPCondition*num_stacks-num_stacks+1
    REPPCondition_dyns(1) = 1;
    if ii == 1
        REPPCondition_idx = 1;
    elseif size(scan,4) == size(scan_T_avg,4)
        REPPCondition_idx = find(ii>=REPPCondition_dyns);
    else
        REPPCondition_idx = find(collapse_partitions(ii,1)>=REPPCondition_dyns);
    end
    REPPCondition_idx = REPPCondition_idx(end);

    disp(['dyn ' num2str(ii) ' - W/F separated reference number ' num2str(REPPCondition_idx)]);
    disp(['    W/F reference indices =  ' num2str(REPPCondition_dyns)]);

    for slice = 1:size(scan,3)

        % Evaluate a single dynamic
        if scanp.prc == 1
            dyn_img = squeeze(scan_T_avg(:,:,slice,:,ii));
        else
            dyn_img = conj(squeeze(scan_T_avg(:,:,slice,:,ii)));
        end
        
        % Scans to base reference image upon
        if scanp.prc == 1
            lib.drift_lib.imgs_ref = mean(squeeze(scan(:,:,slice,1,WFlib_partitions(:,REPPCondition_idx)')),4);
        else
            lib.drift_lib.imgs_ref = conj(mean(squeeze(scan(:,:,slice,1,WFlib_partitions(:,REPPCondition_idx)')),4));
        end
    
        % Update W/F separated reference data onto the parameter structures
        lib.Wlib = Wlib(:,:,slice,REPPCondition_idx);
        lib.Flib = Flib(:,:,slice,REPPCondition_idx);
        lib.dw0lib = dw0lib(:,:,slice,REPPCondition_idx);
        lib.R2starlib = R2starlib(:,:,slice,REPPCondition_idx);
        lib.fatmodel = fatmodel;
        algp.body_mask = body_masks(:,:,slice);

        % Use previously predicted drift field coefficients as initialization for the current dynamic
        disp(['    slice ' num2str(slice) ' out of ' num2str(size(scan,3))]);
        lib.dyn = ii;
        if ii > init_ix
            lib.cinit = squeeze(cfunc(:,slice,ii-step));
        end
        if (isfield(lib,'cinit')) && (ii == WFlib_partitions(1,REPPCondition_idx))
            lib.cinit = zeros(size(cfunc(:,slice,ii-step)));
        end
        if initWithPreviousEstimate && ii == init_ix
            % Optimization of the first dynamic
            [mfunc(:,:,slice,ii),wts(:,slice,ii),Acfunc(:,:,slice,ii),cfunc(:,slice,ii),A,phi(:,:,slice,ii),Ac_phase_accumulation(:,:,slice,ii)] = ...
                OptimizeHIMMBodyNew(dyn_img, algp, scanp, lib,algp.Acswitch,1,0,zeros(scanp.dim));
        elseif initWithPreviousEstimate && ii > init_ix
            % Optimization of remainder dynamics (with nonzero initialization of mfunc and drift)
            [mfunc(:,:,slice,ii),wts(:,slice,ii),Acfunc(:,:,slice,ii),cfunc(:,slice,ii),A,phi(:,:,slice,ii),Ac_phase_accumulation(:,:,slice,ii)] = ...
                OptimizeHIMMBodyNew(dyn_img, algp, scanp, lib,algp.Acswitch,1,0,mfunc(:,:,slice,ii-step));
        else
            % Phi switch off, m switch on (solve for temperature without Tx/Rx gain)
            [mfunc(:,:,slice,ii),wts(:,slice,ii),Acfunc(:,:,slice,ii),cfunc(:,slice,ii),A,phi(:,:,slice,ii),Ac_phase_accumulation(:,:,slice,ii)] = ...
                        OptimizeHIMMBodyNew(dyn_img, algp, scanp, lib,algp.Acswitch,1,0);
            %Optimize_HIMM_only_better(dyn_img, algp, scanp, lib,algp.Acswitch,1,0);

            
        end
    
        % Combine data of current reference to the one of the previous
        mfuncLib(:,:,slice,ii) = mfunc(:,:,slice,ii);
        AcfuncLib(:,:,slice,ii) = Acfunc(:,:,slice,ii);

        REPPCondition_ref_idx = find(ii>=WFSBeginCondition_first);
        if ~isempty(REPPCondition_ref_idx)
            REPPCondition_ref_idx = REPPCondition_ref_idx(end);
        else
            REPPCondition_ref_idx = 0;
        end
        if REPPCondition_ref_idx > 1
            idx_first = WFSBeginCondition_first(REPPCondition_ref_idx)-5;
            idx_last = WFSBeginCondition_first(REPPCondition_ref_idx)-1;
            if idx_first < 1
                idx_first = 1;
            end
            mfunc(:,:,slice,ii) = mfunc(:,:,slice,ii) + mean(mfunc(:,:,slice,idx_first:idx_last),4);
            Acfunc(:,:,slice,ii) = Acfunc(:,:,slice,ii) + mean(Acfunc(:,:,slice,idx_first:idx_last),4);
        end

    end

    % Give visual update during optimization
    % close all; 
    % figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4]);
    % for slice = 1:size(scan,3)
    %     subplot(1,size(scan,3)+1,slice); imagesc(mfunc(:,:,slice,ii)*algp.rad2degC.*muscle_mask(:,:,slice), [-1.5 1.5]); axis image; axis off; colorbar;
    %    %subplot(1,size(scan,3)+1,size(scan,3)+1); plot(1:ii,squeeze(mfunc(110,50,slice,:)*algp.rad2degC)); axis square; ylim([-2.0 2.0]); xlim([1 size(scan,5)]); hold on;
    %     subplot(1,size(scan,3)+1,size(scan,3)+1); plot(1:ii,squeeze(mfunc(35,24,slice,:)*algp.rad2degC)); axis square; ylim([-2.0 2.0]); xlim([1 size(scan,5)]); hold on;
    % end
    % drawnow;

end
toc
%save('optimized_data_pillow_TandChi.mat','mfunc','Acfunc');


%% Reshape A

% reshaped_A = reshape(A(:,size(A,2)), [size(scan,1),size(scan,2),length(slices)]);
% Acfunc = reshape(Acfunc, [size(Acfunc, 1), size(Acfunc, 2), 1,size(Acfunc, 3)]);

%% Quick analysis
% rad2Hz = 1 / (2 * pi);
% figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4]);
% if size(scan, 3) > 1
%     for slice = 1:size(scan,3)
%         subplot(1,size(scan,3)+1,slice); imagesc(Acfunc(:,:,slice,ii)*rad2Hz.*body_mask(:,:,slice), [3 4]); axis image; axis off; colorbar;
%         subplot(1,size(scan,3)+1,size(scan,3)+1); plot(1:ii,squeeze(Acfunc(y_visualize,x_visualize,slice,:)*rad2Hz)); axis square; xlim([1 size(scan,5)]); hold on;
%     end
% else
%     % 2D case where there is no slice
%     subplot(1,2,1); imagesc(Acfunc(:,:,slice,ii)*rad2Hz.*body_mask(:,:), [3 4]); axis image; axis off; colorbar;
%     subplot(1,2,2); plot(1:ii,squeeze(Acfunc(y_visualize,x_visualize,slice,:)*rad2Hz)); axis square; xlim([1 size(scan,5)]); hold on;
% end
% drawnow;
%% 2D plot
figure('Visible', 'off');
p1y = 37; p1x = 24;
p2y = 40; p2x = 37;
subplot(1,2,1); 
dyn_preview_idx = min(250, size(mfunc, 4));
imagesc(mean(squeeze(mfunc(:,:,slice, dyn_preview_idx)),3)*rad2degC.*muscle_mask(:,:,1), [-1.0 2.5]); 
axis image; axis off; colorbar; hold on
plot3(p1x, p1y, 0, 'marker', 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0, 0.4470, 0.7410]);
plot3(p2x, p2y, 0, 'marker', 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.8500, 0.3250, 0.0980]);
subplot(1,2,2); 
plot(squeeze(mfunc(p1y,p1x,slice,:))*rad2degC, 'LineWidth', 2); hold on;
plot(squeeze(mfunc(p2y,p2x,slice,:))*rad2degC, 'LineWidth', 2); hold on;
ylim([-1 2.5]); xlim([1 size(scan,5)]); axis square;
%% Save data as .mat file and as NIfTI files

%% Save data as .mat file and as NIfTI files
project_root = '../../'; 
output_dir      = fullfile(project_root, 'data', '03_processed', 'prfs_thermometry');
for slice = slices

    % Create output directory
    img_name_append = strcat('_bspline90_registered_slice',num2str(slice-1));
    output_dir = ['\Optimization_multiref_poly3\',sequence,img_name_append];
    if ~exist(output_dir, 'dir')    % Check if directory exists and create if not
       mkdir(output_dir)
    end
    
    % Save the optimized temperature and drift together with the body and fat mask as NIfTI files
    niftiwrite(squeeze(mfunc(:,:,slice,:)),[output_dir,'\recon_temperature.nii']);      	   % Reconstructed temperature
    niftiwrite(squeeze(Acfunc(:,:,slice,:)),[output_dir,'\recon_drift.nii']);                  % Reconstructed drift
    niftiwrite(double(squeeze(body_mask(:,:,slice))),[output_dir,'\opt_body_mask.nii']);       % Body mask (first dimension is z-direction)
    niftiwrite(double(squeeze(fat_mask(:,:,slice))),[output_dir,'\opt_fat_mask.nii']);         % Fat mask (first dimension is z-direction)
    niftiwrite(double(squeeze(muscle_mask(:,:,slice))),[output_dir,'\opt_muscle_mask.nii']);   % Muscle mask (first dimension is z-direction)
    
end

% Save the relevant parameters as .mat file in the output directory
img_name_append = strcat('_bspline90_registered_mat');
output_dir = ['\Optimization_multiref_poly3\',sequence,img_name_append];
if ~exist(output_dir, 'dir')    % Check if directory exists and create if not
   mkdir(output_dir)
end
save([output_dir,'\optimized_data.mat'],'mfunc','Acfunc','mfuncLib','AcfuncLib','wts','scanp','lib','rad2degC');




%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
temp_prfs = mfunc .*rad2degC.*muscle_mask;
temp_prfs(temp_prfs == 0) = NaN;
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Define output filenames (use gzip NIfTI extension)
file_temp_prfs       = fullfile(out_dir, 'temp_prfs.nii.gz');       % temperature PRFs (rad2degC applied)
file_Acfunc          = fullfile(out_dir, 'Acfunc.nii.gz');          % reconstructed drift
file_mfunc           = fullfile(out_dir, 'mfunc.nii.gz');           % reconstructed temperature
file_Wlib            = fullfile(out_dir, 'Wlib.nii.gz');            % library W (mean of complex values saved as real/imag channels)
file_Flib            = fullfile(out_dir, 'Flib.nii.gz');            % library F (mean of complex values saved as real/imag channels)
file_R2star          = fullfile(out_dir, 'R2starlib.nii.gz');       % R2* library (mean across 4th dim)
file_off_resonance   = fullfile(out_dir, 'w0lib.nii.gz');           % off-resonance (dw0) library (mean across 4th dim)
file_tes             = fullfile(out_dir, 'tes.nii.gz');             % echo times
file_scan            = fullfile(out_dir, 'scan.nii.gz');            % original scan volume
file_masks           = fullfile(out_dir, 'masks.nii.gz');
% Notes:
% - Only write files for variables that exist in the workspace below.
% - Complex library volumes (Wlib, Flib) are averaged across the 4th dimension,
%   then split into real and imaginary components and stacked along the 4th dim
%   so they can be written as real-valued NIfTI files.

if exist('temp_prfs','var')
    niftiwrite(single(temp_prfs), file_temp_prfs);
end
if exist('Acfunc','var')
    niftiwrite(single(Acfunc), file_Acfunc);
end
if exist('mfunc','var')
    niftiwrite(single(mfunc), file_mfunc);
end
if exist('Wlib','var')
    Wlib_mean = single(mean(Wlib, 4));
    Wlib_stacked = cat(4, real(Wlib_mean), imag(Wlib_mean));
    niftiwrite(Wlib_stacked, file_Wlib);
end
if exist('Flib','var')
    Flib_mean = single(mean(Flib, 4));
    Flib_stacked = cat(4, real(Flib_mean), imag(Flib_mean));
    niftiwrite(Flib_stacked, file_Flib);
end

if exist('R2starlib','var')
    niftiwrite(single(mean(R2starlib, 4)), file_R2star);
end
if exist('dw0lib','var')
    niftiwrite(single(mean(dw0lib, 4)), file_off_resonance);
end

if exist('tes','var')
    niftiwrite(single(tes), file_tes);
end

if exist('scan','var')
    scan_stacked = cat(6, real(scan), imag(scan));
    niftiwrite(single(scan_stacked), file_scan);
end

if exist('body_masks','var') && exist('muscle_mask','var') && exist('fat_mask','var')
    % Stack body, muscle and fat masks along the 4th dimension and write as single
    masks_stack = cat(3, single(body_masks), single(muscle_mask), single(fat_mask));
    niftiwrite(masks_stack, fullfile(out_dir, 'masks.nii.gz'));
    
end
