function s05_t1_thermometry(data_dir, opts)
% S05_T1_THERMOMETRY  Relative T1-weighted signal change dS/S0 per GRE block from multi-echo MP2RAGE data.
%
%   s05_t1_thermometry(data_dir, opts)
%
%   opts.baseline_dynamics   number of initial dynamics averaged into S0 (paper: 75 in vivo, 10 ex vivo)
%   opts.smoothing_sigma_px  Gaussian sigma in pixels on a 5x5 kernel, 0 = off (paper: 1.0 in vivo, 0 ex vivo)
%   opts.apply_correction    divide by the B1/T1 relative sensitivity map from s03 (paper: true in vivo, false ex vivo)
%   opts.s0_fit              'fixed_r2star' (default; used for the paper results) or 'wls' (per-voxel, per-dynamic
%                            weighted log-linear fit of ln S = ln S0 - R2* TE with weights S.^2)
%
% Steps: (1) water-fat separation of the baseline-averaged GRE 2 data gives a shared R2* map and the
% fat/muscle masks; (2) per block and dynamic, S0 = sum(S .* E) / sum(E.^2) with E = exp(-R2* TE)
% ('fixed_r2star'), or S0 and R2* from the weighted log-linear fit ('wls');
% (3) optional Gaussian smoothing within the muscle mask; (4) dS/S0 relative to the baseline mean,
% optionally divided by the relative sensitivity map.
% Step 1 requires the external HIMM / fat-water code (see setup_paths and THIRD_PARTY.md); its results are
% saved as derived/r2star_map.nii.gz and derived/muscle_mask_s05.nii.gz and reused when present.
%
% Input : data_dir/mp2rage_image.nii.gz [x,y,slice,echo,dynamic,block,(magnitude,phase)] + mp2rage_image.json
% Output: data_dir/derived/signal_change_img<block>.nii.gz [x,y,slice,dynamic]

scan = double(niftiread(fullfile(data_dir, 'mp2rage_image.nii.gz')));
info = jsondecode(fileread(fullfile(data_dir, 'mp2rage_image.json')));
TEs = info.imgdef.echo_time.uniq(:) * 1e-3;                     % s
[ny, nx, nz, ne, ndyn, nblocks, ~] = size(scan);
out_dir = fullfile(data_dir, 'derived');
if ~isfolder(out_dir), mkdir(out_dir); end
nbase = min(opts.baseline_dynamics, ndyn);

%% 1. Shared R2* map and muscle mask from GRE 2 (baseline average)
% Saved to derived/ and reused on later runs, so steps 2-4 run without the external water-fat code.
r2_file = fullfile(out_dir, 'r2star_map.nii.gz');
mask_file = fullfile(out_dir, 'muscle_mask_s05.nii.gz');
if isfile(r2_file) && isfile(mask_file)
    R2star = double(niftiread(r2_file));
    muscle = niftiread(mask_file) > 0;
    fprintf('s05: reusing saved R2* map and muscle mask from %s\n', out_dir);
else
    b0 = 3; frequency_system = b0 * 42.5778; prc = +1;
    waterfatparams = struct('water_shift_ppm', -0.22, 'R2star_range', [10 301 301], ...
        'off_resonance_range', [-600 600 301], 'itr', 40, 'lambda', 0.01, 'threshold', 0.01, 'num_fat_peaks', 7);
    ref = scan(:, :, :, :, :, min(2, nblocks), 1) .* exp(1i * scan(:, :, :, :, :, min(2, nblocks), 2));
    R2star = zeros(ny, nx, nz); muscle = false(ny, nx, nz);
    for z = 1:nz
        baseline_mean = mean(squeeze(ref(:, :, z, :, 1:min(75, ndyn))), 4, 'omitnan');
        [W, F, ~, R2star(:, :, z), ~] = perform_water_fat_separation(baseline_mean, info, frequency_system, b0, prc, waterfatparams);
        s = sum(squeeze(sum(abs(ref(:, :, z, :, :)), 4)) ./ max(sum(abs(ref(:, :, z, :, :)), 4), [], 'all'), 3);
        body = imclose(double(s > 0.75 * mean(s, 'all')), strel('disk', 10));
        [~, ~, m] = get_masks(body, W, F, 1, [ny nx], 0.80, 0, 0);
        muscle(:, :, z) = m > 0;
    end
    niftiwrite(single(R2star), fullfile(out_dir, 'r2star_map.nii'), 'Compressed', true);
    niftiwrite(uint8(muscle), fullfile(out_dir, 'muscle_mask_s05.nii'), 'Compressed', true);
end

if opts.apply_correction
    sens = double(niftiread(fullfile(data_dir, 'corrections', 'relative_sensitivity_map.nii.gz')));
end
E = exp(-R2star .* reshape(TEs(1:ne), 1, 1, 1, ne));

for block = 1:nblocks
    %% 2. S0 per dynamic
    mag = abs(scan(:, :, :, :, :, block, 1));
    if isfield(opts, 's0_fit') && strcmp(opts.s0_fit, 'wls')
        S0 = reshape(wls_log_fit(mag, TEs(1:ne)), ny, nx, nz, ndyn);
    else  % shared R2* map, linear least squares
        S0 = reshape(sum(mag .* E, 4) ./ sum(E.^2, 4), ny, nx, nz, ndyn);
    end

    %% 3. Optional Gaussian smoothing restricted to muscle (NaN outside)
    S0(~repmat(muscle, 1, 1, 1, ndyn)) = NaN;
    if opts.smoothing_sigma_px > 0
        S0 = smooth_in_mask(S0, muscle, opts.smoothing_sigma_px);
    end

    %% 4. Relative signal change, optional B1/T1 sensitivity correction
    dS = (S0 - mean(S0(:, :, :, 1:nbase), 4, 'omitnan')) ./ mean(S0(:, :, :, 1:nbase), 4, 'omitnan');
    if opts.apply_correction
        dS = dS ./ sens(:, :, :, block);
    end
    dS(dS == 0 | isinf(dS)) = NaN;
    niftiwrite(single(dS), fullfile(out_dir, sprintf('signal_change_img%d.nii', block)), 'Compressed', true);
end
fprintf('s05: signal change maps written to %s\n', out_dir);
end

function S0 = wls_log_fit(mag, TEs)
% Weighted least squares of ln S = ln S0 - R2* TE over the echo dimension (4), weights w = S.^2
% (the variance of ln S scales with 1/SNR^2). Closed-form 2-parameter solution per voxel and dynamic.
te = reshape(TEs, 1, 1, 1, []);
valid = mag > 0;
w = (mag.^2) .* valid;
y = log(max(mag, realmin));
Sw = sum(w, 4); Swx = sum(w .* te, 4); Swxx = sum(w .* te.^2, 4);
Swy = sum(w .* y, 4); Swxy = sum(w .* te .* y, 4);
denom = Sw .* Swxx - Swx.^2;
S0 = exp((Swxx .* Swy - Swx .* Swxy) ./ denom);
S0(~(denom > 0) | sum(valid, 4) < 2) = NaN;
end

function out = smooth_in_mask(vol, mask, sigma)
% Normalised 5x5 Gaussian convolution that ignores NaN / out-of-mask voxels.
% Identical to smoothdata2(img, 'gaussian', 5, 'omitnan') for sigma = 1 (MATLAB uses sigma = window/5).
[x, y] = meshgrid(-2:2);
k = exp(-(x.^2 + y.^2) / (2 * sigma^2));
out = vol;
for z = 1:size(vol, 3)
    for t = 1:size(vol, 4)
        img = vol(:, :, z, t);
        valid = isfinite(img);
        img(~valid) = 0;
        s = conv2(img, k, 'same') ./ conv2(double(valid), k, 'same');
        s(~mask(:, :, z)) = NaN;
        out(:, :, z, t) = s;
    end
end
end
