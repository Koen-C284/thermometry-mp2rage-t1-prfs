function s03_epg_dictionary(data_dir, dataset)
% S03_EPG_DICTIONARY  EPG lookup dictionary and B1/T1-dependent relative T1 sensitivity maps.
%
%   s03_epg_dictionary(data_dir, 'invivo')
%
% Simulates the MP2RAGE signal S and its T1 derivative dS/dT1 for both GRE blocks over a T1 x B1 grid
% (sequence timing from config/sequence_parameters.json), interpolates (dS/dT1)/S at the baseline T1 and
% B1 of each voxel and normalises it to a tissue mean of 1. Stage s05 divides dS/S0 by this map.
%
% Inputs  (data_dir): t1_map.nii.gz [s], b1_map.nii.gz [% of nominal], registered to the dynamic series
% Outputs (data_dir/corrections): relative_sensitivity_map.nii.gz [x, y, slice, block], dictionaries

repo = fileparts(fileparts(mfilename('fullpath')));
cfg = jsondecode(fileread(fullfile(repo, 'config', 'sequence_parameters.json')));
seq = cfg.(dataset);
grid = cfg.dictionary_grid;
t1 = linspace(grid.T1_ms(1), grid.T1_ms(2), grid.T1_ms(3));
b1_scale = linspace(grid.B1_scale(1), grid.B1_scale(2), grid.B1_scale(3));

params = struct('TR', seq.TR_ms, 'T2', seq.T2_ms, 'PE_steps', seq.PE_steps, ...
    'number_of_reps', seq.number_of_repetitions, 'PD', 1, 'TI1', seq.TI1_ms, 'TI2', seq.TI2_ms, ...
    'MP2RAGETR', seq.MP2RAGE_TR_ms, 'alpha1', seq.alpha1_deg, 'alpha2', seq.alpha2_deg, ...
    'inv_eff', seq.inversion_efficiency, 'phase_increment', seq.rf_phase_increment_deg);

out_dir = fullfile(data_dir, 'corrections');
if ~isfolder(out_dir), mkdir(out_dir); end
S = double(mp2rage_grid_to_nifti(fullfile(out_dir, 'dictionary_signal.nii.gz'), t1, b1_scale, params));
dS = double(mp2rage_grid_to_nifti_dt1(fullfile(out_dir, 'dictionary_dsignal_dT1.nii.gz'), t1, b1_scale, params));

t1_info = niftiinfo(fullfile(data_dir, 't1_map.nii.gz'));
t1_map = double(niftiread(t1_info)) * 1000;             % s -> ms
b1_map = double(niftiread(fullfile(data_dir, 'b1_map.nii.gz')));
assert(isequal(size(t1_map), size(b1_map)), 'T1 and B1 maps must have the same size.');

rel_sens = zeros([size(t1_map, 1), size(t1_map, 2), max(1, size(t1_map, 3)), 2], 'single');
for block = 1:2
    lookup = dS(:, :, block) ./ (S(:, :, block) + eps);   % (dS/dT1) / S over the grid
    m = reshape(interp2(t1, b1_scale * 100, lookup.', t1_map(:), b1_map(:), 'linear'), size(t1_map));
    valid = (t1_map > 0) & isfinite(m);                   % outside the grid -> excluded (grey voxels)
    m = m ./ mean(m(valid));
    m(~valid) = 0;
    rel_sens(:, :, :, block) = single(reshape(m, size(rel_sens, 1:3)));
end

meta = t1_info;
meta.Datatype = 'single';
meta.ImageSize = size(rel_sens);
meta.PixelDimensions = [t1_info.PixelDimensions(1:min(end, 3)), ones(1, 4 - min(numel(t1_info.PixelDimensions), 3))];
meta.Description = 'Normalised relative T1 sensitivity (dS/dT1)/S, block 1 and 2';
niftiwrite(rel_sens, fullfile(out_dir, 'relative_sensitivity_map.nii'), meta, 'Compressed', true);
fprintf('s03: relative sensitivity map written to %s\n', out_dir);
end
