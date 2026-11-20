function s01_compute_maps(raw_dir, out_dir)
% S01_COMPUTE_MAPS  Convert Philips PAR/REC scans to NIfTI and compute the baseline T1 and B1 maps.
%
%   s01_compute_maps(raw_dir, out_dir)
%
% raw_dir must contain ll_data.par/.rec (Look-Locker), afi_data.par/.rec (AFI B1) and mp2rage_data.par/.rec.
% Writes mp2rage_image, t1_map(_image) and b1_map(_image) NIfTI files plus a JSON sidecar that keeps only
% the acquisition parameters used downstream (no patient or examination identifiers are written).
% Requires loadPARREC on the path (see THIRD_PARTY.md).

if ~isfolder(out_dir), mkdir(out_dir); end
[img_ll, info_ll] = load_mri_data(fullfile(raw_dir, 'll_data.par'));
[img_afi, info_afi] = load_mri_data(fullfile(raw_dir, 'afi_data.par'));
[img_mp2rage, info_mp2rage] = load_mri_data(fullfile(raw_dir, 'mp2rage_data.par'));

[t1_map, t1_map_image] = compute_T1(img_ll, info_ll);
[b1_map, b1_map_image] = compute_B1(img_afi, info_afi);

save_nii = @(x, name) niftiwrite(single(x), fullfile(out_dir, [name '.nii']), 'Compressed', true);
save_nii(img_mp2rage, 'mp2rage_image');
save_nii(t1_map, 't1_map');
save_nii(t1_map_image, 't1_map_image');
save_nii(b1_map, 'b1_map');
save_nii(b1_map_image, 'b1_map_image');
writelines(jsonencode(minimal_sidecar(info_mp2rage), 'PrettyPrint', true), fullfile(out_dir, 'mp2rage_image.json'));
fprintf('s01: maps written to %s\n', out_dir);
end

function s = minimal_sidecar(info)
% Keep only the fields read by stages s04/s05 and the dynamic timing used for probe alignment.
d = info.imgdef;
if isfield(d, 'slice_thickness_in_mm'), thick = d.slice_thickness_in_mm; else, thick = d.slice_thickness; end
s.imgdef.echo_time.uniq = d.echo_time.uniq;
s.imgdef.pixel_spacing_x_y.uniq = d.pixel_spacing_x_y.uniq;
s.imgdef.slice_thickness_in_mm.uniq = thick.uniq;
s.imgdef.slice_orientation_TRA_SAG_COR.uniq = d.slice_orientation_TRA_SAG_COR.uniq;
s.pardef.Water_Fat_shift_pixels = info.pardef.Water_Fat_shift_pixels;
s.dynamic_start_times_s = d.dyn_scan_begin_time.uniq;
end
