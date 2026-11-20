function run_exvivo(external_root, recompute_external)
% RUN_EXVIVO  MATLAB stages for the ex vivo phantom dataset (data/exvivo).
%
%   run_exvivo()                         s05 from the shipped intermediates (only this repository's code)
%   run_exvivo(external_root, true)      also recompute s04 (PRFS) and s05 step 1; needs the ISMRM
%                                        Fat-Water Toolbox in <external_root>/fatwater_toolbox
%
% Afterwards run the analysis and figures with:  python python/s06_exvivo.py
% (s01 is not needed: its output mp2rage_image.nii.gz is shipped; s02 and s03 are not used ex vivo.)

if nargin < 1, external_root = ''; end
if nargin < 2, recompute_external = false; end
here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(here);
setup_paths(external_root);

data_dir = fullfile(repo, 'data', 'exvivo');
cfg = jsondecode(fileread(fullfile(repo, 'config', 'exvivo.json')));

if recompute_external
    s04_prfs(data_dir);
    delete(fullfile(data_dir, 'derived', 'r2star_map.nii.gz'));
    delete(fullfile(data_dir, 'derived', 'muscle_mask_s05.nii.gz'));
end
s05_t1_thermometry(data_dir, cfg.s05);
fprintf('\nMATLAB stages done. Next: python python/s06_exvivo.py\n');
end
