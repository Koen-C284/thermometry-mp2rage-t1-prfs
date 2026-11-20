function s04_prfs(data_dir)
% S04_PRFS  PRFS thermometry (HIMM) for one dataset; wrapper around the bundled HIMM driver.
%
%   s04_prfs(data_dir)
%
% Runs matlab/third_party/HIMM/s04_prfs_thermometry.m on data_dir/mp2rage_image.nii.gz (+ .json) and copies
% its results to data_dir/derived/temp_prfs.nii.gz and data_dir/derived/masks.nii.gz.
% Requires the ISMRM Fat-Water Toolbox on the path (see setup_paths and THIRD_PARTY.md).
%
% The driver expects <dir>/02_interim/03_registered/ as input, writes <dir>/03_processed/prfs_thermometry/,
% and also writes figures and copies to the root of the current drive ('\...'). It is therefore run from a
% temporary folder that is mapped to a spare drive letter (Windows), so nothing is written outside that folder.

assert(exist('fw_i2cm1i_3pluspoint_hernando_graphcut', 'file') > 0, ...
    's04_prfs:toolbox', 'ISMRM Fat-Water Toolbox not found; call setup_paths(<external folder>) first.');
work = tempname;
reg = fullfile(work, '02_interim', '03_registered');
mkdir(reg);
copyfile(fullfile(data_dir, 'mp2rage_image.nii.gz'), reg);
copyfile(fullfile(data_dir, 'mp2rage_image.json'), reg);

state = struct('data_dir', data_dir, 'work', work, 'here', pwd, 'matlab_path', path, 'drive', '');
if ispc
    letters = 'ZYXWVUTSRQPONMLKJIHGFE';
    state.drive = letters(find(arrayfun(@(c) ~isfolder([c ':\']), letters), 1));
    [status, msg] = system(sprintf('subst %s: "%s"', state.drive, work));
    assert(status == 0, 's04_prfs:subst', 'Could not map a drive letter to %s: %s', work, msg);
    run_dir = [state.drive ':\'];
else
    run_dir = work;
end
setappdata(0, 's04_prfs_state', state);
cd(run_dir);
try
    run_driver(run_dir);
catch err
    cleanup(getappdata(0, 's04_prfs_state'));
    rethrow(err);
end

state = getappdata(0, 's04_prfs_state');
out = fullfile(state.data_dir, 'derived');
if ~isfolder(out), mkdir(out); end
res = fullfile(state.work, '03_processed', 'prfs_thermometry');
copyfile(fullfile(res, 'temp_prfs.nii.gz'), out);
copyfile(fullfile(res, 'masks.nii.gz'), out);
cleanup(state);
fprintf('s04: PRFS temperature and masks written to %s\n', out);
end

function run_driver(scan_data_dir) %#ok<INUSD> used by the driver script
s04_prfs_thermometry;   % bundled HIMM driver script (clears its workspace except scan_data_dir)
end

function cleanup(state)
cd(state.here);
path(state.matlab_path);   % the driver adds its working folder to the path
close all;
if ~isempty(state.drive)
    system(sprintf('subst %s: /d', state.drive));
end
if isfolder(state.work)
    rmdir(state.work, 's');
end
rmappdata(0, 's04_prfs_state');
end
