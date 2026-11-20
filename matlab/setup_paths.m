function setup_paths(external_root)
% SETUP_PATHS  Add this repository's MATLAB code, the bundled third-party code and external dependencies to the path.
%
%   setup_paths()               uses the folder in the environment variable T1THERMO_EXTERNAL
%   setup_paths(external_root)  uses the given folder
%
% Bundled (matlab/third_party): EPG-X (MIT) and HIMM (CC BY-NC 4.0, see third_party/HIMM/NOTICE.md).
% external_root must contain the code that is not redistributed here (see THIRD_PARTY.md):
%   <external_root>/fatwater_toolbox   ISMRM Fat-Water Toolbox (hernando/*), needed for s04 and s05 step 1
%   <external_root>/parrec             Philips PAR/REC reader (loadPARREC), only needed for stage s01

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'functions'));
addpath(genpath(fullfile(here, 'third_party', 'EPG-X')));
addpath(genpath(fullfile(here, 'third_party', 'HIMM')));

if nargin < 1 || isempty(external_root)
    external_root = getenv('T1THERMO_EXTERNAL');
end
if isempty(external_root)
    return  % stages that only use shipped intermediates still run
end
for sub = {'fatwater_toolbox', 'parrec'}
    p = fullfile(external_root, sub{1});
    if isfolder(p)
        addpath(genpath(p));
    end
end
end
