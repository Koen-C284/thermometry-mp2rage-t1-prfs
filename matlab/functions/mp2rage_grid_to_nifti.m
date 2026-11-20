function vol = mp2rage_grid_to_nifti(outputNii, T1n, b1_scalen, params)
% MP2RAGE_GRID_TO_NIFTI
%   Simulate MP2RAGE signals over a grid of T1 and B1+ scaling values and
%   save only a NIfTI file.
%
%   mp2rage_grid_to_nifti(outputNii, T1n, b1_scalen, params)
%
% Inputs
%   outputNii  : string, output NIfTI filename (e.g., 'mp2rage_dictionary.nii.gz')
%   T1n        : vector of T1 values [ms]
%   b1_scalen  : vector of B1+ scale factors (e.g., 0.6:0.05:1.4)
%   params     : struct with scalar fields (all scalars):
%       .TR              [ms]    GRE TR inside each readout train
%       .T2              [ms]    T2
%       .PE_steps        [#]     # of phase-encode steps per readout
%       .number_of_reps  [#]     # of MP2RAGE repetitions (>=1)
%       .PD              [arb]   proton density scaling (1 recommended)
%       .TI1             [ms]    inversion time to k-space center of S1
%       .TI2             [ms]    inversion time to k-space center of S2
%       .MP2RAGETR       [ms]    repetition time of the full MP2RAGE module
%       .alpha1          [deg]   flip angle for readout 1 (nominal)
%       .alpha2          [deg]   flip angle for readout 2 (nominal)
%       .inv_eff         [0-1]   inversion efficiency (default 0.96)
%       .phase_increment [deg]   RF spoiling phase increment (default 120)
%
% Output NIfTI (float32):
%   size = [numel(T1n)  x  numel(b1_scalen)  x  2]
%   volume(:,:,1) = S1 at k-space center (final repetition)
%   volume(:,:,2) = S2 at k-space center (final repetition)
%
% Requirements:
%   - EPG_GRE2, EPG_GRE3 functions available on the MATLAB path.
%   - RF_phase_cycle available; if not, a fallback is implemented below.
%
% Notes:
%   - We sample the complex GRE signal at the k-space center line.
%   - If timing constraints are violated, we store NaN at that grid point.
%   - For multiple reps, we propagate the state (finit) through TC to reach
%     steady-state before the final S1/S2 readouts.

tic;

% -------- Parameters & defaults --------
mustBeVector(T1n); mustBeVector(b1_scalen);
assert(isstruct(params), 'params must be a struct of scalar fields.');

TR            = fetchfield(params, 'TR, required');
T2            = fetchfield(params, 'T2, required');
PE_steps      = fetchfield(params, 'PE_steps, required');
number_of_reps= fetchfield(params, 'number_of_reps, required');
PD            = fetchfield(params, 'PD, required');
TI1           = fetchfield(params, 'TI1, required');
TI2           = fetchfield(params, 'TI2, required');
MP2RAGETR     = fetchfield(params, 'MP2RAGETR, required');
alpha1         = fetchfield(params, 'alpha1, required');
alpha2         = fetchfield(params, 'alpha2, required');
inv_eff       = fetchfield(params, 'inv_eff', 0.96);
phase_increment = fetchfield(params, 'phase_increment', 120);

% Safety checks
assert(isscalar(TR) && isscalar(T2) && isscalar(PE_steps) && isscalar(number_of_reps) ...
    && isscalar(PD) && isscalar(TI1) && isscalar(TI2) && isscalar(MP2RAGETR) ...
    && isscalar(alpha1) && isscalar(alpha2), ...
    'All parameters except T1n and b1_scalen must be scalars.');
assert(PE_steps >= 2 && mod(PE_steps,2)==0, 'PE_steps must be an even integer.');
assert(number_of_reps >= 1, 'number_of_reps must be >= 1');

% RF phase cycle (Philips typical spoiling uses ~120 deg increment)
if exist('RF_phase_cycle','file')
    mp2rage_phase = RF_phase_cycle(PE_steps, phase_increment).';
else
    % fallback: 120 deg incremental spoiling (most common)
    mp2rage_phase = deg2rad(mod(cumsum(0:PE_steps-1)*phase_increment, 360)).';
end

kspace_centre = PE_steps/2 + 1;

% Delays per MP2RAGE (assuming k-space center is in the middle of the train)
prep1.flip   = pi * inv_eff;  % inversion
prep1.t_delay= TI1 - TR*(PE_steps/2);  % "TA"
prep2.flip   = 0;
prep2.t_delay= TI2 - (TI1 + PE_steps*TR); % "TB"
TC           = MP2RAGETR - (TI2 + TR*PE_steps/2); % "TC"

% Output volume: [T1 x B1 x 2]
nT1  = numel(T1n);
nB1  = numel(b1_scalen);

% Preallocate 1D buffers (for parfor)
S1vec = nan(nT1*nB1, 1, 'single');
S2vec = nan(nT1*nB1, 1, 'single');

% Precompute scaled flip arrays indices
alpha1_nom = alpha1*pi/180;
alpha2_nom = alpha2*pi/180;

% Parallelize over grid
parfor idx = 1:(nT1*nB1)

    [iT, jB] = ind2sub([nT1, nB1], idx);
    T1       = T1n(iT);
    b1_scale = b1_scalen(jB);

    % Feasibility: ensure timing constraints
    feasible = ((TI2-TI1)>PE_steps*TR) && ...
               (MP2RAGETR>(TI2+TR*(PE_steps-(kspace_centre-1))));

    if ~feasible
        S1vec(idx) = single(NaN);
        S2vec(idx) = single(NaN);
        continue;
    end

    alfa1_scale = repmat(alpha1_nom * b1_scale, [PE_steps, 1]); % radians
    alfa2_scale = repmat(alpha2_nom * b1_scale, [PE_steps, 1]); % radians

    % -------- First repetition --------
    [F0_1, ~, ~, F_1] = EPG_GRE2(alfa1_scale, mp2rage_phase, TR, T1, T2, ...
        'kmax', inf, 'prep', prep1, 'zinit', PD);
    finit = F_1(:, end);

    % k-space center signal of first readout
    S1 = F0_1(kspace_centre);

    [F0_2, ~, ~, F_2] = EPG_GRE3(alfa2_scale, mp2rage_phase, TR, T1, T2, ...
        'kmax', inf, 'prep', prep2, 'finit', finit, 'zinit', PD);
    finit = F_2(:, end);

    % Relax through TC to next rep
    if number_of_reps > 1
        ridx = length(finit);
        E1_3 = exp(-TC/T1);
        E2_3 = exp(-TC/T2);
        b3   = zeros(ridx,1); b3(3) = PD*(1 - E1_3);
        E_3  = spdiags(repmat([E2_3 E2_3 E1_3],[1 ridx])', 0, ridx, ridx);
        finit = E_3 * finit + b3;

        % Additional reps
        for rep = 2:number_of_reps
            [F0_1, ~, ~, F_1] = EPG_GRE3(alfa1_scale, mp2rage_phase, TR, T1, T2, ...
                'kmax', inf, 'prep', prep1, 'finit', finit, 'zinit', PD);
            finit = F_1(:, end);
            S1    = F0_1(kspace_centre); % overwrite with latest (steady)

            [F0_2, ~, ~, F_2] = EPG_GRE3(alfa2_scale, mp2rage_phase, TR, T1, T2, ...
                'kmax', inf, 'prep', prep2, 'finit', finit, 'zinit', PD);
            finit = F_2(:, end);

            % relax TC again except after last rep (harmless if applied)
            ridx = length(finit);
            E1_3 = exp(-TC/T1);
            E2_3 = exp(-TC/T2);
            b3   = zeros(ridx,1); b3(3) = PD*(1 - E1_3);
            E_3  = spdiags(repmat([E2_3 E2_3 E1_3],[1 ridx])', 0, ridx, ridx);
            finit = E_3 * finit + b3;
        end
    end

    % k-space center signal of second readout after final rep
    S2 = F0_2(kspace_centre);

    S1vec(idx) = single(abs(S1));
    S2vec(idx) = single(abs(S2));
end

% Alternative axis layout: [nT1 x nB1 x 2]
S1t1b1 = reshape(S1vec, [nT1, nB1]);  % transpose to get [nT1 x nB1]
S2t1b1 = reshape(S2vec, [nT1, nB1]);
vol = cat(3, S1t1b1, S2t1b1);

% --------- Save NIfTI only ----------
% NOTE: these are not spatial axes; we store identity spacing and annotate.
try
    meta = struct();
    meta.Datatype    = 'single';
    meta.Description = sprintf('MP2RAGE S1/S2 @k0 over T1 x B1 grid | TR=%.4f TI1=%.4f TI2=%.4f FA1=%.1f FA2=%.1f PE=%d Reps=%d', ...
        TR, TI1, TI2, alpha1, alpha2, PE_steps, number_of_reps);
    meta.PixelDimensions = [1 1 1];
    niftiwrite(vol, outputNii, meta, 'Compressed', endsWith(lower(outputNii), '.gz'));
catch
    % Minimal write if meta not accepted in your MATLAB version
    niftiwrite(vol, outputNii);
end

fprintf('Saved NIfTI: %s  [size %dx%dx%d]\n', outputNii, size(vol));

toc;
end

% ------------- Helpers ----------------
function val = fetchfield(S, name, default)
parts = strsplit(name, ',');
key   = strtrim(parts{1});
if isfield(S, key)
    val = S.(key);
    return;
end
if nargin == 3
    val = default;
else
    error('Missing required field params.%s', key);
end
end

function mustBeVector(x)
if ~(isvector(x) && ~isempty(x))
    error('Input must be a non-empty vector.');
end
end
