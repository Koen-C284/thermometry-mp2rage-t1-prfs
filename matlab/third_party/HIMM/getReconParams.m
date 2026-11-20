function [algoParams] = getReconParams(R2star_range,fm_range,itr,lambda,th,fatmodel)
% Define the reconstruction parameters and store them in a struct
% The struct is passed onto the ISMRM toolbox for WF reconstruction

% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;

if ~exist('fatmodel','var')
%     fwater = 4.7;
%     fwater = 0;
%     fatppms = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]-4.7;% 9-peak model from welch MEP
%     fatalphas = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; % Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
%     fatalphas = fatalphas./sum(fatalphas);
    algoParams.species(2).name = 'fat';
    algoParams.species(2).frequency = [-3.8 -3.4 -2.6 -1.94 -0.39 0.6];
    algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

else
    algoParams.species(2).name = 'fat';
    algoParams.species(2).frequency = fatmodel.frequency;
    algoParams.species(2).relAmps = fatmodel.relAmps;
end

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [R2star_range(1) R2star_range(2)]; % Range of R2* values
algoParams.NUM_R2STARS = R2star_range(3); % Number of R2* values for quantization
algoParams.range_fm = [fm_range(1) fm_range(2)]; % Range of field map values
algoParams.NUM_FMS = fm_range(3); % Number of field map values to discretize
algoParams.NUM_ITERS = itr; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % KC changed to 1 % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = lambda; % Regularization parameter
algoParams.LMAP_EXTRA = 0; % Changed to zero KC (lambda); % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
algoParams.THRESHOLD = th;

end

