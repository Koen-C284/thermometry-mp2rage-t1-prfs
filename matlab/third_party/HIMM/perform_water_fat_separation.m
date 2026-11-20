function [Wlib,Flib,dw0lib,R2starlib,fatmodel] = perform_water_fat_separation(imgs,GRE_info,frequency_system,b0,prc,waterfatparams)
% Apply water/fat separation to baseline data 

% Selection of fat model and construction of algorithm parameters
fatmodel = get_fatmodel(waterfatparams.num_fat_peaks,prc,waterfatparams.water_shift_ppm);
algoParams = getReconParams(waterfatparams.R2star_range,waterfatparams.off_resonance_range,waterfatparams.itr,waterfatparams.lambda,waterfatparams.threshold,fatmodel);

% WF separation is done for each dynamic individually with image size of [ nx | ny | 1 | ncoils | nTE ]
for ii = 1:size(imgs,4)
    % Image data parameters  
    imDataParams = getImDataParams(GRE_info,imgs(:,:,:,ii),frequency_system,b0,prc);
    %algoParams.fieldmap = dw0lib(:,:,ii)/(2*pi);
    %algoParams.r2starmap = R2starlib(:,:,ii);

    % Perform Mixed-Magnitude method initialized with the Graph Cut method
    outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams,algoParams);
    %outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit(imDataParams,algoParams);
    % In perform_water_fat_separation.m:
    if isfield(waterfatparams, 'initial_fm') && isfield(waterfatparams, 'initial_r2')
        algoParams.fieldmap = waterfatparams.initial_fm;
        algoParams.r2starmap = waterfatparams.initial_r2;
        algoParams.NUM_ITERS = 15; 
    end
    % Combine data (water, fat, off-resonance and R2*) of individual dynamics
    Wlib(:,:,ii) = outParams.species(1).amps;
    Flib(:,:,ii) = outParams.species(2).amps;
    dw0lib(:,:,ii) = outParams.fieldmap*2*pi;
    R2starlib(:,:,ii) = outParams.r2starmap;
    
end

end
