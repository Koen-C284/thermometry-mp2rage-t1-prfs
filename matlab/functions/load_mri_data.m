% ========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [img, info] = load_mri_data(filePath)
    % LOAD_MRI_DATA Standardized loader for PAR/REC, XML/REC, and MAT
    % Usage: [img, info] = load_mri_data('path/to/file')

    [folder, name, ext] = fileparts(filePath);
    fullBase = fullfile(folder, name);

    switch lower(ext)
        case {'.par', '.rec'}
            % Couple loadPARREC.m
            fprintf('Loading Philips PAR/REC: %s\n', name);
            [img, info] = loadPARREC([fullBase, '.par']);

        case {'.xml'}
            % Couple loadXMLREC.m
            fprintf('Loading Philips XML/REC: %s\n', name);
            [img, info] = loadXMLREC([fullBase, '.xml']);

        case '.mat'
            fprintf('Loading MAT file: %s\n', name);
            temp = load(filePath);
            vars = fields(temp);
            img = temp.(vars{1}); 
            info = temp;
            
        otherwise
            % Make sure an error is displayed when it breaks
            error('Unsupported file extension: %s', ext);
    end
    
end