% =========================================================================
% HELPER FUNCTIONS: Actual Flip Angle (AFI) B1 mapping
%
% Extracts and masks the Philips AFI B1 map (central slice when 3D/multi-slice).
% =========================================================================

function [b1_map, b1_map_image] = compute_B1(img, info) %#ok<INUSD>

slice_idx = min(2, size(img, 3));
b1_map = img(:, :, slice_idx, 1, 1, 1, 1, 3);
image = img(:, :, slice_idx, 1, 1, 1, 1, 1);
% create mask where signal intensity > 10% of max (avoid division by zero)
signal = abs(image);
threshold = 0.10 * max(signal(:));
mask = signal > threshold;

b1_map = squeeze(mask .* b1_map);
b1_map_image = squeeze(image);
fprintf("Done reshaping AFI B1 maps ... \n")

end