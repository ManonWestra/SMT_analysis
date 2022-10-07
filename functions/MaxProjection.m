function [MAX_TIRF_image] = MaxProjection(filename, imsize, frames)
% Create max projection in Matlab instead of ImageJ
% for widefield/tirf images from ONI
% Manon Westra, January 2022

% read image sequence
TIRF_image = zeros([imsize frames],'uint16');
for frame = 1:frames
    [TIRF_image(:,:,frame)] = imread(filename,frame);
end

% max intensity
MAX_TIRF_image = max(TIRF_image, [], 3);