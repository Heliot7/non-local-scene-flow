%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Builds a pyramid for both input images with a given number of levels.
% INPUT PARAMETERS:
% - Input images at frames 't' and 't+1' for all cameras
% - Number of levels of our pyramid
% - Scale factor of downsampling

function [t0ImgPyramid, t1ImgPyramid] = buildPyramid(images, numLevels, factor)

    % Define pyramid at frames 't' and 't+1'
    t0ImgPyramid = cell(numLevels, length(images(1,:)));
    t1ImgPyramid = cell(numLevels, length(images(1,:)));
    
    % Downsampling from finest level upwards to fill the whole pyramid
    scale = 1.0;
    for LEVEL = numLevels:-1:1
        for CAMERA = 1:length(images(1,:))
        t0ImgPyramid{LEVEL, CAMERA} = imresize(images{1, CAMERA}, scale);
        % imshow(t0ImgPyramid{LEVEL, CAMERA});
        t1ImgPyramid{LEVEL, CAMERA} = imresize(images{2, CAMERA}, scale);   
        % imshow(t1ImgPyramid{LEVEL, CAMERA});
        end
        scale = scale * factor;
    end
end