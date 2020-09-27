%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Builds all projection matrices per camera and multires. pyramid level.
% INPUT PARAMETERS:
% - Information about all camera parameters
% - Number of pyramid levels
% - Factor of level reduction

function [projMatrices, focals, principals] = buildProjMatrices(infoCameras, numLevels, factor)

    % Initialize proj. matrices, focal lengths and principal points.
    projMatrices = cell(numLevels, length(infoCameras));
    focals = cell(numLevels, length(infoCameras));
    principals = cell(numLevels, length(infoCameras));
    scale = 1.0;
    % Proceeds to fill the pyramid for all views and levels.
    for LVLi = numLevels:-1:1
        for CAMi = 1:length(infoCameras)
            projMatrices{LVLi, CAMi} = infoCameras(CAMi).buildProjMatrix(scale);
            focals{LVLi, CAMi} = infoCameras(CAMi).focalLength * scale;
            principals{LVLi, CAMi} = infoCameras(CAMi).principalPoint * scale;
        end
        scale = scale * factor;
    end

end
