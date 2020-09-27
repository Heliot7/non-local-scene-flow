%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of a 3D Scene Flow technique given a sequence of images.
% ====================================================================

% Cleans MATLAB output screen.
clc;
% Closes all plots, figures or similars that are currently open.
close all;
% Removes current data structures
clear input;

% INPUT PARAMETERS:
input = InputParameters;

if(~input.onSequence)
    path = 'data/';
else
    path = '/work/panareda/data/';
end

if(input.isRealScenario)
    % Gathers all necessary data for 3x4 projection matrices (for each camera).
    infoCameras = CameraParameters.empty(input.numCameras, 0);
    for CAMi = 1:input.numCameras
        infoCameras(CAMi) = CameraParameters.getCameraParameters([path input.dataset '/ProjMat'], CAMi);
    end
else
    % Builds a 3D artificial scene for testing purposes.
    [images, infoCameras] = artificialScene();
end

% Reading images on the fly from current video sequence.
currentFrame = 0;
goSequence = true;
while(goSequence)
    
    if(input.isRealScenario)
        % Get new image frames to process.
        images = getFrameImages([path input.dataset '/' input.dataset '_t'], input.numCameras, currentFrame, input.defaultImgSize);
    end

    % Pre-timing
    t0 = cputime;
    
    % Call 3D scene flow method.  
    % - Matlab call:
%     profile on;
    [Z, u, v, w] = nonLocalSceneFlow(images, infoCameras);
%     profile off;
%     profile viewer;
%     keyboard;
    % - C++ call:
%     [Z, u, v, w] = mexNonLocalSceneFlow(images, infoCameras);
    
    % Post-timing
    runtime = cputime - t0;
    fprintf('Frame %d - time: %f sec\n', currentFrame, runtime);
    
    % Show results in a 3D space.
    showResults(images, infoCameras, Z, u, v, w, runtime, currentFrame);

    % Forcing the 'while' loop to finish after just 1 iteration.
    currentFrame = currentFrame + 1;
    goSequence = (input.onSequence && ...
        exist([path input.dataset '/' input.dataset '_t' num2str(currentFrame+1) '_0' input.imgExt], 'file'));
end
