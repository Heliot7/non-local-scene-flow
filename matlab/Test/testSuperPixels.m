% Clean stuff
clc;
close all;
clear all;

% Files path
path = 'data/';

% Input parameters
input = InputParameters;

% Camera parameters
infoCameras = CameraParameters.empty(input.numCameras, 0);
for CAMi = 1:input.numCameras
    infoCameras(CAMi) = CameraParameters.getCameraParameters([path input.dataset '/ProjMat'], CAMi);
end

% Get images
images = getFrameImages([path input.dataset '/' input.dataset '_t'], input.numCameras, 0, input.defaultImgSize);
% Create multi-resolution pyramid for all camera input images.
[t0ImgPyramid, t1ImgPyramid] = buildPyramid(images, input.numLevels, input.factor);

lvl = input.numLevels;
% lvl = input.initLevel;

img = t0ImgPyramid{lvl, 1};
imgRows = length(t0ImgPyramid{lvl, 1}(:,1,1));
imgCols = length(t0ImgPyramid{lvl, 1}(1,:,1));

fprintf('\nResolution image: %dx%d\n', imgRows, imgCols);

% Depth (Z) initialisation.
Z = zeros(imgRows,imgCols);
u = zeros(imgRows,imgCols);
v = zeros(imgRows,imgCols);
w = zeros(imgRows,imgCols);

% (1.1) KMeans segmentation:
t0 = cputime;
segL = mexKMeans9D(img, Z, u, v, w, input);
KMeansTime = cputime - t0;
numSegL = max(segL(:));
fprintf('KMeans %d segments: %f sec\n', numSegL, KMeansTime);
% - Results:
viewSegmentation(img, segL);

% (1.2) Connected components with minimum size:
t0 = cputime;
segL = mexConnectedSegments(img, segL, Z, u, v, w);
ConnectedTime = cputime - t0;
numSegL = max(segL(:));
fprintf('Connected %d segments: %f sec\n', numSegL, ConnectedTime);
% - Results:
viewSegmentation(img, segL);
% viewSegmentContour(img, segL);

% (1.3) Merge similar segments in colour, space and motion (threshold):
t0 = cputime;
segL = mexMergedSegments(img, segL, Z, u, v, w, input.mergedThreshold);
MergedTime = cputime - t0;
numSegL = max(segL(:));
fprintf('Merged %d segments: %f sec\n', numSegL, MergedTime);
% - Results:
viewSegmentation(img, segL);
% viewSegmentContour(img, segL);

fprintf('TOTAL: %f sec\n', ...
    KMeansTime + ConnectedTime + MergedTime);

