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

imgRows = length(t0ImgPyramid{input.initLevel,1}(:,1,1));
imgCols = length(t0ImgPyramid{input.initLevel,1}(1,:,1));

fprintf('\nResolution image: %dx%d\n', imgRows, imgCols);

% Depth (Z) initialisation.
Z = zeros(imgRows,imgCols);
% if(input.onZinit)
%     Z = dpth(:,:,1);
% end
% 3D-flow (u,v,w) initialisation.
u = zeros(imgRows,imgCols);
v = zeros(imgRows,imgCols);
w = zeros(imgRows,imgCols);

% (1.1) KMeans segmentation:
t0 = cputime;
segL = mexKMeans9D(t0ImgPyramid{input.initLevel, 1}, Z, u, v, w, input);
KMeansTime = cputime - t0;
numSegL = max(segL(:));
fprintf('KMeans %d segments: %f sec\n', numSegL, KMeansTime);
% - Results:
% viewSegmentation(t0ImgPyramid{input.initLevel,1}, segL);

% (1.2) Connected components with minimum size:
t0 = cputime;
segL = mexConnectedSegments(t0ImgPyramid{input.initLevel, 1}, segL, Z, u, v, w);
ConnectedTime = cputime - t0;
numSegL = max(segL(:));
fprintf('Connected %d segments: %f sec\n', numSegL, ConnectedTime);
% - Results:
% viewSegmentation(t0ImgPyramid{input.initLevel, 1}, segL);
% viewSegmentContour(t0ImgPyramid{input.initLevel,1}, segL);

% (1.3) Merge similar segments in colour, space and motion (threshold):
t0 = cputime;
segL = mexMergedSegments(t0ImgPyramid{input.initLevel, 1}, segL, Z, u, v, w, input.mergedThreshold);
MergedTime = cputime - t0;
numSegL = max(segL(:));
fprintf('Merged %d segments: %f sec\n', numSegL, MergedTime);
% - Results:
viewSegmentation(t0ImgPyramid{input.initLevel,1}, segL);
% viewSegmentContour(t0ImgPyramid{input.initLevel,1}, segL);

% (1.4) Proceed again from 1.1 to 1.3 with right image:
empty = zeros(imgRows, imgCols);
t0 = cputime;
segR = mexKMeans9D(t0ImgPyramid{input.initLevel, 2}, Z, u, v, w, input);
segR = mexConnectedSegments(t0ImgPyramid{input.initLevel, 2}, segR, empty, empty, empty, empty);
segR = mexMergedSegments(t0ImgPyramid{input.initLevel, 2}, segR, empty, empty, empty, empty, input.mergedThreshold);
numSegR = max(segR(:));
OtherSegmentationTime = cputime - t0;
fprintf('Again, all segmentation steps for right image: %f sec\n', OtherSegmentationTime);
% - Results:
viewSegmentation(t0ImgPyramid{input.initLevel, 2}, segR);

% (2.1) Initial disparity estimation:
[projMatrices, focals, principals] = buildProjMatrices(infoCameras, input.numLevels, input.factor);
baseline = uint8(abs(infoCameras(1).COP(1) - infoCameras(2).COP(1)) / projMatrices{input.initLevel,1}(3,3)+1);
t0 = cputime;
[dispL, scoresL, dispR, scoresR] = mexComputeLocalDisparities(segL, numSegL, segR, numSegR, ...
    t0ImgPyramid{input.initLevel, 1:2}, baseline, input);
DisparityTime = cputime - t0;
fprintf('Initial disparity: %f sec\n', DisparityTime);
% - Results:
figure; imagesc(dispL); colormap(gray); figure; imagesc(dispR); colormap(gray);

% (2.2) Cross Validation:
t0 = cputime;
[disparity stereoConsistency] = mexCrossValidation(dispL, dispR);
CrossValidationTime = cputime - t0;
fprintf('Cross-validation: %f sec\n', CrossValidationTime);
% - Results:
% figure; imagesc(stereoConsistency); colormap(gray);
figure; imagesc(disparity); colormap(gray);

% (2.3) Plane fitting
t0 = cputime;
[disparity planes] = mexComputeSegmentPlanes (segL, numSegL, disparity, scoresL, stereoConsistency, input);
PlaneFittingTime = cputime - t0;
fprintf('Plane Fitting: %f sec\n', PlaneFittingTime);
% - Results:
figure; imagesc(disparity); colormap(gray);

% (2.4) Fill inconsistent pixels coherent depths
t0 = cputime;
disparity = mexConsistentFill(disparity);
ConsistentFill = cputime - t0;
fprintf('Consistent Filling: %f sec\n', ConsistentFill);
% - Results:
figure; imagesc(disparity); colormap(gray);

% (2.5) Compute similarities
t0 = cputime;
similarities = mexSegmentSimilarities(segL, numSegL, planes, disparity, stereoConsistency, input, input.initLevel);
SimilaritiesTime = cputime - t0;
fprintf('NonLocal Similarities: %f sec\n', SimilaritiesTime);
% - Results:
% minSimilarities = min(similarities,[],3); figure; imagesc(minSimilarities); colormap(gray);

fprintf('TOTAL: %f sec\n', ...
    KMeansTime + ConnectedTime + MergedTime + OtherSegmentationTime + ...
    DisparityTime + CrossValidationTime + ConsistentFill + PlaneFittingTime + SimilaritiesTime);

