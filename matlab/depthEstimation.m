%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns the depth a 2 given rectified images.
% INPUT PARAMETERS:
% - left and right rectified images
% - depth and 3D flow current approximations
% - baseline of camera set
% - focal length
% - factor of image downsampling

function [newZ segmentPenalties] = depthEstimation(imgL, imgR, Z, u, v, w, baseline, focals, factor, LEVEL)

    input = InputParameters;
    imgRows = length(imgL(:,1,1));
    imgCols = length(imgR(1,:,1));
    
    % Prepare segmented image.
    segL = mexKMeans9D(imgL, Z, u, v, w, input);
%     viewSegmentation(t0ImgPyramid{LEVEL,1}, segL);
    segL = mexConnectedSegments(imgL, segL, Z, u, v, w);    
%     viewSegmentation(t0ImgPyramid{LEVEL,1}, segL);
    segL = mexMergedSegments(imgL, segL, Z, u, v, w, input.mergedThreshold);
    
    empty = zeros(imgRows, imgCols);
    segR = mexKMeans9D(imgR, Z, u, v, w, input);
    segR = mexConnectedSegments(imgR, segR, empty, empty, empty, empty);
    segR = mexMergedSegments(imgR, segR, empty, empty, empty, empty, input.mergedThreshold);
%     viewSegmentation(t0ImgPyramid{LEVEL,1}, segL);

    % Initial disparity estimation.
    numSegL = max(segL(:)); numSegR = max(segR(:));
    baselineFactor = baseline / factor;
    [dispL, scoresL, dispR] = mexComputeLocalDisparities(segL, numSegL, segR, numSegR, imgL, imgR, baselineFactor, input);
%     figure; imagesc(dispL); colormap(gray);

%     keyboard;

    % Cross-Validation
    [disparity stereoConsistency] = mexCrossValidation(dispL, dispR);
%     figure; imagesc(disparity); colormap(gray);

    % Plane fitting
    [disparity planes] = mexComputeSegmentPlanes (segL, numSegL, disparity, scoresL, stereoConsistency, input);
%     figure; imagesc(disparity); colormap(gray);

%     keyboard;

    % Inconsistencies
    disparity = mexConsistentFill(disparity);
%     figure; imagesc(disparity); colormap(gray);

    % Plane Similarity estimation
    segmentPenalties = mexSegmentSimilarities(segL, numSegL, planes, disparity, stereoConsistency, input, LEVEL);
%     minSimilarities = min(segmentPenalties,[],3); figure; imagesc(minSimilarities); colormap(gray);

    % Assign the disparities as new estimated depth Z.
    newZ = getZfromD(disparity, focals, baseline);
    % Apply median filter to get rid of possible remaining outliers.
    newZ = mexMedianFilter(newZ, 3);
%     figure; imagesc(newZ); colorbar;    
    
end