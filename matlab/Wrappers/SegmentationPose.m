function SegmentationPose(image, numSegments)

    Z = zeros(imgRows,imgCols);
    u = zeros(imgRows,imgCols);
    v = zeros(imgRows,imgCols);
    w = zeros(imgRows,imgCols);
    input.numK = numSegments;
    input.onConnectedComponents = true;
    input.numIterKMeans = 5;

    segmentation = mexKMeans9D(image, Z, u, v, w, input);
    viewSegmentation(image, segmentation);
    segmentation = mexConnectedSegments(image, segmentation);
    viewSegmentation(image, segmentation);
    viewSegmentContour(image, segmentation);
end
