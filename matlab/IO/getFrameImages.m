%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores current and following images of a given frame for all cameras.
% INPUT PARAMETERS:
% - Path of the folder with the stored images
% - Number of cameras
% - Number of frame we are currently processing

function images = getFrameImages(imgPath, numCameras, currentFrame, defaultImgSize)

    input = InputParameters;

    if(input.onSequence)
        if(currentFrame < 10)
            sFrame0 = [num2str(0) num2str(0) num2str(0) num2str(currentFrame)];
            if(currentFrame+1 < 10)
                sFrame1 = [num2str(0) num2str(0) num2str(0) num2str(currentFrame+1)];
            else
                sFrame1 = [num2str(0) num2str(0) num2str(currentFrame+1)];
            end
        elseif(currentFrame < 100)
            sFrame0 = [num2str(0) num2str(0) num2str(currentFrame)];
            if(currentFrame+1 < 100)
                sFrame1 = [num2str(0) num2str(0) num2str(currentFrame+1)];
            else
                sFrame1 = [num2str(0) num2str(currentFrame+1)];
            end
        else
            sFrame0 = [num2str(0) num2str(currentFrame)];
            if(currentFrame+1 < 1000)
                sFrame1 = [num2str(0) num2str(currentFrame+1)];
            else
                sFrame1 = num2str(currentFrame+1);
            end
        end
    else
        sFrame0 = num2str(currentFrame);
        sFrame1 = num2str(currentFrame+1);
    end
    images = cell(2, numCameras);
    % Get images for frames 't' and 't+1'.
    for CAMERA = 1:numCameras
        % Load image file and store it into the data structure.
        filename1 = [imgPath sFrame0 '_' num2str(CAMERA-1) input.imgExt];
        filename2 = [imgPath sFrame1 '_' num2str(CAMERA-1) input.imgExt];
        % If default image size is negative, no resizing
        if defaultImgSize(1) < 0 || defaultImgSize(2) < 0
            images{1, CAMERA} = imread(filename1);
            images{2, CAMERA} = imread(filename2);
        else % Otherwise, resize input image
            images{1, CAMERA} =  imresize(imread(filename1), defaultImgSize);
            images{2, CAMERA} =  imresize(imread(filename2), defaultImgSize);
        end
        % In case the image is 16 bit encoded, downgrade to 8-bit image.
        if(isa(images{1, CAMERA}, 'uint16'))
            images{1, CAMERA} = uint8(255*double(images{1, CAMERA})/65535);
            images{2, CAMERA} = uint8(255*double(images{2, CAMERA})/65535);
        end
    end
end