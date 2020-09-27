%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualises the segmentation with average of segment colour value.
% INPUT PARAMETERS:
% - 3D colour image
% - Segmentation.

function viewSegmentation(img, segmentation)

    colourSegment = zeros(max(segmentation(:)),3);

    segmentation = cat(3, segmentation, segmentation, segmentation);
    if(length(img(1,1,:)) == 1)
        img = cat(3, img, img, img);
    end
    
    for SEGMENT = 1:length(colourSegment)
        
        colours = img(segmentation == SEGMENT);
        colours = reshape(colours,length(colours)/3,3);
        % fprintf('Lenfth segment %d: %d\n', SEGMENT, length(colours(:,1)));
        colourSegment(SEGMENT,:) = sum(colours, 1) ./ length(colours(:,1));

    end

    imgRows = length(img(:,1,1));
    imgCols = length(img(1,:,1));
    for ROW = 1:imgRows
        for COL = 1:imgCols
            img(ROW,COL,:) = colourSegment(segmentation(ROW,COL),:);
        end
    end
    
    figure; 
    image(uint8(round(img)));
    
end


