%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualises the segmentation by highlighting the countours.
% INPUT PARAMETERS:
% - 3D colour image
% - Segmentation.

function viewSegmentContour(img, segmentation)

    if(length(img(1,1,:)) == 1)
        img = cat(3, img, img, img);
    end
    
    imgRows = length(img(:,1,1));
    imgCols = length(img(1,:,1));
    for ROW = 1:imgRows
        for COL = 1:imgCols
            if(isBorderPxl(segmentation,ROW,COL))
                img(ROW,COL,:) = [0 0 0];                
            end
        end
    end
    
    figure; 
    image(uint8(round(img)));
    
end

function isBorder = isBorderPxl(seg, ROW, COL)

    if(ROW > 1)
        isBorderTop = seg(ROW,COL) ~= seg(ROW-1,COL);
    else
        isBorderTop = false;
    end
    
    if(ROW < length(seg(:,1)))
        isBorderBottom = seg(ROW,COL) ~= seg(ROW+1,COL);
    else
        isBorderBottom = false;
    end

    if(COL > 1)
        isBorderLeft = seg(ROW,COL) ~= seg(ROW,COL-1);
    else
        isBorderLeft = false;
    end
    
    if(COL < length(seg(1,:)))
        isBorderRight = seg(ROW,COL) ~= seg(ROW,COL+1);
    else
        isBorderRight = false;
    end 
    
 isBorder = isBorderTop || isBorderBottom || isBorderLeft || isBorderRight;

end
