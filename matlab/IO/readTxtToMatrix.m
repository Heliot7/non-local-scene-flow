%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reads from a text file to a matrix.
% INPUT PARAMETERS:
% - Path
% - Image size

function M = readTxtToMatrix(path, imgSize)

    % Define matrix
    M = zeros(imgSize(1), imgSize(2));

    numRow = 1;
    fid = fopen(path);
    % Iterative until all rows have been read and added to the matrix.
    while(~feof(fid))
        row = textscan(fid, '%n', imgSize(2));
        % Last read line will be empty, ignore it!
        if(numRow <= imgSize(1))
            M(numRow,:) = row{1}';
            numRow = numRow + 1;
        end
    end
    fclose(fid);

end

