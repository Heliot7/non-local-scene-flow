%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores from a matrix to a txt file in a given path.
% INPUT PARAMETERS:
% - Path
% - Matrix

function storeMatrixToTxt(path, M)

    dlmwrite(path, M, 'delimiter', '\t');

end

