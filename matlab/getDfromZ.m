%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes the disparity from a given 3D Z value matrix.
% INPUT PARAMETERS:
% - Z values
% - focal length 
% - baseline between both cameras
% Assumption #1: both cameras are on the same XY-XZ plane
% Assumption #2: pixel sizes in x-y are equal

function d = getDfromZ(Z, focalLength, baseline)

    d = (focalLength*baseline)./Z;

end

