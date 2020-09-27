%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes the Z value in 3D space from disparity map.
% INPUT PARAMETERS:
% - disparity
% - focal length 
% - baseline between both cameras
% Assumption #1: both cameras are on the same XY-XZ plane
% Assumption #2: pixel sizes in x-y are equal

function Z = getZfromD(d, focalLength, baseline)

    Z = (focalLength*baseline)./d;

end

