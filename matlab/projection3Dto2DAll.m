%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates the projection of all 3D points from a given image.
% INPUT PARAMETERS:
% - Projection matrix of the camera that took the image
% - All image 3D points P

function p = projection3Dto2DAll(projMatrix, P)

    P = [P; ones(1,length(P(1,:)))];

    den = (projMatrix(3,:)*P);
    p = (projMatrix(1:2,:)*P);
    p(1,:) = p(1,:) ./ den;
    p(2,:) = p(2,:) ./ den;
    
    % Round at 4th decimal for precision after projection
    p = round(p*10000)/10000;
    
end
