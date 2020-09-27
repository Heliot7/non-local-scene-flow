%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Class: Gathers all camera parameters.

classdef CameraParameters
    
    properties
        C;
        R;
        T;
        focalLength;
        principalPoint;
        COP;
    end
    
    methods(Static)
        % Get camera info from files for a specific camera.
        function params = getCameraParameters(filePath, CAMi)
            params = CameraParameters;
            
            % Intrinsic parameters.
            fid = fopen([filePath '/CMat' num2str(CAMi-1) '.txt']);
            params.C = cell2mat(reshape(textscan(fid,'%f %f %f %f %f %f %f %f %f', 'delimiter', '\n'), 3, 3)');
            params.focalLength = [params.C(1,1), params.C(2,2)];
            params.principalPoint = [params.C(1,3), params.C(2,3)];
            fclose(fid);
            % Extrinsic parameters (rotations).
            fid = fopen([filePath '/RMat' num2str(CAMi-1) '.txt']);
            params.R = cell2mat(reshape(textscan(fid,'%f %f %f %f %f %f %f %f %f', 'delimiter', '\n'), 3, 3)');
            fclose(fid);
            % Extrinsic parameters (translation).
            fid = fopen([filePath '/TMat' num2str(CAMi-1) '.txt']);
            params.T = cell2mat(textscan(fid,'%f %f %f', 'delimiter', '\n')');
            fclose(fid);
            
            % Calculates the centre of projection for the current camera.
            params = params.calculateCOP;
        end
        % Returns a built projection matrix with given parameters.
    end
    methods
        % Current multiresolution pyramid scale factor is added.
        function projMatrix = buildProjMatrix(params, scale)
            scaledC = params.C;
            scaledC(3,3) = 1/scale;
            projMatrix = [scaledC, zeros(3,1)]*[params.R params.T; 0 0 0 1];
        end
        % Calculates Centre of Projection for a given view.
        function cam = calculateCOP(cam)
            % Assign COP for each view, moving from camera to world coordinates.
            % We assume that: Camera point = (0,0,0) in camera coordinates.
            % PointWorld = (PointCamera - Translation)* Rotation^-1
            cam.COP = cam.R \ -cam.T;
        end
        
    end
    
end