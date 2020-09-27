%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
% Source code provided by Esther Horbert         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creates a frienly virtual environment to test Scene Flow algorithms.

function [images, infoCameras] = artificialScene()

    input = InputParameters;

    % Step 1: Build an artificial 3D scene:
    % =====================================
    
    sizeSceneCoord = 256;
    % Note: -1 means transparent.
    scene1 = -ones([sizeSceneCoord, sizeSceneCoord, sizeSceneCoord]); 
    scene2 = -ones([sizeSceneCoord, sizeSceneCoord, sizeSceneCoord]);

    % Movement of scene between the two frames.
    % Note: [0; 0; 0] = No movement.
    displacement = [1; 0; 0]; 
    % displacement = [1; 0; 0]; % Move down.
    % displacement = [0; 1; 0]; % Move to the right.
    % displacement = [0; 0; -5]; % Move towards the camera.

    
    % Plane drawn at the back of the scene.
    % Random matrix.
    % Note: same configuration at every run, due to same chosen seed.
    RandStream.setGlobalStream( RandStream('swb2712', 'Seed', 3) );
    scene1(:, :, sizeSceneCoord) = round( rand(sizeSceneCoord)*255 );
%     scene1(:,:,sizeSceneCoord) = 0;
%     scene1(15:end-15,15:end-15,sizeSceneCoord) = 255;
    
    % Check display of random coloured plane.
%     if (input.onDisplay)
%         figure;
%         imagesc(scene1(:, :, sizeSceneCoord));
%     end
    % Check display contour of 3d scene.
%     if (input.onDisplay)
%         figure;
%         displaycontour3d(scene1);
%         title('scene1 view');
%     end

    % Assigns displaced values of scene1 into scene2.
    % It takes into account 3D scene boundaries.
    scene2(max(1,1+displacement(1)):min(end,end+displacement(1)),max(1,1+displacement(2)):min(end,end+displacement(2)),max(1,1+displacement(3)):min(end,end+displacement(3))) = ...
    scene1(max(1,1-displacement(1)):min(end,end-displacement(1)),max(1,1-displacement(2)):min(end,end-displacement(2)),max(1,1-displacement(3)):min(end,end-displacement(3)));

    % Step 2: Take artificial images:
    % ===============================
    
    % Size of the artificial image.
    imgSize = [100, 100]; 

    % Camera parameters.
    %camPos1 = [7; 7; 0];
    %camPos2 = [8; 7; 0];
    camPos1 = [128; 128; 0];
    camPos2 = [129; 128; 0];
    % camPos1 = [0.5*size(scene1,2); 0.5*size(scene1,1); -0.75*size(scene1,3)];
    % camPos2 = [0.75*size(scene1,2); 0.5*size(scene1,1); -0.75*size(scene1,3)];
    sizeSceneCoord1 = size(scene1, 3);
    sizeSceneCoord2 = sizeSceneCoord1;
    cam1 = getCamParams(camPos1, sizeSceneCoord1, imgSize);
    cam2 = getCamParams(camPos2, sizeSceneCoord2, imgSize);
    % Gathers information.
    infoCameras = CameraParameters.empty(input.numCameras, 0);
    camParams1 = CameraParameters;
    camParams1.C = cam1.K;
    camParams1.R = cam1.R;
    camParams1.T = cam1.t;
    camParams1.focalLength = [sizeSceneCoord1 sizeSceneCoord1];
    camParams1.principalPoint = imgSize / 2;
    camParams1.COP = cam1.pos;
    infoCameras(1) = camParams1;
    camParams2 = CameraParameters;
    camParams2.C = cam2.K;
    camParams2.R = cam2.R;
    camParams2.T = cam2.t;
    camParams2.focalLength = [sizeSceneCoord2 sizeSceneCoord2];
    camParams2.principalPoint = imgSize / 2;
    camParams2.COP = cam2.pos;
    infoCameras(2) = camParams2;
    

    numFrames = 2;
    images = cell(numFrames, input.numCameras);
    % Proceed to obtain images from 3D world.
    [pim11 dim] = SceneToImg(scene1, cam1, imgSize);
    pim11 = uint8(pim11);
    images{1,1} = pim11;
    [pim12 dim] = SceneToImg(scene1, cam2, imgSize);
    pim12 = uint8(pim12);
    images{1,2} = pim12;
    [pim21 dim] = SceneToImg(scene2, cam1, imgSize);
    pim21 = uint8(pim21);
    images{2,1} = pim21;
    [pim22 dim] = SceneToImg(scene2, cam2, imgSize);
    pim22 = uint8(pim22);
    images{2,2} = pim22;

    % Check produced images.
    if (0)
        figure;
        subplot(2,2,1);
        imagesc(pim11); title('frame 1, camera 1');
        
        subplot(2,2,2);
        imagesc(pim12); title('frame 1, camera 2');
        
        subplot(2,2,3);
        imagesc(pim21); title('frame 2, camera 1');
        
        subplot(2,2,4);
        imagesc(pim22); title('frame 2, camera 2');
    end
    
end

% Gets a 
function [pim dim] = SceneToImg(scene, cam, imgSize)

    lsf2world = eye(4);
    [pim dim] = lsf2depth(scene, cam, imgSize, lsf2world);
    
end

% Note: This only works if the camera isn't rotated!
function [valueImg depthImg] = lsf2depth(scene, cam, imgSize, lsf2world)

    % Init and extreme values.
    depthImg = ones(imgSize(1),imgSize(2)) * -10;
    valueImg = zeros(imgSize(1),imgSize(2));

    % Perform raytracing for each pixel towards scene from camera pos.
    rayDir = [-floor(imgSize(2)/2); -floor(imgSize(1)/2); cam.lengthScene; 0];  

    % For each camera pixel calculate the depth and the intensity values.
    for ind1 = 1:imgSize(1) % Y
        for ind2 = 1:imgSize(2) % X
            
            % Camera coordinates
            pixp = [ind2; ind1; 0; 1] + rayDir; 
            % Camera center is at origin in camera coordinates.
            x = pixp - [0; 0; 0; 1]; 
            % Every "step" has length one in z-direction
            x = x / x(3); 
            % From camera to world coordinates
            x = cam.Rt4 \ x; 

            % This is an upper bound, not computing the distance between
            % camera and position of lsf here.
            for z = 1:5*size(scene,3) 
                
                xzw = cam.pos + z * x;
                xz = lsf2world \ xzw;
                xz = xz/xz(4);

                % This is discrete case, should maybe use interpolation.
                i1 = round(xz(1)); 
                i2 = round(xz(2));
                i3 = round(xz(3));

                if ( i1 > 0 && i1 <= size(scene,2) && i2 > 0 && i2 <= size(scene,1) ...
                        && i3 > 0 && i3 <= size(scene,3) )
                    % Not transparent case (should enter here whole time).
                    if (scene(i2,i1,i3) >= 0)
                        % Depth of intsersected voxel in scene.
                        depthImg(ind1,ind2) = norm(xzw-cam.pos);
                        % Intersected voxel value.
                        valueImg(ind1,ind2) = scene(i2,i1,i3);
                        break;
                    end
                end
            end
        end
    end

    % Treatment of out of bounds values.
    valueImg(valueImg < 0) = 0;
    valueImg(valueImg > 255) = 255;

end

% Returns the rotation matrix given the angles of rotation for 3D axis.
function R = getRotation(alphaX, alphaY, alphaZ)

    Rx = [1 0 0; 0 cos(alphaX) -sin(alphaX); 0 sin(alphaX) cos(alphaX)];
    Ry = [cos(alphaY) 0 sin(alphaY); 0 1 0; -sin(alphaY) 0 cos(alphaY)];
    Rz = [cos(alphaZ) -sin(alphaZ) 0; sin(alphaZ) cos(alphaZ) 0;0 0 1];
    R = Rx*Ry*Rz;

end

function cam = getCamParams(camPos, sizeSceneCoord, imgSize, alphaX, alphaY, alphaZ)

    if (nargin <= 3)
        % No rotation.
        alphaX = 0;
        alphaY = 0;
        alphaZ = 0;
    else
        % Angles in degrees.
        alphaX = alphaX * pi / 180;
        alphaY = alphaY * pi / 180;
        alphaZ = alphaZ * pi / 180;
    end
    
    % 3x3 Intrinsic matrix.
    cam.K = [sizeSceneCoord 0 imgSize(1)/2; 0 sizeSceneCoord imgSize(2)/2; 0 0 1];
    % 3x3 Rotation matrix.
    cam.R = getRotation(alphaX,alphaY,alphaZ);
    % 3x1 Translation matrix.
    cam.t = -cam.R*camPos;
    % 4x4 Extrinsic matrix.
    cam.Rt4 = [[cam.R cam.t]; 0 0 0 1];
    % COP = -cam.R \ cam.t;
    cam.pos = [camPos; 1]; 
    % lengthScene = cam.K(1,1) = focalLength_xy;
    cam.lengthScene = sizeSceneCoord;
end

% Shows the countour of a given discretized 3D scene with a 3D matrix.
% Note 1: >= 0 means solid.
% Note 2: < 0 means transparent.
% Note 3: This function will not consider any voxel colours, just surface.
function displaycontour3d(img, color)

    if (nargin < 2)
        color = 'red';
    end
    
    imgSize = size(img);
    p = patch(isosurface(img,-1));
    set(p,'FaceColor',color,'EdgeColor','none');
    daspect([1 1 1]);
    % Make surface smoother
    isonormals(img,p);
    lighting gouraud;
    % ... otherwise it doesn't look 3D.
    camlight; 
    view([-10,-50,-10]);
    % Another light source.
    camlight('headlight'); 
    hold on;
    plot3([0 0 imgSize(2) imgSize(2) 0 ],[0 imgSize(1) imgSize(1) 0 0],[0 0 0 0 0],'b');
    plot3([0 0 imgSize(2) imgSize(2) 0 ],[0 imgSize(1) imgSize(1) 0 0],[imgSize(3) imgSize(3) imgSize(3) imgSize(3) imgSize(3) ],'b');
    plot3([0 0],[0 0],[0 imgSize(3)],'b'); 
    plot3([imgSize(2) imgSize(2)],[0 0],[0 imgSize(3)],'b');
    plot3([imgSize(2) imgSize(2)],[imgSize(1) imgSize(1)],[0 imgSize(3)],'b');
    plot3([0 0],[imgSize(1) imgSize(1)],[0 imgSize(3)],'b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
end
