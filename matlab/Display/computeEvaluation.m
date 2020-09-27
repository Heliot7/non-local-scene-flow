%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluates output variables with ground truth.
% INPUT PARAMETERS:
% - estimated Z, u, v, w

function eval = computeEvaluation(X, Y, Z, u, v, w, uProj, vProj, infoCameras, params)

    input = InputParameters;
    
    if(~input.isRealScenario)
        input.dataset = 'Plane';
    end
    
    % Adapt results (which are possibly at a lower resolution to real one.
    X = imresize(X, params.imgSize);
    Y = imresize(Y, params.imgSize);
    Z = imresize(Z, params.imgSize);
    u = imresize(u, params.imgSize);
    v = imresize(v, params.imgSize);
    w = imresize(w, params.imgSize);
    uProj = imresize(uProj, params.imgSize);
    vProj = imresize(vProj, params.imgSize);
    
    % Get the ground truth from the dataset folder.
    if(~input.onSequence)
        path = ['data/' input.dataset '/GroundTruth'];
    else
        path = ['/work/panareda/data/' input.dataset '/GroundTruth'];
    end
    if(params.frame < 10)
            sFrame = [num2str(0) num2str(0) num2str(0) num2str(params.frame)];
    elseif(params.frame < 100)
            sFrame = [num2str(0) num2str(0) num2str(params.frame)];
    else
            sFrame = [num2str(0) num2str(params.frame)];
    end
    if(exist([path '/u'], 'dir'))
        Z_ground = readTxtToMatrix([path '/Z/Z_' sFrame '.txt'], params.imgSize);
        u_ground = readTxtToMatrix([path '/u/u_' sFrame '.txt'], params.imgSize);
        v_ground = readTxtToMatrix([path '/v/v_' sFrame '.txt'], params.imgSize);
        w_ground = readTxtToMatrix([path '/w/w_' sFrame '.txt'], params.imgSize);

        % Calculate the X,Y amd uProj and vProj values for ground truth.
        [projMatrices, focals, principals] = buildProjMatrices(infoCameras, input.numLevels, input.factor);
        projMatrix = projMatrices{input.numLevels, 1};
        
        % - Calculate X and Y values per pixel in the 3D space.
        [X_ground, Y_ground] = mexGetXYfromZ(Z, focals{input.numLevels,1}, principals{input.numLevels,1}, infoCameras(1).COP);

        % - Calculate u and v projection to reference view from 3D flow.
        XYZ1_ground = [(X_ground(:)+u_ground(:))'; (Y_ground(:)+v_ground(:))'; (Z_ground(:)+w_ground(:))'];
        XYZ0_ground = [(X_ground(:))'; (Y_ground(:))'; (Z_ground(:))'];
        pt1 = projection3Dto2DAll(projMatrix, XYZ1_ground);
        pt0 = projection3Dto2DAll(projMatrix, XYZ0_ground);
        flow = pt1 - pt0;
        uProj_ground = reshape(flow(1,:), params.imgSize(1), params.imgSize(2));
        vProj_ground = reshape(flow(2,:), params.imgSize(1), params.imgSize(2));
        
        
        % - Calculate disparity map from Z values in 3D space.
        % Ground truth:
        d_ground = getDfromZ(Z_ground, focals{input.numLevels, 1}(1), abs(infoCameras(1).COP(1) - infoCameras(2).COP(1)));
        % Output:
        d = getDfromZ(Z, focals{input.numLevels, 1}(1), abs(infoCameras(1).COP(1) - infoCameras(2).COP(1)));
        
        numTotal = (params.imgSize(1)*params.imgSize(2));
        
        % > NRMS: Normalized Root Mean Square Error for 3D points (%).
        rmseXYZ = sqrt(sum(sum((X - X_ground).^2 + (Y - Y_ground).^2 + (Z - Z_ground).^2)) / numTotal);
        XYZ = [X_ground(:), Y_ground(:), Z_ground(:)];
        normXYZ = sqrt(sum(XYZ.^2,2));
        eval.NRMS_XYZ = rmseXYZ / (max(normXYZ) - min(normXYZ)) * 100;
        
        % > NRMS: Normalized Root Mean Square Error for 3D flow vectors (%).
        rmseUVW = sqrt(sum(sum((u - u_ground).^2 + (v - v_ground).^2 + (w - w_ground).^2)) / numTotal);
        UVW = [u_ground(:), v_ground(:), w_ground(:)];
        normUVW = sqrt(sum(UVW.^2,2));        
        normalize = max(normUVW) - min(normUVW);
        if(normalize < 1)
            normalize = 1;
        end
        eval.NRMS_UVW = rmseUVW / normalize * 100;
        
        % > AAE: Average Angular Error for projected 2D-uv flow (deg).
        w_estim = [uProj(:), vProj(:), ones(length(vProj(:)),1)];
        normEstim = sqrt(sum(w_estim.^2,2));
        w_truth = [uProj_ground(:), vProj_ground(:), ones(length(vProj_ground(:)),1)];
        normTruth = sqrt(sum(w_truth.^2,2));
        eval.AAE_uv = mean(radtodeg( acos( sum(w_estim.*w_truth, 2) ./ (normEstim.*normTruth) ) ));
        
        % > AAE: Average Angular Error for 3D flow (deg).
        w_estim3 = [u(:), v(:), w(:), ones(length(w(:)),1)];
        normEstim3 = sqrt(sum(w_estim3.^2,2));
        w_truth3 = [u_ground(:), v_ground(:), w_ground(:), ones(length(w_ground(:)),1)];
        normTruth3 = sqrt(sum(w_truth3.^2,2));
        eval.AAE_UVW = mean(radtodeg( acos( sum(w_estim3.*w_truth3, 2) ./ (normEstim3.*normTruth3) ) ));

        % > AEP_uv: Average End-point Error for image optical flow (mm).
        eval.AEP_uv = sqrt(sum(sum((uProj_ground - uProj).^2 + (vProj_ground - vProj).^2 )) / numTotal);
        
        % > AEP_d: Average End-point Error for image disparity (mm).
        eval.AEP_d = sqrt(sum(sum((d_ground - d).^2)) / numTotal);
    else
        disp(' ');
        disp('No ground truth available for evaluation');
        eval.NRMS_XYZ = -99;
        eval.NRMS_UVW = -99;
        eval.AAE_UVW = -99;
        eval.AAE_uv = -99;
        eval.AEP_uv = -99;
        eval.AEP_d = -99;
    end

end


