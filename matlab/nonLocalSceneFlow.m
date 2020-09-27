%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D Scene Flow method from Basha et al (2010) publication.
% INPUT PARAMETERS:
% - Images for frames 't' and 't+1'
% - Information involving all cameras
% OUTPUT PARAMETERS:
% - Dense disparity map (depth Z) w.r.t reference view
% - 3D disparity map w.r.t. reference view (u,v,w components)

function [Z, u, v, w] = nonLocalSceneFlow(images, infoCameras)

    % Input parameters
    input = InputParameters;

    % Create multi-resolution pyramid for all camera input images.
    [t0ImgPyramid, t1ImgPyramid] = buildPyramid(images, input.numLevels, input.factor);
    
    % Create pyramid for all project matrices.
    [projMatrices, focals, principals] = buildProjMatrices(infoCameras, input.numLevels, input.factor);
    
    imgRows = length(t0ImgPyramid{input.initLevel,1}(:,1,1));
    imgCols = length(t0ImgPyramid{input.initLevel,1}(1,:,1));
    imgChannels = length(t0ImgPyramid{input.initLevel,1}(1,1,:));

    % Depth (Z) initialisation.
    Z = zeros(imgRows,imgCols);
    Z(:) = deal(input.Zinit);
    % 3D-flow (u,v,w) initialisation.
    u = zeros(imgRows,imgCols);
    v = zeros(imgRows,imgCols);
    w = zeros(imgRows,imgCols);
   
    %%% DISPLAY %%%%
    % Test Z results in outer loop.
    if(input.onFigure)
        figLvl = figure;
        set(figLvl, 'Position', [50 250 1000 800 ] );
        set(figLvl,'Name','Z [Level]');
        figOut = figure;
        set(figOut, 'Position', [50 250 1000 800 ] );
        set(figOut,'Name','Z [Outer Loop]');
        figIn = figure;
        set(figIn, 'Position', [50 250 1000 800 ] );
        set(figIn,'Name','dZ [Inner Loop]');
    end
    if(input.onResidual)
        fResAcc = figure;
        posResidual = 1;
        listResidual = zeros(1, input.maxIn*input.maxOut*(input.initLevel-input.endLevel+1));
        posResidualAcc = 1;
    end
    %%% End DISPLAY %%%
    
    % Go along all pyramid levels.
    for LEVEL = input.initLevel:input.endLevel
        % Define some variables depending on pyramid level.
        dZ = zeros(imgRows, imgCols);       % Increment of Z
        du = zeros(imgRows, imgCols);       % Increment of u
        dv = zeros(imgRows, imgCols);       % Increment of v
        dw = zeros(imgRows, imgCols);       % Increment of w
        DdZ = zeros(imgRows, imgCols);      % Delta variation of Z incr
        Ddu = zeros(imgRows, imgCols);      % Delta variation of u incr
        Ddv = zeros(imgRows, imgCols);      % Delta variation of v incr
        Ddw = zeros(imgRows, imgCols);      % Delta variation of w incr
    
        %%% MESSAGE %%%
        printText(true, '\nPyramid level: %d, Img size: %dx%d\n', [LEVEL, imgRows, imgCols])
        %%% End MESSAGE %%%
        
%         figure; imagesc(Z); % caxis([300 800]);

        %%% SEGMENTATION BASED APPROACH %%%
        if(input.onSegmentationBased)
            
            [Z, segmentPenalties] = depthEstimation(t0ImgPyramid{LEVEL, 1:2}, Z, u, v, w, ...
                abs(infoCameras(1).COP(1) - infoCameras(2).COP(1)), focals{LEVEL, 1}(1), projMatrices{LEVEL,1}(3,3), LEVEL);
            
        elseif(input.onZinit && LEVEL == input.initLevel)
            
            empty = zeros(length(t0ImgPyramid{input.endLevel,1}(:,1,1)),length(t0ImgPyramid{input.endLevel,1}(1,:,1)));
            depth = depthEstimation(t0ImgPyramid{input.numLevels, 1:2}, empty, empty, empty, empty, ...
                abs(infoCameras(1).COP(1) - infoCameras(2).COP(1)), focals{input.numLevels, 1}(1), 1, LEVEL);
            % Adapt estimation to current downsampling factor
            Z = imresize(depth, input.factor^(abs(input.numLevels - LEVEL) ));
            segmentPenalties = ones(imgRows, imgCols, 4);
            
        else
            segmentPenalties = ones(imgRows, imgCols, 4);
        end
        %%% End SEGMENTATION BASED APPROACH %%%

        % 2 Fixed-point iterations (Outer & Inner loops).        
        % >>> OUTER LOOP <<<
        itOut = 1;
        while( itOut <= input.maxOut && input.epsOut < max([norm(dZ) norm(du) norm(dv) norm(dw)]) )    

            %%% MESSAGE %%%
            printText(true, '- OUTER -> %d\n', itOut);
            %%% End MESSAGE %%%
            
            % Calculate X and Y values from current unknown Z.
            [X, Y] = mexGetXYfromZ(Z, focals{LEVEL,1}, principals{LEVEL,1}, infoCameras(1).COP);
            
            % Compute projection maps w.r.t. reference view pixels.
            % t0ProjImgs => [CAM, ROW, COL, (x,y)]
            [t0ProjImgs, t1ProjImgs] = mexBuildProjectionMaps(projMatrices(LEVEL,:), X, Y, Z, u, v, w);
            
            % Computation of occlusion maps w.r.t. reference view pixels.
            [t0occlusionMaps, t1occlusionMaps] = mexBuildOcclusionMaps(t0ProjImgs, t1ProjImgs, ...
                X, Y, Z, u, v, w, {infoCameras(:).COP});
            occlusionBool = mexBuildOcclusionBinary(t0ProjImgs, t1ProjImgs, t0occlusionMaps, t1occlusionMaps);

            % Calculate warped images from all views to reference view.
            [t0WarpedImgs, t1WarpedImgs] = mexBuildWarpedImages(t0ImgPyramid(LEVEL,:), ...
                t1ImgPyramid(LEVEL,:), t0ProjImgs, t1ProjImgs, t0occlusionMaps, t1occlusionMaps);
            
            % Smoothing of warped images (done in Basha's code)
%             smoothFilter = fspecial('gaussian');
%             for CAMi = 1:input.numCameras                
%                 t0WarpedImgs(CAMi,:,:) = imfilter(t0WarpedImgs(CAMi,:,:), smoothFilter, 'same');
%                 t1WarpedImgs(CAMi,:,:) = imfilter(t1WarpedImgs(CAMi,:,:), smoothFilter, 'same');
%             end

            % Calculate gradient of warped images.
            t0WarpedImgs_x = zeros(input.numCameras, imgRows, imgCols, imgChannels);
            t0WarpedImgs_y = zeros(input.numCameras, imgRows, imgCols, imgChannels);
            t1WarpedImgs_x = zeros(input.numCameras, imgRows, imgCols, imgChannels);
            t1WarpedImgs_y = zeros(input.numCameras, imgRows, imgCols, imgChannels);
            for CAMi = 1:input.numCameras
                for CHNi = 1:imgChannels
                    warpImg = reshape(t0WarpedImgs(CAMi,:,:,CHNi), imgRows, imgCols);
                    [t0WarpedImgs_x(CAMi,:,:,CHNi), t0WarpedImgs_y(CAMi,:,:,CHNi)] = mexGradient(warpImg);
                    warpImg = reshape(t1WarpedImgs(CAMi,:,:,CHNi), imgRows, imgCols);
                    [t1WarpedImgs_x(CAMi,:,:,CHNi), t1WarpedImgs_y(CAMi,:,:,CHNi)] = mexGradient(warpImg);
                end
            end
            
            % Calculate unknown derivatives
            [Z_x, Z_y] = mexGradient(Z);
            [u_x, u_y] = mexGradient(u);
            [v_x, v_y] = mexGradient(v);
            [w_x, w_y] = mexGradient(w);
            
            % Calculate partial derivative_unknowns for all values I_i.
            % Note: Per each pixel (at each view), we store the Z partial
            % derivative for frame t and all 4 unknowns for frame t+1.
            [t0I_unknown, t1I_unknown] = mexComputePartialsUnknowns(projMatrices(LEVEL,:), ...
                focals{LEVEL,1}, principals{LEVEL,1}, infoCameras(1).COP, ...
                Z, Z_x, Z_y, u, u_x, u_y, v, v_x, v_y, w, w_x, w_y, ...
                t0WarpedImgs_x, t0WarpedImgs_y, t1WarpedImgs_x, t1WarpedImgs_y);
            
            % >>> INNER LOOP <<<
            itIn = 1;
            % Initalise increments to 0.
            dZ(:) = deal(0.0); du(:) = deal(0.0); dv(:) = deal(0.0); dw(:) = deal(0.0);
            while( itIn <= input.maxIn && input.epsIn < max([norm(DdZ) norm(Ddu) norm(Ddv) norm(Ddw)]) )
                
                %%% MESSAGES %%%
                printText(true, '- - INNER -> %d\n', itIn);
                %%% End MESSAGES %%%g 
                
                % Calculate increment derivatives
                [dZ_x, dZ_y] = mexGradient(dZ);
                [du_x, du_y] = mexGradient(du);
                [dv_x, dv_y] = mexGradient(dv);
                [dw_x, dw_y] = mexGradient(dw);
                
                % Define occlusions for the equations with new increments.
%                 [X, Y] = mexGetXYfromZ(Z+dZ, focals{LEVEL,1}, principals{LEVEL,1}, infoCameras(1).COP);
%                 [t0ProjImgs, t1ProjImgs] = mexBuildProjectionMaps(projMatrices(LEVEL,:), X, Y, Z, u, v, w);
%                 occlusionBool = mexBuildOcclusionBinary(t0ProjImgs, t1ProjImgs, t0occlusionMaps, t1occlusionMaps);

                % Calculate the coefficient of diffusivity at smooth equation part.
                psiSmooth = mexBuildDiffusivitySmooth(Z_x, Z_y, u_x, u_y, v_x, v_y, w_x, w_y, ...
                    dZ_x, dZ_y, du_x, du_y, dv_x, dv_y, dw_x, dw_y, input);

                % Compute Robustness functions Psi'() with fixed increments for data term.
                psiData = mexBuildRobustnessData(dZ, du, dv, dw, t0WarpedImgs, ...
                    t1WarpedImgs, t0I_unknown, t1I_unknown, input);
            
%                 profile on;
                
                % Build Ax = B (with fixed Psi' with new updated incr. in every loop).
%                 for i = 1:10
                [matrixA, matrixB] = mexBuildMatricesAb(Z, u, v, w, t0WarpedImgs, t1WarpedImgs, ...
                    t0I_unknown, t1I_unknown, psiData, psiSmooth, segmentPenalties, occlusionBool, LEVEL, input);
%                 end

%                 profile off;
%                 profile viewer;
%                 keyboard;

%                 for i = 1:10
                % Iterative solver for all pixels.
                [NEWdZ, NEWdu, NEWdv, NEWdw, evolResidual] = ...
                    mexIterativeSolver(matrixA, matrixB, psiSmooth, segmentPenalties, LEVEL, input);
%                 end

                % Show residual evolution over the iterations at this loop.
                if(input.onResidual)
                    evolResidual = evolResidual(evolResidual(:)~=-1);
                    listResidual(posResidual) = evolResidual(length(evolResidual));
                    figure(fResAcc);
                    hold on;
                    % axis([0, input.maxIter, 0, 1]);
                    plot(posResidualAcc:posResidualAcc+length(evolResidual)-1, evolResidual);
                    if(itIn == 1)
                        txtAcc = [ num2str(evolResidual(1)) ' % -> ' num2str(listResidual(posResidual)) ' %']; %...
                            %' (' num2str(LEVEL) ':' num2str(itOut) ')' ];
                        text(posResidualAcc, evolResidual(1), txtAcc);
                    end
                    hold off;
                    posResidual = posResidual + 1;
                    posResidualAcc = posResidualAcc + length(evolResidual);
                end
               
                % Update values of increments.
                [DdZ, Ddu, Ddv, Ddw] = deal(abs(dZ - NEWdZ), abs(du - NEWdu), abs(dv - NEWdv), abs(dw - NEWdw));
                [dZ ,du, dv, dw] = deal(NEWdZ, NEWdu, NEWdv, NEWdw);
                
                %%% DISPLAY %%% 
                % Test dZ results in inner loop.
                if(input.onFigure)
                    figure(figIn);
                    subplot(4,5,itIn);
                    imagesc(dZ);
                    title(sprintf('itOut = %d, itIn = %d', itOut, itIn));
                    axis equal; colorbar; drawnow;
                end
                %%% End DISPLAY %%% 
                
                itIn = itIn + 1;                
            end
            
            % Update values of unknowns
            [Z, u, v, w] = deal(Z+dZ, u+du, v+dv, w+dw);            

            %%% DISPLAY %%% 
            % Test Z results in outer loop.
            if(input.onFigure)
                figure(figOut);
                subplot(4,5,itOut);
                imagesc(Z);
                title(sprintf('itOut = %d',itOut));
                axis equal; colorbar; drawnow;
            end
            %%% End DISPLAY %%%
            
            itOut = itOut + 1;
        end
        
        %%% DISPLAY %%%
        if(input.onFigure)
            figure(figLvl);
            subplot(4, 5, LEVEL);
            imagesc(Z);
            title(sprintf('LEVEL = %d', LEVEL));
            axis equal; colorbar; drawnow;
        end
        %%% End DISPLAY %%%
        
        % Upsample solutions of current); pyramid level to init finer one.
        if(LEVEL < input.endLevel)
            imgRows = length(t0ImgPyramid{LEVEL+1,1}(:,1,1));
            imgCols = length(t0ImgPyramid{LEVEL+1,1}(1,:,1));
            Z = imresize(Z, [imgRows imgCols]);
            u = imresize(u, [imgRows imgCols]);
            v = imresize(v, [imgRows imgCols]);
            w = imresize(w, [imgRows imgCols]);
        end
    end
    
end

