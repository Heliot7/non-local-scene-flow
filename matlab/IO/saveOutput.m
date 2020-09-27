%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saves program output information.
% INPUT PARAMETERS:
% - dataset name
% - unknowns output: Z, u, v, w
% - evaluation results
% - output figures

function saveOutput(Z, u, v, w, figs, eval, params)

    input = InputParameters;

    if(~input.isRealScenario)
        dataset = 'Plane';
    else
        dataset = input.dataset;
    end
    
    % If not done before, creates a new directory of the dataset.
    pathInit = '/work/panareda/results/';
    if(size(dir([pathInit dataset]), 1) == 0)
        mkdir(pathInit, dataset);
    end
    % Then makes a new folder for the current run with its exec time.
    % If already created, exits because there is already info.
    if(params.frame == 0)
        timeStamp = datestr(now, 'yyyy.mm.dd.HH.MM');
    else
        listing = dir([pathInit dataset '/']);
        timeStamp = listing(length(listing)).name;
    end
    if(size(dir([pathInit dataset '/' timeStamp]), 1) == 0 || params.frame > 0)
        
        if(params.frame == 0)
            mkdir([pathInit dataset '/'], timeStamp);
        end
        % Within the output folder:
        if(params.frame == 0)
            % - Folder with current source code
            % >>> matlab files
            mkdir([pathInit dataset '/' timeStamp], 'src/matlab');
            copyfile('matlab', [pathInit dataset '/' timeStamp '/src/matlab']);
            % >>> MEX-file connectors
            mkdir([pathInit dataset '/' timeStamp], 'src/MEX');
            copyfile('MEX-files/src/MEXs/mex*', [pathInit dataset '/' timeStamp '/src/MEX']);
            % >>> C++ files
            mkdir([pathInit dataset '/' timeStamp], 'src/C++');
            copyfile('MEX-files/src/libs/libSceneFlow', [pathInit dataset '/' timeStamp '/src/C++']);
        end
        
        % - Figures
        % >>> .fig and .png of depth and 3D flow (u,v,w separated and together)
        if(params.frame == 0)
            mkdir([pathInit dataset '/' timeStamp], 'figs');
        end
        figFolder = [pathInit dataset '/' timeStamp '/figs'];
        if(params.frame < 10)
            sFrame = [num2str(0) num2str(0) num2str(0) num2str(params.frame)];
        elseif(params.frame < 100)
            sFrame = [num2str(0) num2str(0) num2str(params.frame)];
        else
            sFrame = [num2str(0) num2str(params.frame)];
        end
        saveas(figs.outZ, [figFolder '/outZ_' sFrame '.fig']);
        saveas(figs.outZ, [figFolder '/outZ_' sFrame '.png']);
        saveas(figs.outU, [figFolder '/outU_' sFrame '.fig']);
        saveas(figs.outU, [figFolder '/outU_' sFrame '.png']);
        saveas(figs.outV, [figFolder '/outV_' sFrame '.fig']);
        saveas(figs.outV, [figFolder '/outV_' sFrame '.png']);
        saveas(figs.outW, [figFolder '/outW_' sFrame '.fig']);
        saveas(figs.outW, [figFolder '/outW_' sFrame '.png']);
        saveas(figs.outUVW, [figFolder '/outUVW_' sFrame '.fig']);
        saveas(figs.outProjFlow, [figFolder '/outUVW_' sFrame '.png']);
        if(~input.onDisplay)
            close all;
        end

        % - Txt with output
        dlmwrite([ pathInit dataset '/' timeStamp '/Z_' sFrame '.txt'], Z, 'delimiter', '\t');    
        dlmwrite([ pathInit dataset '/' timeStamp '/u_' sFrame '.txt'], u, 'delimiter', '\t');    
        dlmwrite([ pathInit dataset '/' timeStamp '/v_' sFrame '.txt'], v, 'delimiter', '\t');    
        dlmwrite([ pathInit dataset '/' timeStamp '/w_' sFrame '.txt'], w, 'delimiter', '\t');

        % - Txt with parameters
        fid = fopen([pathInit dataset '/' timeStamp '/parameters.txt'], 'w');
        fprintf(fid, '[Scene]\n');
        fprintf(fid, 'numCameras = %d\n', input.numCameras);
        fprintf(fid, 'imgSize = [%d, %d]\n', params.imgSize);
        fprintf(fid, 'Zinit = %d\n', input.onZinit);
        fprintf(fid, '[Pyramid]\n');
        fprintf(fid, 'numLevels = %d\n', input.numLevels);
        fprintf(fid, 'initLevel = %d\n', input.initLevel);
        fprintf(fid, 'endLevel = %d\n', input.endLevel);
        fprintf(fid, 'factor = %f\n', input.factor);
        fprintf(fid, '[Loops]\n');
        fprintf(fid, 'maxOut = %d\n', input.maxOut);
        fprintf(fid, 'maxIn = %d\n', input.maxIn);
        fprintf(fid, '[SOLVER]\n');
        fprintf(fid, 'itSolver = %s\n', input.itSolver);
        fprintf(fid, 'preconditioner = %s\n', input.preconditioner);
        fprintf(fid, 'maxIter = %d\n', input.maxIter);
        fprintf(fid, 'errorSolver = %f\n', input.errorSolver);
        fprintf(fid, '[Weights]\n');
        fprintf(fid, 'alphaZ = %f\n', input.alphaZ);
        fprintf(fid, 'alphaUV = %f\n', input.alphaUV);
        fprintf(fid, 'alphaW = %f\n', input.alphaW);
        fprintf(fid, 'muZ = %f\n', input.muZ);
        fprintf(fid, 'muUVW = %f\n', input.muUVW);
        if(input.onSegmentationBased)
            fprintf(fid, 'numK = %f\n', input.numK);
            fprintf(fid, 'numIterKMeans = %f\n', input.numIterKMeans);
            fprintf(fid, 'mergedThreshold = %f\n', input.mergedThreshold);
            fprintf(fid, 'maxDisparity = %f\n', input.maxDisparity);
            fprintf(fid, 'maxIterPlane = %f\n', input.maxIterPlane);
            fprintf(fid, 'pImpact = %f\n', input.pImpact);
            fprintf(fid, 'aggregatedWindowSize = %f\n', input.aggregatedWindowSize);
            fprintf(fid, 'dispLambda = %f\n', input.dispLambda);
        end
        fclose(fid);

        % - Txt with output evaluation
        fid = fopen([pathInit dataset '/' timeStamp '/evaluation_' sFrame '.txt'], 'w');
        % >>> Total time
        fprintf(fid, 'TOTAL Execution Time = %f sec\n', params.time);
        % >>> Evalutation parameters
        fprintf(fid, '\n[3D ERROR]\n');
        fprintf(fid, 'NRMS_XYZ = %f %%\n', eval.NRMS_XYZ);
        fprintf(fid, 'NRMS_UVW = %f %%\n', eval.NRMS_UVW);
        fprintf(fid, 'AAE_UVW = %f deg\n', eval.AAE_UVW);
        fprintf(fid, '\n[2D ERROR]\n');
        fprintf(fid, 'AEP_uv = %f mm\n', eval.AEP_uv);
        fprintf(fid, 'AEP_d = %f mm\n', eval.AEP_d);
        fprintf(fid, 'AAE_uv = %f deg\n', eval.AAE_uv);
        fclose(fid);
    else
        disp('Output not stored. There are already previous results!');
    end
    
end
