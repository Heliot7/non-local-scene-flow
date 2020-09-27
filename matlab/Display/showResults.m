%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show results of output 3D points and disparities.
% INPUT PARAMETERS:
% - 4 method unknowns: Z, u, v and w.

function figs = showResults(images, infoCameras, Z, u, v, w, runtime, frame)

    input = InputParameters;
    
    imgRows = length(Z(:,1));
    imgCols = length(Z(1,:));
    step = imgRows/(imgRows/4);
    
    [projMatrices, focals, principals] = buildProjMatrices(infoCameras, input.numLevels, input.factor);
    projMatrix = projMatrices{input.endLevel, 1};
    % - Calculate X and Y values per pixel in the 3D space.
    [X, Y] = mexGetXYfromZ(Z, focals{input.endLevel,1}, principals{input.endLevel,1}, infoCameras(1).COP);
    % - Calculate u and v projection to reference view from 3D flow.
    uProj = zeros(imgRows,imgCols);
    vProj = zeros(imgRows,imgCols);
    for r = 1:imgRows
        for c = 1:imgCols
            pt1 = mexProjection3Dto2D(projMatrix, [X(r,c)+u(r,c); Y(r,c)+v(r,c); Z(r,c)+w(r,c)]);
            pt0 = mexProjection3Dto2D(projMatrix, [X(r,c); Y(r,c); Z(r,c)]);
            flow = pt1 - pt0;
            uProj(r,c) = flow(1);
            vProj(r,c) = flow(2);
        end
    end
    
    show = 'on';
    if(~input.onDisplay)
        show = 'off';
    end
    
    % Show unknown Z.
    figs.outZ = figure('visible', show);
    set(figs.outZ, 'Name', 'Dense disparity map');
    [xAxis, yAxis] = meshgrid(1:imgCols, 1:imgRows);
    tri = delaunay(xAxis, yAxis);
    trisurf(tri, xAxis, yAxis, Z(:), 'EdgeColor','none');
    title(sprintf('%s [Dense disparity map]', input.dataset));
    colorbar;
    axis off;
    campos([imgCols/2, imgRows/2, min(Z(:))]);
    camtarget([imgCols/2, imgRows/2, max(Z(:))]);
    
    % Show unknown u.
    figs.outU = figure('visible', 'off');
    set(figs.outU, 'Name', 'u - horizontal flow');
    imagesc(u);
    title(sprintf('%s [u - horizontal flow]', input.dataset));
    colorbar;
    axis image;

    % Show unknown v.
    figs.outV = figure('visible', 'off');
    set(figs.outV, 'Name', 'v - vertical flow');
    imagesc(v);
    title(sprintf('%s [v - horizontal flow]', input.dataset));
    colorbar;
    axis image;

    % Show unknown w.
    figs.outW = figure('visible', 'off');
    set(figs.outW, 'Name', 'w - depth flow');
    imagesc(w);
    title(sprintf('%s [w - horizontal flow]', input.dataset));
    colorbar;
    axis image;

    % Show 3D flow.
    figs.outUVW = figure('visible', 'off');
    set(figs.outUVW, 'Name', '3D flow');
    img = imresize(images{1}, input.factor^(input.numLevels - input.endLevel));
    hold on;
    surface('XData',[0 imgCols; 0 imgCols],'YData',[0 0; imgRows imgRows], ...
        'ZData',[max(Z(:)) max(Z(:)); max(Z(:)) max(Z(:))],'CData',img, ...
        'FaceColor','texturemap','EdgeColor','none');
    colormap gray;
    axis image;
    campos([imgCols/2, imgRows/2, min(Z(:))]);
    camtarget([imgCols/2, imgRows/2, max(Z(:))]);
    camroll(0);
    Zstep = Z(1:step:imgRows, 1:step:imgCols);
    Ustep = u(1:step:imgRows, 1:step:imgCols);
    Vstep = v(1:step:imgRows, 1:step:imgCols);
    Wstep = w(1:step:imgRows, 1:step:imgCols);
    quiver3(1:step:imgCols, 1:step:imgRows, Zstep, Ustep, Vstep, Wstep, 'color','r');
	hold off;

    % Show projected flow onto the reference view.
    figs.outProjFlow = figure('visible', show);
    set(figs.outProjFlow, 'Name', '2D flow projected');
    imagesc( imresize(images{1}, input.factor^(input.numLevels - input.endLevel)));
    axis image;
    colormap gray;
    hold on;
    vProjStep = uProj(1:step:imgRows, 1:step:imgCols);
    uProjStep = vProj(1:step:imgRows, 1:step:imgCols);
    quiver(1:step:imgCols, 1:step:imgRows, vProjStep, uProjStep, 'color', 'r');
    hold off;

    % Show 3D points.
    figs.out3D = figure('visible', 'off');
    set(figs.out3D, 'Name', '3D points');
    set(figs.out3D, 'Position', [50 250 1000 800 ] );
	plot3(X(1:step:imgRows*imgCols), Y(1:step:imgRows*imgCols), Z(1:step:imgRows*imgCols),'.','color','k');
    axis equal;
    campos([imgCols/2, imgRows/2, min(Z(:))]);
    camtarget([imgCols/2, imgRows/2, max(Z(:))]);
    camroll(180);

    % Save results in a folder.
    params.time = runtime;
    params.imgSize = [length(images{1,1}(:,1)) length(images{1,1}(1,:))];
    params.frame = frame;
    disp('... evaluating result...');
    eval = computeEvaluation(X, Y, Z, u, v, w, uProj, vProj, infoCameras, params);
    fprintf('\n[3D ERROR]\n');
    fprintf('NRMS_XYZ = %f %%\n', eval.NRMS_XYZ);
    fprintf('NRMS_UVW = %f %%\n', eval.NRMS_UVW);
    fprintf('AAE_UVW = %f deg\n', eval.AAE_UVW);
    fprintf('\n[2D ERROR]\n');
    fprintf('AEP_uv = %f mm\n', eval.AEP_uv);
    fprintf('AEP_d = %f mm\n', eval.AEP_d);
    fprintf('AAE_uv = %f deg\n', eval.AAE_uv);
    if(input.onStore)
        disp('... saving output...');
        saveOutput(Z, u, v, w, figs, eval, params);
    end
end

