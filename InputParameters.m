%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              .:: Master thesis::.              %
% Title: 3D Scene Flow with a rigid motion prior %
% Author: Pau Panareda Busto                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Class file with all input parameters.

classdef InputParameters
    
    properties
        
        % Global Paremeters:
        % - read whole sequence? Or just first pair of frames...
        onSequence = false;
        % - is it real or artificial scenario?
        isRealScenario = true;
        %  Name of the input dataset.
        dataset = 'Daimler'; imgExt = '.pgm';
        % - Number of cameras used in the sequence.
        numCameras = 2;
        % - Fixes input images into a default size.
        defaultImgSize = [-1 -1];
        % - Stores execution information in a 'results' subfolder.
        onStore = true;
        % - Visualizes the resulting figures.
        onDisplay = true;
        % - Shows partial results of execution.
        onPrint = false;
        onFigure = false;
        onResidual = false;
        % 3D Scene Flow Parameters:
        % - Multiresolution pyramid.
        numLevels = 10;      % Number of levesl
        initLevel = 1;      % Coarsest level
        endLevel = 10;       % Finest level
        factor = 0.75;         % Downsampling factor (0.93)
        % - 2 Fixed point iterations.
        maxOut = 5;             % Maximum #loops in outer iteration (2)
        maxIn = 5;              % Maximum #loops in inner iteration (15)
        epsOut = -0.001;         % epsilon (cutoff value) to exit outer loops (0.01)
        epsIn = -0.001;          % epsilon (cutoff value) to exit outer loops (0.01)
        % - Initialisation of unknowns.
        % TODO: 1000 seems to work fine for first pyramid level. We need to
        % find a proper initialization of Z value and a proper upsampling for
        % subsequent levels.
        Zinit = 50;           % Z initialisation.
        onZinit = false;        % Disparity algorithim inits Z or not.
        % - Occlusions parameters.
        onOcclusions = true;
        % - Weights for divergence coefficients.
        alphaZ = 2;
        alphaUV = 30;
        alphaW = 30;
        mu = 1.0;
        muZ = 1.0;
        muUVW = 1.0;
        % - iterative solver parameters.
        itSolver = 'SOR';      % Type of iterative solver
        preconditioner = 'NO'; % Type of preconditioner (SSOR, ILU, ILUT, ILUM)
        ILUp = 0;
        maxIter = 20;               % Stablish max. #iterations
        weightSOR = 1.8;           % Weight factor
        errorSolver = 1e-2;         % Maximum relative residual
        % - Epsilon Psi robust functions
        epsSmooth = 0.0001;         % Epsilon for diffusivity
        epsData = 0.0001;           % Epsilon for data robust function

        %%% SEGMENTATION BASED %%%
        % - Segmentation:
        onSegmentationBased = true;
        numK = 100;
        numIterKMeans = 5;
        onConnectedComponents = true; % Segmentation with superpixels
        mergedThreshold = 0.1; % 10% of overall difference
        % - Disparity:
        typeMatchCost = 'SAD';
        maxDisparity = 0.20;  % 20-25% of pixel width)
        % - Plane fitting:
        maxIterPlane = 5;
        % - Penalties among segments:
        pImpact = 0.0;
        onImpact = true;
        % Must be half size -> (w*2+1) to ensure odd number.
        aggregatedWindowSize = 0.05;
        dispLambda = 0.01;
        %%% End SEGMENTATION BASED %%%
    end
        
end
