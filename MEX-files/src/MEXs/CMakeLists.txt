# Define some Matlab variables
set(MATLAB_MEX_LIBS libmex.so libmx.so libmat.so)
set(MATLAB_MEXFILE_EXT .mexa64)
set(MEX_FILES mexNonLocalSceneFlow mexGetXYfromZ mexBuildProjectionMaps mexBuildOcclusionMaps mexBuildOcclusionBinary mexBuildWarpedImages mexComputePartialsUnknowns mexBuildDiffusivitySmooth mexBuildRobustnessData mexBuildMatricesAb mexIterativeSolver mexProjection3Dto2D mexKMeans9D mexConnectedSegments mexMergedSegments mexComputeLocalDisparities mexCrossValidation mexComputeSegmentPlanes mexConsistentFill mexSegmentSimilarities mexMedianFilter mexGradient)

# Locate where Matlab includes are
include_directories(${MATLAB_ROOT}/extern/include)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/libs)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/structs)

foreach(mex_file ${MEX_FILES})

  # Append files to the coming library
  add_library(${mex_file} SHARED ${mex_file}.cpp)

  # Also link Matlab-Mex libraries  
  foreach(mex_lib ${MATLAB_MEX_LIBS})
    target_link_libraries(${mex_file} ${MATLAB_ROOT}/bin/glnxa64/${mex_lib})
  endforeach(mex_lib)
  target_link_libraries(${mex_file} NonLocalSceneFlow) 

  # Define the future mex file format
  set_target_properties(${mex_file} PROPERTIES PREFIX "" LINKER_LANGUAGE CXX)
  set_target_properties(${mex_file} PROPERTIES SUFFIX ${MATLAB_MEXFILE_EXT})

endforeach(mex_file)
