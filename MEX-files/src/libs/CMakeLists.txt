# Locate where Matlab includes are
include_directories(${MATLAB_ROOT}/extern/include)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/libs)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/structs)

add_library(IO SHARED libIO/lodepng.cpp libIO/pgm.cpp libIO/file.cpp libIO/viewMatlab.cpp)

add_library(Maths SHARED libMath/matrix.cpp libMath/metric.cpp)
target_link_libraries(Maths Structs)

add_library(ImgProc SHARED libImgProc/differentiation.cpp libImgProc/filters.cpp libImgProc/segmentation.cpp libImgProc/disparity.cpp libImgProc/colour.cpp)
target_link_libraries(ImgProc Maths)

add_library(NonLocalSceneFlow SHARED libSceneFlow/nonLocalSceneFlow.cpp)
target_link_libraries(NonLocalSceneFlow ImgProc IO)
target_link_libraries(NonLocalSceneFlow ${MATLAB_ROOT}/bin/glnxa64/libeng.so)
