# Locate where we are
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/structs)

add_library(Structs SHARED structs.cpp)

