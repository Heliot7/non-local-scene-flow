# Version used in CMake
# Tells older version that cannot compile and newer ones that some functionalities are not available
cmake_minimum_required (VERSION 2.6)

# Name given to the project itself
# Used for internal purposes such as reference within CMakeList files (i.e. ${TEST_SOURCE_DIR})
project(NONLOCAL-SCENEFLOW-MEX)

# Define some Matlab variables
set(MATLAB_ROOT /opt/matlab)

# Added directories to this root project file
add_subdirectory(structs)
add_subdirectory(libs)
add_subdirectory(MEXs)

# Add those library directories which will be needed to compile
include_directories(${MATLAB_ROOT}/extern/include)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/libs)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/MEXs)
include_directories(${NONLOCAL-SCENEFLOW-MEX_SOURCE_DIR}/structs)

# Link these libraries
link_directories(${NONLOCAL-SCENEFLOW-MEX_BINARY_DIR}/libs)
link_directories(${NONLOCAL-SCENEFLOW-MEX_BINARY_DIR}/MEXs)
link_directories(${NONLOCAL-SCENEFLOW-MEX_BINARY_DIR}/structs)

# Configures the final executable which will trigger the added file
add_executable(NonLocalSceneFlow-Mex main.cpp)
# Links the libraries to this executable
target_link_libraries(NonLocalSceneFlow-Mex NonLocalSceneFlow IO mexNonLocalSceneFlow)

## Compiler flags
#if(CMAKE_COMPILER_IS_GNUCXX)
#  set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
#  set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
#endif()

# DEBUG MODE
set(CMAKE_BUILD_TYPE Debug)
