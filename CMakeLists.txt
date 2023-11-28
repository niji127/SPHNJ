cmake_minimum_required(VERSION 3.26)
string(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
project(${CURRENT_FOLDER} LANGUAGES CUDA CXX)

find_package(CUDA REQUIRED)

include(CheckLanguage)
check_language(CUDA)

# file
file(GLOB_RECURSE CUDA_FILES "src/*.cu")
file(GLOB_RECURSE CPP_FILES "src/*.cpp")

# Set compiler options
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_90 -w")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS} -Xptxas -O3,-dlcm=cg")

# set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -Xptxas -v")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -use_fast_math")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -Xcompiler /O2")

# include
include_directories(${PROJECT_SOURCE_DIR}/deps/Eigen/eigen-3.4.0)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CUDA_INCLUDE_DIRS})

# Add the executable
add_executable(${CURRENT_FOLDER} ${CPP_FILES} ${CUDA_FILES})

# add_executable(${CURRENT_FOLDER} src/test.cu)

# set executable path
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)
target_link_libraries(${CURRENT_FOLDER} Thrust)

# Set CUDA specific properties
set_target_properties(${CURRENT_FOLDER} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CXX_STANDARD 17
)