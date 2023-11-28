cmake_minimum_required(VERSION 3.16)
string(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
project(${CURRENT_FOLDER} LANGUAGES CXX)

# set
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} /Od /wd4819 /wd4068 /favor:INTEL64")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} /O2 /wd4819 /wd4068 /favor:INTEL64")

# file
file(GLOB_RECURSE CPP_FILES "particle_generate.cpp")

# include
include_directories(D:/Code/Eigen/eigen-3.4.0) # # Eigen directory

# target
add_executable(${PROJECT_NAME} ${CPP_FILES})