cmake_minimum_required(VERSION 3.10)

# Project name
project(AllanDeviationLab)

# Enable export of compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(AllanDeviationLab src/main.cpp)

