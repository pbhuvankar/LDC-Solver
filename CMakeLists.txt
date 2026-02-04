cmake_minimum_required(VERSION 3.12)
project(LDC_Solver VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect macOS and set up GCC
if(APPLE)
    message(STATUS "macOS detected - looking for Homebrew GCC")
    
    # Try to find g++-15, g++-14, g++-13, etc.
    find_program(GCC_EXECUTABLE 
        NAMES g++-15 g++-14 g++-13 g++-12 g++-11
        PATHS /usr/local/bin /opt/homebrew/bin
    )
    
    if(GCC_EXECUTABLE)
        set(CMAKE_CXX_COMPILER ${GCC_EXECUTABLE})
        message(STATUS "Found GCC: ${GCC_EXECUTABLE}")
    else()
        message(WARNING "Homebrew GCC not found. Install with: brew install gcc")
        message(WARNING "OpenMP may not be available with Apple Clang")
    endif()
endif()

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Options
option(BUILD_OPENMP "Build with OpenMP support" ON)
option(BUILD_SERIAL "Build serial version" ON)
option(BUILD_TEST "Build test version" ON)
option(ENABLE_PROFILING "Enable detailed profiling" OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

if(ENABLE_PROFILING)
    add_compile_definitions(ENABLE_PROFILING)
endif()

# Find OpenMP
if(BUILD_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found - building parallel version")
    else()
        message(WARNING "OpenMP not found - only serial version will be built")
        if(APPLE)
            message(WARNING "On macOS, make sure you have: brew install gcc")
            message(WARNING "Or set CXX environment variable: export CXX=g++-15")
        endif()
        set(BUILD_OPENMP OFF)
    endif()
endif()

# Add subdirectory for C++ source
add_subdirectory(src/cpp/core)

# Print configuration summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "LDC Solver Configuration Summary")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OpenMP support: ${BUILD_OPENMP}")
message(STATUS "Serial version: ${BUILD_SERIAL}")
message(STATUS "Test version: ${BUILD_TEST}")
message(STATUS "Profiling: ${ENABLE_PROFILING}")
if(APPLE)
    message(STATUS "Platform: macOS")
endif()
message(STATUS "========================================")
message(STATUS "")
