#!/bin/bash
# Build script for LDC Solver

set -e  # Exit on error

echo "======================================"
echo "Building LDC Solver"
echo "======================================"

# Detect macOS and set up GCC
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - looking for Homebrew GCC..."
    
    # Try to find g++-15, g++-14, g++-13, etc.
    for version in 15 14 13 12 11; do
        if command -v g++-$version &> /dev/null; then
            export CXX=g++-$version
            export CC=gcc-$version
            echo "Found and using: $CXX"
            break
        fi
    done
    
    if [ -z "$CXX" ]; then
        echo "WARNING: Homebrew GCC not found!"
        echo "Install with: brew install gcc"
        echo "Falling back to system compiler (OpenMP may not work)"
    fi
fi

# Parse arguments
BUILD_TYPE="Release"
BUILD_OPENMP=ON
BUILD_SERIAL=ON
CLEAN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --no-openmp)
            BUILD_OPENMP=OFF
            shift
            ;;
        --no-serial)
            BUILD_SERIAL=OFF
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --debug       Build in debug mode (default: Release)"
            echo "  --no-openmp   Don't build OpenMP version"
            echo "  --no-serial   Don't build serial version"
            echo "  --clean       Clean build directory before building"
            echo "  --help        Show this help message"
            echo ""
            echo "macOS users:"
            echo "  Requires Homebrew GCC for OpenMP support"
            echo "  Install with: brew install gcc"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="$PROJECT_DIR/build"

# Clean if requested
if [ "$CLEAN" = true ]; then
    echo "Cleaning build directory..."
    rm -rf "$BUILD_DIR"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Run CMake
echo ""
echo "Running CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DBUILD_OPENMP="$BUILD_OPENMP" \
    -DBUILD_SERIAL="$BUILD_SERIAL"

# Build
echo ""
echo "Building..."
make -j$(nproc)

echo ""
echo "======================================"
echo "Build complete!"
echo "======================================"
echo "Executables are in: $BUILD_DIR/bin/"
ls -lh "$BUILD_DIR/bin/" 2>/dev/null || echo "No executables built"
echo ""
