# macOS Build Instructions

## TL;DR for macOS Users

```bash
# 1. Install Homebrew GCC (if not already installed)
brew install gcc

# 2. Build with automatic GCC detection
make

# 3. Run
cd examples/basic_cavity
export OMP_NUM_THREADS=8
../../build/bin/ldc_openmp
```

That's it! The build system now automatically detects and uses Homebrew GCC.

---

## Why macOS is Different

Apple's default `clang` compiler **does not support OpenMP**. To get OpenMP support on macOS, you need to:

1. Install Homebrew GCC: `brew install gcc`
2. Use `g++-15` (or `g++-14`, `g++-13`, etc.) instead of `clang`

## What We Fixed

### Automatic GCC Detection

The CMake build system now:
1. Detects when running on macOS
2. Searches for Homebrew GCC (`g++-15`, `g++-14`, etc.)
3. Automatically uses it if found
4. Shows helpful error messages if not found

### Your Compiler Setup

You mentioned you have `g++-15` working. Perfect! The build system will find and use it automatically.

```bash
# This now works without any extra configuration:
make

# Under the hood, it's doing:
# export CXX=g++-15
# cmake .. -DCMAKE_CXX_COMPILER=g++-15
```

## Build Options

### Option 1: Let It Auto-Detect (Recommended)

```bash
make
```

The build script will:
- Detect macOS
- Find your Homebrew GCC (`g++-15`)
- Configure CMake to use it
- Build both serial and OpenMP versions

### Option 2: Manually Specify Compiler

```bash
# Set compiler explicitly
export CXX=g++-15
export CC=gcc-15

# Then build
make
```

### Option 3: Direct CMake

```bash
mkdir build && cd build

# Explicitly tell CMake to use g++-15
cmake .. -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_C_COMPILER=gcc-15

# Build
make -j$(sysctl -n hw.ncpu)
```

### Option 4: Use Your Original Makefile

Your original Makefile (with `CXX = g++-15`) works great! You can use it directly:

```bash
cd src/cpp/core

# Copy your Makefile
cp /path/to/your/Makefile .

# Build OpenMP version
make openmp

# Run
./ldc
```

## Verifying Your Build

### Check Compiler

```bash
# See which compiler CMake chose
cmake .. 2>&1 | grep "C++ compiler:"

# Should show something like:
# -- C++ compiler: /opt/homebrew/bin/g++-15
```

### Check OpenMP

```bash
# Look for this message during build:
# -- OpenMP found - building parallel version

# Check the executables
ls -lh build/bin/
# Should see both:
# ldc_serial
# ldc_openmp
```

### Test OpenMP

```bash
cd examples/basic_cavity

# Test with different thread counts
for threads in 1 2 4 8; do
    echo "Testing with $threads threads..."
    export OMP_NUM_THREADS=$threads
    time ../../build/bin/ldc_openmp
done
```

## Troubleshooting

### "OpenMP not found" Warning

**Problem:**
```
CMake Warning at CMakeLists.txt:XX (message):
  OpenMP not found - only serial version will be built
```

**Solution:**

1. **Install Homebrew GCC:**
   ```bash
   brew install gcc
   ```

2. **Check installation:**
   ```bash
   which g++-15
   # Should show: /opt/homebrew/bin/g++-15 (or /usr/local/bin/g++-15)
   ```

3. **Clean and rebuild:**
   ```bash
   make clean
   make
   ```

### Compiler Not Found

**If you see: "Homebrew GCC not found"**

```bash
# Install it
brew install gcc

# Check what version was installed
ls /opt/homebrew/bin/g++-*

# You should see: g++-15 (or g++-14, etc.)
```

### Apple Clang Instead of GCC

**If CMake is using Apple Clang:**

```bash
# Check current compiler
which g++
# If it shows /usr/bin/g++, that's Apple Clang (symlink)

# Force CMake to use Homebrew GCC
rm -rf build
export CXX=/opt/homebrew/bin/g++-15
export CC=/opt/homebrew/bin/gcc-15
make
```

### Wrong GCC Version

**If you have multiple GCC versions:**

```bash
# See all installed versions
ls /opt/homebrew/bin/g++-*

# Pick your preferred version
export CXX=/opt/homebrew/bin/g++-15
make clean
make
```

## Performance Notes for macOS

### Apple Silicon (M1/M2/M3)

If you're on Apple Silicon:

```bash
# Use all performance cores
export OMP_NUM_THREADS=8  # Adjust based on your Mac

# Set thread affinity for better performance
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Run
../../build/bin/ldc_openmp
```

### Intel Mac

```bash
# Use all cores
export OMP_NUM_THREADS=$(sysctl -n hw.ncpu)

# Run
../../build/bin/ldc_openmp
```

## Comparison: CMake vs Your Makefile

### Your Original Makefile

**Pros:**
- Simple and direct
- Explicitly uses `g++-15`
- Easy to understand
- Has useful targets (`openmp`, `debug`, `profile`)

**Cons:**
- Manual compilation only
- No automatic dependency tracking
- Files must be in same directory

### Our CMake Build

**Pros:**
- Organized directory structure
- Automatic dependency management
- Builds multiple versions simultaneously
- Cross-platform (Linux, macOS, Windows)
- Automatic GCC detection on macOS

**Cons:**
- More complex
- Requires CMake installation

### Best of Both Worlds

You can keep **both**:

1. **Use CMake** for the full project structure:
   ```bash
   make
   ```

2. **Use your Makefile** for quick iterations during development:
   ```bash
   cd src/cpp/core
   make openmp
   ```

## macOS-Specific CMake Configuration

If you want more control, create a `build-macos.sh` script:

```bash
#!/bin/bash
# build-macos.sh - macOS-specific build

# Find Homebrew GCC
export CXX=$(which g++-15 || which g++-14 || which g++-13)
export CC=$(echo $CXX | sed 's/g++/gcc/')

echo "Using compiler: $CXX"

# Create build directory
mkdir -p build
cd build

# Configure with macOS-specific settings
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_C_COMPILER=$CC \
    -DBUILD_OPENMP=ON \
    -DBUILD_SERIAL=ON

# Build using all cores
make -j$(sysctl -n hw.ncpu)

echo ""
echo "Build complete!"
echo "Executables in: build/bin/"
ls -lh bin/
```

Make it executable and run:
```bash
chmod +x build-macos.sh
./build-macos.sh
```

## Summary

✅ **The build system now automatically detects and uses your `g++-15`**  
✅ **No manual configuration needed**  
✅ **Just run `make` and it works**  
✅ **Both serial and OpenMP versions built**  

Your original Makefile approach was correct - we've now integrated that wisdom into the CMake build system so it "just works" on macOS!

## Quick Reference

```bash
# First time setup
brew install gcc

# Build (auto-detects GCC)
make

# Run
cd examples/basic_cavity
export OMP_NUM_THREADS=8
../../build/bin/ldc_openmp

# Clean rebuild
make clean
make

# Manual compiler selection (if needed)
export CXX=g++-15
make
```

Happy computing! 🚀
