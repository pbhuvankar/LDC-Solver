# 2D Lid-Driven Cavity Flow Solver

A high-performance 2D incompressible Navier-Stokes solver for the classic lid-driven cavity benchmark problem, featuring both Python prototype and optimized C++ implementations with parallel computing support.

## Features

- **Dual Implementation**: Python prototype for rapid development and production C++ code
- **Multiple Numerical Schemes**: Explicit and semi-implicit time-stepping
- **Parallel Computing**: OpenMP support for multi-core CPUs
- **Validated Results**: Benchmarked against Ghia et al. (1982) for Re = 100, 400, 1000
- **Visualization Tools**: Comprehensive plotting for velocity fields, streamlines, vorticity
- **CMake Build System**: Modern, cross-platform build configuration

## Project Structure

```
LDC-Solver/
├── CMakeLists.txt              # Root CMake configuration
├── README.md                   # This file
├── src/
│   ├── cpp/core/               # C++ implementation
│   │   ├── LDC.hpp            # Main header
│   │   ├── main.cpp           # Entry point
│   │   ├── BoundaryConditions.cpp
│   │   ├── Convection.cpp
│   │   ├── Diffusion.cpp
│   │   ├── Initialize.cpp
│   │   ├── Projection.cpp
│   │   └── Solver.cpp
|   |   |__ tests.cpp
│   └── python/
│       ├── prototype/          # Python prototype
│       └── visualization/      # Plotting tools
│           ├── visualize_ldc.py
│           └── benchmark_ldc.py
|           |__ visualize_tests.py
├── benchmarks/
│   └── data/                   # Reference data and input files
├── examples/
│   └── basic_cavity/           # Example simulations
├── scripts/
│   └── build.sh               # Build automation
├── docs/                       # Documentation (TODO)
├── build/                      # CMake build directory (generated)
├── output/                     # Simulation results (generated)
├── Profile/                    # Performance profiling (generated)
└── tests/                      # Test of convection, diffusion & pressure solver
```

## Dependencies

### C++ Build
- CMake >= 3.12
- C++11 compatible compiler (GCC, Clang, MSVC)
- OpenMP (optional, for parallel version)

### Python Visualization
- Python >= 3.7
- NumPy
- Matplotlib
- SciPy

Install Python dependencies:
```bash
pip install numpy matplotlib scipy
```

## Building

### Quick Start

```bash
# Build both serial and OpenMP versions (Release mode)
./scripts/build.sh

# macOS users: Install Homebrew GCC first if needed
# brew install gcc
# The build system auto-detects and uses g++-15, g++-14, etc.

# Build in debug mode
./scripts/build.sh --debug

# Clean build
./scripts/build.sh --clean

# Build only OpenMP version
./scripts/build.sh --no-serial
```

### Manual CMake Build

```bash
mkdir build && cd build

# Configure
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_OPENMP=ON

# Build
make -j$(nproc)

# Executables will be in build/bin/
```

### Build Options

- `CMAKE_BUILD_TYPE`: `Release` (default) or `Debug`
- `BUILD_OPENMP`: Build OpenMP parallel version (default: ON)
- `BUILD_SERIAL`: Build serial version (default: ON)
- `BUILD_TEST`: Build test version (default: ON)
- `ENABLE_PROFILING`: Enable detailed timing profiling (default: OFF)

## Usage

### Running the C++ Solver

1. **Prepare input file** (see `examples/basic_cavity/input` for template):
   ```
   Nx = 64          # Grid points in X
   Ny = 64          # Grid points in Y
   Xlen = 1.0       # Domain length in X
   Ylen = 1.0       # Domain length in Y
   Tend = 50.0      # End time
   mu1 = 0.01       # Kinematic viscosity (1/Re)
   Semi_implicit = T # Implicit diffusion
   which_solver = 1  # 1=RedBlackSOR, 2=JOR
   ...
   ```

2. **Run the simulation**:
   ```bash
   # Serial version
   cd examples/basic_cavity
   ../../build/bin/ldc_serial
   
   # OpenMP version (set threads)
   cd examples/basic_cavity
   export OMP_NUM_THREADS=4
   ../../build/bin/ldc_openmp

   # test version 
   cd test/ 
   ../build/bin/ldc_test
   ```

3. **Output files**:
   - `Visualize/output_{x,y,u,v,p}.csv` - Field data
   - `Visualize/out.m` - MATLAB script
   - `Profile/Time_Breakdown*.txt` - Performance data

### Visualization

```bash
# Generate all plots
cd examples/basic_cavity
python ../../src/python/visualization/visualize_ldc.py

# Compare with benchmark data
python ../../src/python/visualization/benchmark_ldc.py

# Generate test plots for convection, diffusion with manufactured solution: 
# U = xy(1-x)(1-y), V = 2xy(1-x)(1-y),  Poisson test with: Nabla^2(p) = -8pi^2 cos(2pi x)cos(2pi y)
cd test/
python ../src/python/visualization/visualize_tests.py
```

Generated plots:
- `velocity_contours.png` - U and V velocity fields
- `streamlines.png` - Flow streamlines
- `vorticity.png` - Vorticity contours
- `pressure.png` - Pressure field
- `centerline_profiles.png` - Velocity profiles
- `summary_plot.png` - 4-panel overview
- `centerline_velocity_ReXXX.png` - Benchmark comparison

## Benchmark Cases

### Reynolds Number = 100
```bash
# Set mu1 = 0.01 in input file (Re = 1/mu1 = 100)
Nx = 64
Ny = 64
mu1 = 0.01
Tend = 50.0
```

### Reynolds Number = 400
```bash
# Set mu1 = 0.0025 in input file (Re = 400)
Nx = 128
Ny = 128
mu1 = 0.0025
Tend = 100.0
```

### Reynolds Number = 1000
```bash
# Set mu1 = 0.001 in input file (Re = 1000)
Nx = 256
Ny = 256
mu1 = 0.001
Tend = 200.0
```

## Algorithm Overview

### Numerical Method
- **Discretization**: Finite volume method on staggered grid
- **Convection**: 2nd-order upwind interpolation
- **Diffusion**: 2nd-order central differencing (explicit or implicit)
- **Pressure-Velocity Coupling**: Fractional step (projection) method
- **Pressure Solver**: Red-Black SOR or Jacobi over-relaxation
- **Time Stepping**: Adaptive based on CFL condition

### Parallelization (OpenMP)
Parallelized operations:
- Boundary condition updates
- Convection term calculation
- Diffusion coefficient setup
- Pressure Poisson equation (Red-Black coloring)
- Projection step
- Residual calculation with reduction

## Performance

### Typical Performance (Example: 256×256 grid, Re=1000)
### 5000 time steps 
| Version | Threads | Time (s) | Speedup |
|---------|---------|----------|---------|
| Serial  | 1       |   296    | 1.0×   |
| OpenMP  | 2       |   161    | 1.83×  |
| OpenMP  | 4       |   108    | 2.73×  |
| OpenMP  | 8       |    81    | 3.62x  |

*Note: Run your own benchmarks and update this table*

### Profiling

Enable detailed profiling:
```bash
cmake .. -DENABLE_PROFILING=ON
```

Profile output saved to `Profile/Time_Breakdown_<threads>_<gridsize>.txt`

## Python Prototype

The Python prototype is **complete** and provides:
- Pure Python + NumPy implementation using Chorin's projection method
- Explicit time-stepping with uniform staggered grid
- Validated against Ghia et al. (1982) benchmark at Re=1000
- Simple, readable code perfect for learning and algorithm development
- See `src/python/prototype/README.md` for details

**Quick start:**
```bash
cd src/python/prototype
python LDC.py
```

Generates flow visualizations and comparison with benchmark data in ~30-60 seconds.

See `docs/python_vs_cpp_comparison.md` for detailed comparison with C++ version

## Troubleshooting

**Build fails with "Permission denied":**
```bash
chmod +x scripts/build.sh
make
```
*Note: If you cloned from Git, this shouldn't be necessary as Git preserves executable permissions.*

**macOS: "OpenMP not found" warning:**
```bash
# Install Homebrew GCC
brew install gcc

# Rebuild
make clean
make
```
See `docs/macOS_build_guide.md` for detailed macOS instructions.

**OpenMP not found (Linux):**
- Install: `sudo apt-get install libomp-dev` (Ubuntu/Debian)
- Or build serial only: `make serial`

**Python errors:**
```bash
pip install numpy matplotlib scipy
```

**Simulation diverges:**
- Reduce CFL: set `cfl = 0.3` in input file
- Use semi-implicit: set `Semi_implicit = T`
- Reduce timestep: set `Max_dt = 0.001`

## References

- Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. *Journal of Computational Physics*, 48(3), 387-411.

## License

MIT

## Contributing

-

## Authors

Pramod Bhuvankar

## Acknowledgments

- Benchmark data from Ghia et al. (1982)
# LDC-Solver
