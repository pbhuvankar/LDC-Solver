# Quick Start Guide

## Build and Run in 5 Minutes

### 1. Build the solver

```bash
cd LDC-Solver
make
# or equivalently: ./scripts/build.sh
```

This creates:
- `build/bin/ldc_serial` - Serial version
- `build/bin/ldc_openmp` - Parallel OpenMP version

### 2. Run a simulation

```bash
cd examples/basic_cavity

# Run serial version
../../build/bin/ldc_serial

# Or run OpenMP version with 4 threads
export OMP_NUM_THREADS=4
../../build/bin/ldc_openmp
```

### 3. Visualize results

```bash
# Generate all plots
python ../../src/python/visualization/visualize_ldc.py

# Compare with benchmark (enter Re when prompted)
python ../../src/python/visualization/benchmark_ldc.py
```

Your plots will be in `Visualize/` directory!

## Understanding the Input File

The `input` file controls simulation parameters:

```
Nx = 64              # Grid resolution in X (try: 32, 64, 128, 256)
Ny = 64              # Grid resolution in Y
Xlen = 1.0           # Domain size X
Ylen = 1.0           # Domain size Y
mesh_ref = 0.5       # Grid stretching (0=uniform, >0=refined near walls)

Tend = 50.0          # Simulation end time
Max_dt = 0.01        # Maximum timestep
cfl = 0.5            # CFL number (stability: typically 0.3-0.8)

mu1 = 0.01           # Kinematic viscosity = 1/Re
                     # Re=100  -> mu1=0.01
                     # Re=400  -> mu1=0.0025  
                     # Re=1000 -> mu1=0.001

Beta = 1.7           # SOR relaxation parameter (1.0-2.0, typically 1.7-1.9)
maxit = 10000        # Max iterations for pressure solver

which_solver = 1     # 1 = Red-Black SOR (recommended)
                     # 2 = Jacobi

Semi_implicit = T    # T = Implicit diffusion (stable, larger dt)
                     # F = Explicit diffusion (requires small dt)

Breakdown_Time = T   # T = Print timing breakdown
```

## Typical Workflow

### 1. Low Reynolds Number (Re=100) - Quick Test
```bash
# Edit input file:
Nx = 64
Ny = 64
mu1 = 0.01
Tend = 50.0

# Run (takes ~30 seconds)
../../build/bin/ldc_openmp

# Visualize and benchmark
python ../../src/python/visualization/visualize_ldc.py
python ../../src/python/visualization/benchmark_ldc.py  # enter: 100
```

### 2. Medium Reynolds Number (Re=400)
```bash
# Edit input file:
Nx = 128
Ny = 128
mu1 = 0.0025
Tend = 100.0

# Run (takes ~5 minutes)
export OMP_NUM_THREADS=8
../../build/bin/ldc_openmp

# Visualize
python ../../src/python/visualization/visualize_ldc.py
python ../../src/python/visualization/benchmark_ldc.py  # enter: 400
```

### 3. High Reynolds Number (Re=1000) - Production Run
```bash
# Edit input file:
Nx = 256
Ny = 256
mu1 = 0.001
Tend = 200.0

# Run (takes ~30-60 minutes)
export OMP_NUM_THREADS=16
../../build/bin/ldc_openmp

# Visualize
python ../../src/python/visualization/visualize_ldc.py
python ../../src/python/visualization/benchmark_ldc.py  # enter: 1000
```

## Performance Tips

1. **Use OpenMP version**: 3-8x faster than serial
   ```bash
   export OMP_NUM_THREADS=$(nproc)  # Use all cores
   ```

2. **Use Semi-implicit**: More stable, allows larger timesteps
   ```
   Semi_implicit = T
   ```

3. **Tune Beta**: Optimal SOR parameter (problem-dependent)
   ```
   Beta = 1.7   # Good starting point
   Beta = 1.8   # Try if convergence is slow
   ```

4. **Grid refinement near walls**: Better boundary layer resolution
   ```
   mesh_ref = 0.5   # Moderate refinement
   mesh_ref = 1.0   # Strong refinement
   ```

## Output Files

After running, you'll find:

```
Visualize/
├── output_x.csv          # X coordinates
├── output_y.csv          # Y coordinates
├── output_u.csv          # U velocity field
├── output_v.csv          # V velocity field
├── output_p.csv          # Pressure field
├── out.m                 # MATLAB script
└── *.png                 # Generated plots

Profile/
└── Time_Breakdown_*.txt  # Performance profiling
```

## Troubleshooting

**Problem: Solver diverges**
- Reduce `cfl` (try 0.3)
- Use `Semi_implicit = T`
- Reduce `Max_dt`

**Problem: Convergence too slow**
- Increase `Beta` (try 1.8 or 1.9)
- Use `which_solver = 1` (Red-Black SOR)
- Increase `maxit` if hitting iteration limit

**Problem: Results don't match benchmark**
- Increase grid resolution (Nx, Ny)
- Run longer (increase Tend)
- Check Re is correct: Re = 1/mu1

**Problem: OpenMP not working**
- Check CMake found OpenMP: `cmake .. | grep OpenMP`
- Rebuild: `make clean && make`
- Set threads: `export OMP_NUM_THREADS=4`

## Next Steps

- Run benchmark cases (Re=100, 400, 1000)
- Experiment with grid sizes and parameters
- Profile performance: check `Profile/Time_Breakdown*.txt`
- Develop Python prototype in `src/python/prototype/`
- Document performance model

Happy computing! 🚀
