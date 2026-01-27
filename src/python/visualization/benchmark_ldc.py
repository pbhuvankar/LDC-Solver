import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

def get_benchmark_data(Re):
    """
    Get benchmark data for vertical centerline u-velocity
    Data from Ghia, Ghia & Shin (1982)
    
    Parameters:
    -----------
    Re : int
        Reynolds number (100, 400, or 1000)
    
    Returns:
    --------
    y_vertical : np.array
        Y coordinates for vertical centerline
    u_vertical : np.array
        U velocity at vertical centerline
    """
    # Y coordinates for vertical centerline (same for all Re)
    y_vertical = np.array([
        0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 
        0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 
        0.9688, 0.9766, 1.0000
    ])
    
    if Re == 100:
        u_vertical = np.array([
            0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
            -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151,
            0.68717, 0.73722, 0.78871, 0.84123, 1.00000
        ])        
        
    elif Re == 400:
        u_vertical = np.array([
            0.00000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299,
            -0.32726, -0.17119, -0.11477, 0.02135, 0.16256, 0.29093,
            0.55892, 0.61756, 0.68439, 0.75837, 1.00000
        ])        
        
    elif Re == 1000:
        u_vertical = np.array([
            0.00000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289,
            -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304,
            0.46604, 0.51117, 0.57492, 0.65928, 1.00000
        ])
    else:
        raise ValueError(f"Reynolds number {Re} not supported. Use 100, 400, or 1000.")
    
    return y_vertical, u_vertical


def read_csv_data(visualize_dir='Visualize'):
    """
    Read coordinate and velocity data from CSV files
    
    Parameters:
    -----------
    visualize_dir : str
        Directory containing the CSV files
    
    Returns:
    --------
    x : np.array
        X coordinates
    y : np.array
        Y coordinates
    u : np.array
        U velocity field (2D array)
    """
    # Read coordinate files
    x_file = os.path.join(visualize_dir, 'output_x.csv')
    y_file = os.path.join(visualize_dir, 'output_y.csv')
    u_file = os.path.join(visualize_dir, 'output_u.csv')
    
    # Check if files exist
    for file in [x_file, y_file, u_file]:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File not found: {file}")
    
    # Read X coordinates (1D array)
    x = np.loadtxt(x_file, delimiter=',')
    
    # Read Y coordinates (1D array)
    y = np.loadtxt(y_file, delimiter=',')
    
    # Read U velocity field (2D array: rows=y, cols=x)
    u = np.loadtxt(u_file, delimiter=',')
    
    print(f"Data loaded successfully:")
    print(f"  X grid: {len(x)} points, range [{x.min():.4f}, {x.max():.4f}]")
    print(f"  Y grid: {len(y)} points, range [{y.min():.4f}, {y.max():.4f}]")
    print(f"  U velocity field shape: {u.shape}")
    
    return x, y, u


def extract_centerline_velocity(x, y, u, x_center=0.5):
    """
    Extract vertical u-velocity profile at specified x location
    Uses linear interpolation if x_center doesn't align with grid
    
    Parameters:
    -----------
    x : np.array
        X coordinates (1D)
    y : np.array
        Y coordinates (1D)
    u : np.array
        U velocity field (2D: rows=y, cols=x)
    x_center : float
        X location for vertical profile (default: 0.5)
    
    Returns:
    --------
    y_profile : np.array
        Y coordinates for profile
    u_profile : np.array
        U velocity along vertical centerline
    """
    # Find the closest x indices to x_center
    idx = np.argmin(np.abs(x - x_center))
    
    # Check if we need interpolation
    if np.abs(x[idx] - x_center) < 1e-10:
        # Exact match - no interpolation needed
        print(f"Exact match at x = {x[idx]:.6f}")
        u_profile = u[:, idx]
    else:
        # Need interpolation
        print(f"Interpolating between x = {x[idx-1]:.6f} and x = {x[idx]:.6f}")
        
        # Get velocities at neighboring x locations
        u_left = u[:, idx-1] if idx > 0 else u[:, 0]
        u_right = u[:, idx+1] if idx < len(x)-1 else u[:, -1]
        x_left = x[idx-1] if idx > 0 else x[0]
        x_right = x[idx+1] if idx < len(x)-1 else x[-1]
        
        # Linear interpolation
        weight = (x_center - x_left) / (x_right - x_left)
        u_profile = u_left + weight * (u_right - u_left)
    
    return y, u_profile


def calculate_error(y_sim, u_sim, y_bench, u_bench):
    """
    Calculate error metrics between simulation and benchmark
    
    Parameters:
    -----------
    y_sim : np.array
        Simulation y coordinates
    u_sim : np.array
        Simulation u velocities
    y_bench : np.array
        Benchmark y coordinates
    u_bench : np.array
        Benchmark u velocities
    
    Returns:
    --------
    l2_error : float
        L2 norm error
    max_error : float
        Maximum absolute error
    """
    # Interpolate simulation data at benchmark locations
    interp_func = interp1d(y_sim, u_sim, kind='cubic', fill_value='extrapolate')
    u_sim_interp = interp_func(y_bench)
    
    # Calculate errors
    l2_error = np.sqrt(np.mean((u_sim_interp - u_bench)**2))
    max_error = np.max(np.abs(u_sim_interp - u_bench))
    
    return l2_error, max_error


def plot_centerline_comparison(y_sim, u_sim, y_bench, u_bench, Re, 
                               save_fig=True, output_dir='Visualize'):
    """
    Plot comparison between simulation and benchmark data
    
    Parameters:
    -----------
    y_sim : np.array
        Simulation y coordinates
    u_sim : np.array
        Simulation u velocities
    y_bench : np.array
        Benchmark y coordinates
    u_bench : np.array
        Benchmark u velocities
    Re : int
        Reynolds number
    save_fig : bool
        Whether to save the figure
    output_dir : str
        Directory to save the figure
    """
    # Calculate error metrics
    l2_error, max_error = calculate_error(y_sim, u_sim, y_bench, u_bench)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Plot simulation results
    ax.plot(y_sim, u_sim, 'b-', linewidth=2, label='Present Simulation')
    
    # Plot benchmark data
    ax.plot(y_bench, u_bench, 'ro', markersize=8, markerfacecolor='none', 
            markeredgewidth=2, label='Ghia et al. (1982)')
    
    # Formatting
    ax.set_ylabel('U Velocity', fontsize=14, fontweight='bold')
    ax.set_xlabel('Y Position', fontsize=14, fontweight='bold')
    ax.set_title(f'Vertical Centerline U-Velocity Profile (Re = {Re})', 
                 fontsize=16, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(fontsize=12, loc='best')
    
    # Add error text box
    textstr = f'L2 Error: {l2_error:.6f}\nMax Error: {max_error:.6f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
    
    # Set axis limits
    ax.set_xlim([0, 1])
    ax.set_ylim([min(u_sim.min(), u_bench.min()) - 0.1, 1.1])
    
    plt.tight_layout()
    
    # Save figure
    if save_fig:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        figname = os.path.join(output_dir, f'centerline_velocity_Re{Re}.png')
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved: {figname}")
    
    plt.show()
    
    # Print error metrics
    print(f"\n{'='*60}")
    print(f"Error Analysis (Re = {Re}):")
    print(f"{'='*60}")
    print(f"L2 Norm Error:     {l2_error:.6e}")
    print(f"Max Absolute Error: {max_error:.6e}")
    print(f"{'='*60}\n")


def main():
    """
    Main function to process and plot centerline velocity
    """
    # Set Reynolds number (change this as needed)
    Re = input("Enter Re ")#1000  # Options: 100, 400, 1000
    Re = float(Re)
    
    print(f"\n{'='*60}")
    print(f"Lid-Driven Cavity: Centerline Velocity Analysis (Re = {Re})")
    print(f"{'='*60}\n")
    
    # Read simulation data
    try:
        x, y, u = read_csv_data(visualize_dir='Visualize')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("\nMake sure you have run the simulation and generated CSV files.")
        return
    
    # Extract centerline velocity
    print(f"\nExtracting vertical centerline at x = 0.5...")
    y_sim, u_sim = extract_centerline_velocity(x, y, u, x_center=0.5)
    
    # Get benchmark data
    print(f"Loading benchmark data for Re = {Re}...")
    y_bench, u_bench = get_benchmark_data(Re)
    
    # Plot comparison
    print(f"Creating comparison plot...")
    plot_centerline_comparison(y_sim, u_sim, y_bench, u_bench, Re, 
                               save_fig=True, output_dir='Visualize')
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
