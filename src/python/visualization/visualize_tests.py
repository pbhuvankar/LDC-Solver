#!/usr/bin/env python3
"""
Visualization script for LDC numerical scheme validation
Plots centerline comparisons between numerical and analytical solutions
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

def read_centerline_data(filename):
    """Read centerline data from validation output file"""
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found!")
        return None, None
    
    # Read vertical centerline data
    vertical_data = []
    horizontal_data = []
    current_section = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                if 'VERTICAL_CENTERLINE' in line:
                    current_section = 'vertical'
                elif 'HORIZONTAL_CENTERLINE' in line:
                    current_section = 'horizontal'
                continue
            
            if line and current_section:
                values = [float(x) for x in line.split()]
                if current_section == 'vertical':
                    vertical_data.append(values)
                elif current_section == 'horizontal':
                    horizontal_data.append(values)
    
    vertical_data = np.array(vertical_data) if vertical_data else None
    horizontal_data = np.array(horizontal_data) if horizontal_data else None
    
    return vertical_data, horizontal_data

def plot_convection_validation():
    """Plot convection scheme validation"""
    print("Plotting convection validation...")
    vertical, horizontal = read_centerline_data('validation_convection.dat')
    
    if vertical is None or horizontal is None:
        print("  Error: Could not read convection data")
        return
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    # Vertical centerline - U momentum
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(vertical[:, 1], vertical[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(vertical)//15))
    ax1.plot(vertical[:, 1], vertical[:, 3], 'r--', linewidth=2, label='Analytical')
    ax1.set_xlabel('y')
    ax1.set_ylabel('du/dt · Δt (U-momentum)')
    ax1.set_title('Convection: U-momentum (Vertical Centerline)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Horizontal centerline - U momentum
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(horizontal[:, 0], horizontal[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(horizontal)//15))
    ax2.plot(horizontal[:, 0], horizontal[:, 3], 'r--', linewidth=2, label='Analytical')
    ax2.set_xlabel('x')
    ax2.set_ylabel('du/dt · Δt (U-momentum)')
    ax2.set_title('Convection: U-momentum (Horizontal Centerline)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Vertical centerline - V momentum
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(vertical[:, 1], vertical[:, 4], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(vertical)//15))
    ax3.plot(vertical[:, 1], vertical[:, 5], 'r--', linewidth=2, label='Analytical')
    ax3.set_xlabel('y')
    ax3.set_ylabel('dv/dt · Δt (V-momentum)')
    ax3.set_title('Convection: V-momentum (Vertical Centerline)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Horizontal centerline - V momentum
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(horizontal[:, 0], horizontal[:, 4], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(horizontal)//15))
    ax4.plot(horizontal[:, 0], horizontal[:, 5], 'r--', linewidth=2, label='Analytical')
    ax4.set_xlabel('x')
    ax4.set_ylabel('dv/dt · Δt (V-momentum)')
    ax4.set_title('Convection: V-momentum (Horizontal Centerline)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Error plots - U momentum
    ax5 = fig.add_subplot(gs[2, 0])
    error_u_vert = np.abs(vertical[:, 2] - vertical[:, 3])
    error_u_horiz = np.abs(horizontal[:, 2] - horizontal[:, 3])
    ax5.semilogy(vertical[:, 1], error_u_vert, 'b-o', linewidth=2, markersize=4, 
                 label='Vertical centerline', markevery=max(1, len(vertical)//15))
    ax5.semilogy(horizontal[:, 0], error_u_horiz, 'g-s', linewidth=2, markersize=4, 
                 label='Horizontal centerline', markevery=max(1, len(horizontal)//15))
    ax5.set_xlabel('Position')
    ax5.set_ylabel('Absolute Error')
    ax5.set_title('Convection: U-momentum Error')
    ax5.legend()
    ax5.grid(True, alpha=0.3, which='both')
    
    # Error plots - V momentum
    ax6 = fig.add_subplot(gs[2, 1])
    error_v_vert = np.abs(vertical[:, 4] - vertical[:, 5])
    error_v_horiz = np.abs(horizontal[:, 4] - horizontal[:, 5])
    ax6.semilogy(vertical[:, 1], error_v_vert, 'b-o', linewidth=2, markersize=4, 
                 label='Vertical centerline', markevery=max(1, len(vertical)//15))
    ax6.semilogy(horizontal[:, 0], error_v_horiz, 'g-s', linewidth=2, markersize=4, 
                 label='Horizontal centerline', markevery=max(1, len(horizontal)//15))
    ax6.set_xlabel('Position')
    ax6.set_ylabel('Absolute Error')
    ax6.set_title('Convection: V-momentum Error')
    ax6.legend()
    ax6.grid(True, alpha=0.3, which='both')
    
    plt.suptitle('Convection Scheme Validation: U = x(1-x)y(1-y) & V = x(1-x)y(1-y)', fontsize=14, fontweight='bold', y=0.995)
    plt.savefig('test_convection.png', dpi=300, bbox_inches='tight')
    print("  Saved: test_convection.png")
    plt.close()

def plot_diffusion_validation():
    """Plot diffusion scheme validation"""
    print("Plotting diffusion validation...")
    vertical, horizontal = read_centerline_data('validation_diffusion.dat')
    
    if vertical is None or horizontal is None:
        print("  Error: Could not read diffusion data")
        return
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    # Vertical centerline - U momentum
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(vertical[:, 1], vertical[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(vertical)//15))
    ax1.plot(vertical[:, 1], vertical[:, 3], 'r--', linewidth=2, label='Analytical')
    ax1.set_xlabel('y')
    ax1.set_ylabel('∂²u/∂x² + ∂²u/∂y² (scaled)')
    ax1.set_title('Diffusion: U-momentum (Vertical Centerline)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Horizontal centerline - U momentum
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(horizontal[:, 0], horizontal[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(horizontal)//15))
    ax2.plot(horizontal[:, 0], horizontal[:, 3], 'r--', linewidth=2, label='Analytical')
    ax2.set_xlabel('x')
    ax2.set_ylabel('∂²u/∂x² + ∂²u/∂y² (scaled)')
    ax2.set_title('Diffusion: U-momentum (Horizontal Centerline)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Vertical centerline - V momentum
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(vertical[:, 1], vertical[:, 4], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(vertical)//15))
    ax3.plot(vertical[:, 1], vertical[:, 5], 'r--', linewidth=2, label='Analytical')
    ax3.set_xlabel('y')
    ax3.set_ylabel('∂²v/∂x² + ∂²v/∂y² (scaled)')
    ax3.set_title('Diffusion: V-momentum (Vertical Centerline)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Horizontal centerline - V momentum
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(horizontal[:, 0], horizontal[:, 4], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(horizontal)//15))
    ax4.plot(horizontal[:, 0], horizontal[:, 5], 'r--', linewidth=2, label='Analytical')
    ax4.set_xlabel('x')
    ax4.set_ylabel('∂²v/∂x² + ∂²v/∂y² (scaled)')
    ax4.set_title('Diffusion: V-momentum (Horizontal Centerline)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Error plots - U momentum
    ax5 = fig.add_subplot(gs[2, 0])
    error_u_vert = np.abs(vertical[:, 2] - vertical[:, 3])
    error_u_horiz = np.abs(horizontal[:, 2] - horizontal[:, 3])
    ax5.semilogy(vertical[:, 1], error_u_vert, 'b-o', linewidth=2, markersize=4, 
                 label='Vertical centerline', markevery=max(1, len(vertical)//15))
    ax5.semilogy(horizontal[:, 0], error_u_horiz, 'g-s', linewidth=2, markersize=4, 
                 label='Horizontal centerline', markevery=max(1, len(horizontal)//15))
    ax5.set_xlabel('Position')
    ax5.set_ylabel('Absolute Error')
    ax5.set_title('Diffusion: U-momentum Error')
    ax5.legend()
    ax5.grid(True, alpha=0.3, which='both')
    
    # Error plots - V momentum
    ax6 = fig.add_subplot(gs[2, 1])
    error_v_vert = np.abs(vertical[:, 4] - vertical[:, 5])
    error_v_horiz = np.abs(horizontal[:, 4] - horizontal[:, 5])
    ax6.semilogy(vertical[:, 1], error_v_vert, 'b-o', linewidth=2, markersize=4, 
                 label='Vertical centerline', markevery=max(1, len(vertical)//15))
    ax6.semilogy(horizontal[:, 0], error_v_horiz, 'g-s', linewidth=2, markersize=4, 
                 label='Horizontal centerline', markevery=max(1, len(horizontal)//15))
    ax6.set_xlabel('Position')
    ax6.set_ylabel('Absolute Error')
    ax6.set_title('Diffusion: V-momentum Error')
    ax6.legend()
    ax6.grid(True, alpha=0.3, which='both')
    
    plt.suptitle('Diffusion Scheme Validation: U = x(1-x)y(1-y) & V = x(1-x)y(1-y)', fontsize=14, fontweight='bold', y=0.995)
    plt.savefig('test_diffusion.png', dpi=300, bbox_inches='tight')
    print("  Saved: test_diffusion.png")
    plt.close()

def plot_poisson_validation():
    """Plot Poisson solver validation"""
    print("Plotting Poisson validation...")
    vertical, horizontal = read_centerline_data('validation_poisson.dat')
    
    if vertical is None or horizontal is None:
        print("  Error: Could not read Poisson data")
        return
    
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # Vertical centerline
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(vertical[:, 1], vertical[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(vertical)//15))
    ax1.plot(vertical[:, 1], vertical[:, 3], 'r--', linewidth=2, label='Analytical')
    ax1.set_xlabel('y')
    ax1.set_ylabel('Pressure')
    ax1.set_title('Poisson Solution (Vertical Centerline)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Horizontal centerline
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(horizontal[:, 0], horizontal[:, 2], 'b-o', linewidth=2, markersize=4, 
             label='Numerical', markevery=max(1, len(horizontal)//15))
    ax2.plot(horizontal[:, 0], horizontal[:, 3], 'r--', linewidth=2, label='Analytical')
    ax2.set_xlabel('x')
    ax2.set_ylabel('Pressure')
    ax2.set_title('Poisson Solution (Horizontal Centerline)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Error plot - Vertical centerline
    ax3 = fig.add_subplot(gs[1, 0])
    error_vert = np.abs(vertical[:, 2] - vertical[:, 3])
    ax3.semilogy(vertical[:, 1], error_vert, 'b-o', linewidth=2, markersize=5, 
                 label='Absolute Error', markevery=max(1, len(vertical)//15))
    ax3.set_xlabel('y')
    ax3.set_ylabel('Absolute Error')
    ax3.set_title('Poisson Error (Vertical Centerline)')
    ax3.legend()
    ax3.grid(True, alpha=0.3, which='both')
    
    # Error plot - Horizontal centerline
    ax4 = fig.add_subplot(gs[1, 1])
    error_horiz = np.abs(horizontal[:, 2] - horizontal[:, 3])
    ax4.semilogy(horizontal[:, 0], error_horiz, 'b-o', linewidth=2, markersize=5, 
                 label='Absolute Error', markevery=max(1, len(horizontal)//15))
    ax4.set_xlabel('x')
    ax4.set_ylabel('Absolute Error')
    ax4.set_title('Poisson Error (Horizontal Centerline)')
    ax4.legend()
    ax4.grid(True, alpha=0.3, which='both')
    
    plt.suptitle('Poisson Solver Validation: Nabla^2(p) = -8pi^2 cos(2pi x)cos(2pi y)', fontsize=14, fontweight='bold', y=0.995)
    plt.savefig('test_poisson.png', dpi=300, bbox_inches='tight')
    print("  Saved: test_poisson.png")
    plt.close()

def create_summary_plot():
    """Create a summary comparison plot"""
    print("Creating summary plot...")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Convection
    vertical_conv, horizontal_conv = read_centerline_data('validation_convection.dat')
    if vertical_conv is not None:
        ax = axes[0]
        ax.plot(vertical_conv[:, 1], vertical_conv[:, 2], 'b-', linewidth=2, 
                label='U-mom (numerical)', alpha=0.7)
        ax.plot(vertical_conv[:, 1], vertical_conv[:, 3], 'b--', linewidth=2, 
                label='U-mom (analytical)', alpha=0.7)
        ax.plot(vertical_conv[:, 1], vertical_conv[:, 4], 'r-', linewidth=2, 
                label='V-mom (numerical)', alpha=0.7)
        ax.plot(vertical_conv[:, 1], vertical_conv[:, 5], 'r--', linewidth=2, 
                label='V-mom (analytical)', alpha=0.7)
        ax.set_xlabel('y (Vertical Centerline)')
        ax.set_ylabel('Convection Terms')
        ax.set_title('Convection')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Diffusion
    vertical_diff, horizontal_diff = read_centerline_data('validation_diffusion.dat')
    if vertical_diff is not None:
        ax = axes[1]
        ax.plot(vertical_diff[:, 1], vertical_diff[:, 2], 'b-', linewidth=2, 
                label='U-mom (numerical)', alpha=0.7)
        ax.plot(vertical_diff[:, 1], vertical_diff[:, 3], 'b--', linewidth=2, 
                label='U-mom (analytical)', alpha=0.7)
        ax.plot(vertical_diff[:, 1], vertical_diff[:, 4], 'r-', linewidth=2, 
                label='V-mom (numerical)', alpha=0.7)
        ax.plot(vertical_diff[:, 1], vertical_diff[:, 5], 'r--', linewidth=2, 
                label='V-mom (analytical)', alpha=0.7)
        ax.set_xlabel('y (Vertical Centerline)')
        ax.set_ylabel('Diffusion Terms')
        ax.set_title('Diffusion')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Poisson
    vertical_pois, horizontal_pois = read_centerline_data('validation_poisson.dat')
    if vertical_pois is not None:
        ax = axes[2]
        ax.plot(vertical_pois[:, 1], vertical_pois[:, 2], 'b-', linewidth=2, 
                label='Numerical', alpha=0.7)
        ax.plot(vertical_pois[:, 1], vertical_pois[:, 3], 'r--', linewidth=2, 
                label='Analytical', alpha=0.7)
        ax.set_xlabel('y (Vertical Centerline)')
        ax.set_ylabel('Pressure')
        ax.set_title('Poisson Solver')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Numerical Scheme Validation Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('test_summary.png', dpi=300, bbox_inches='tight')
    print("  Saved: test_summary.png")
    plt.close()

def main():
    """Main function to generate all validation plots"""
    print("\n" + "="*60)
    print("  LDC Numerical Scheme Validation Visualization")
    print("="*60 + "\n")
    
    # Check if data files exist
    files = ['validation_convection.dat', 'validation_diffusion.dat', 'validation_poisson.dat']
    missing = [f for f in files if not os.path.exists(f)]
    
    if missing:
        print("Error: The following data files are missing:")
        for f in missing:
            print(f"  - {f}")
        print("\nPlease run the validation test program first.")
        sys.exit(1)
    
    # Generate all plots
    plot_convection_validation()
    plot_diffusion_validation()
    plot_poisson_validation()
    create_summary_plot()
    
    print("\n" + "="*60)
    print("  Visualization Complete!")
    print("="*60)
    print("\nGenerated files:")
    print("  - test_convection.png")
    print("  - test_diffusion.png")
    print("  - test_poisson.png")
    print("  - test_summary.png")
    print()

if __name__ == "__main__":
    main()
