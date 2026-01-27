import numpy as np
import matplotlib.pyplot as plt

class LDCVisualizer:
    """Visualization class for lid-driven cavity results"""
    
    def __init__(self, prefix='Visualize/output'):
        """Load data from CSV files"""
        print("Loading data...")
        self.x = np.loadtxt(f'{prefix}_x.csv', delimiter=',')
        self.y = np.loadtxt(f'{prefix}_y.csv', delimiter=',')
        self.u = np.loadtxt(f'{prefix}_u.csv', delimiter=',')
        self.v = np.loadtxt(f'{prefix}_v.csv', delimiter=',')
        self.p = np.loadtxt(f'{prefix}_p.csv', delimiter=',')
        
        # Create meshgrid
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        # Calculate derived quantities
        self.vel_mag = np.sqrt(self.u**2 + self.v**2)
        self.calculate_vorticity()
        
        print(f"Grid size: {self.u.shape[0]} x {self.u.shape[1]}")
        print(f"Data loaded successfully\n")
    
    def calculate_vorticity(self):
        """Calculate vorticity (omega_z = dv/dx - du/dy)"""
        ny, nx = self.u.shape
        self.vorticity = np.zeros_like(self.u)
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                dvdx = (self.v[j, i+1] - self.v[j, i-1]) / (self.x[i+1] - self.x[i-1])
                dudy = (self.u[j+1, i] - self.u[j-1, i]) / (self.y[j+1] - self.y[j-1])
                self.vorticity[j, i] = dvdx - dudy
    
    def plot_velocity_contours(self):
        """Plot U and V velocity contours"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # U velocity
        cf1 = axes[0].contourf(self.X, self.Y, self.u, levels=20, cmap='RdBu_r')
        axes[0].set_xlabel('X')
        axes[0].set_ylabel('Y')
        axes[0].set_title('U-Velocity', fontweight='bold')
        axes[0].set_aspect('equal')
        plt.colorbar(cf1, ax=axes[0])
        
        # V velocity
        cf2 = axes[1].contourf(self.X, self.Y, self.v, levels=20, cmap='RdBu_r')
        axes[1].set_xlabel('X')
        axes[1].set_ylabel('Y')
        axes[1].set_title('V-Velocity', fontweight='bold')
        axes[1].set_aspect('equal')
        plt.colorbar(cf2, ax=axes[1])
        
        plt.tight_layout()
        plt.savefig('Visualize/velocity_contours.png', dpi=300, bbox_inches='tight')
        print("Saved: velocity_contours.png")
        plt.close()
    
    def plot_streamlines(self):
        """Plot streamlines - interpolate to uniform grid if needed"""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Check if grid is uniform
        dx_uniform = np.allclose(np.diff(self.x), np.diff(self.x)[0])
        dy_uniform = np.allclose(np.diff(self.y), np.diff(self.y)[0])
        
        if dx_uniform and dy_uniform:
            # Use original grid
            ax.contourf(self.X, self.Y, self.vel_mag, levels=20, cmap='YlOrRd', alpha=0.6)
            ax.streamplot(self.X, self.Y, self.u, self.v, color='black', linewidth=1, density=2)
        else:
            # Interpolate to uniform grid
            from scipy.interpolate import griddata
            nx, ny = len(self.x), len(self.y)
            x_uniform = np.linspace(self.x.min(), self.x.max(), nx)
            y_uniform = np.linspace(self.y.min(), self.y.max(), ny)
            X_uniform, Y_uniform = np.meshgrid(x_uniform, y_uniform)
            
            points = np.column_stack([self.X.ravel(), self.Y.ravel()])
            u_uniform = griddata(points, self.u.ravel(), (X_uniform, Y_uniform), method='linear')
            v_uniform = griddata(points, self.v.ravel(), (X_uniform, Y_uniform), method='linear')
            vel_mag_uniform = griddata(points, self.vel_mag.ravel(), (X_uniform, Y_uniform), method='linear')
            
            ax.contourf(X_uniform, Y_uniform, vel_mag_uniform, levels=20, cmap='YlOrRd', alpha=0.6)
            ax.streamplot(X_uniform, Y_uniform, u_uniform, v_uniform, color='black', linewidth=1, density=2)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Streamlines', fontweight='bold')
        ax.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('Visualize/streamlines.png', dpi=300, bbox_inches='tight')
        print("Saved: streamlines.png")
        plt.close()
    
    def plot_vorticity(self):
        """Plot vorticity contours"""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        vort_max = max(abs(self.vorticity.min()), abs(self.vorticity.max()))
        levels = np.linspace(-vort_max, vort_max, 20)
        
        cf = ax.contourf(self.X, self.Y, self.vorticity, levels=levels, cmap='RdBu_r')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Vorticity', fontweight='bold')
        ax.set_aspect('equal')
        plt.colorbar(cf, ax=ax)
        
        plt.tight_layout()
        plt.savefig('Visualize/vorticity.png', dpi=300, bbox_inches='tight')
        print("Saved: vorticity.png")
        plt.close()
    
    def plot_pressure(self):
        """Plot pressure contours"""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        cf = ax.contourf(self.X, self.Y, self.p, levels=20, cmap='viridis')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Pressure', fontweight='bold')
        ax.set_aspect('equal')
        plt.colorbar(cf, ax=ax)
        
        plt.tight_layout()
        plt.savefig('Visualize/pressure.png', dpi=300, bbox_inches='tight')
        print("Saved: pressure.png")
        plt.close()
    
    def plot_centerline_profiles(self):
        """Plot velocity profiles along centerlines"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        nx, ny = len(self.x), len(self.y)
        i_center = nx // 2
        j_center = ny // 2
        
        # U velocity along vertical centerline
        u_vertical = self.u[:, i_center]
        axes[0].plot(u_vertical, self.y, 'b-', linewidth=2)
        axes[0].set_xlabel('U velocity')
        axes[0].set_ylabel('Y')
        axes[0].set_title('U-velocity along vertical centerline', fontweight='bold')
        axes[0].grid(True, alpha=0.3)
        
        # V velocity along horizontal centerline
        v_horizontal = self.v[j_center, :]
        axes[1].plot(self.x, v_horizontal, 'b-', linewidth=2)
        axes[1].set_xlabel('X')
        axes[1].set_ylabel('V velocity')
        axes[1].set_title('V-velocity along horizontal centerline', fontweight='bold')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('Visualize/centerline_profiles.png', dpi=300, bbox_inches='tight')
        print("Saved: centerline_profiles.png")
        plt.close()
    
    def plot_vector_field(self, skip=3):
        """Plot velocity vector field"""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        ax.contourf(self.X, self.Y, self.vel_mag, levels=20, cmap='YlOrRd', alpha=0.5)
        ax.quiver(self.X[::skip, ::skip], self.Y[::skip, ::skip], 
                 self.u[::skip, ::skip], self.v[::skip, ::skip])
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Velocity Vector Field', fontweight='bold')
        ax.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('Visualize/vector_field.png', dpi=300, bbox_inches='tight')
        print("Saved: vector_field.png")
        plt.close()
    
    def create_summary_plot(self):
        """Create a 4-panel summary figure"""
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # Check if grid is uniform for streamlines
        dx_uniform = np.allclose(np.diff(self.x), np.diff(self.x)[0])
        dy_uniform = np.allclose(np.diff(self.y), np.diff(self.y)[0])
        
        # Streamlines
        if dx_uniform and dy_uniform:
            axes[0,0].streamplot(self.X, self.Y, self.u, self.v, 
                                color='black', linewidth=1, density=2)
        else:
            # Use interpolation for non-uniform grids
            from scipy.interpolate import griddata
            nx, ny = len(self.x), len(self.y)
            x_uniform = np.linspace(self.x.min(), self.x.max(), nx)
            y_uniform = np.linspace(self.y.min(), self.y.max(), ny)
            X_uniform, Y_uniform = np.meshgrid(x_uniform, y_uniform)
            
            points = np.column_stack([self.X.ravel(), self.Y.ravel()])
            u_uniform = griddata(points, self.u.ravel(), (X_uniform, Y_uniform), method='linear')
            v_uniform = griddata(points, self.v.ravel(), (X_uniform, Y_uniform), method='linear')
            
            axes[0,0].streamplot(X_uniform, Y_uniform, u_uniform, v_uniform, 
                                color='black', linewidth=1, density=2)
        
        axes[0,0].set_xlabel('X')
        axes[0,0].set_ylabel('Y')
        axes[0,0].set_title('Streamlines', fontweight='bold')
        axes[0,0].set_aspect('equal')
        
        # Vorticity
        vort_max = max(abs(self.vorticity.min()), abs(self.vorticity.max()))
        levels_v = np.linspace(-vort_max, vort_max, 15)
        cf2 = axes[0,1].contourf(self.X, self.Y, self.vorticity, levels=levels_v, cmap='RdBu_r')
        axes[0,1].set_xlabel('X')
        axes[0,1].set_ylabel('Y')
        axes[0,1].set_title('Vorticity', fontweight='bold')
        axes[0,1].set_aspect('equal')
        plt.colorbar(cf2, ax=axes[0,1])
        
        # Velocity magnitude
        cf3 = axes[1,0].contourf(self.X, self.Y, self.vel_mag, levels=15, cmap='YlOrRd')
        axes[1,0].set_xlabel('X')
        axes[1,0].set_ylabel('Y')
        axes[1,0].set_title('Velocity Magnitude', fontweight='bold')
        axes[1,0].set_aspect('equal')
        plt.colorbar(cf3, ax=axes[1,0])
        
        # Pressure
        cf4 = axes[1,1].contourf(self.X, self.Y, self.p, levels=15, cmap='viridis')
        axes[1,1].set_xlabel('X')
        axes[1,1].set_ylabel('Y')
        axes[1,1].set_title('Pressure', fontweight='bold')
        axes[1,1].set_aspect('equal')
        plt.colorbar(cf4, ax=axes[1,1])
        
        plt.tight_layout()
        plt.savefig('Visualize/summary_plot.png', dpi=300, bbox_inches='tight')
        print("Saved: summary_plot.png")
        plt.close()
    
    def generate_all_plots(self):
        """Generate all visualization plots"""
        print("Generating plots...\n")
        self.plot_velocity_contours()
        self.plot_streamlines()
        self.plot_vorticity()
        self.plot_pressure()
        self.plot_centerline_profiles()
        self.plot_vector_field()
        self.create_summary_plot()
        print("\nAll plots generated successfully!")


if __name__ == "__main__":
    viz = LDCVisualizer()
    viz.generate_all_plots()