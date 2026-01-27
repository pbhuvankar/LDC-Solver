#include "LDC.hpp"

// Projection step - correct velocities with pressure gradient - OPENMP OPTIMIZED
void projection(Grid& grid, Momentum& mom) 
{
    const int Nxt= grid.Nxt;
    const int is= grid.is;
    const int ie= grid.ie;
    const int ieu= grid.ieu;
    const int js= grid.js;
    const int je= grid.je;
    const int jev= grid.jev;
    const double dt= mom.dt;
    
    // Update U velocity - PARALLELIZED
    // Independent updates, perfect for parallelization
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i= is; i<= ieu; i++) 
    {
        for (int j= js; j<= je; j++) 
        {
            int ij= ij_k(i, j, Nxt);
            mom.u[ij]= mom.u0[ij] + 
                       (2.0 * dt * (mom.p[ij] - mom.p[ij_k(i+1, j, Nxt)]) / 
                        grid.dxh[i]);
        }
    }
    
    // Update V velocity - PARALLELIZED
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j= js; j<= jev; j++) 
    {
        for (int i= is; i<= ie; i++) 
        {
            int ij= ij_k(i, j, Nxt);
            mom.v[ij]= mom.v0[ij] + 
                       (2.0 * dt * (mom.p[ij] - mom.p[ij_k(i, j+1, Nxt)]) / 
                        grid.dyh[j]);
        }
    }
}

// Calculate adaptive time step based on CFL condition - OPENMP OPTIMIZED
void timestep(Grid& grid, Momentum& mom) 
{
    const int is= grid.is;
    const int ie= grid.ie;
    const int js= grid.js;
    const int je= grid.je;
    const int Nt= grid.Nt;
    
    // Find minimum grid spacing - sequential (small overhead)
    double mindx= *std::min_element(grid.dx.begin() + is, grid.dx.begin() + ie + 1);
    double mindy= *std::min_element(grid.dy.begin() + js, grid.dy.begin() + je + 1);
    mindx= std::min(mindx, mindy);
    
    // Find maximum velocity magnitude - PARALLELIZED WITH REDUCTION
    // reduction(max:maxu) ensures thread-safe maximum finding
    double maxu= 0.0;
    #pragma omp parallel for reduction(max:maxu) schedule(static)
    for (int i= 0; i< Nt; i++) 
    {
        double vel_mag= std::sqrt(mom.u[i]*mom.u[i] + mom.v[i]*mom.v[i]);
        maxu= std::max(maxu, vel_mag);
    }
    
    // CFL-based time step
    double dtc= mindx / (maxu + 1e-16);  // Add small epsilon to avoid division by zero
    
    // Diffusion-based time step
    double dtn= 100000.0;
    if (!grid.Semi_implicit) 
    {
        dtn= (mindx * mindx) / mom.mu1;
    }
    
    // Choose time step
    mom.dt= grid.cfl * std::min(dtc, dtn);
    mom.dt= std::min(mom.dt, grid.Max_dt);
}
