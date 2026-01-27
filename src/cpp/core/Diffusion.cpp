#include "LDC.hpp"
// Central differencing scheme
// U-momentum diffusion
void diffusionU(std::vector<double>& du1, int flg, Grid& grid, Momentum& mom, Solver& solver) 
{
    const int Nxt= grid.Nxt;
    const int Nyt= grid.Nyt;
    const int is= grid.is;
    const int ieu= grid.ieu;
    const int js= grid.js;
    const int je= grid.je;
    const double dt= mom.dt;
    const double mu1= mom.mu1;

    double time1, time2;
    #ifdef _OPENMP
    time1 = omp_get_wtime();
    #endif
    
    std::fill(solver.A1.begin(), solver.A1.end(), 0.0);
    std::fill(solver.A2.begin(), solver.A2.end(), 0.0);
    std::fill(solver.A3.begin(), solver.A3.end(), 0.0);
    std::fill(solver.A4.begin(), solver.A4.end(), 0.0);
    std::fill(solver.A5.begin(), solver.A5.end(), 0.0);
    std::fill(solver.A6.begin(), solver.A6.end(), 0.0);
    
    // Calculate diffusion coefficients
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i= is; i<= ieu; i++) 
    {
        for (int j= js; j<= je; j++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            solver.A1[ij]= dt * mu1 / (grid.dxh[i] * grid.dx[i]);
            solver.A2[ij]= dt * mu1 / (grid.dxh[i] * grid.dx[i+1]);
            solver.A3[ij]= mu1 * dt / (grid.dyh[j-1] * grid.dy[j]);
            solver.A4[ij]= mu1 * dt / (grid.dyh[j] * grid.dy[j]);
            
            solver.A5[ij]= -(solver.A1[ij] + solver.A2[ij] + 
                             solver.A3[ij] + solver.A4[ij]);
            
            if (flg == 0) 
            {
                // Explicit diffusion calculation
                du1[ij]= du1[ij] + 
                         solver.A1[ij] * mom.u[ij_k(i-1, j, Nxt)] + 
                         solver.A2[ij] * mom.u[ij_k(i+1, j, Nxt)] + 
                         solver.A3[ij] * mom.u[ij_k(i, j-1, Nxt)] + 
                         solver.A4[ij] * mom.u[ij_k(i, j+1, Nxt)] + 
                         solver.A5[ij] * mom.u[ij];
            }
            
            if (flg == 1) 
            {
                // Implicit setup
                solver.A5[ij]= 1.0 - solver.A5[ij];
                solver.A6[ij]= mom.u[ij] + du1[ij];
            }
        }
    }
    
    if (flg == 1) 
    {
        // Boundary conditions for implicit solver
        for (int j= 1; j<= Nyt; j++) 
        {
            solver.A1[ij_k(is, j, Nxt)]= 0.0;
            solver.A2[ij_k(ieu, j, Nxt)]= 0.0;
        }
        
        for (int i= 1; i<= Nxt; i++) 
        {
            // Bottom wall
            int ij_js= ij_k(i, js, Nxt);
            solver.A5[ij_js]= solver.A5[ij_js] + solver.A3[ij_js];
            solver.A6[ij_js]= solver.A6[ij_js] + (2.0 * 0.0 * solver.A3[ij_js]);
            solver.A3[ij_js]= 0.0;
            
            // Top wall (velocity= 1.0---Hardcoded)
            int ij_je= ij_k(i, je, Nxt);
            solver.A5[ij_je]= solver.A5[ij_je] + solver.A4[ij_je];
            solver.A6[ij_je]= solver.A6[ij_je] + (2.0 * 1.0 * solver.A4[ij_je]);
            solver.A4[ij_je]= 0.0;
        }
    }
    
    mom.ind[0]= is;
    mom.ind[1]= ieu;
    mom.ind[2]= js;
    mom.ind[3]= je;

    #ifdef _OPENMP
    time2 = omp_get_wtime();
    #endif

    grid.diff_time += time2 - time1;
}

// V-momentum diffusion - OPENMP OPTIMIZED
void diffusionV(std::vector<double>& dv1, int flg, Grid& grid, Momentum& mom, Solver& solver) 
{
    const int Nxt= grid.Nxt;
    const int Nyt= grid.Nyt;
    const int is= grid.is;
    const int ie= grid.ie;
    const int js= grid.js;
    const int jev= grid.jev;
    const double dt= mom.dt;
    const double mu1= mom.mu1;

    double time1, time2;

    #ifdef _OPENMP
    time1 = omp_get_wtime();
    #endif
    
    std::fill(solver.A1.begin(), solver.A1.end(), 0.0);
    std::fill(solver.A2.begin(), solver.A2.end(), 0.0);
    std::fill(solver.A3.begin(), solver.A3.end(), 0.0);
    std::fill(solver.A4.begin(), solver.A4.end(), 0.0);
    std::fill(solver.A5.begin(), solver.A5.end(), 0.0);
    std::fill(solver.A6.begin(), solver.A6.end(), 0.0);
    
    // Calculate diffusion coefficients - PARALLELIZED
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j= js; j<= jev; j++) 
    {
        for (int i= is; i<= ie; i++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            solver.A3[ij]= dt * mu1 / (grid.dyh[j] * grid.dy[j]);
            solver.A4[ij]= dt * mu1 / (grid.dyh[j] * grid.dy[j+1]);
            solver.A1[ij]= mu1 * dt / (grid.dxh[i-1] * grid.dx[i]);
            solver.A2[ij]= mu1 * dt / (grid.dxh[i] * grid.dx[i]);
            
            solver.A5[ij]= -(solver.A1[ij] + solver.A2[ij] + 
                             solver.A3[ij] + solver.A4[ij]);
            
            if (flg == 0) 
            {
                // Explicit diffusion calculation
                dv1[ij]= dv1[ij] + 
                         solver.A1[ij] * mom.v[ij_k(i-1, j, Nxt)] + 
                         solver.A2[ij] * mom.v[ij_k(i+1, j, Nxt)] + 
                         solver.A3[ij] * mom.v[ij_k(i, j-1, Nxt)] + 
                         solver.A4[ij] * mom.v[ij_k(i, j+1, Nxt)] + 
                         solver.A5[ij] * mom.v[ij];
            }
            
            if (flg == 1) // Implicit diffusion calculation
            {
                // Implicit equation system diagonal coefficient & rhs setup
                solver.A5[ij]= 1.0 - solver.A5[ij];
                solver.A6[ij]= mom.v[ij] + dv1[ij];
            }
        }
    }
    
    if (flg == 1) 
    {
        // Boundary conditions for implicit solver        
        for (int i= 1; i<= Nxt; i++) 
        {
            solver.A3[ij_k(i, js, Nxt)]= 0.0;
            solver.A4[ij_k(i, jev, Nxt)]= 0.0;
        }
        
        for (int j= 1; j<= Nyt; j++) 
        {
            // Left wall----zero tangential velocity hard-coded
            int ij_is= ij_k(is, j, Nxt);
            solver.A5[ij_is]= solver.A5[ij_is] + solver.A1[ij_is];
            solver.A6[ij_is]= solver.A6[ij_is] + (2.0 * 0.0 * solver.A1[ij_is]);
            solver.A1[ij_is]= 0.0;
            
            // Right wall----zero tangential velocity hard-coded
            int ij_ie= ij_k(ie, j, Nxt);
            solver.A5[ij_ie]= solver.A5[ij_ie] + solver.A2[ij_ie];
            solver.A6[ij_ie]= solver.A6[ij_ie] + (2.0 * 0.0 * solver.A2[ij_ie]);
            solver.A2[ij_ie]= 0.0;
        }
    }
    
    mom.ind[0]= is;
    mom.ind[1]= ie;
    mom.ind[2]= js;
    mom.ind[3]= jev;

    #ifdef _OPENMP
    time2 = omp_get_wtime();
    #endif

    grid.diff_time += time2 - time1;
}
