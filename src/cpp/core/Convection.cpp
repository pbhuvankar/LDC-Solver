#include "LDC.hpp"

// Convection terms using 2nd order upwind interpolation
void convection(std::vector<double>& du1, std::vector<double>& dv1, 
                //const 
                Grid& grid, const Momentum& mom) {
    const int Nxt= grid.Nxt;
    const int is= grid.is;
    const int ie= grid.ie;
    const int ieu= grid.ieu;
    const int js= grid.js;
    const int je= grid.je;
    const int jev= grid.jev;
    const double dt= mom.dt;

    double time1, time2;
    #ifdef _OPENMP
    time1 = omp_get_wtime();
    #endif
    
    std::vector<double> flux1(grid.Nt, 0.0);
    std::vector<double> flux2(grid.Nt, 0.0);
    
    std::fill(du1.begin(), du1.end(), 0.0);
    std::fill(dv1.begin(), dv1.end(), 0.0);
    
    // U-momentum convection - PARALLELIZED
    // merged for-loop parallelization
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i= is; i<= ie; i++) 
    {
        for (int j= js-1; j<= je; j++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            // X-direction flux
            if (mom.u[ij] + mom.u[ij_k(i-1, j, Nxt)] > 0.0) 
            {
                flux1[ij]= interp(grid.xh[i], grid.xh[i-1], grid.xh[i-2],
                                   mom.u[ij], mom.u[ij_k(i-1, j, Nxt)], 
                                   mom.u[ij_k(i-2, j, Nxt)], grid.x[i]);
            } 
            else 
            {
                flux1[ij]= interp(grid.xh[i-1], grid.xh[i], grid.xh[i+1],
                                   mom.u[ij_k(i-1, j, Nxt)], mom.u[ij], 
                                   mom.u[ij_k(i+1, j, Nxt)], grid.x[i]);
            }
            
            // Y-direction flux
            if (mom.v[ij] + mom.v[ij_k(i+1, j, Nxt)] > 0.0) 
            {
                flux2[ij]= interp(grid.y[j], grid.y[j-1], grid.y[j+1],
                                   mom.u[ij], mom.u[ij_k(i, j-1, Nxt)], 
                                   mom.u[ij_k(i, j+1, Nxt)], grid.yh[j]);
            } 
            else 
            {
                flux2[ij]= interp(grid.y[j], grid.y[j+1], grid.y[j+2],
                                   mom.u[ij], mom.u[ij_k(i, j+1, Nxt)], 
                                   mom.u[ij_k(i, j+2, Nxt)], grid.yh[j]);
            }
        }
    }
    
    // Calculate U-momentum convection terms - PARALLELIZED
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i= is; i<= ieu; i++) 
    {
        for (int j= js; j<= je; j++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            du1[ij]= ((flux1[ij_k(i+1, j, Nxt)] * flux1[ij_k(i+1, j, Nxt)]) - 
                       (flux1[ij] * flux1[ij])) / grid.dxh[i] +
                      ((0.5 * flux2[ij] * (mom.v[ij] + mom.v[ij_k(i+1, j, Nxt)])) - 
                       (0.5 * flux2[ij_k(i, j-1, Nxt)] * 
                        (mom.v[ij_k(i, j-1, Nxt)] + mom.v[ij_k(i+1, j-1, Nxt)]))) / 
                       grid.dy[j];
        }
    }
    
    std::fill(flux1.begin(), flux1.end(), 0.0);
    std::fill(flux2.begin(), flux2.end(), 0.0);
    
    // V-momentum convection 
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j= js; j<= je; j++) 
    {
        for (int i= is-1; i<= ie; i++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            // Y-direction flux
            if (mom.v[ij] + mom.v[ij_k(i, j-1, Nxt)] > 0.0) 
            {
                flux1[ij]= interp(grid.yh[j], grid.yh[j-1], grid.yh[j-2],
                                   mom.v[ij], mom.v[ij_k(i, j-1, Nxt)], 
                                   mom.v[ij_k(i, j-2, Nxt)], grid.y[j]);
            } 
            else 
            {
                flux1[ij]= interp(grid.yh[j-1], grid.yh[j], grid.yh[j+1],
                                   mom.v[ij_k(i, j-1, Nxt)], mom.v[ij], 
                                   mom.v[ij_k(i, j+1, Nxt)], grid.y[j]);
            }
            
            // X-direction flux
            if (mom.u[ij] + mom.u[ij_k(i, j+1, Nxt)] > 0.0) 
            {
                flux2[ij]= interp(grid.x[i], grid.x[i-1], grid.x[i+1],
                                   mom.v[ij], mom.v[ij_k(i-1, j, Nxt)], 
                                   mom.v[ij_k(i+1, j, Nxt)], grid.xh[i]);
            } 
            else 
            {
                flux2[ij]= interp(grid.x[i], grid.x[i+1], grid.x[i+2],
                                   mom.v[ij], mom.v[ij_k(i+1, j, Nxt)], 
                                   mom.v[ij_k(i+2, j, Nxt)], grid.xh[i]);
            }
        }
    }
    
    // Calculate V-momentum convection terms - PARALLELIZED
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j= js; j<= jev; j++) 
    {
        for (int i= is; i<= ie; i++) 
        {
            int ij= ij_k(i, j, Nxt);
            
            dv1[ij]= ((flux1[ij_k(i, j+1, Nxt)] * flux1[ij_k(i, j+1, Nxt)]) - 
                       (flux1[ij] * flux1[ij])) / grid.dyh[j] +
                      ((0.5 * flux2[ij] * (mom.u[ij] + mom.u[ij_k(i, j+1, Nxt)])) - 
                       (0.5 * flux2[ij_k(i-1, j, Nxt)] * 
                        (mom.u[ij_k(i-1, j, Nxt)] + mom.u[ij_k(i-1, j+1, Nxt)]))) / 
                       grid.dx[i];
        }
    }
    
    //  = (du/dt) * delta_t
    #pragma omp parallel for schedule(static)
    for (size_t i= 0; i< du1.size(); i++) 
    {
        du1[i] *= -dt;
        dv1[i] *= -dt;
    }
    #ifdef _OPENMP
    time2 = omp_get_wtime();
    #endif
    grid.conv_time += time2 - time1;
}
