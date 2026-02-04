#include "LDC.hpp"

double err_calc(std::vector<double>& , std::vector<double>& , Grid&, const double, const int, const int,
    const int, const int);
void reset_val(Momentum& );
void output_centerline_data(const Grid& grid, const Momentum& mom, const Solver& solver);

int main() 
{
    Grid grid;
    Momentum mom;
    Solver solver;
    
    //Initialize simulation
    read_params(grid, mom, solver);
    initialize(grid, mom, solver);

    mom.mu1 = 1.0;
    grid.beta = 1.7;
    solver.maxit = 10000000;

    double cu, cv, diffu, diffv, pe;
    
    std::cout<<"========================================\n";
    std::cout<<"  Component Tests\n";
    std::cout<<"========================================\n";
    std::cout<<"Grid: "<<grid.Nx<<" x "<<grid.Ny<<"\n";
    std::cout<<"----------------------------------------\n\n";

    mom.dt = 1.0;
    
    // ========== CONVECTION TEST ==========
    std::cout<<"1. Testing Convection Scheme...\n";
    reset_val(mom);        

    for (int i=grid.is-2; i<=grid.ie+2; i++)
    {
        for (int j=grid.js-2; j<=grid.je+2; j++)
        {
            int ij= ij_k(i, j, grid.Nxt);
            mom.u[ij] = grid.xh[i]*(1.0 - grid.xh[i])*grid.y[j]*(1.0 - grid.y[j]); 
            mom.v[ij] = 2.0*grid.x[i]*(1.0 - grid.x[i])*grid.yh[j]*(1.0 - grid.yh[j]);    
            
            mom.u0[ij] =  -2.0*mom.u[ij]*( (1.0 - grid.xh[i])*(1.0 - grid.y[j])*(grid.y[j] + 2.0*grid.xh[i]) - 
                grid.xh[i]*grid.y[j]*(3.0 - grid.y[j] - 2.0*grid.xh[i]) );
            
            mom.v0[ij] =  -2.0*mom.v[ij]*( (1.0 - grid.x[i])*(1.0 - grid.yh[j])*(grid.yh[j] + 2.0*grid.x[i]) - 
                grid.x[i]*grid.yh[j]*(3.0 - grid.yh[j] - 2.0*grid.x[i]) );
        }
    }

    convection(mom.du, mom.dv, grid, mom);

    cu = err_calc(mom.du, mom.u0, grid, mom.dt, grid.is, grid.ie-1, grid.js, grid.je);
    cv = err_calc(mom.dv, mom.v0, grid, mom.dt, grid.is, grid.ie, grid.js, grid.je-1);

    std::cout<<"  U-momentum L1 error: "<<std::scientific<<std::setprecision(6)<<cu<<"\n";
    std::cout<<"  V-momentum L1 error: "<<cv<<"\n\n";

    //Save convection data for visualization
    {
        std::ofstream outfile("validation_convection.dat");
        outfile<<std::setprecision(16);
        outfile<<"# Convection Test - Centerline Data\n";
        outfile<<"# x, y, u_numerical, u_analytical, v_numerical, v_analytical\n";
        
        // Vertical centerline (constant i)
        int i_center = grid.is + grid.Nx/2;
        outfile<<"# VERTICAL_CENTERLINE (x="<<grid.x[i_center]<<")\n";
        for (int j=grid.js; j<=grid.je; j++)
        {
            int ij = ij_k(i_center, j, grid.Nxt);
            outfile<<grid.x[i_center]<<" "<<grid.y[j]<<" "
                   <<mom.du[ij]<<" "<<mom.u0[ij]<<" "
                   <<mom.dv[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        
        // Horizontal centerline (constant j)
        int j_center = grid.js + grid.Ny/2;
        outfile<<"\n# HORIZONTAL_CENTERLINE (y="<<grid.y[j_center]<<")\n";
        for (int i=grid.is; i<=grid.ie; i++)
        {
            int ij = ij_k(i, j_center, grid.Nxt);
            outfile<<grid.x[i]<<" "<<grid.y[j_center]<<" "
                   <<mom.du[ij]<<" "<<mom.u0[ij]<<" "
                   <<mom.dv[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        outfile.close();
    }

    //========== DIFFUSION TEST ==========
    std::cout<<"2. Testing Diffusion Scheme...\n";
    reset_val(mom);

    // Set up initial field for diffusion test
    for (int i=grid.is-2; i<=grid.ie+2; i++)
    {
        for (int j=grid.js-2; j<=grid.je+2; j++)
        {
            int ij= ij_k(i, j, grid.Nxt);
            mom.u[ij] = grid.xh[i]*(1.0 - grid.xh[i])*grid.y[j]*(1.0 - grid.y[j]); 
            mom.v[ij] = 2.0*grid.x[i]*(1.0 - grid.x[i])*grid.yh[j]*(1.0 - grid.yh[j]);    
        }
    }

    diffusionU(mom.du, grid.Sflag, grid, mom, solver);
    diffusionV(mom.dv, grid.Sflag, grid, mom, solver);

    for (int i=grid.is; i<=grid.ie; i++)
    {
        for (int j=grid.js; j<=grid.je; j++)
        {
            int ij= ij_k(i, j, grid.Nxt);                                
            mom.u0[ij] =  -2.0*(grid.y[j]*(1.0-grid.y[j]) + grid.xh[i]*(1.0-grid.xh[i]));
            mom.v0[ij] =  -4.0*(grid.yh[j]*(1.0-grid.yh[j]) + grid.x[i]*(1.0-grid.x[i]));                   
        }
    }
    
    diffu = err_calc(mom.du, mom.u0, grid, mom.dt, grid.is, grid.ie-1, grid.js, grid.je);
    diffv = err_calc(mom.dv, mom.v0, grid, mom.dt, grid.is, grid.ie, grid.js, grid.je-1);

    std::cout<<"  U-momentum L1 error: "<<std::scientific<<std::setprecision(6)<<diffu<<"\n";
    std::cout<<"  V-momentum L1 error: "<<diffv<<"\n\n";

    //Save diffusion data for visualization
    {
        std::ofstream outfile("validation_diffusion.dat");
        outfile<<std::setprecision(16);
        outfile<<"# Diffusion Test - Centerline Data\n";
        outfile<<"# x, y, u_numerical, u_analytical, v_numerical, v_analytical\n";
        
        // Vertical centerline
        int i_center = grid.is + grid.Nx/2;
        outfile<<"# VERTICAL_CENTERLINE (x="<<grid.x[i_center]<<")\n";
        for (int j=grid.js; j<=grid.je-1; j++)
        {
            int ij = ij_k(i_center, j, grid.Nxt);
            outfile<<grid.x[i_center]<<" "<<grid.y[j]<<" "
                   <<mom.du[ij]<<" "<<mom.u0[ij]<<" "
                   <<mom.dv[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        
        //Horizontal centerline
        int j_center = grid.js + grid.Ny/2;
        outfile<<"\n# HORIZONTAL_CENTERLINE (y="<<grid.y[j_center]<<")\n";
        for (int i=grid.is; i<=grid.ie-1; i++)
        {
            int ij = ij_k(i, j_center, grid.Nxt);
            outfile<<grid.x[i]<<" "<<grid.y[j_center]<<" "
                   <<mom.du[ij]<<" "<<mom.u0[ij]<<" "
                   <<mom.dv[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        outfile.close();
    }

    //========== POISSON EQUATION TEST (FIXED) ==========
    std::cout<<"3. Testing Pressure Poisson Solver...\n";
    reset_val(mom);

    SetPressurePoi(grid, mom, solver);

    const double pi = 3.14159265358979323846;
    
    // FIX: Use FIXED physical location for rooting
    const double x_root = grid.x[grid.is+4];//0.5;
    const double y_root = grid.y[grid.js+4];;//0.5;
    double rootval = cos(2.0*pi*x_root)*cos(2.0*pi*y_root);
    
    for (int i=grid.is; i<=grid.ie; i++)
    {
        for (int j=grid.js; j<=grid.je; j++)
        {
            int ij= ij_k(i, j, grid.Nxt);                                
            solver.A6[ij] = 8.0*pi*pi*cos(2.0*pi*grid.x[i])*cos(2.0*pi*grid.y[j]);
            mom.v0[ij] = cos(2.0*pi*grid.x[i])*cos(2.0*pi*grid.y[j]) - rootval;
        }
    }
    
    int root_ij_k = ij_k(grid.is+4, grid.js+4, grid.Nxt);
    solver.A6[root_ij_k] = 0.0;

    double res;
    int ct;
    RedBlackSOR(mom.p, mom.ind, res, ct, grid, solver);
    
    pe = err_calc(mom.p, mom.v0, grid, mom.dt, grid.is, grid.ie, grid.js, grid.je);

    std::cout<<"  Iterations: "<<ct<<"\n";
    std::cout<<"  Final residual: "<<std::scientific<<std::setprecision(6)<<res<<"\n";
    std::cout<<"  L1 error: "<<pe<<"\n\n";

    // Save Poisson data for visualization
    {
        std::ofstream outfile("validation_poisson.dat");
        outfile<<std::setprecision(16);
        outfile<<"# Poisson Test - Centerline Data\n";
        outfile<<"# x, y, p_numerical, p_analytical\n";
        
        // Vertical centerline
        int i_center = grid.is + grid.Nx/2;
        outfile<<"# VERTICAL_CENTERLINE (x="<<grid.x[i_center]<<")\n";
        for (int j=grid.js; j<=grid.je; j++)
        {
            int ij = ij_k(i_center, j, grid.Nxt);
            outfile<<grid.x[i_center]<<" "<<grid.y[j]<<" "
                   <<mom.p[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        
        // Horizontal centerline
        int j_center = grid.js + grid.Ny/2;
        outfile<<"\n# HORIZONTAL_CENTERLINE (y="<<grid.y[j_center]<<")\n";
        for (int i=grid.is; i<=grid.ie; i++)
        {
            int ij = ij_k(i, j_center, grid.Nxt);
            outfile<<grid.x[i]<<" "<<grid.y[j_center]<<" "
                   <<mom.p[ij]<<" "<<mom.v0[ij]<<"\n";
        }
        outfile.close();
    }

    // ========== SUMMARY ==========
    std::cout<<"========================================\n";
    std::cout<<"           VALIDATION SUMMARY\n";
    std::cout<<"========================================\n";
    std::cout<<std::left<<std::setw(25)<<"Test"<<std::setw(15)<<"L1 Error\n";
    std::cout<<"----------------------------------------\n";
    std::cout<<std::setw(25)<<"Convection U"<<std::scientific<<std::setprecision(4)<<cu<<"\n";
    std::cout<<std::setw(25)<<"Convection V"<<cv<<"\n";
    std::cout<<std::setw(25)<<"Diffusion U"<<diffu<<"\n";
    std::cout<<std::setw(25)<<"Diffusion V"<<diffv<<"\n";
    std::cout<<std::setw(25)<<"Poisson Pressure"<<pe<<"\n";
    std::cout<<"========================================\n\n";
    
    std::cout<<"Output files created:\n";
    std::cout<<"  - validation_convection.dat\n";
    std::cout<<"  - validation_diffusion.dat\n";
    std::cout<<"  - validation_poisson.dat\n\n";
    std::cout<<"Run 'python3 visualize_validation.py' to generate plots.\n";

    return 0;
}

double err_calc(std::vector<double>& phi, std::vector<double>& phi_a, Grid& grid, const double dt, const int i1,
    const int i2, const int j1, const int j2)
{
    double err = 0.0;
    for (int i=i1; i<=i2; i++)
    {
        for (int j=j1; j<=j2; j++)
        {
            int ij= ij_k(i, j, grid.Nxt);
            err += std::abs((phi[ij]) - phi_a[ij]);
        }
    }
    err /= static_cast<double>((i2-i1+1)*(j2-j1+1));
    return err;  //L1 norm
}

void reset_val(Momentum& mom)
{
    std::fill(mom.du.begin(), mom.du.end(), 0.0);
    std::fill(mom.dv.begin(), mom.dv.end(), 0.0); 
    std::fill(mom.u0.begin(), mom.u0.end(), 0.0); 
    std::fill(mom.v0.begin(), mom.v0.end(), 0.0); 
    std::fill(mom.p.begin(), mom.p.end(), 0.0); 
}
