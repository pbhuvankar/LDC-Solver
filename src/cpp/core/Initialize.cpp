#include "LDC.hpp"
#include<sstream>

void read_params(Grid& grid, Momentum& mom, Solver& solver) 
{
    // Reading input parameters
    std::ifstream infile("input");
    if (!infile.is_open()) 
    {
        std::cerr<< "Error: Cannot open input file"<< std::endl;
        exit(1);
    }
    
    std::string line;
    while (std::getline(infile, line)) 
    {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '/' || line[0] == '&' || line[0] == '*') 
            continue;
            
        std::istringstream iss(line);
        std::string key;
        
        if (std::getline(iss, key, '=')) 
        {            
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            
            if (key == "Nx") iss >> grid.Nx;
            else if (key == "Ny") iss >> grid.Ny;
            else if (key == "Xlen") iss >> grid.Xlen;
            else if (key == "Ylen") iss >> grid.Ylen;
            else if (key == "mesh_ref") iss >> grid.mesh_ref;
            else if (key == "cfl") iss >> grid.cfl;
            else if (key == "Tend") iss >> mom.Tend;
            else if (key == "Max_dt") iss >> grid.Max_dt;
            else if (key == "mu1") iss >> mom.mu1;
            else if (key == "Beta") iss >> grid.beta;
            else if (key == "maxit") iss >> solver.maxit;
            else if (key == "t_scheme") iss >> solver.t_scheme;
            else if (key == "which_solver") iss >> solver.which_solver;
            else if (key == "Semi_implicit") 
            {
                char val;
                iss >> val;
                grid.Semi_implicit= (val == 'T' || val == 't');
            }
            else if (key == "Breakdown_Time") 
            {
                char val;
                iss >> val;
                grid.Breakdown_Time= (val == 'T' || val == 't');
            }
        }
    }
    infile.close();
}
void initialize(Grid& grid, Momentum& mom, Solver& solver) 
{    
    // Set grid indices
    grid.is= 3;
    grid.ie= grid.Nx + 2;
    grid.Nxt= grid.Nx + 5;
    
    grid.js= 3;
    grid.je= grid.Ny + 2;
    grid.Nyt= grid.Ny + 5;
    
    grid.Nt= grid.Nxt * grid.Nyt;
    
    grid.ieu= grid.ie - 1;
    grid.jev= grid.je - 1;
    
    grid.Sflag= 0;
    if (grid.Semi_implicit) grid.Sflag= 1;
    
    // Allocating arrays
    mom.u.resize(grid.Nt, 0.0);
    mom.v.resize(grid.Nt, 0.0);
    mom.p.resize(grid.Nt, 0.0);
    mom.u0.resize(grid.Nt, 0.0);
    mom.v0.resize(grid.Nt, 0.0);
    mom.p0.resize(grid.Nt, 0.0);
    mom.du.resize(grid.Nt, 0.0);
    mom.dv.resize(grid.Nt, 0.0);
    mom.mask.resize(grid.Nt, 0.0);

    mom.u_old.resize(grid.Nt, 0.0);
    mom.v_old.resize(grid.Nt, 0.0);
    mom.p_old.resize(grid.Nt, 0.0);
    
    grid.x.resize(grid.Nxt, 0.0);
    grid.y.resize(grid.Nyt, 0.0);
    grid.xh.resize(grid.Nxt, 0.0);
    grid.yh.resize(grid.Nyt, 0.0);
    grid.dx.resize(grid.Nxt, 0.0);
    grid.dy.resize(grid.Nyt, 0.0);
    grid.dxh.resize(grid.Nxt, 0.0);
    grid.dyh.resize(grid.Nyt, 0.0);
    
    solver.A1.resize(grid.Nt, 0.0);
    solver.A2.resize(grid.Nt, 0.0);
    solver.A3.resize(grid.Nt, 0.0);
    solver.A4.resize(grid.Nt, 0.0);
    solver.A5.resize(grid.Nt, 0.0);
    solver.A6.resize(grid.Nt, 0.0);
    
    
    //Creating grid for u-velocity (horizontal face locations). Uses staggered grid.
    for (int i= grid.is-1; i<= grid.ie; i++) 
    {
        double ep= static_cast<double>(i-2) / static_cast<double>(grid.ie-2);
        grid.xh[i]= (ep + (grid.mesh_ref * ep * (ep - 0.5) * (1.0 - ep))) * grid.Xlen;
    }
    grid.xh[grid.ie+1]=  2.0 * grid.xh[grid.ie] - grid.xh[grid.ie-1];
    grid.xh[grid.ie+2]=  2.0 * grid.xh[grid.ie] - grid.xh[grid.ie-2];
    grid.xh[grid.is-2]=  2.0 * grid.xh[grid.is-1] - grid.xh[grid.is];
    
    //Creating grid for v-velocity (vertical face locations). Uses staggered grid.
    for (int j= grid.js-1; j<= grid.je; j++) 
    {
        double ep= static_cast<double>(j-2) / static_cast<double>(grid.je-2);
        grid.yh[j]= (ep + (grid.mesh_ref * ep * (ep - 0.5) * (1.0 - ep))) * grid.Ylen;
    }
    grid.yh[grid.je+1]=  2.0 * grid.yh[grid.je] - grid.yh[grid.je-1];
    grid.yh[grid.je+2]=  2.0 * grid.yh[grid.je] - grid.yh[grid.je-2];
    grid.yh[grid.js-2]=  2.0 * grid.yh[grid.js-1] - grid.yh[grid.js];
    
    //Cell centers
    for (int i= grid.is; i<= grid.ie; i++) 
    {
        grid.x[i]=  0.5 * (grid.xh[i] + grid.xh[i-1]);
    }
    grid.x[grid.is-1]= -grid.x[grid.is];
    grid.x[grid.is-2]= -grid.x[grid.is+1];
    grid.x[grid.ie+1]=  2.0 * grid.Xlen - grid.x[grid.ie];
    grid.x[grid.ie+2]=  2.0 * grid.Xlen - grid.x[grid.ie-1];
    
    for (int j= grid.js; j<= grid.je; j++) 
    {
        grid.y[j]= 0.5 * (grid.yh[j] + grid.yh[j-1]);
    }
    grid.y[grid.js-1]= -grid.y[grid.js];
    grid.y[grid.js-2]= -grid.y[grid.js+1];
    grid.y[grid.je+1]=  2.0 * grid.Ylen - grid.y[grid.je];
    grid.y[grid.je+2]=  2.0 * grid.Ylen - grid.y[grid.je-1];
    
    //Cell widths computed
    for (int i= grid.is; i<= grid.ie; i++) 
    {
        grid.dx[i]= grid.xh[i] - grid.xh[i-1];
    }
    grid.dx[grid.is-1]=  grid.dx[grid.is];
    grid.dx[grid.is-2]=  grid.dx[grid.is+1];
    grid.dx[grid.ie+1]=  grid.dx[grid.ie];
    grid.dx[grid.ie+2]=  grid.dx[grid.ie-1];
    
    //Cell heights computed
    for (int j= grid.js; j<= grid.je; j++) 
    {
        grid.dy[j]=  grid.yh[j] - grid.yh[j-1];
    }
    grid.dy[grid.js-1]=  grid.dy[grid.js];
    grid.dy[grid.js-2]=  grid.dy[grid.js+1];
    grid.dy[grid.je+1]=  grid.dy[grid.je];
    grid.dy[grid.je+2]=  grid.dy[grid.je-1];
    
    //Widths for U-velocity control volumes
    for (int i= grid.is-2; i<= grid.ie+1; i++) 
    {
        grid.dxh[i]=  grid.x[i+1] - grid.x[i];
    }
    grid.dxh[grid.ie + 2]=  grid.dxh[grid.ie+1];
    
    //Heights for V-velocity control volumes
    for (int j= grid.js - 2; j<= grid.je+1; j++) 
    {
        grid.dyh[j]= grid.y[j+1] - grid.y[j];
    }
    grid.dyh[grid.je+2]=  grid.dyh[grid.je+1];
    
    //Initializing ghost mask for pressure Poisson BC
    for (int i= grid.is; i<= grid.ie; i++) 
    {
        for (int j= grid.js; j<= grid.je; j++) 
        {
            mom.mask[ij_k(i, j, grid.Nxt)]= 1.0;
        }
    }
    
    grid.iEnd= 100000000; //Large value. Adjustable.
    mom.time= 0.0;
    grid.solve_time= 0.0;
    grid.conv_time= 0.0;
    grid.diff_time= 0.0;
}
