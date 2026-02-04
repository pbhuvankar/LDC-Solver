#include "LDC.hpp"
#include <sys/stat.h>
#include <sys/types.h>

//Create directory for outputting data
void createDirectory(const std::string& dirname) {
    #ifdef _WIN32
        _mkdir(dirname.c_str());
    #else
        mkdir(dirname.c_str(), 0755);
    #endif
}

//Set up pressure Poisson equation coefficients
void SetPressurePoi(Grid& grid, Momentum& mom, Solver& solver) 
{
    const int Nxt=grid.Nxt;
    const int is=grid.is;
    const int ie=grid.ie;
    const int js=grid.js;
    const int je=grid.je;
    const double dt=mom.dt;
    
    // Initialize coefficient arrays
    std::fill(solver.A1.begin(), solver.A1.end(), 0.0);
    std::fill(solver.A2.begin(), solver.A2.end(), 0.0);
    std::fill(solver.A3.begin(), solver.A3.end(), 0.0);
    std::fill(solver.A4.begin(), solver.A4.end(), 0.0);
    std::fill(solver.A5.begin(), solver.A5.end(), 0.0);
    std::fill(solver.A6.begin(), solver.A6.end(), 0.0);
    
    //Calculate coefficients: 
    //Format:  A5 * p_i,j = A1 * p_i-1,j + A2 * p_i+1,j + A3 * p_i,j-1 + A4 * p_i,j+1 + A6
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i=is; i<=ie; i++) 
    {
        for (int j=js; j<=je; j++) 
        {
            int ij=ij_k(i, j, Nxt);
            
            solver.A1[ij]=mom.mask[ij_k(i-1, j, Nxt)] * (dt / (grid.dxh[i-1] * grid.dx[i]));
            solver.A2[ij]=mom.mask[ij_k(i+1, j, Nxt)] * (dt / (grid.dxh[i] * grid.dx[i]));
            solver.A3[ij]=mom.mask[ij_k(i, j-1, Nxt)] * (dt / (grid.dyh[j-1] * grid.dy[j]));
            solver.A4[ij]=mom.mask[ij_k(i, j+1, Nxt)] * (dt / (grid.dyh[j] * grid.dy[j]));
            
            solver.A6[ij]=-((mom.u0[ij] - mom.u0[ij_k(i-1, j, Nxt)]) / grid.dx[i]) - 
                            ((mom.v0[ij] - mom.v0[ij_k(i, j-1, Nxt)]) / grid.dy[j]);
            
            solver.A5[ij]=solver.A1[ij] + solver.A2[ij] + solver.A3[ij] + solver.A4[ij];
        }
    }
    
    //Rooting pressure (else rank deficient)
    int root_ij_k=ij_k(is+4, js+4, Nxt);
    solver.A1[root_ij_k]=0.0;
    solver.A2[root_ij_k]=0.0;
    solver.A3[root_ij_k]=0.0;
    solver.A4[root_ij_k]=0.0;
    solver.A5[root_ij_k]=1.0;//solver.A5[ij_k(is+2, js+2, Nxt)];
    solver.A6[root_ij_k]=0.0;
    
    mom.ind[0]=is;
    mom.ind[1]=ie;
    mom.ind[2]=js;
    mom.ind[3]=je;
}

//Successive Over-Relaxation (GS) Red-Black Type (For parallel implementation)
//Red-Black ordering 
void RedBlackSOR(std::vector<double>& p, const int ind[4], double& res, 
                 int& count, Grid& grid, const Solver& solver) 
{
    const int Nxt=grid.Nxt;
    const int Nx=grid.Nx;
    const int Ny=grid.Ny;
    const double beta=grid.beta;
    const int is1=ind[0];
    const int ie1=ind[1];
    const int js1=ind[2];
    const int je1=ind[3];
    
    count=0;
    res=100.0;

    double time1, time2;
    #ifdef _OPENMP
    time1 = omp_get_wtime();
    #endif
    
    //Iterate until convergence
    while (res > 1e-5 && count< solver.maxit) 
    {
        count++;
        
        
        //Points with (i+j) even are updated together
        //Points with (i+j) odd are updated together
        for (int color= 0; color< 2; color++) 
        {
            //All points of the same color can be updated in parallel
            #pragma omp parallel for collapse(2) schedule(static)
            for (int i=is1; i<=ie1; i++) 
            {
                for (int j=js1; j<=je1; j++) 
                {
                    if ((i + j) % 2 == color) 
                    {
                        int ij=ij_k(i, j, Nxt);
                        
                        p[ij]=beta * ((solver.A1[ij] * p[ij_k(i-1, j, Nxt)] + 
                                         solver.A2[ij] * p[ij_k(i+1, j, Nxt)] + 
                                         solver.A3[ij] * p[ij_k(i, j-1, Nxt)] + 
                                         solver.A4[ij] * p[ij_k(i, j+1, Nxt)] + 
                                         solver.A6[ij]) / solver.A5[ij]) + 
                                (1.0 - beta) * p[ij];
                    }
                }
            }
        }
        
        //Residue
        res=0.0;
        #pragma omp parallel for collapse(2) reduction(+:res) schedule(static)
        for (int i=is1; i<=ie1; i++) 
        {
            for (int j=js1; j<=je1; j++) 
            {
                int ij=ij_k(i, j, Nxt);
                
                res +=std::abs(solver.A1[ij] * p[ij_k(i-1, j, Nxt)] + 
                               solver.A2[ij] * p[ij_k(i+1, j, Nxt)] + 
                               solver.A3[ij] * p[ij_k(i, j-1, Nxt)] + 
                               solver.A4[ij] * p[ij_k(i, j+1, Nxt)] + 
                               solver.A6[ij] - 
                               solver.A5[ij] * p[ij]);
            }
        }
        res /=static_cast<double>(Nx * Ny);
    }
    #ifdef _OPENMP
    time2 = omp_get_wtime();
    #endif
    grid.solve_time += time2 - time1;    
}

//Standard Jacobi Over-Relaxation (JOR)
void JORSolver(std::vector<double>& p, const int ind[4], double& res, 
               int& count, Grid& grid, const Solver& solver) 
{
    const int Nxt=grid.Nxt;
    const int Nx=grid.Nx;
    const int Ny=grid.Ny;
    const double beta=grid.beta;
    const int is1=ind[0];
    const int ie1=ind[1];
    const int js1=ind[2];
    const int je1=ind[3];
    
    std::vector<double> p0(grid.Nt);
    count=0;
    res=100.0;
    int s_time, e_time;
    #ifdef _OPENMP
        s_time = omp_get_wtime();
    #endif
    
    while (res > 1e-5 && count< solver.maxit) 
    {
        count++;

        double s_time, e_time;        

        p0=p;
                
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i=is1; i<=ie1; i++) 
        {
            for (int j=js1; j<=je1; j++) 
            {
                int ij=ij_k(i, j, Nxt);
                
                p[ij]=beta * ((solver.A1[ij] * p0[ij-1] + 
                                 solver.A2[ij] * p0[ij+1] + 
                                 solver.A3[ij] * p0[ij-Nxt] + 
                                 solver.A4[ij] * p0[ij+Nxt] + 
                                 solver.A6[ij]) / solver.A5[ij]) + 
                        (1.0 - beta) * p0[ij];
            }
        }
        
        

        //Residual
        res=0.0;
        #pragma omp parallel for collapse(2) reduction(+:res) schedule(static)
        for (int i=is1; i<=ie1; i++) 
        {
            for (int j=js1; j<=je1; j++) 
            {
                int ij=ij_k(i, j, Nxt);
                
                res +=std::abs(solver.A1[ij] * p[ij-1] + 
                               solver.A2[ij] * p[ij+1] + 
                               solver.A3[ij] * p[ij-Nxt] + 
                               solver.A4[ij] * p[ij+Nxt] + 
                               solver.A6[ij] - 
                               solver.A5[ij] * p[ij]);
            }
        }
        res /=static_cast<double>(Nx * Ny);

    }
    #ifdef _OPENMP
        e_time = omp_get_wtime();
    #endif
        grid.solve_time += (e_time - s_time);
}

//Profiling
void TimeAnalysis(const Solver& solver, const Grid& grid)
{
    // Create Profile directory
    createDirectory("Profile");

    int N = grid.Nth;
    //#ifdef _OPENMP        
    //N = omp_get_num_threads();    
    //#else 
    //N = 1;
    //#endif
    
    std::ofstream tfile("Profile/Time_Breakdown"+std::to_string(N)+"_"+
    std::to_string(grid.Nx)+".txt");

    tfile<<"---------------Configuration: "<<grid.Nx<<"x"<<grid.Ny<<"-----------------\n";
    tfile<<"Total time steps= "<<grid.itime<<"\n";
    tfile<<"Total run time= "<<grid.prog_time<<"\n";
    tfile<<"Diffusion time= "<<grid.diff_time<<"\n";
    tfile<<"Convection time= "<<grid.conv_time<<"\n";
    tfile<<"Solver time= "<<grid.solve_time<<"\n";

    tfile.close();
    
    std::cout << "Profiling data saved to Profile/Time_Breakdown.txt" << std::endl;
}

//Output for flow visualization
void OutputCSV(const Grid& grid, const Momentum& mom) 
{
    const int Nxt=grid.Nxt;
    const int is=grid.is;
    const int ie=grid.ie;
    const int js=grid.js;
    const int je=grid.je;
    
    // Create Visualize directory
    createDirectory("Visualize");
    
    std::ofstream xfile("Visualize/output_x.csv");
    std::ofstream yfile("Visualize/output_y.csv");
    
    for (int i=is; i<=ie; i++) 
    {
        xfile<< grid.x[i];
        if (i< ie) xfile<< ",";
    }
    xfile.close();
    
    for (int j=js; j<=je; j++) 
    {
        yfile<< grid.y[j];
        if (j< je) yfile<< ",";
    }
    yfile.close();
    
    std::ofstream ufile("Visualize/output_u.csv");
    std::ofstream vfile("Visualize/output_v.csv");
    std::ofstream pfile("Visualize/output_p.csv");
    
    for (int j=js; j<=je; j++) 
    {
        for (int i=is; i<=ie; i++) 
        {
            int ij=ij_k(i, j, Nxt);
            
            double u_center=0.5 * (mom.u[ij] + mom.u[ij_k(i-1, j, Nxt)]);
            double v_center=0.5 * (mom.v[ij] + mom.v[ij_k(i, j-1, Nxt)]);
            
            ufile<< u_center;
            vfile<< v_center;
            pfile<< mom.p[ij];
            
            if (i< ie) 
            {
                ufile<< ",";
                vfile<< ",";
                pfile<< ",";
            }
        }
        if (j< je) 
        {
            ufile<< "\n";
            vfile<< "\n";
            pfile<< "\n";
        }
    }
    
    ufile.close();
    vfile.close();
    pfile.close();    
    std::cout<< "CSV output files created in Visualize/ directory" << std::endl;
}

//For Matlab visualization
void Output(const Grid& grid, const Momentum& mom) 
{
    const int Nxt=grid.Nxt;
    const int is=grid.is;
    const int ie=grid.ie;
    const int js=grid.js;
    const int je=grid.je;
    
    // Create Visualize directory
    createDirectory("Visualize");
    
    std::ofstream outfile("Visualize/out.m");
    outfile<< std::setprecision(16);
    
    outfile<< "clear all\n";
    outfile<< "is="<< is<< "; js="<< js<< ";\n";
    outfile<< "ie="<< ie<< "; je="<< je<< ";\n";
    
    for (int i=is; i<=ie; i++) 
    {
        outfile<< "x("<< i-is+1<< ")="<< grid.x[i]<< ";\n";
    }
    
    for (int i=js; i<=je; i++) 
    {
        outfile<< "y("<< i-js+1<< ")="<< grid.y[i]<< ";\n";
    }
    
    for (int i=is; i<=ie; i++) 
    {
        for (int j=js; j<=je; j++) 
        {
            int ii=i - is + 1;
            int jj=j - js + 1;
            int ij=ij_k(i, j, Nxt);
            
            outfile<< "p("<< ii<< ","<< jj<< ")="<< mom.p[ij]<< ";\n";
            outfile<< "u("<< ii<< ","<< jj<< ")=" 
                   << 0.5 * (mom.u[ij] + mom.u[ij_k(i-1, j, Nxt)])<< ";\n";
            outfile<< "v("<< ii<< ","<< jj<< ")=" 
                   << 0.5 * (mom.v[ij] + mom.v[ij_k(i, j-1, Nxt)])<< ";\n";
        }
    }
    
    outfile<< "[X,Y]=meshgrid(x,y);\n";
    outfile<< "figure(1);contourf(Y,X,p)\n";
    outfile<< "figure(2);contourf(Y,X,u)\n";
    outfile<< "figure(3);contourf(Y,X,v)\n";
    
    outfile.close();
    
    std::cout << "MATLAB visualization script saved to Visualize/out.m" << std::endl;
}
