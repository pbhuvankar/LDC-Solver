#include "LDC.hpp"

int main() {
    Grid grid;
    Momentum mom;
    Solver solver;
    
    //Initialize simulation
    read_params(grid, mom, solver);

    initialize(grid, mom, solver);
    
    //Print OpenMP information
    grid.Nth = 1;
    #ifdef _OPENMP
    int num_threads = 1;

    #pragma omp parallel
    {
        #pragma omp master
        num_threads = omp_get_num_threads();
    }    
    std::cout << "OpenMP ENABLED with " << num_threads << " threads" << std::endl;
    std::cout << "=================================================" << std::endl;

    grid.Nth = num_threads;
    #else
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Running SEQUENTIAL version" << std::endl;
    std::cout << "Compile with 'make openmp' for parallel" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    #endif
    
    //Set boundary conditions
    double u_w = 1.0;
    UVBoundaryCond(mom.u, mom.v, grid, u_w);
    
    //Calculate time step
    timestep(grid, mom);
    
    grid.itime=0;
    mom.time=0.0;

    double start_time, end_time, prog_start_time, prog_end_time;
    
    #ifdef _OPENMP
    start_time = omp_get_wtime();
    prog_start_time = start_time;
    #endif    

    //Main time loop
    while (mom.time<=mom.Tend && grid.itime<=grid.iEnd) 
    {
        grid.itime++;
        mom.time +=mom.dt;

        double t_count = 0;
        RecordOldFields(grid, mom);

        while(t_count < solver.t_scheme)
        {
            t_count++;
                
            //Reset dU & dV
            std::fill(mom.du.begin(), mom.du.end(), 0.0);
            std::fill(mom.dv.begin(), mom.dv.end(), 0.0);        
        
            //Convection
            convection(mom.du, mom.dv, grid, mom);        
            
            //Diffusion for U-velocity
            diffusionU(mom.du, grid.Sflag, grid, mom, solver);        
            
            if (grid.Semi_implicit) {
                double res;
                int ct;
                //JORSolver(mom.u0, mom.ind, res, ct, grid, solver);
                RedBlackSOR(mom.u0, mom.ind, res, ct, grid, solver);
            }
        
            //Diffusion for V-velocity
            diffusionV(mom.dv, grid.Sflag, grid, mom, solver);        
        
            if (grid.Semi_implicit) {
                double res;
                int ct;
                //JORSolver(mom.v0, mom.ind, res, ct, grid, solver);
                RedBlackSOR(mom.v0, mom.ind, res, ct, grid, solver);
            }
        

            //Update intermediate velocities
            if (!grid.Semi_implicit) {
                for (size_t i=0; i< mom.u.size(); i++) {
                    mom.u0[i]=mom.u[i] + mom.du[i];
                    mom.v0[i]=mom.v[i] + mom.dv[i];
                }
            }
        
            //Apply boundary conditions to intermediate velocities
            double u_w = 1.0;
            UVBoundaryCond(mom.u0, mom.v0, grid, u_w);
        
            //Set up Pressure Poisson equation
            SetPressurePoi(grid, mom, solver);
                
            double res;
            int ct;        

            //Solve Pressure Poisson equation
            if(solver.which_solver==1){
              RedBlackSOR(mom.p, mom.ind, res, ct, grid, solver);    
            }
            else
            {
              JORSolver(mom.p, mom.ind, res, ct, grid, solver);
            }
        

            //Print output to screen every 50 time-steps
            if ((grid.itime +1) % 50 ==0) {
                
                #ifdef _OPENMP
                end_time = omp_get_wtime();
                #endif

                std::cout<< "I="<< grid.itime<< " dt="<< mom.dt 
                         << " time="<< mom.time<< std::endl;
                std::cout<< "Pressure residue: "<< res<< std::endl;
                std::cout<< "Pressure iterations: "<< ct<< std::endl;
                std::cout<< "===========================================" 
                         << "==========================================="<< std::endl;
                std::cout<< "Wtime= "<<end_time-start_time<<"\n";                     

                #ifdef _OPENMP
                start_time = omp_get_wtime();   
                #endif      
            }
        
            //Projection: Correct U, V from P field
            projection(grid, mom);
                
            // pply boundary conditions U & V
            UVBoundaryCond(mom.u, mom.v, grid, u_w);        
            
            //Output results right before final time
            if (mom.time + mom.dt > mom.Tend) {

                #ifdef _OPENMP
                prog_end_time = omp_get_wtime();
                #endif
                grid.prog_time = prog_end_time - prog_start_time;

                Output(grid, mom);
                OutputCSV(grid, mom);            
                TimeAnalysis(solver, grid);
            }
        
            //Calculate time step
            timestep(grid, mom); 
        }
        if(solver.t_scheme==2)
        {
            LeapFrog(grid, mom);
        }       
    }

    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Simulation completed successfully!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
