#ifndef LDC_HPP
#define LDC_HPP

#include<vector>
#include<string>
#include<cmath>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<iomanip>

#ifdef _OPENMP
#include<omp.h>
#endif

struct Grid //Grid Parameters
{
    int Nx, Ny;                             //Number of gird points
    int is, ie, js, je;                     //start/end index of cell center--Pressure nodes
    int ieu, jev;                           //start/end index of--X & Y velocity nodes
    int Nxt, Nyt, Nt;                       //Total number of points in X direction, Y direction & in total
    int Ng, itime, iOut, iEnd;              //Ghost points, time iteration index
    double Xlen, Ylen, beta, cfl, Max_dt, mesh_ref; //Geometrical & Time stepping parameters, Mesh refinement parameter

    int Nth;
    
    std::vector<double> x, y, dx, dy, dxh, dyh;
    std::vector<double> xh, yh;
    
    int BC[4];
    int Sflag;                              
    bool Semi_implicit, Breakdown_Time;     //=True for Semi-Implicit (Implicit diffusion)
    double solve_time, conv_time, diff_time, prog_time;    //Profiling variables
};

struct Momentum 
{
    std::vector<double> u, v, p;            //u, v, p flattened vectors
    std::vector<double> u0, v0, p0;
    std::vector<double> du, dv, mask;
    
    double dt, time;
    double mu1, Tend;
    int ind[4];
};


struct Solver 
{
    std::vector<double> A1, A2, A3, A4, A5, A6; //Store coefficients of linear equations with inhomogeneous spacing
    int which_solver;                           // 1--SOR (Red/Black), 2--JOR
    int maxit;
};

// Function declarations
void initialize(Grid& grid, Momentum& mom, Solver& solver); //Initialize the grid, fields, and solver

void timestep(Grid& grid, Momentum& mom);                   //Computes time-stepping based on CFL criterion & scheme style

void UVBoundaryCond(std::vector<double>& u, std::vector<double>& v, const Grid& grid); 

void convection(std::vector<double>& du1, std::vector<double>& dv1, 
    //const 
    Grid& grid, const Momentum& mom);                       //Computes descretized convection terms for U & V

void diffusionU(std::vector<double>& du1, int flg, Grid& grid, Momentum& mom, Solver& solver); //Implicit diffusion if required

void diffusionV(std::vector<double>& dv1, int flg, Grid& grid, Momentum& mom, Solver& solver); // " "

void SetPressurePoi(Grid& grid, Momentum& mom, Solver& solver);                                 //Sets up coefficients for pressure equation

void JORSolver(std::vector<double>& p, const int ind[4], double& res, int& count, Grid& grid, const Solver& solver); //JOR--Traditional, Optional

//SOR-Red-Black solver--- Preferred choice for simplicity & convergence---Potential improvements--CG
void RedBlackSOR(std::vector<double>& p, const int ind[4], double& res, int& count, Grid& grid, const Solver& solver);

void projection(Grid& grid, Momentum& mom);                  //Computes corrected velocity
void TimeAnalysis(const Solver& solver, const Grid& grid);   //Manual profiling, OpenMP only!

void Output(const Grid& grid, const Momentum& mom);          //Outputs matlab file

void OutputCSV(const Grid& grid, const Momentum& mom);       //Outputs python-readable csv file

// Interpolation Function
inline double interp(double x1, double x2, double x3, double T1, double T2, 
                     double T3, double x4) 
{
    return (T1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))) + 
           (T2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))) + 
           (T3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)));
}

// 1D to 2D index mapping: (i,j)--->k
inline int ij_k(int i, int j, int Nxt) {
    return i + (j-1)*Nxt;
}

#endif
