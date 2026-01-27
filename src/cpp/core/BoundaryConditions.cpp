#include "LDC.hpp"

// Boundary conditions for U and V velocities - OPENMP OPTIMIZED
void UVBoundaryCond(std::vector<double>& u, std::vector<double>& v, const Grid& grid) {
    const int Nxt= grid.Nxt;
    const int is= grid.is;
    const int ie= grid.ie;
    const int js= grid.js;
    const int je= grid.je;
    
    // Left and right boundary conditions - PARALLELIZED
    // Independent updates for each j-level
    #pragma omp parallel for schedule(static)
    for (int j= js; j<= je; j++) {
        // Left boundary
        u[ij_k(is-1, j, Nxt)]= 0.0;
        u[ij_k(is-2, j, Nxt)]= -u[ij_k(is, j, Nxt)];
        v[ij_k(is-1, j, Nxt)]= -v[ij_k(is, j, Nxt)];
        v[ij_k(is-2, j, Nxt)]= -v[ij_k(is+1, j, Nxt)];
        
        // Right boundary
        u[ij_k(ie, j, Nxt)]= 0.0;
        u[ij_k(ie+1, j, Nxt)]= -u[ij_k(ie-1, j, Nxt)];
        v[ij_k(ie+1, j, Nxt)]= -v[ij_k(ie, j, Nxt)];
        v[ij_k(ie+2, j, Nxt)]= -v[ij_k(ie-1, j, Nxt)];
    }
    
    // Bottom and top boundary conditions - PARALLELIZED
    // Independent updates for each i-level
    #pragma omp parallel for schedule(static)
    for (int i= is; i<= ie; i++) {
        // Bottom boundary
        v[ij_k(i, js-1, Nxt)]= 0.0;
        v[ij_k(i, js-2, Nxt)]= -v[ij_k(i, js, Nxt)];
        u[ij_k(i, js-1, Nxt)]= -u[ij_k(i, js, Nxt)];
        u[ij_k(i, js-2, Nxt)]= -u[ij_k(i, js+1, Nxt)];
        
        // Top wall - lid-driven cavity with velocity= 1.0
        v[ij_k(i, je, Nxt)]= 0.0;
        v[ij_k(i, je+1, Nxt)]= -v[ij_k(i, je-1, Nxt)];
        u[ij_k(i, je+1, Nxt)]= 2.0*1.0 - u[ij_k(i, je, Nxt)];
        u[ij_k(i, je+2, Nxt)]= 2.0*1.0 - u[ij_k(i, je-1, Nxt)];
    }
}
