//
//  linearTools.cpp
//  project1
//
//  Created by Curtis Rau on 1/28/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>                 // Remove after debugging
#include "linearTools.hpp"

namespace linTools {
    
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    void triDiagMatSolve(int N, double* a, double* b, double* c, double* g, double* x)
    /* Solves a linear system of the form:
     | b0  c0                    | |  x0  |   |  g0  |
     | a0  b1  c1                | |  x1  |   |  g1  |
     |     a1  b2   c2           | |  x2  |   |  g2  |
     |         **   **   **      | |  **  | = |  **  |
     |            aN-3 bN-2 cN-2 | | xN-2 |   | gN-2 |
     |                 aN-2 bN-1 | | xN-1 |   | gN-1 |
     */
    {

        // Forward substitution
        for (int i=1; i<N; i++) {
            b[i] = b[i] - (a[i-1]*c[i-1])/b[i-1];
            g[i] = g[i] - (a[i-1]*g[i-1])/b[i-1];
        }
        
        x[N-1] = g[N-1] / b[N-1];
        
        // Backwards substitution
        for (int i = N-2; i >= 0; i--) {
            g[i] = g[i] - c[i]*g[i+1]/b[i+1];
            x[i] = g[i]/b[i];
        }
    }
    
    
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    void poisonSolver1D(int N, double L, double v0, double vN, double* f, double* v)
    /*
     Solves the Poison Equation of the form:
    
        -v'' = f
     
     with boundary conditions
        v_0 = v0;
        v_N = vN;
    */
    {
        // Boundary Conditions
        v[0]   = v0;
        v[N-1] = vN;
        
        double h2 = L*L / ((N-1)*(N-1));                            // The step size squared;
        
        // Setup source vector "f"
        f[0] = f[0] * h2 + v0;
        for (int i = 1; i < N-3; i++) {
            f[i] = f[i] * h2;
        }
        f[N-3] = f[N-3] * h2 + vN;

        double* b = new double[N-2];
        for (int i = 0; i < N-2; i++) {
            b[i] = (i+2.0)/(i+1.0);
        }
        
        // Forward substitution
        for (int i=1; i<N-2; i++) {
            f[i] = f[i] + f[i-1]/b[i-1];
        }
        
        v[N-2] = f[N-3] / b[N-3];
        
        // Backwards substitution
        for (int i = N-3; i > 0; i--) {
            f[i-1] = f[i-1] + f[i] / b[i];
            v[i]   = f[i-1] / b[i-1];
        }
        
        delete [] b;
    }
    
    
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

    void poisonSolverDirichletBC1D(int N, double L, double* f, double* v)
    /*
         Solves the Poison Equation of the form:
         
         -v'' = f
         
         with boundary conditions
         v_0 = 0;
         v_N = 0;
     */
    {
        // Boundary Conditions
        v[0]   = 0.0;
        v[N-1] = 0.0;
        
        double h2 = L*L / ((N-1)*(N-1));                            // The step size squared;
        
        // Forward substitution
        for (int i=1; i<N-2; i++) {
            f[i] = f[i] + f[i-1] * i / (i+1.0);
        }
        
        v[N-2] = f[N-3] * h2 * (N-2.0) / (N-1.0);
        
        // Backwards substitution
        for (int i = N-3; i > 0; i--) {
            f[i-1] = f[i-1] + f[i] * (i+1.0) / (i+2.0);
            v[i]   = f[i-1] * h2 * i / (i+1.0);
        }
    }
}