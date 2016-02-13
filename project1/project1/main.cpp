//  main.cpp
//  computationalPhysics1
//
//  Created by Curtis Rau on 1/15/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>
#include <new>          // Is this necessary?
#include <math.h>       /* exp squt */
#include "linearTools.hpp"

// Declare Functions Here
double source (double x) {
    return 10.0 * exp (-10.0 * x);
};


int main(int argc, const char* argv[]) {
    // Using the flag "-plot" plots the results in the terminal.
    
        // BEGINNING CODE COMMON TO BOTH METHODS (DO NOT COMMENT OUT)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
    const int N = std::stoi(argv[1]);          // takes the first argument, converts it to an integer, and stores it as "n".
                                                    // N is the number of points, including end points
    double L = 1.0;                                 // The length of the sample;
    double h = L / (N-1);                           // The step size;
    
// Setup solution vector "u"
    double* u = new double[N];

    
        // CODE FOR GENERAL SOLVER
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
////-------------------------------------------------------------------------------------------------------------------------------------------------------------// Setup source vector "f"
//    double* f = new double[N-2];
//    f[0] = -source(1*h) * (h*h) - 0;
//    for (int i = 1; i < N-3; i++) {
//        f[i] = -source((i+1) * h) * (h*h);
//    }
//    f[N-3] = -source((N-2) * h) * (h*h) - 0;
//    
//    const clock_t begin_time = std::clock();
//    
//// Boundary conditions:
//    u[0] = 0.0;                   //The solution at x=0
//    u[N-1] = 0.0;                 //The solution at x=N*h
//    
//// Create the vectors representing the "A" matrix
//    double* a = new double[N-3];                              // Could replace array with constant, but leaving it general for now
//    for (int i = 0; i < N-3; i++) {
//        a[i] = 1.0;
//    }
//    
//    double* b = new double[N-2];                              // There are N+1 points, but remove end points, so N-1 points to solve for.
//    for (int i = 0; i < N-2; i++) {
//        b[i] = -2.0;
//    }
//    
//    double* c = new double[N-3];                              // Could replace array with constant, but leaving it general for now
//    for (int i = 0; i < N-3; i++) {
//        c[i] = 1.0;
//    }
//
//    linTools::triDiagMatSolve(N, a, b, c, f, u);
//    std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
//
//    delete [] a;
//    delete [] b;
//    delete [] c;
//    delete [] f;
    
    
        // CODE FOR "RACECAR" (300,5e-5; 1000,8.5e-5)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//// Setup source vector "f"
//    double* f = new double[N-2];
//    f[0] = source(1*h);
//    for (int i = 1; i < N-3; i++) {
//        f[i] = source((i+1) * h);
//    }
//    f[N-3] = source((N-2) * h);
//    
//    const clock_t begin_time = std::clock();
//    linTools::poisonSolver1D(N, L, 0.0, 0.0, f, u);
//    std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
//    
//    delete [] f;

    
        // CODE FOR "SUPER RACECAR" (300,4.8e-5; 1000,8.3e-5; 4000; 0.000233)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Setup source vector "f"
    double* f = new double[N-2];
    for (int i = 0; i < N-2; i++) {
        f[i] = source((i+1) * h);
    }
    
    const clock_t begin_time = std::clock();
    linTools::poisonSolverDirichletBC1D(N, L, f, u);
    std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
    
    delete [] f;
    
        // CODE FOR RUDIMENTARY PLOTTER
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
    if (strncmp(argv[0], "-plot", 5)) {
        for (int i = 0; i < N; i = i+1) {
            double num = 1000 * u[i];
            //std:: cout << "Here's that val:" << num;                    // For debugging purposes
            for (int i = -1; i < num; i++) {
                std::cout << "-";
            }
            std::cout << "\r";
        }
    }
    
    
    delete [] u;
    
    return 0;
}