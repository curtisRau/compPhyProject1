//  main.cpp
//  computationalPhysics1
//
//  Created by Curtis Rau on 1/15/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>             // For working with terminal output: std::cout?
#include <fstream>              // for working with files.
#include <new>                  // Is this necessary?
#include <math.h>               /* exp squt */
#include "linearTools.hpp"

// Declare Functions Here
double source (double x) {
    return 10.0 * exp (-10.0 * x);
};


int main(int argc, const char* argv[]) {
    // The first argument should be N -- the number of points, including endpoints.
    // Using the flag "-plot" plots the results in the terminal.
    
    const clock_t;
    
    bool PLOT = false;
    bool SAVE = false;
    const char* filename;
    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i], "-plot")==0) {PLOT = true;}
        
        if (strcmp(argv[i], "-save")==0) {                  // If the "-save" flag is passed to main()
            filename = argv[i+1];                           // The filename should be the next argument passed to main()
        }
    }
    
    // Check if filename is valid
    if (!filename) {                                // If the filename = NULL...
        std::cout << "No path to file specified.  Was expecting '-save path_to_file' /r";
    } else {                                        // If the filename != NULL...
        std::ofstream outputFile;
        outputFile.open(filename);                  // Try opening the file
        if (!outputFile.is_open()) {                // If the file doesn't open...
            std::cout << "Could not open file.  Was expecting '-save path_to_file' /r";
        } else {                                    // If the file can be opened...
            SAVE = true;
        }
    }
    
    const unsigned char METHOD = 3;
    
        // BEGINNING CODE COMMON TO ALL METHODS (DO NOT COMMENT OUT)
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    const int N = 100000000;//std::stoi(argv[1]);          // takes the first argument, converts it to an integer, and stores it as "n".
                                                    // N is the number of points, including end points
    double L = 1.0;                                 // The length of the sample;
    
// Setup solution vector "u"
    double* u = new double[N];

    
        // CODE FOR GENERAL SOLVER
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    if (METHOD == 1) {
        // Setup source vector "f"
        
        double h = L / (N+1);                           // The step size;
        
        double* f = new double[N];
        f[0] = -source(1*h) * (h*h) - 0;
        for (int i = 1; i < N-1; i++) {
            f[i] = -source((i+1) * h) * (h*h);
        }
        f[N-1] = -source(N * h) * (h*h) - 0;
        
        clock_t begin_time = std::clock();
        
        // Create the vectors representing the "A" matrix
        double* a = new double[N-1];                              // Could replace array with constant, but leaving it general for now
        for (int i = 0; i < N-1; i++) {
            a[i] = 1.0;
        }
        
        double* b = new double[N];                              // There are N+1 points, but remove end points, so N-1 points to solve for.
        for (int i = 0; i < N; i++) {
            b[i] = -2.0;
        }
        
        double* c = new double[N-1];                              // Could replace array with constant, but leaving it general for now
        for (int i = 0; i < N-1; i++) {
            c[i] = 1.0;
        }
        
        linTools::triDiagMatSolve(N, a, b, c, f, u);
        std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
        
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] f;
    }
    
    
        // CODE FOR "RACECAR"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    if (METHOD == 2) {
        double h = L / (N+1);                           // The step size;
        // Setup source vector "f"
        double* f = new double[N-2];
        f[0] = source(1*h);
        for (int i = 1; i < N-3; i++) {
            f[i] = source((i+1) * h);
        }
        f[N-3] = source((N-2) * h);
        
        clock_t begin_time = std::clock();
        linTools::poisonSolver1D(N, L, 0.0, 0.0, f, u);
        std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
        
        delete [] f;
    }

    
        // CODE FOR "SUPER RACECAR"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    if (METHOD == 3) {
        double h = L / (N+1);                           // The step size;
        // Setup source vector "f"
        double* f = new double[N-2];
        for (int i = 0; i < N-2; i++) {
            f[i] = source((i+1) * h);
        }
        
        clock_t begin_time = std::clock();
        linTools::poisonSolverDirichletBC1D(N, L, f, u);
        std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\r";
        
        delete [] f;
    }
    
    
        // CODE FOR RUDIMENTARY PLOTTER
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    if (false) {              // If the flag "-plot" is passed to main() and is in the argument vector.

        for (int i = 0; i < N; i = i+1) {           // Plot the result.  The x-scale is controlled by "i = i + 1"
            double num = 1000 * u[i];               // The y-scale is controlled by "1000 * u[i]"
            for (int i = -1; i < num; i++) {
                std::cout << "-";
            }
            std::cout << "\r";
        }
    }
    
    // CODE FOR WRITING OUTPUT TO A .CSV FILE
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    if (SAVE) {              // If the flag "-save" is passed to main() and is in the argument vector.
        std::ofstream outputFile;
        outputFile.open(filename, std::ios::out | std::ios::trunc);         // Open a file for output and overwrite current content if it exists.
            
        if (outputFile.is_open()) {                                         // If the file is open...
            outputFile << u[0];
            for (int i=1; i<N; i++) {
                outputFile << "," << u[i];
            }
        } else {
            std::cout << "File '" << filename << "' did not open /r";
        }
            
        outputFile.close();
    }
    
    
    // GENERAL CODE
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    delete [] u;
    
    return 0;
}