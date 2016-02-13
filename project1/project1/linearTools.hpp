//
//  linearTools.hpp
//  project1
//
//  Created by Curtis Rau on 1/28/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#ifndef linearTools_hpp
#define linearTools_hpp

#include <stdio.h>

#endif /* linearTools_hpp */

// Declare the public class "linTools" which will have functions for solving linear equations
namespace linTools {
    void triDiagMatSolve (int, double*, double*, double*, double*, double*);
    void poisonSolver1D(int N, double L, double v0, double vN, double* f, double* v);
    void poisonSolverDirichletBC1D(int N, double L, double* f, double* v);
};


