//
//  genProgTest.cpp
//  244GeneticProgramming
//  21297899
//  Created by Miriam Lennig on 2016-03-28.
//  Copyright © 2016 Miriam Lennig. All rights reserved.
//

#include <ctime>
#include <iostream>
#include "genProg.hpp"
using namespace std;

int main() {
    // Generate and evolve a CubicFit object and a LangFit object. Run through 10 different sets of 30 points for CubicFit
    
    srand(unsigned(time(NULL)));
    
    // The following 5 lines are parameter values to use
    const double modelRange = 1000;     // Uniform in [-modelRange, modelRange]
    const unsigned numPoints = 30;    // # of (x, y) points for cubic MSE
    const double indivRange = 1000;         // Uniform in [-indivRange, indivRange]
    const unsigned populSize = 1000;  // Population size S
    const unsigned maxIter = 1000;    // Maximum # of iterations
    const bool verbose = true;    // Print intermediate results
    const bool silent = false;
    
    // Generate a random Langermann function (select the a & c parameters randomly in [-modelRange, modelRange])
    LangFit modelL(modelRange);
    cout << "\nHere are the random parameters that were generated for the Langermann function:\n";
    modelL.print();
    cout << "\nShowing intermediate results on every iteration, run genetic algorithm to find the (x, y) that minimizes the Langermann function:\n";
    Evolution eL(modelL, indivRange, populSize, maxIter, verbose);
    
    // Generate a random set of points for cubic MSE fit
    CubicFit modelC(modelRange, numPoints);
    cout << "Here is an example set of points for the MSE cubic fit:\n";
    modelC.print();
    cout << "\nShowing intermediate results on every iteration, run genetic algorithm to find the cubic coefficients that have the lowest MSE.\n";
    Evolution eC(modelC, indivRange, populSize, maxIter, verbose);
    
    cout << "\nWithout displaying intermediate values, apply genetic algorithm to minimize MSE for 10 different sets of 30 points for cubic fit:\n";    // The set of (x, y) points are not displayed in order to save space. To print the individual sets, model.print() can be added to the beginning of the loop.
    for(int i = 0; i < 10; i++){
        CubicFit model(modelRange, numPoints);
        cout << "Set #" << i << " of " << numPoints << " random (x, y) points converged after ";
        Evolution e(model, indivRange, populSize, maxIter, silent);
        e.print();
    }
    cout << "\nThe 10 different sets of 30 points for the cubic function have been fit.\n";
}
