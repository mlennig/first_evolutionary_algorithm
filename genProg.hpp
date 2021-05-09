//
//  genProg.hpp
//  244GeneticProgramming
//  21297899
//  Created by Miriam Lennig on 2016-03-28.
//  Copyright © 2016 Miriam Lennig. All rights reserved.
//

#ifndef genProg_hpp
#define genProg_hpp

#include <iostream>
#include <vector>
#include <list>
using namespace std;

class Model{
// Abstract base class
public:
    virtual unsigned numFeatures() const = 0;                               // Pure virtual function
    virtual double computeFitness(vector<double>&) const = 0;   // Pure virtual function
};

class LangFit : public Model{
// Langermann Model
public:
    LangFit(double); // Randomizes constants
    void print() const;
    double computeFitness(vector<double>&) const;   // Returns fitness
    unsigned numFeatures() const;       // Returns # of features
private:
    double c[5];
    double a[5][2];
    unsigned nFeatures = 2;
    
};

class CubicFit : public Model{
// Cubic Model
public:
    CubicFit(double, int n);
    void print() const;
    double computeFitness(vector<double>&) const;
    unsigned numFeatures() const;
private:
    vector <double> x, y;   // The set of points that the cubic function has to fit
    unsigned nFeatures = 4;
};

class Individual{
// This is a universal individual class which can work with any function that maps an n-dimensional vector to a real number
// It has a pointer to a Model that is shared with many other individuals
public:
    void print() const;                                             // Prints the n-dimensional feature vector that characterizes this individual
    double fitness() const;                                      // Returns real valued fitness
    Individual& procreate(const Individual&) const;   // Performs crossover & mutation to yield a child
    Individual(Model*, double);                                // Randomly assigns values in [-range, range] to the feature vector
    Individual(Model*, vector<double>&);                // Assigns predetermined values to the feature vector
private:
    Model* model;                                                 // Contains the function to optimize; shared with many other individuals
    vector<double> x;                                           // Feature vector
    double storedFitness;                                      // Store fitness to avoid recomputing every time
    void mutate();                                                // Called by procreate
};

class Evolution{
// Runs the main loop of the genetic algorithm for a Model
public:
    Evolution(Model&, double, unsigned, unsigned, bool verbose);
    void print();
private:
    Model& model;   // Reference to an object of a derived class of Model
    unsigned popSize;
    double range;
    list<Individual> population;    // Container for population
    list<Individual> children;        // Container for offspring
    bool stopCriterion();            
    unsigned numIterations = 0;  // Iteration counter
    unsigned maxIterations;
    Individual& chooseParent();   // Tournament selection from population
    void cull();    // Eliminate undesirables
};

#endif /* genProg_hpp */
