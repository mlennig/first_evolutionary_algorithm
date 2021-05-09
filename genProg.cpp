//
//  genProg.cpp
//  244GeneticProgramming
//  21297899
//  Created by Miriam Lennig on 2016-03-28.
//  Copyright © 2016 Miriam Lennig. All rights reserved.
//

#include "genProg.hpp"
#include <cmath>
#include <vector>
using namespace std;

const unsigned long maxRand = 2147483647;

Individual::Individual(Model* m, double range) : model(m){
    // Fill feature vector x with random doubles in [-range, range] and compute fitness
    unsigned nfeat = model->numFeatures();
    for(int i = 0; i < nfeat; i++)
        x.push_back(2 * range * rand() / maxRand - range);
    storedFitness = model->computeFitness(x);
}

Individual::Individual(Model* m, vector<double>& xvalues) : model(m){
    // Fill feature vector x with values from vector passed to it and compute fitness
    unsigned nfeat = model->numFeatures();
    for(int i = 0; i < nfeat; i++)
        x.push_back(xvalues[i]);
    storedFitness = model->computeFitness(x);
}

void Individual::print() const{
    // Prints the n-dimensional feature vector that characterizes this individual
    unsigned long nfeat = x.size();
    for(int i = 0; i < nfeat; i ++)
        cout << x[i] << "\t\t";
    cout << endl;
}

double Individual::fitness() const{
    return storedFitness;
}

void Individual::mutate(){
    // Mutates the individual in place using real valued mutation
    // http://www.geatbx.com/docu/algindex-03.html#P550_28854
    unsigned nfeat = model->numFeatures();  // Get number of features
    double mutateProb = 1./nfeat;   // Probability of mutation is 1/nfeat
    // For explanation of domain, r, k, u, s, a, see the web reference above
    double domain = 100;
    double r = 0.1 * domain; //  range = r in [0.1, 10^-6]
    unsigned k = 15;  // precision = k in {4, ... ,20}
    double u, s, a;
    for(int i = 0; i < nfeat; i++)
        if(rand() < mutateProb * maxRand){
            u = double(rand()) / maxRand;
            s = 2. * rand() / maxRand - 1;
            a = pow(2, -u * k);
            x[i] +=  s * r * a;
        }
    storedFitness = model->computeFitness(x);   // Recompute fitness of mutated child
}

Individual& Individual::procreate(const Individual& mate) const {
    // Create a new Individual using intermediate recombination
    // http://www.geatbx.com/docu/algindex-03.html#P550_28854
    unsigned nfeat = model->numFeatures();
    double d = 0.25;
    double a;
    vector<double> z;
    
    for(int i = 0; i < nfeat; i++){
        a = (1 + 2 * d) * rand() / maxRand - d;    // a in [-d, 1 + d]
        z.push_back(a * x[i] + (1 - a) * mate.x[i]);
    }
    
    Individual* pChild = new Individual (model, z);
    pChild->mutate();
    return *pChild;
}

LangFit::LangFit(double range){
    // Assign random values to matrices c, a
    for(int i = 0; i < 5; i++){
        c[i] = 2 * range *double(rand())/maxRand - range;
        for(int j = 0; j < 2; j++)
            a[i][j] = 2 * range *double(rand())/maxRand - range;
    }
}

unsigned LangFit::numFeatures() const {
    return nFeatures;
}

void LangFit::print() const {
    //  Print constants of the LangFit object
    cout << "Constants of the LangFit model object: " << endl;
    cout << "The c matrix is:\n";
    for (int i = 0; i < 5; i++)
        cout << c[i] << "\t\t";
    cout << endl << endl;
    cout << "The a matrix is:\n";
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 2; j++)
            cout << a[i][j] << "\t\t";
        cout << endl;
    }
    cout << endl;
}

double LangFit::computeFitness(vector<double>& x) const{
    // Generates the value of the Langermann function
    const double pi = 3.1415926535;
    const double e = 2.71828182846;
    double f = 0;
    for(int i = 0; i < 5; i++){
        double secondSum = 0;
        for(int j = 0; j < 2; j++)
            secondSum += pow((x[j] - a[i][j]), 2);
        double thirdSum = 0;
        for(int j = 0; j < 2; j++)
            thirdSum += pow((x[j] - a[i][j]), 2);
        f += c[i] * pow(e, -secondSum / pi) * cos(pi * thirdSum);
    }
    return f;
}

CubicFit::CubicFit(double range, int n){
    //  Initialize n random (x, y) points
    for(int i = 0; i < n; i++){
        x.push_back(2 * range * rand() / maxRand - range);
        y.push_back(2 * range * rand() / maxRand - range);
    }
}

void CubicFit::print() const{
    cout << "Printing (x, y) coordinate pairs for MSE cubic fit:\n";
    unsigned long n = x.size();
    for (int i = 0; i < n; i++)
        cout << "(" << x[i] << ", " << y[i] << ")\n";
    cout << endl;
}

double CubicFit::computeFitness(vector<double>& coeffs) const {
    //  Fitness function is the MSE fit to the point cloud
    double mse = 0;
    unsigned long n = x.size();
    double yp;
    for(int i = 0; i < n; i ++){
        yp = coeffs[0] + coeffs[1]*x[i] + coeffs[2]*pow(x[i], 2) + coeffs[3]*pow(x[i], 3);
        mse += pow(yp - y[i], 2);
    }
    mse /= n;
    return mse;
}

unsigned CubicFit::numFeatures() const {
    return nFeatures;
}

bool cmp(const Individual& i1, const Individual& i2){
    return i1.fitness() <= i2.fitness();
}


Individual& Evolution::chooseParent(){
    // Use tournament selection to choose one parent from the population pool
    list<Individual>::iterator it1 = population.begin();
    list<Individual>::iterator it2 = population.begin();
    list<Individual>::iterator it3 = population.begin();
    
    //  Select 3 random candidates from population
    unsigned long index = rand() % popSize;    // Index number of item in population list
    advance(it1, index);    // Increments it1 by index
    index = rand() % popSize;   // Re-randomize index number
    advance(it2, index);
    index = rand() % popSize;
    advance(it3, index);
    
    //  Perform tournament selection to choose 1 parent
    if ((it1->fitness() < it2->fitness()) && (it1->fitness() < it3->fitness()))
        return *it1;
    else if ((it2->fitness() < it1->fitness()) && (it2->fitness() < it3->fitness()))
        return *it2;
    else
        return *it3;
}

void Evolution::cull(){
    // Delete all but the first popSize individuals from population
    list<Individual>::iterator iter = population.begin();
    advance(iter, popSize);     // iter is now pointing to where we want the end to be, to return to original population size
    population.erase(iter, population.end());   // Erase all individuals after popSize
}

bool Evolution::stopCriterion(){
    // Stop when the ratio of the fitness of worst and best individuals is equal to 1.0 in single precision
    // or if maximum # of iterations have been reached using a ratio
    float ratio;
    if(numIterations > maxIterations){
        cout << "Max # of iterations has been reached\n\n";
        return true;
    }
    double bestFitness = population.front().fitness();
    double worstFitness = population.back().fitness();
    
    if (worstFitness != 0)
        ratio = bestFitness / worstFitness;
    else if (bestFitness != 0)
        ratio = worstFitness / bestFitness;
    else if (worstFitness == bestFitness)
        ratio = 1;
    else
        ratio = 10;
    
    if (ratio == 1) {
        return true;
    }
    return false;
}

void Evolution::print(){
    // Displays intermediate results
    cout << "Number of iterations: " << numIterations << endl;
    cout << "The best individual has fitness: " << population.front().fitness() << " and the following features:\n";
    population.front().print();
    cout << "The worst individual has fitness: " << population.back().fitness() << " and the following features:\n";
    population.back().print();
    cout << endl;
}

Evolution::Evolution(Model& m, double r, unsigned pSize, unsigned maxIter, bool verbose) : model(m), maxIterations(maxIter), popSize(pSize), range(r) {
    // Runs main genetic algorithm
    // Create popSize random individuals & put them into the population
    for(int i = 0; i < popSize; i++)
        population.push_back(Individual(&model, range));
    
    
    population.sort(cmp);   // Sort from fittest to least fit
    
    if(verbose)
        print();
    
    while(!stopCriterion()){    // Main loop
        // Procreation loop
        while(children.size() < 10 * popSize){
            // Choose 2 parents from the population by tournament selection
            Individual& parent1 = chooseParent();
            Individual& parent2 = chooseParent();
            
            // Make 2 children from chosen parents, push them into the children pool
            children.push_back(parent1.procreate(parent2));
            children.push_back(parent1.procreate(parent2));
        }
        
        population.merge(children, cmp);
        population.sort(cmp);   // Sort from lowest fitness to highest fitness
        cull();     // Erase all individuals after popSize one
        numIterations++;    // Increment iteration counter
        if(verbose)
            print();    // Print intermediate results
        
    }
    
    
}

