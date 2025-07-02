/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Population__
#define __Population__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "individual.h" 
#include "random.h"

using namespace std;
using namespace arma;

class Population {
private:
  int _M ;          // size of the population
  double _p ;       // exponent to select individual
  int _generation ; // number of generations
  string _data_file; // path to province data file
  int _N_migr; // migration frequency (every N generations)

public:
  Random _rnd;                // Random number generator
  field<Individual> _population; // vector to store individuals

  void initialize();              // Initialize the population
  int selection();                // Select an individual based on fitness
  void sort(int gen, bool output = true);             // Sort the population based on fitness  
  vec compute_fitness();          // Compute fitness for each individual
  void evolve();                  // Evolve the population
  int get_generation() { return _generation; }; // Get the number of generations
  int get_M() { return _M; }; // Get the size of the population
  double get_average_best_half(); // Get the average fitness of the best half of the population
  Individual get_best();          // Get the best individual in the population
  double get_best_fitness();      // Get the fitness of the best individual
  void replace_worst(const Individual& immigrant); // Replace worst with immigrant
  int get_N_migr() const { return _N_migr; } // Migration frequency (every N generations)
};

#endif // __Population__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/