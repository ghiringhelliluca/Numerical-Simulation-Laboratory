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
  double _p ;      // exponent to select individual
  int _generation ; // number of generations
  int _sim_type; // type of simulation

public:
  Random _rnd;                // Random number generator
  field <Individual> _population; // vector to store individuals

  void initialize();              // Initialize the population
  int selection();                // Select an individual based
  void sort(int gen);                    // Sort the population based on fitness  
  vec compute_fitness();          // Compute fitness for each individual
  void evolve();                  // Evolve the population
  int get_generation() { return _generation; }; // Get the number of generations
  int get_M() { return _M; }; // Get the size of the population
  double get_average_best_half(); // Get the average fitness of the best half of the population

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
