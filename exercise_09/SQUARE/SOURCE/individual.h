/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Individual__
#define __Individual__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Individual {

private:
  int _ndim ; // Dimensionality of the system depends on the problem circle or square
  double _r = 1.0 ; // Radius of the circuference
  double _l = 1.0 ; // Length of the square
  int _ncity = 34; // Number of cities

public: 

  field <vec> _individual;  // field of vectors for the individual (group of cities)

  void initialize(Random &_rnd, int sim_type); // Initialize the individual
  bool check();                          // Check if the individual is valid
  double loss(int sim_type);                         // Calculate the loss function (distance from first to first city)
  double loss_circle();                         // Calculate the loss function for the circle
  double loss_square();                         // Calculate the loss function for the square 
  void permutation(Random &_rnd);                    // Random permutation of the cities
  void permutation_pair(Random &_rnd); // Random permutation of a pair of cities
  void mutation_swap_blocks(Random &_rnd);           // Swap two blocks of cities
  void mutation_shift(Random &_rnd);       // shift m city +n positions
  void mutation_inversion(Random &_rnd);   // Inversion of m cities
  pair<Individual, Individual> crossover_pair(const Individual& partner, Random &_rnd); // Crossover between two individuals
  void print_path(); // Print the path of the individual
  void save_path(ofstream& path);       // save path of cities 
  void save_individual();             // save dictionary of cities with ID and position
  int get_ncity() { return _ncity; }; // Get the number of cities

};

#endif // __Individual__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
