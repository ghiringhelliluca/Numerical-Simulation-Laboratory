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
  int _ndim = 2; // Always 2D for geographic coordinates (lon, lat)
  int _ncity; // Number of provinces (read from file)
  field<vec> _base_individual; // Base coordinates of provinces

public: 
  // field<vec> _base_individual; // Base coordinates of provinces
  field<vec> _individual;  // field of vectors for the individual (group of provinces)

  void initialize(Random &_rnd, const string& filename); // Initialize from file
  void initialize_from_base(Random &_rnd); // Initialize from base coordinates
  bool check();                          // Check if the individual is valid
  double loss();                         // Calculate the loss function (total distance)
  void permutation(Random &_rnd);        // Random permutation of the provinces
  void mutation_swap_blocks(Random &_rnd);           // Swap two blocks of provinces
  void mutation_shift(Random &_rnd);     // shift m provinces +n positions
  void mutation_inversion(Random &_rnd); // Inversion of m provinces
  pair<Individual, Individual> crossover_pair(const Individual& partner, Random &_rnd); // Crossover between two individuals
  void print_path(); // Print the path of the individual
  void save_path(ofstream& path);       // save path of provinces 
  void save_individual();             // save dictionary of provinces with ID and position
  int get_ncity() { return _ncity; }; // Get the number of provinces
  
  // MPI communication methods
  vector<int> encode() const;
  void decode(const vector<int>& data);
  void set_base_individual(const field<vec>& base) { _base_individual = base; }
  field<vec> get_base_individual() const { return _base_individual; }
  
private:
  double GetXCoordinate(int province_id) const;
  double GetYCoordinate(int province_id) const;
  double haversine_distance(double lon1, double lat1, double lon2, double lat2) const;
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