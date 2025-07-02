/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef HYDROGEN_ATOM_H
#define HYDROGEN_ATOM_H

#include <vector>
#include <string>
#include <tuple>
#include <array>
# include <cmath>
#include "random.h"

using namespace std;

struct BlockResult {
    vector<double> mean;
    vector<double> error;
};

struct SimulationParams {
    int n_eq;           
    int M;              
    int N;              
    double sigma_100;   
    double delta_100;   
    double sigma_210;   
    double delta_210;   
    double x_start;     
    double y_start;     
    double z_start_100; 
    double z_start_210; 
};

inline double distance_from_origin_squared(double x, double y, double z) {
    return x * x + y * y + z * z;
}

inline double distance_from_origin(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

inline double psi100_density(double x, double y, double z) {
    return exp(-2.0 * sqrt(x * x + y * y + z * z));
}

inline double psi210_density(double x, double y, double z) {
    return z * z * exp(-sqrt(x * x + y * y + z * z));
}

double statistical_error(double average, double average_squared, int n);
BlockResult blocking_method(const vector<double>& data, int N, int L);

tuple<vector<double>, vector<double>> metropolis_sampling_optimized(
    const SimulationParams& params,
    Random& rnd, 
    double (*target_density)(double, double, double),
    int steps, 
    double sigma, 
    double delta, 
    double x_start, 
    double y_start, 
    double z_start, 
    const string& output_file_unif, 
    const string& output_file_gauss,
    ostream& out,
    bool save_coordinates = false
);

Random initialize_random_generator();
SimulationParams read_parameters(const string& filename);
void save_equilibration_data(const vector<double>& data_unif, 
                           const vector<double>& data_gauss, 
                           const string& filename);
void save_block_results(const vector<double>& means,
                       const vector<double>& errors,
                       int L, const string& filename);

#endif // HYDROGEN_ATOM_H