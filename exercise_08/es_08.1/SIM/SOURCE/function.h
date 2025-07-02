/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Function__
#define __Function__

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include "random.h"

using namespace std;

class QuantumSimulation {
private:
    // Configuration parameters
    double temperature;
    double delta;
    int nblocks;
    int nsteps;
    int equilibration_blocks;
    double mu;
    double sigma;
    
    // Simulation state
    double x_position;
    int total_accepted;
    int total_attempted;
    
    // Block averaging variables
    vector<double> block_energies;  // Store energy for each block
    double sum_prog;               // Progressive sum
    double sum2_prog;              // Progressive sum of squares
    
    // Private methods
    double psi_component(double x, double mu_val, double sigma_val) const;
    void reset_accumulators();
    void write_output_headers() const;
    double calculate_error(double sum_avg, double sum2_avg, int n) const;
    
public:
    // Initialization methods
    void initialize_from_file();
    Random setup_random_generator() const;
    
    // Physics methods
    double calculate_energy(double x_val, double mu_val, double sigma_val) const;
    void metropolis_step(Random &rnd, int &accepted);
    
    // Block averaging methods
    void reset_for_new_block();
    void add_block_energy(double energy);
    void calculate_progressive_averages(int block_num, double &prog_avg, double &prog_error);
    
    // Simulation methods
    tuple<double, double, double> run_simulation(Random &rnd);
    void equilibrate(Random &rnd);
    
    // Output methods
    void save_results(int block, double block_energy, double prog_avg, double prog_error, double acceptance_rate) const;
    void save_final_summary(double final_energy, double final_error, double final_acceptance) const;
    void clear_output_files() const;
    
    // Getters
    double get_temperature() const { return temperature; }
    double get_delta() const { return delta; }
    double get_mu() const { return mu; }
    double get_sigma() const { return sigma; }
    double get_position() const { return x_position; }
    int get_nblocks() const { return nblocks; }
    int get_nsteps() const { return nsteps; }
    int get_equilibration_blocks() const { return equilibration_blocks; }
    double get_acceptance_rate() const { 
        return (total_attempted > 0) ? (double)total_accepted / total_attempted : 0.0; 
    }
    
    // Setters
    void set_position(double x) { x_position = x; }
    void set_parameters(double mu_val, double sigma_val) { mu = mu_val; sigma = sigma_val; }
    
    void print_parameters() const;
};

#endif