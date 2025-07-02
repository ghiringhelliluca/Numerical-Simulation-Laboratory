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
#include "random.h"

using namespace std;

class QuantumSimulation {
private:
    // Configuration parameters
    double temperature;
    double delta;
    int nblocks;
    int nblocks_final;
    int nsteps;
    int equilibration_blocks;
    double mu;
    double sigma;
    double alpha;
    
    // Simulation state
    double x_position;
    int total_accepted;
    int total_attempted;
    
    // Block averaging variables
    double block_average;
    double global_average;
    double global_average_squared;
    
    // SA parameters
    double step_mu;
    double step_sigma;
    int max_sa_steps;
    
    double best_energy;
    double best_mu;
    double best_sigma;
    double best_error;
    double best_acceptance;
    
    int accepted_moves;
    int total_moves;

    int steps_without_improvement;
    int max_steps_without_improvement;
    
    double psi_component(double x, double mu_val, double sigma_val) const;
    void reset_block_accumulators();
    void write_output_headers() const;
    void write_sa_header() const;
    double calculate_error(double acc, double acc2, int blocks) const;
    
public:

    void initialize_from_file();
    void initialize_best_parameters();
    
    double calculate_energy(double x_val, double mu_val, double sigma_val) const;
    void metropolis_step(Random &rnd, int &accepted);
    void measure_energy();
    
    double process_block_averages(int block_number);
    void reset_for_new_temperature();
    
    tuple<double, double, double> run_energy_calculation(Random &rnd, int sa_iteration = 0);
    tuple<double, double, double> run_energy_calculation_no_save(Random &rnd, int sa_iteration = 0);
    tuple<double, double, double> run_final_energy_calculation(Random &rnd);
    tuple<double, double> propose_new_parameters(Random &rnd) const;
    bool metropolis_accept_sa(double E_old, double E_new, Random &rnd) const;
    void update_parameters(double new_mu, double new_sigma);
    void cool_temperature();
    
    void update_best_parameters(double energy, double error, double acceptance);
    void restore_best_parameters();

    void save_trajectory_step(double energy, double error, double acceptance_rate, double best_energy) const;
    void save_position() const;
    void clear_position_file() const;
    
    double get_temperature() const { return temperature; }
    double get_mu() const { return mu; }
    double get_sigma() const { return sigma; }
    double get_position() const { return x_position; }
    double get_best_energy() const { return best_energy; }
    double get_best_mu() const { return best_mu; }
    double get_best_sigma() const { return best_sigma; }
    double get_best_error() const { return best_error; }
    double get_best_acceptance() const { return best_acceptance; }
    int get_sa_iteration_count() const { return max_sa_steps; }
    int get_steps_without_improvement() const { return steps_without_improvement; }
    int get_max_steps_without_improvement() const { return max_steps_without_improvement; }

    void print_current_state() const;
    void print_best_state() const;
    Random setup_random_generator() const;
};

#endif