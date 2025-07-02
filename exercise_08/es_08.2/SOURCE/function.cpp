/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <iomanip>
#include <limits>
#include "function.h"

using namespace std;

// Initialize from input file
void QuantumSimulation::initialize_from_file() {
    write_output_headers();
    write_sa_header();
    
    ifstream input("../INPUT/input.dat");  
    ofstream output("../OUTPUT/output.dat");

    string property;
    while (!input.eof()) {
        input >> property;
        if (property == "TEMP") {
            input >> temperature;
            output << "TEMPERATURE= " << temperature << endl;
        } else if (property == "DELTA") {
            input >> delta;
            output << "DELTA= " << delta << endl;
        } else if (property == "NBLOCKS") {
            input >> nblocks;
            output << "NBLOCKS= " << nblocks << endl;
        } else if (property == "NSTEPS") {
            input >> nsteps;
            output << "NSTEPS= " << nsteps << endl;
        } else if (property == "NEQUILIBRATION") {
            input >> equilibration_blocks;
            output << "NEQUILIBRATION= " << equilibration_blocks << endl;
        } else if (property == "MU") {
            input >> mu;
            output << "MU= " << mu << endl; 
        } else if (property == "SIGMA") {
            input >> sigma;
            output << "SIGMA= " << sigma << endl;
        } else if (property == "ALPHA") {
            input >> alpha;
            output << "ALPHA= " << alpha << endl;
        } else if (property == "STEP_MU") {
            input >> step_mu;
            output << "STEP_MU= " << step_mu << endl;
        } else if (property == "STEP_SIGMA") {
            input >> step_sigma;
            output << "STEP_SIGMA= " << step_sigma << endl;
        } else if (property == "NBLOCKS_FINAL") {
            input >> nblocks_final;
            output << "NBLOCKS_FINAL= " << nblocks_final << endl;
        } else if (property == "MAX_SA_STEPS") {
            input >> max_sa_steps;
            output << "MAX_SA_STEPS= " << max_sa_steps << endl;
        }  else if (property == "MAX_STEPS_NO_IMPROVE") {
            input >> max_steps_without_improvement;
            output << "MAX_STEPS_NO_IMPROVE= " << max_steps_without_improvement << endl;
        } else if (property == "ENDINPUT") {
            output << "Reading input completed!" << endl;
            break;
        } else if (!property.empty()) {
            cerr << "WARNING: Unknown input property: " << property << endl;
        }
    }
    
    input.close();
    output.close();
    return;
}

void QuantumSimulation::initialize_best_parameters() {
    best_energy = numeric_limits<double>::max();
    best_mu = mu;
    best_sigma = sigma;
    best_error = 0.0;
    best_acceptance = 0.0;
    steps_without_improvement = 0;

}

double QuantumSimulation::psi_component(double x, double mu_val, double sigma_val) const {
    return exp(-pow((x - mu_val), 2) / (2 * pow(sigma_val, 2)));
}

void QuantumSimulation::reset_block_accumulators() {
    block_average = 0.0;
    global_average = 0.0;
    global_average_squared = 0.0;
}

void QuantumSimulation::write_output_headers() const {
    ofstream acceptance_file("../OUTPUT/acceptance.dat");
    acceptance_file << "#   N_BLOCK:  ACCEPTANCE:" << endl;
    acceptance_file.close();
    
    ofstream energy_file("../OUTPUT/energy.dat");
    energy_file << "#     BLOCK:   ACTUAL_EN:     EN_AVE:       ERROR:" << endl;
    energy_file.close();
}

void QuantumSimulation::write_sa_header() const {
    ofstream sa_file("../OUTPUT/sa_trajectory.dat");
    sa_file << "T:    MU:       SIGMA:     ENERGY:      EN_ERR:    MEAN_ACCEPTANCE:   BEST_ENERGY:"  << endl;
    sa_file.close();
}

double QuantumSimulation::calculate_error(double acc, double acc2, int blocks) const {
    if (blocks <= 1) return 0.0;
    return sqrt(fabs(acc2/double(blocks) - pow(acc/double(blocks), 2))/double(blocks));
}

double QuantumSimulation::calculate_energy(double x, double mu_val, double sigma_val) const {
    double sigma2 = sigma_val * sigma_val;
    double sigma4 = sigma2 * sigma2;
    
    double psi1 = psi_component(x, mu_val, sigma_val);
    double psi2 = psi_component(x, -mu_val, sigma_val);
    double psi_total = psi1 + psi2;
    
    double d2psi1 = ((pow(x - mu_val, 2) / sigma4) - (1.0 / sigma2)) * psi1;
    double d2psi2 = ((pow(x + mu_val, 2) / sigma4) - (1.0 / sigma2)) * psi2;
    
    double kinetic = -0.5 * (d2psi1 + d2psi2) / psi_total;
    double potential = pow(x, 4) - 5.0/2.0 * pow(x, 2);
    
    return kinetic + potential;
}

void QuantumSimulation::metropolis_step(Random &rnd, int &accepted) {
    double x_new = x_position + rnd.Rannyu(-delta, delta);
    double p_old = pow(psi_component(x_position, mu, sigma) + psi_component(x_position, -mu, sigma), 2);
    double p_new = pow(psi_component(x_new, mu, sigma) + psi_component(x_new, -mu, sigma), 2);
    double acceptance_probability = min(1.0, p_new / p_old);
    
    if (rnd.Rannyu() < acceptance_probability) {
        x_position = x_new;
        accepted++;
    }
}

void QuantumSimulation::measure_energy() {
    double energy = calculate_energy(x_position, mu, sigma);
    block_average += energy;
}

// Block averaging methods
double QuantumSimulation::process_block_averages(int block_number) {
    double average = block_average / double(nsteps);
    global_average += average;
    global_average_squared += average * average;
    
    // Write energy data
    ofstream energy_file("../OUTPUT/energy.dat", ios::app);
    energy_file << setw(12) << block_number
                << setw(12) << average
                << setw(12) << global_average/double(block_number)
                << setw(12) << calculate_error(global_average, global_average_squared, block_number) << endl;
    energy_file.close();
    
    // Calculate and write acceptance rate
    double acceptance_rate = (total_attempted > 0) ? double(total_accepted)/double(total_attempted) : 0.0;
    ofstream acceptance_file("../OUTPUT/acceptance.dat", ios::app);
    acceptance_file << setw(12) << block_number << setw(12) << acceptance_rate << endl;
    acceptance_file.close();
    
    // Reset block accumulator
    block_average = 0.0;
    
    return acceptance_rate;
}

void QuantumSimulation::reset_for_new_temperature() {
    write_output_headers();
    
    ofstream output("../OUTPUT/output.dat", ios::app);
    output << "Blocking method temperature " << temperature << " completed." << endl << endl;
    output.close();
}

// Simulated Annealing methods 
tuple<double, double, double> QuantumSimulation::run_energy_calculation_no_save(Random &rnd, int sa_iteration) {
    int accepted = 0;
    
    // Equilibration (only for first SA iteration)
    if (sa_iteration == 0) {
        for (int i = 0; i < equilibration_blocks; ++i) {
            for (int j = 0; j < nsteps; ++j) {
                metropolis_step(rnd, accepted);
            }
        }
    }
    
    reset_block_accumulators();
    total_accepted = 0;
    total_attempted = 0;
    double final_acceptance_rate = 0.0;
    
    // Run simulation blocks without saving positions
    for (int block = 0; block < nblocks; ++block) {
        accepted = 0;
        for (int step = 0; step < nsteps; ++step) {
            metropolis_step(rnd, accepted);
            measure_energy();
            total_attempted++;
        } 
        total_accepted += accepted;
        
        double average = block_average / double(nsteps);
        global_average += average;
        global_average_squared += average * average;
        final_acceptance_rate = (total_attempted > 0) ? double(total_accepted)/double(total_attempted) : 0.0;
        block_average = 0.0;
    }
    
    double final_energy = global_average / double(nblocks);
    double final_error = calculate_error(global_average, global_average_squared, nblocks);
    
    return make_tuple(final_energy, final_error, final_acceptance_rate);
}

// FINAL calculation WITH save for optimal parameters
tuple<double, double, double> QuantumSimulation::run_final_energy_calculation(Random &rnd) {
    int accepted = 0;
    
    // Clear previous position file
    clear_position_file();
    
    // Equilibration
    for (int i = 0; i < equilibration_blocks; ++i) {
        for (int j = 0; j < nsteps; ++j) {
            metropolis_step(rnd, accepted);
        }
    }
    
    reset_block_accumulators();
    total_accepted = 0;
    total_attempted = 0;
    double final_acceptance_rate = 0.0;
    
    // Run simulation blocks and save positions
    for (int block = 0; block < nblocks_final; ++block) {
        accepted = 0;
        for (int step = 0; step < nsteps; ++step) {
            metropolis_step(rnd, accepted);
            measure_energy();
            save_position(); // Save every position for final histogram
            total_attempted++;
        } 
        total_accepted += accepted;
        final_acceptance_rate = process_block_averages(block + 1);
    }
    
    double final_energy = global_average / double(nblocks_final);
    double final_error = calculate_error(global_average, global_average_squared, nblocks_final);
    
    return make_tuple(final_energy, final_error, final_acceptance_rate);
}

tuple<double, double> QuantumSimulation::propose_new_parameters(Random &rnd) const {
    double mu_new =abs( mu + rnd.Rannyu(-step_mu * temperature, step_mu * temperature));
    double sigma_new = sigma + rnd.Rannyu(-step_sigma * temperature, step_sigma * temperature);
    return make_tuple(mu_new, sigma_new);
}

bool QuantumSimulation::metropolis_accept_sa(double E_old, double E_new, Random &rnd) const {
    double dE = E_new - E_old;
    if (dE < 0.0) {
        return true;
    } else {
        return rnd.Rannyu() < exp(-dE / temperature);
    }
}

void QuantumSimulation::update_parameters(double new_mu, double new_sigma) {
    mu = new_mu;
    sigma = new_sigma;
}

void QuantumSimulation::cool_temperature() {
    temperature *= alpha;
}


void QuantumSimulation::update_best_parameters(double energy, double error, double acceptance) {
    if (energy < best_energy) {
        best_energy = energy;
        best_mu = mu;
        best_sigma = sigma;
        best_error = error;
        best_acceptance = acceptance;
        steps_without_improvement = 0;  // Reset counter

        ofstream output("../OUTPUT/output.dat", ios::app);
        output << "*** NEW BEST ENERGY FOUND: " << best_energy << " ± " << best_error << " ***" << endl;
        output << "    Best parameters: mu=" << best_mu << ", sigma=" << best_sigma << endl;
        output << "    Acceptance rate: " << best_acceptance << endl;
        output.close();
    } else {
        steps_without_improvement++;  // Increment counter
    }
}

void QuantumSimulation::restore_best_parameters() {
    mu = best_mu;
    sigma = best_sigma;
}

// Output methods
void QuantumSimulation::save_trajectory_step(double energy, double error, double acceptance_rate, double best_energy) const {
    ofstream sa_file("../OUTPUT/sa_trajectory.dat", ios::app);
    sa_file << temperature
            << setw(12) << mu
            << setw(12) << sigma
            << setw(12) << energy
            << setw(12) << error
            << setw(12) << acceptance_rate 
            << setw(12) << best_energy << endl;
    sa_file.close();
}

void QuantumSimulation::save_position() const {
    static bool first_call = true;
    ofstream position_file("../OUTPUT/x.dat", first_call ? ios::out : ios::app);
    position_file << x_position << endl;
    position_file.close();
    first_call = false;
}

void QuantumSimulation::clear_position_file() const {
    ofstream position_file("../OUTPUT/x.dat");
    position_file.close();
}


void QuantumSimulation::print_current_state() const {
    cout << "Current state:" << endl;
    cout << "  Temperature: " << temperature << endl;
    cout << "  Mu: " << mu << ", Sigma: " << sigma << endl;
    cout << "  Position: " << x_position << endl;
}

void QuantumSimulation::print_best_state() const {
    cout << "Best state found:" << endl;
    cout << "  Energy: " << best_energy << " ± " << best_error << endl;
    cout << "  Mu: " << best_mu << ", Sigma: " << best_sigma << endl;
    cout << "  Acceptance: " << best_acceptance << endl;
}

Random QuantumSimulation::setup_random_generator() const {
    Random rnd;
    int seed[4];
    int p1, p2;
    
    ifstream primes_file("../INPUT/Primes");
    if (primes_file.is_open()) {
        primes_file >> p1 >> p2;
        primes_file.close();
    } else {
        cerr << "ERROR: Unable to open ../INPUT/Primes" << endl;
    }
    
    ifstream seed_file("../INPUT/seed.in");
    if (seed_file.is_open()) {
        string property;
        while (!seed_file.eof()) {
            seed_file >> property;
            if (property == "RANDOMSEED") {
                seed_file >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
                break;
            }
        }
        seed_file.close();
    } else {
        cerr << "ERROR: Unable to open ../INPUT/seed.in" << endl;
    }
    
    return rnd;
}

