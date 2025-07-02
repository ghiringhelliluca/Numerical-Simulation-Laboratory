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
#include <numeric>
#include <iomanip>
#include "function.h"

using namespace std;

double QuantumSimulation::psi_component(double x, double mu_val, double sigma_val) const {
    return exp(-pow((x - mu_val), 2) / (2 * pow(sigma_val, 2)));
}

void QuantumSimulation::reset_accumulators() {
    block_energies.clear();
    sum_prog = 0.0;
    sum2_prog = 0.0;
}

void QuantumSimulation::write_output_headers() const {
    ofstream results("../OUTPUT/energy.dat");
    results << left << setw(8) << "Block" 
            << setw(15) << "Block_Energy" 
            << setw(18) << "Progressive_Avg" 
            << setw(18) << "Progressive_Err" 
            << setw(15) << "Acceptance" << endl;
    results.close();
}

double QuantumSimulation::calculate_error(double sum_avg, double sum2_avg, int n) const {
    if (n <= 1) return 0.0;
    return sqrt((sum2_avg - sum_avg*sum_avg) / (n - 1));
}

void QuantumSimulation::initialize_from_file() {
    write_output_headers();
    
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) {
        cerr << "Error: Cannot open input.dat" << endl;
        return;
    }
    
    string property;
    while (!input.eof()) {
        input >> property;
        if (property == "TEMP") {
            input >> temperature;
        } else if (property == "DELTA") {
            input >> delta;
        } else if (property == "NBLOCKS") {
            input >> nblocks;
        } else if (property == "NSTEPS") {
            input >> nsteps;
        } else if (property == "NEQUILIBRATION") {
            input >> equilibration_blocks;
        } else if (property == "MU") {
            input >> mu;
        } else if (property == "SIGMA") {
            input >> sigma;
        } else if (property == "ENDINPUT") {
            break;
        } else if (!property.empty()) {
            cerr << "PROBLEM: unknown input " << property << endl;
        }
    }
    input.close();
}

Random QuantumSimulation::setup_random_generator() const {
    Random rnd;
    int seed[4];
    int p1, p2;
    
    ifstream Primes("../INPUT/Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();
    
    ifstream input("../INPUT/seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }
    
    return rnd;
}

double QuantumSimulation::calculate_energy(double x_val, double mu_val, double sigma_val) const {
    double sigma2 = sigma_val * sigma_val;
    double sigma4 = sigma2 * sigma2;

    double psi1 = psi_component(x_val, mu_val, sigma_val);
    double psi2 = psi_component(x_val, -mu_val, sigma_val);
    double psi_total = psi1 + psi2;

    // Check for numerical issues
    if (psi_total < 1e-30) {
        cerr << "Warning: Very small wavefunction at x = " << x_val << endl;
        return 0.0;
    }

    double d2psi1 = ((pow(x_val - mu_val, 2) / sigma4) - (1.0 / sigma2)) * psi1;
    double d2psi2 = ((pow(x_val + mu_val, 2) / sigma4) - (1.0 / sigma2)) * psi2;

    double kin = -0.5 * (d2psi1 + d2psi2) / psi_total;
    double pot = pow(x_val, 4) - 5.0/2.0 * pow(x_val, 2);

    return kin + pot;
}

void QuantumSimulation::metropolis_step(Random &rnd, int &accepted) {
    double x_new = x_position + rnd.Rannyu(-delta, delta);
    double p_old = pow(psi_component(x_position, mu, sigma) + psi_component(x_position, -mu, sigma), 2);
    double p_new = pow(psi_component(x_new, mu, sigma) + psi_component(x_new, -mu, sigma), 2);
    
    double alpha = min(1.0, p_new / p_old);

    total_attempted++;
    if (rnd.Rannyu() < alpha) {
        x_position = x_new;
        accepted++;
        total_accepted++;
    }
}

void QuantumSimulation::reset_for_new_block() {
    // This method can be used for any per-block initialization if needed
}

void QuantumSimulation::add_block_energy(double energy) {
    block_energies.push_back(energy);
}

void QuantumSimulation::calculate_progressive_averages(int block_num, double &prog_avg, double &prog_error) {
    if (block_energies.empty() || block_num <= 0) {
        prog_avg = 0.0;
        prog_error = 0.0;
        return;
    }
    
    // Calculate sum and sum of squares up to current block
    sum_prog = 0.0;
    sum2_prog = 0.0;
    
    for (int i = 0; i < block_num; i++) {
        sum_prog += block_energies[i];
        sum2_prog += block_energies[i] * block_energies[i];
    }
    
    prog_avg = sum_prog / block_num;
    double prog_avg2 = sum2_prog / block_num;
    prog_error = calculate_error(prog_avg, prog_avg2, block_num);
}

tuple<double, double, double> QuantumSimulation::run_simulation(Random &rnd) {
    x_position = 0.0;
    total_accepted = 0;
    total_attempted = 0;
    reset_accumulators();
    
    equilibrate(rnd);
    
    // Main simulation loop
    for (int block = 1; block <= nblocks; block++) {
        reset_for_new_block();
        int block_accepted = 0;
        double block_energy_sum = 0.0;
        
        // Single block measurement
        for (int step = 0; step < nsteps; step++) {
            metropolis_step(rnd, block_accepted);
            double energy = calculate_energy(x_position, mu, sigma);
            block_energy_sum += energy;
        }
        
        // Calculate block average
        double block_energy = block_energy_sum / nsteps;
        add_block_energy(block_energy);
        
        // Calculate progressive averages and error
        double prog_avg, prog_error;
        calculate_progressive_averages(block, prog_avg, prog_error);
        
        // Current acceptance rate
        double acceptance_rate = get_acceptance_rate();
        
        // Save results
        save_results(block, block_energy, prog_avg, prog_error, acceptance_rate);
        
    }
    
    // Final results
    double final_prog_avg, final_prog_error;
    calculate_progressive_averages(nblocks, final_prog_avg, final_prog_error);
    double final_acceptance = get_acceptance_rate();
    
    save_final_summary(final_prog_avg, final_prog_error, final_acceptance);
    
    return make_tuple(final_prog_avg, final_prog_error, final_acceptance);
}

void QuantumSimulation::equilibrate(Random &rnd) {
    // Save current counters
    int saved_accepted = total_accepted;
    int saved_attempted = total_attempted;
    
    for (int i = 0; i < equilibration_blocks; i++) {
        for (int j = 0; j < nsteps; j++) {
            int temp_accepted = 0;
            metropolis_step(rnd, temp_accepted);
        }
    }
    
    // Reset counters for actual measurement
    total_accepted = saved_accepted;
    total_attempted = saved_attempted;
}

void QuantumSimulation::save_results(int block, double block_energy, double prog_avg, double prog_error, double acceptance_rate) const {
    ofstream results("../OUTPUT/energy.dat", ios::app);
    results << left << setw(8) << block 
            << setw(15) << fixed << setprecision(8) << block_energy 
            << setw(18) << prog_avg 
            << setw(18) << prog_error 
            << setw(15) << setprecision(4) << acceptance_rate << endl;
    results.close();
}

void QuantumSimulation::save_final_summary(double final_energy, double final_error, double final_acceptance) const {
    ofstream summary("../OUTPUT/summary.dat", ios::app);
    summary << "\n=== SIMULATION RESULTS ===" << endl;
    summary << "Final Energy: " << fixed << setprecision(8) << final_energy << " ± " << final_error << endl;
    summary << "Final Acceptance Rate: " << setprecision(4) << final_acceptance << endl;
    summary << "Total blocks: " << nblocks << endl;
    summary << "Steps per block: " << nsteps << endl;
    summary << "Equilibration blocks: " << equilibration_blocks << endl;
    summary.close();
}

void QuantumSimulation::clear_output_files() const {
    ofstream("../OUTPUT/energy.dat").close();
    ofstream("../OUTPUT/summary.dat").close();
}


void QuantumSimulation::print_parameters() const {
    ofstream summary("../OUTPUT/summary.dat");
    summary << "=== SIMULATION PARAMETERS ===" << endl;
    summary << "Temperature: " << temperature << endl;
    summary << "Delta: " << delta << endl;
    summary << "Mu: " << mu << endl;
    summary << "Sigma: " << sigma << endl;
    summary << "N Blocks: " << nblocks << endl;
    summary << "N Steps: " << nsteps << endl;
    summary << "Equilibration Blocks: " << equilibration_blocks << endl;
    summary.close();
    

}