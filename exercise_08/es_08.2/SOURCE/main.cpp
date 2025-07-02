/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "function.h"
#include <tuple>

using namespace std;

int main() {
    
    // Create and initialize the quantum simulation
    QuantumSimulation simulation;
    simulation.initialize_from_file();
    simulation.initialize_best_parameters();

    Random rnd = simulation.setup_random_generator();
    
    ofstream output("../OUTPUT/output.dat", ios::app);
    output << "\n=== Starting Simulated Annealing optimization ===" << endl;
    output << "Initial parameters:" << endl;
    output << "  Temperature: " << simulation.get_temperature() << endl;
    output << "  Mu: " << simulation.get_mu() << ", Sigma: " << simulation.get_sigma() << endl;
    output.close();
    
    // Initial energy 
    auto initial_result = simulation.run_energy_calculation_no_save(rnd, 0);
    double current_energy = get<0>(initial_result);
    double current_error = get<1>(initial_result);
    double current_acceptance = get<2>(initial_result);
    
    // Update best parameters with initial values
    simulation.update_best_parameters(current_energy, current_error, current_acceptance);
    
    output.open("../OUTPUT/output.dat", ios::app);
    output << "Initial energy: " << current_energy << " ± " << current_error << endl;
    output.close();
    
    // Simulated Annealing
    int sa_iteration = 0;
    int max_iterations = simulation.get_sa_iteration_count();
    
    while (sa_iteration < max_iterations) {
        

        simulation.save_trajectory_step(current_energy, current_error, current_acceptance, simulation.get_best_energy());
    
        sa_iteration++;
        
        // Controllo early stopping
        if (simulation.get_steps_without_improvement() >= simulation.get_max_steps_without_improvement()) {
            ofstream output("../OUTPUT/output.dat", ios::app);
            output << "\n=== EARLY STOPPING ===" << endl;
            output << "Stopped after " << simulation.get_steps_without_improvement() 
                << " iterations without improvement (iteration " << sa_iteration << ")" << endl;
            output << "Best energy found: " << simulation.get_best_energy() << " ± " << simulation.get_best_error() << endl;
            output.close();
            break;
        }
        
        //  new parameters
        auto new_parameters = simulation.propose_new_parameters(rnd);
        double mu_new = get<0>(new_parameters);
        double sigma_new = get<1>(new_parameters);
        
        // Skip if sigma becomes negative
        if (sigma_new < 0.0) {
            ofstream output("../OUTPUT/output.dat", ios::app);
            output << "Iteration " << sa_iteration << ": Skipping - negative sigma proposed" << endl;
            output.close();
            simulation.cool_temperature();
            continue;
        }
        
        // Store current parameters for potential rollback
        double mu_old = simulation.get_mu();
        double sigma_old = simulation.get_sigma();
        
        // Update parameters temporarily
        simulation.update_parameters(mu_new, sigma_new);
        
        // Calculate energy with new parameters
        auto new_result = simulation.run_energy_calculation_no_save(rnd, sa_iteration);
        double new_energy = get<0>(new_result);
        double new_error = get<1>(new_result);
        double new_acceptance = get<2>(new_result);
        
        double energy_difference = new_energy - current_energy;
        
        ofstream output("../OUTPUT/output.dat", ios::app);
        
        // Metropolis acceptance for SA
        if (simulation.metropolis_accept_sa(current_energy, new_energy, rnd)) {
            // Accept new parameters
            current_energy = new_energy;
            current_error = new_error;
            current_acceptance = new_acceptance;
            
            // Check if this is the best energy found so far
            simulation.update_best_parameters(current_energy, current_error, current_acceptance);
            
            output << "Iteration " << sa_iteration << ": ACCEPTED" << endl;
            output << "  New energy: " << current_energy << " ± " << current_error << endl;
            output << "  New parameters: mu=" << mu_new << ", sigma=" << sigma_new << endl;
            output << "  Temperature: " << simulation.get_temperature() << endl;
        } else {
            // Reject: restore old parameters
            simulation.update_parameters(mu_old, sigma_old);
            
            output << "Iteration " << sa_iteration << ": REJECTED" << endl;
            output << "  Energy difference: " << energy_difference << endl;
            output << "  Temperature: " << simulation.get_temperature() << endl;
        }
        output.close();
        
        simulation.cool_temperature();
        
        // Print progress every 10 iterations
        if (sa_iteration % 10 == 0) {
            ofstream progress_output("../OUTPUT/output.dat", ios::app);
            progress_output << "\nProgress update (iteration " << sa_iteration << "):" << endl;
            progress_output << "  Temperature: " << simulation.get_temperature() << endl;
            progress_output << "  Current: mu=" << simulation.get_mu() << ", sigma=" << simulation.get_sigma() << endl;
            progress_output << "  Current energy: " << current_energy << " ± " << current_error << endl;
            progress_output << "  Best energy so far: " << simulation.get_best_energy() << " ± " << simulation.get_best_error() << endl;
            progress_output << "  Best parameters: mu=" << simulation.get_best_mu() << ", sigma=" << simulation.get_best_sigma() << endl;
            progress_output.close();
        }
    }
    
    // SA completed - now run final calculation with optimal parameters
    ofstream final_output("../OUTPUT/output.dat", ios::app);
    final_output << "\n=== Simulated Annealing completed! ===" << endl;
    final_output << "Total iterations: " << sa_iteration << endl;
    final_output << "Best energy found: " << simulation.get_best_energy() << " ± " << simulation.get_best_error() << endl;
    final_output << "Best parameters: mu=" << simulation.get_best_mu() << ", sigma=" << simulation.get_best_sigma() << endl;
    final_output << "\n=== Running final calculation with optimal parameters ===" << endl;
    final_output.close();
    
    simulation.restore_best_parameters();
    
    // Reset output files for final calculation
    simulation.reset_for_new_temperature();
    
    // Run final energy calculation with optimal parameters and save positions
    auto final_result = simulation.run_final_energy_calculation(rnd);
    double final_energy = get<0>(final_result);
    double final_error = get<1>(final_result);
    double final_acceptance = get<2>(final_result);
    
    final_output.open("../OUTPUT/output.dat", ios::app);
    final_output << "\n=== FINAL RESULTS ===" << endl;
    final_output << "Optimal parameters found by SA:" << endl;
    final_output << "  Mu: " << simulation.get_best_mu() << endl;
    final_output << "  Sigma: " << simulation.get_best_sigma() << endl;
    final_output << "Final energy calculation with optimal parameters:" << endl;
    final_output << "  Energy: " << final_energy << " ± " << final_error << endl;
    final_output << "  Acceptance rate: " << final_acceptance << endl;
    final_output.close();
    
    cout << "\n=== SIMULATION COMPLETED ===" << endl;
    cout << "Optimal parameters: mu=" << simulation.get_best_mu() << ", sigma=" << simulation.get_best_sigma() << endl;
    cout << "Final energy: " << final_energy << " ± " << final_error << endl;
    
    return 0;
}