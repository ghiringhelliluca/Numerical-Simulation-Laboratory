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
#include <vector>
#include <tuple>
#include <numeric>
#include "function.h"

using namespace std;

int main() {
    
    // simulation object
    QuantumSimulation simulation;
    
    // Initialize parameters from input files
    simulation.initialize_from_file();
    
    Random rnd = simulation.setup_random_generator();
    simulation.print_parameters();
    
    // Run the simulation
    ofstream summary("../OUTPUT/summary.dat", ios::app);    
    summary << "\nStarting quantum simulation..." << endl;
    auto results = simulation.run_simulation(rnd);
    
    double final_energy = get<0>(results);
    double final_error = get<1>(results);
    double final_acceptance = get<2>(results);

    summary << "\nSimulation completed successfully!" << endl;

    summary << "=== QUANTUM SIMULATION RESULTS ===" << endl;
    summary << "Final Energy: " << final_energy << " ± " << final_error << endl;
    summary << "Final Acceptance Rate: " << final_acceptance << endl;
    summary.close();
    
    
    return 0;
}

















