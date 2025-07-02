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
#include <vector>
#include <tuple>
#include <fstream>
#include <iomanip>
#include "hydrogen_atom.h"

using namespace std;

int main(int argc, char* argv[]) {
    
    Random rnd = initialize_random_generator();

    // Load parameters
    string params_file = "../INPUT/input.dat";
    SimulationParams params = read_parameters(params_file);

    ofstream out("../OUTPUT/output.dat");

    out << fixed << setprecision(6);

    out << "Parameters loaded from: " << params_file << '\n';
    out << "- Equilibration steps: " << params.n_eq << '\n';
    out << "- Sampling steps: " << params.M << '\n';
    out << "- Number of blocks: " << params.N << '\n';
    out << "- Block size: " << params.M / params.N << "\n\n";

    int L = params.M / params.N;

    out << "====== SAMPLING ======\n\n";

    // ========== 1S STATE SAMPLING ==========
    out << "1s state sampling:\n";
    auto sampling_100 = metropolis_sampling_optimized(
        params, rnd, psi100_density, params.M, 
        params.sigma_100, params.delta_100, 
        params.x_start, params.y_start, params.z_start_100,
        "../OUTPUT/xyz_100_unif.dat", "../OUTPUT/xyz_100_gauss.dat", 
        out, true
    );
    
    vector<double> r_100_unif, r_100_gauss;
    tie(r_100_unif, r_100_gauss) = sampling_100;

    // Block analysis for 1s state
    BlockResult result_100_unif = blocking_method(r_100_unif, params.N, L);
    BlockResult result_100_gauss = blocking_method(r_100_gauss, params.N, L);

    save_block_results(result_100_unif.mean, result_100_unif.error, L, 
                      "../OUTPUT/output_100_unif.dat");
    save_block_results(result_100_gauss.mean, result_100_gauss.error, L, 
                      "../OUTPUT/output_100_gauss.dat");

    out << '\n';

    // ========== 2P STATE SAMPLING ==========
    out << "2p state sampling:\n";
    auto sampling_210 = metropolis_sampling_optimized(
        params, rnd, psi210_density, params.M, 
        params.sigma_210, params.delta_210, 
        params.x_start, params.y_start, params.z_start_210,
        "../OUTPUT/xyz_210_unif.dat", "../OUTPUT/xyz_210_gauss.dat", 
        out, true
    );
    
    vector<double> r_210_unif, r_210_gauss;
    tie(r_210_unif, r_210_gauss) = sampling_210;

    // Block analysis for 2p state
    BlockResult result_210_unif = blocking_method(r_210_unif, params.N, L);
    BlockResult result_210_gauss = blocking_method(r_210_gauss, params.N, L);

    save_block_results(result_210_unif.mean, result_210_unif.error, L, 
                      "../OUTPUT/output_210_unif.dat");
    save_block_results(result_210_gauss.mean, result_210_gauss.error, L, 
                      "../OUTPUT/output_210_gauss.dat");

    out << "\nSampling phase completed!\n\n";

    // ========== RESULTS ==========
    out << "====== FINAL RESULTS ======\n";
    out << "1s state final estimates:\n";
    out << "- Uniform proposal: <r> = " << result_100_unif.mean.back() 
        << " +- " << result_100_unif.error.back() << '\n';
    out << "- Gaussian proposal: <r> = " << result_100_gauss.mean.back() 
        << " +- " << result_100_gauss.error.back() << '\n';
    
    out << "\n2p state final estimates:\n";
    out << "- Uniform proposal: <r> = " << result_210_unif.mean.back() 
        << " +- " << result_210_unif.error.back() << '\n';
    out << "- Gaussian proposal: <r> = " << result_210_gauss.mean.back() 
        << " +- " << result_210_gauss.error.back() << '\n';

    out << "Simulation completed\n";

    out.close();
    return 0;
}