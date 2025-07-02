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

    string params_file = "../INPUT/input.dat";
    SimulationParams params = read_parameters(params_file);

    ofstream out("../OUTPUT/output.dat");
    
    out << fixed << setprecision(6);

    out << "Parameters loaded from: " << params_file << '\n';
    out << "- Equilibration steps: " << params.n_eq << '\n';


    // ====== EQUILIBRATION ======
    out << "====== EQUILIBRATION ======" << endl << endl;

    // Equilibration for 1s state (near the origin)
    out << "1s state equilibration (starting close to origin):" << endl;
    auto eq_close_100 = metropolis_sampling_optimized(params, rnd, psi100_density, params.n_eq,
                                           params.sigma_100, params.delta_100,
                                           params.x_start, params.y_start, params.z_start_100,
                                           "../OUTPUT/xyz_eq_close100_unif.dat", "../OUTPUT/xyz_eq_close100_gauss.dat", out, false);
    
    vector<double> eq_close_unif_100, eq_close_gauss_100;
    tie(eq_close_unif_100, eq_close_gauss_100) = eq_close_100;
    save_equilibration_data(eq_close_unif_100, eq_close_gauss_100, "../OUTPUT/r_eq_close100.dat");

    out << endl;

    // Equilibration for 1s state (starting far from origin)
    out << "1s state equilibration (starting far from origin):" << endl;
    auto eq_far_100 = metropolis_sampling_optimized(params, rnd, psi100_density, params.n_eq,
                                        params.sigma_100, params.delta_100,
                                        200.0, 200.0, 200.0,
                                        "../OUTPUT/xyz_eq_far100_unif.dat", "../OUTPUT/xyz_eq_far100_gauss.dat", out, false);
    
    vector<double> eq_far_unif_100, eq_far_gauss_100;
    tie(eq_far_unif_100, eq_far_gauss_100) = eq_far_100;
    save_equilibration_data(eq_far_unif_100, eq_far_gauss_100, "../OUTPUT/r_eq_far100.dat");

    out << endl;

    // Equilibration for 2p state (starting close to origin)
    out << "2p state equilibration (starting close to origin):" << endl;
    auto eq_close_210 = metropolis_sampling_optimized(params, rnd, psi210_density, params.n_eq,
                                        params.sigma_210, params.delta_210,
                                        params.x_start, params.y_start, params.z_start_210,
                                          "../OUTPUT/xyz_eq_close210_unif.dat", "../OUTPUT/xyz_eq_close210_gauss.dat", out, false);
    
    vector<double> eq_close_unif_210, eq_close_gauss_210;
    tie(eq_close_unif_210, eq_close_gauss_210) = eq_close_210;
    save_equilibration_data(eq_close_unif_210, eq_close_gauss_210, "../OUTPUT/r_eq_close210.dat");

    out << endl;

    // Equilibration for 2p state (starting far from origin)
    out << "2p state equilibration (starting far from origin):" << endl;
    auto eq_far_210 = metropolis_sampling_optimized(params, rnd, psi210_density, params.n_eq,
                                        params.sigma_210, params.delta_210,
                                        200.0, 200.0, 200.0, 
                                        "../OUTPUT/xyz_eq_far210_unif.dat", "../OUTPUT/xyz_eq_far210_gauss.dat", out,false);
    
    vector<double> eq_far_unif_210, eq_far_gauss_210;
    tie(eq_far_unif_210, eq_far_gauss_210) = eq_far_210;
    save_equilibration_data(eq_far_unif_210, eq_far_gauss_210, "../OUTPUT/r_eq_far210.dat");

    out << endl;
    out << "Equilibration phase completed!" << endl << endl;


    out.close();

    return 0;
}