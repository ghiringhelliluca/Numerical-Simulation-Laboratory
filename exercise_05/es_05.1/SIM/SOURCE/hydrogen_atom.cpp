/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "hydrogen_atom.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
// #include <algorithm>

using namespace std;

double statistical_error(double average, double average_squared, int n) {
    if (n == 0) {
        return 0.0;
    }
    return sqrt((average_squared - average * average) / n);
}


BlockResult blocking_method(const vector<double>& data, int N, int L) {
    vector<double> block_averages;
    block_averages.reserve(N);
    vector<double> block_averages_squared(N);
   
    
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int start_idx = i * L;
        for (int j = 0; j < L; j++) {
            sum += data[start_idx + j];
        }
        block_averages.push_back(sum / L);
        block_averages_squared[i] = block_averages[i] * block_averages[i];

    }

    vector<double> cumulative_means;
    vector<double> cumulative_errors;
    cumulative_means.reserve(N);
    cumulative_errors.reserve(N);


    
    for (int i = 0; i < N; i++) {
        double running_sum = 0.0;
        double running_sum_squared = 0.0;

        for (int j = 0; j <= i; j++) {
            running_sum += block_averages[j];
            running_sum_squared += block_averages_squared[j];
        }

        // running_sum += block_averages[i];
        // running_sum_squared += block_averages[i] * block_averages[i];
        
        double mean = running_sum / (i + 1);
        double mean_of_squares = running_sum_squared / (i + 1);
        double error = statistical_error(mean, mean_of_squares, i);
        
        cumulative_means.push_back(mean);
        cumulative_errors.push_back(error);
    }

    return {cumulative_means, cumulative_errors};
}

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
    bool save_coordinates
) {
    vector<double> r_uniform, r_gaussian;
    r_uniform.reserve(steps );
    r_gaussian.reserve(steps);
    
    // Starting position
    double x = x_start, y = y_start, z = z_start;
    // r_uniform.push_back(distance_from_origin(x, y, z));
    // r_gaussian.push_back(distance_from_origin(x, y, z));

    ofstream out_unif, out_gauss;
    if (save_coordinates) {
        out_unif.open(output_file_unif);
        out_gauss.open(output_file_gauss);
        if (!out_unif.is_open() || !out_gauss.is_open()) {
            cerr << "Error: Unable to open output files" << endl;
            exit(1);
        }
    }

    // ========== UNIFORM PROPOSAL ==========
    int accepted_uniform = 0;
    
    // Equilibration phase
    double current_density = target_density(x, y, z);
    for (int i = 0; i < params.n_eq; i++) {
        double x_new = x + rnd.Rannyu(-delta, delta);
        double y_new = y + rnd.Rannyu(-delta, delta);
        double z_new = z + rnd.Rannyu(-delta, delta);

        double new_density = target_density(x_new, y_new, z_new);
        
        if (rnd.Rannyu() < new_density / current_density) {
            x = x_new;
            y = y_new;
            z = z_new;
            current_density = new_density;
        }
    }
    
    r_uniform.push_back(distance_from_origin(x, y, z));
    for (int i = 0; i < steps; i++) {
        if (save_coordinates) {
            out_unif << x << ' ' << y << ' ' << z << '\n';
        }

        double x_new = x + rnd.Rannyu(-delta, delta);
        double y_new = y + rnd.Rannyu(-delta, delta);
        double z_new = z + rnd.Rannyu(-delta, delta);

        double new_density = target_density(x_new, y_new, z_new);
        
        if (rnd.Rannyu() < new_density / current_density) {
            x = x_new;
            y = y_new;
            z = z_new;
            current_density = new_density;
            accepted_uniform++;
        }
        
        r_uniform.push_back(distance_from_origin(x, y, z));
    }

    out << "- Uniform proposal acceptance rate: " 
        << static_cast<double>(accepted_uniform) / steps << endl;

    // ========== GAUSSIAN PROPOSAL ==========
    int accepted_gaussian = 0;
    x = x_start; y = y_start; z = z_start;
    current_density = target_density(x, y, z);

    // Equilibration phase
    for (int i = 0; i < params.n_eq; i++) {
        double x_new = x + rnd.Gauss(0, sigma);
        double y_new = y + rnd.Gauss(0, sigma);
        double z_new = z + rnd.Gauss(0, sigma);
        
        double new_density = target_density(x_new, y_new, z_new);
        
        if (rnd.Rannyu() < new_density / current_density) {
            x = x_new;
            y = y_new;
            z = z_new;
            current_density = new_density;
        }
    }
    r_gaussian.push_back(distance_from_origin(x, y, z));
    for (int i = 0; i < steps; i++) {
        if (save_coordinates) {
            out_gauss << x << ' ' << y << ' ' << z << '\n';
        }

        double x_new = x + rnd.Gauss(0, sigma);
        double y_new = y + rnd.Gauss(0, sigma);
        double z_new = z + rnd.Gauss(0, sigma);

        double new_density = target_density(x_new, y_new, z_new);
        
        if (rnd.Rannyu() < new_density / current_density) {
            x = x_new;
            y = y_new;
            z = z_new;
            current_density = new_density;
            accepted_gaussian++;
        }
        
        r_gaussian.push_back(distance_from_origin(x, y, z));
    }

    out << "- Gaussian proposal acceptance rate: " 
        << static_cast<double>(accepted_gaussian) / steps << endl;

    if (save_coordinates) {
        out_unif.close();
        out_gauss.close();
    }

    return make_tuple(move(r_uniform), move(r_gaussian));
}

Random initialize_random_generator() {
    Random rnd;
    int seed[4];
    int p1, p2;



    ifstream primes_file("../INPUT/Primes");
    if (primes_file.is_open()) {
        primes_file >> p1 >> p2;
        primes_file.close();
    } else {
        cerr << "ERROR: Unable to open Primes file" << endl;
        exit(1);
    }

    ifstream seed_file("../INPUT/seed.in");
    string property;
    if (seed_file.is_open()) {
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
        cerr << "ERROR: Unable to open seed.in file" << endl;
        exit(1);
    }

    return rnd;
}


SimulationParams read_parameters(const string& filename) {
    SimulationParams params = {}; // Initialize all to zero
    ifstream input(filename);
    
    if (!input.is_open()) {
        cerr << "ERROR: Unable to open parameters file: " << filename << endl;
        exit(1);
    }

    // cout << "Reading parameters from: " << filename << endl;
    
    string property;
    while (!input.eof()) {
        input >> property;
        
        // Skip comments
        if (property[0] == '#' || property.empty()) {
            string line;
            getline(input, line); // Skip rest of line
            continue;
        }
        
        if (property == "N_EQUILIBRATION") {
            input >> params.n_eq;
            // cout << "N_EQUILIBRATION = " << params.n_eq << endl;
        } else if (property == "N_SAMPLING") {
            input >> params.M;
            // cout << "N_SAMPLING = " << params.M << endl;
        } else if (property == "N_BLOCKS") {
            input >> params.N;
            // cout << "N_BLOCKS = " << params.N << endl;
        } else if (property == "SIGMA_100") {
            input >> params.sigma_100;
            // cout << "SIGMA_100 = " << params.sigma_100 << endl;
        } else if (property == "DELTA_100") {
            input >> params.delta_100;
            // cout << "DELTA_100 = " << params.delta_100 << endl;
        } else if (property == "SIGMA_210") {
            input >> params.sigma_210;
            // cout << "SIGMA_210 = " << params.sigma_210 << endl;
        } else if (property == "DELTA_210") {
            input >> params.delta_210;
            // cout << "DELTA_210 = " << params.delta_210 << endl;
        } else if (property == "X_START") {
            input >> params.x_start;
            // cout << "X_START = " << params.x_start << endl;
        } else if (property == "Y_START") {
            input >> params.y_start;
            // cout << "Y_START = " << params.y_start << endl;
        } else if (property == "Z_START_100") {
            input >> params.z_start_100;
            // cout << "Z_START_100 = " << params.z_start_100 << endl;
        } else if (property == "Z_START_210") {
            input >> params.z_start_210;
            // cout << "Z_START_210 = " << params.z_start_210 << endl;
        } else if (property == "ENDINPUT") {
            // cout << "Reading input completed!" << endl;
            break;
        } else {
            cout << "Unknown parameter " << property << endl;
        }
    }
    
    input.close();
    return params;
}


void save_equilibration_data(const vector<double>& data_unif, 
                           const vector<double>& data_gauss, 
                           const string& filename) {
    ofstream output_file(filename);
    if (!output_file.is_open()) {
        cerr << "ERROR: Unable to open output file: " << filename << endl;
        return;
    }
    
    size_t min_size = min(data_unif.size(), data_gauss.size());
    for (size_t i = 0; i < min_size; i++) {
        output_file << i << ' ' << data_unif[i] << ' ' << data_gauss[i] << '\n';
    }
    
    output_file.close();
}

void save_block_results(const vector<double>& means,
                       const vector<double>& errors,
                       int L, const string& filename) {
    ofstream output_file(filename);
    if (!output_file.is_open()) {
        cerr << "ERROR: Unable to open output file: " << filename << endl;
        return;
    }
    
    for (size_t i = 0; i < means.size(); i++) {
        output_file << i * L << ' ' << means[i] << ' ' << errors[i] << '\n';
    }
    
    output_file.close();
}