/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include <mpi.h>
#include "population.h"

using namespace std;
using namespace arma;

void Population::initialize() {
    int rank, size;
    // assign rank to each MPI process, to know which process is executing
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Number of processes (continents)


    // Initialize random number generator with different seeds for each process
    int p1, p2;
    ifstream Primes("../INPUT/Primes");
    Primes >> p1 >> p2;
    Primes.close();
    
    int seed[4];
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    
    // Modify seed for each MPI process
    for (int i = 0; i < 4; i++) {
        seed[i] += rank * 1000;
    }
    
    _rnd.SetRandom(seed, p1, p2);
    Seed.close();


    // Read input parameters
    ifstream input("../INPUT/input.dat");
    
    string property;
    while (!input.eof()) {
        input >> property;
        if (property == "DATA_FILE") {
            input >> _data_file;
        } else if (property == "M") {
            input >> _M;
        } else if (property == "p") {
            input >> _p;
        } else if (property == "GENERATION") {
            input >> _generation;
        } else if (property == "N_MIGR") {
            input >> _N_migr;
        } else if (property == "ENDINPUT") {
            break;
        } else {
            cerr << "PROBLEM: unknown input " << property << endl;
        }
    }
    input.close();

    // Create base individual from province data
    Individual base_indiv;
    base_indiv.initialize(_rnd, _data_file);

    // Save province coordinates (only rank 0)
    if (rank == 0) {
        ofstream pos_provinces("../OUTPUT/provinces.dat");
        pos_provinces << "PROVINCE_ID      LONGITUDE      LATITUDE" << endl;
        for (int i = 0; i < base_indiv.get_ncity(); i++) {
            pos_provinces << base_indiv._individual(i)(0) << "             " 
                         << base_indiv._individual(i)(1) << "             " 
                         << base_indiv._individual(i)(2) << endl;
        }
        pos_provinces.close();
    }

    // Initialize population
    _population.set_size(_M);
    for (int i = 0; i < _M; i++) {
        _population(i).set_base_individual(base_indiv.get_base_individual());
        _population(i).initialize_from_base(_rnd);
        _population(i).permutation(_rnd);  // different permutation

        if (!_population(i).check()) {
            cerr << "Invalid individual generated!" << endl;
            exit(1);
        }
    }

    if (rank == 0) {
      ofstream coutf;
      string output_filename = "../OUTPUT/output.dat";
      coutf.open(output_filename);

      coutf << "DATA FILE= " << _data_file << endl;
      coutf << "NUMBER OF INDIVIDUALS= " << _M << endl;
      coutf << "EXP FOR SELECTION= " << _p << endl;
      coutf << "GENERATIONS= " << _generation << endl;
      if (size > 1) {
          coutf << "MIGRATION EVERY " << _N_migr << " GENERATIONS" << endl;
      } else {
          coutf << "NO MIGRATION (SINGLE PROCESS)" << endl;
      }
      coutf << "Reading input completed!" << endl;

      coutf.close();
    }


    // Initialize output files
    string path_filename = "../OUTPUT/best_path_rank" + to_string(rank) + ".dat";
    ofstream path(path_filename);
    path << "Best path for every generation on rank " << rank << endl;
    path.close();

    string loss_filename = "../OUTPUT/loss_rank" + to_string(rank) + ".dat";
    ofstream loss_progress(loss_filename);
    loss_progress << "Generation Loss_Best Loss_Average_Half" << endl;
    loss_progress.close();
}

int Population::selection() {
    double r = _rnd.Rannyu();
    int i = int(pow(r, _p) * _M);
    return i;
}

void Population::evolve() {
    field<Individual> new_population(_M);

    // Generate new individuals
    for (int i = 0; i < _M; i += 2) {
        int i1 = selection();
        int i2 = selection();

        Individual child1, child2;

        // CROSSOVER
        if (_rnd.Rannyu() < 0.8) {
            tie(child1, child2) = _population(i1).crossover_pair(_population(i2), _rnd);
        } else {
            child1 = _population(i1);
            child2 = _population(i2);
        }

        // MUTATION for child1
        if (_rnd.Rannyu() < 0.1) {
            child1.mutation_swap_blocks(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child1.mutation_shift(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child1.mutation_inversion(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child1.permutation(_rnd);
        }

        // MUTATION for child2
        if (_rnd.Rannyu() < 0.1) {
            child2.mutation_swap_blocks(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child2.mutation_shift(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child2.mutation_inversion(_rnd);
        }
        if (_rnd.Rannyu() < 0.1) {
            child2.permutation(_rnd);
        }

        // CHECK
        if (!child1.check() || !child2.check()) {
            cerr << "Invalid individual generated!" << endl;
            exit(1);
        }

        new_population(i) = child1;
        if (i + 1 < _M) new_population(i + 1) = child2;
    }

    _population = new_population;
}

vec Population::compute_fitness() {
    vec fitness(_M);
    for (int i = 0; i < _M; i++) {
        fitness(i) = _population(i).loss();
    }
    return fitness;
}

void Population::sort(int gen, bool output) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vec fitness = compute_fitness();
    uvec sorted_indices = sort_index(fitness);

    field<Individual> sorted_pop(_M);
    for (int i = 0; i < _M; ++i) {
        sorted_pop(i) = _population(sorted_indices(i));
    }

    _population = sorted_pop;

    // Save best path
    if (output){
      string path_filename = "../OUTPUT/best_path_rank" + to_string(rank) + ".dat";
      ofstream path(path_filename, ios::app);
      _population(0).save_path(path);
      path.close();

      // Calculate and save statistics
      double best_loss = _population(0).loss();
      double avg_loss = get_average_best_half();

      string loss_filename = "../OUTPUT/loss_rank" + to_string(rank) + ".dat";
      ofstream loss_file(loss_filename, ios::app);
      // cout << best_loss << endl;
      loss_file << gen << " " << best_loss << " " << avg_loss << endl;
      loss_file.close();
    }
}

double Population::get_average_best_half() {
    double sum = 0.0;
    for (int i = 0; i < _M / 2; ++i) {
        sum += _population(i).loss();
    }
    return sum / (_M / 2);
}

Individual Population::get_best() {
    return _population(0);
}

double Population::get_best_fitness() {
    return _population(0).loss();
}

void Population::replace_worst(const Individual& immigrant) {
    _population(_M - 1) = immigrant;
    // Note: We don't call sort here to avoid output issues during migration
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/