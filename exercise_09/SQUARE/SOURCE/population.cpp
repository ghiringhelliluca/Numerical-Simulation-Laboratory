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
#include "population.h"

using namespace std;
using namespace arma;


void Population::initialize() {

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
  Seed.close();



  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "GA ON A CIRCLE"  << endl;
      else if(_sim_type == 1) coutf << "GA INSIDE A SQUARE"         << endl;
    } else if( property == "M" ){
      input >> _M;
      coutf << "NUMBER OF INDIVIDUAL= " << _M << endl;
    } else if( property == "p" ){
      input >> _p;
      coutf << "EXP FOR SELECTION= " << _p << endl;
    } else if( property == "GENERATION" ){
      input >> _generation;
      coutf << "GENERATION= " << _generation << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else{ 
      cerr << property << endl;
      cerr << "PROBLEM: unknown input" << endl;
    }
  }
  input.close();


  Individual base_indiv;
  base_indiv.initialize(_rnd, _sim_type);  // base individual

  ofstream pos_individual("../OUTPUT/cities.dat");
  if(_sim_type == 0 ){
      pos_individual << "CITY_ID      ANGLE" << endl;
    for(int i =0; i < base_indiv.get_ncity(); i++){
      pos_individual << base_indiv._individual(i)(0) << "             " << base_indiv._individual(i)(1) << endl;
    }
  } else if(_sim_type == 1){
    pos_individual << "CITY_ID      X            Y" << endl;
    for(int i =0; i < base_indiv.get_ncity(); i++){
      pos_individual << base_indiv._individual(i)(0) << "             " << base_indiv._individual(i)(1) << "             " << base_indiv._individual(i)(2) << endl;
    }
  } else {
    cerr << "PROBLEM: unknown simulation type" << endl;
    exit(EXIT_FAILURE);
  }


  _population.set_size(get_M());
  for (int i = 0; i < get_M(); i++) {
    _population(i) = base_indiv;
    _population(i).permutation(_rnd);  // different permutation

    if (!_population(i).check()) {
      cerr << "Invalid individual generated!" << endl;
      exit(1);
    }  
  }

  coutf << "Population initialized." << endl;
  coutf.close();
  pos_individual.close();

  ofstream path("../OUTPUT/best_path.dat");
  path << "Best path for every generation" << endl;
  path.close();

  ofstream loss_progress("../OUTPUT/loss.dat");
  loss_progress << "Generation Loss_Best Loss_Average_Half" << endl;
  loss_progress.close();

  return;
}


int Population :: selection(){
  double r = _rnd.Rannyu();
  int i = int(pow(r,_p) * _M);
  return i ;
}



void Population::evolve() {

  field<Individual> new_population(_M);

  // // copy the best half of the population
  // for (int i = 0; i < _M / 2; i++) {
  //   new_population(i) = _population(i); 
  // }

  // generate new individuals
  for (int i = 0; i < _M; i += 2) {   // _M / 2
    int i1 = selection();
    int i2 = selection();

    Individual child1, child2;

    // --- CROSSOVER ---
    if (_rnd.Rannyu() < 0.8) { // > 50% crossover
      tie(child1, child2) = _population(i1).crossover_pair(_population(i2), _rnd); // two children
    } else {
      child1 = _population(i1);
      child2 = _population(i2);
    }

    // --- MUTATION ---
    if (_rnd.Rannyu() < 0.1){
      child1.mutation_swap_blocks(_rnd);
    }
    
    if (_rnd.Rannyu() < 0.1){
       child1.mutation_shift(_rnd);
    }
    if (_rnd.Rannyu() < 0.1) {
      child1.mutation_inversion(_rnd);
    }
    if (_rnd.Rannyu() < 0.1) {
      child1.permutation_pair(_rnd);
    }


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
      child2.permutation_pair(_rnd);
    }

    // --- CHECK ---
    if (!child1.check() || !child2.check()) {
      cerr << "Invalid individual generated!" << endl;
      exit(1);
    }

    new_population(i) = child1;
    if (i + 1 < _M) new_population(i + 1) = child2;
  }

  _population = new_population;

  return;
}


vec Population::compute_fitness() {
  vec fitness(_M);
  for (int i = 0; i < _M; i++) {
    fitness(i) = _population(i).loss(_sim_type); 
  }
  return fitness;
}


void Population::sort(int gen) {

  vec fitness = compute_fitness();
  uvec sorted_indices = sort_index(fitness); // uvec = unsigned vector, increasingly sorted indices, the best (low fitness) is first

  field<Individual> sorted_pop(_M);
  for (int i = 0; i < _M; ++i) {
    sorted_pop(i) = _population(sorted_indices(i));
  }

  _population = sorted_pop; 

  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Population sorted." << endl;
  coutf.close();

  ofstream path;
  path.open("../OUTPUT/best_path.dat",ios::app);
  _population(0).save_path(path);
  path.close();

  double best_loss = _population(0).loss(_sim_type); 
  double avg_loss = get_average_best_half(); 

  for(int i = 0; i < _M; ++i) {
    _population(i).loss(_sim_type);
  }

  coutf.open("../OUTPUT/loss.dat",ios::app);
  coutf << gen << " " << best_loss << " " << avg_loss << endl;
  coutf.close();

}


double Population::get_average_best_half() {
  double sum = 0.0;
  for (int i = 0; i < _M / 2; ++i) {
    sum += _population(i).loss(_sim_type);
  }
  return sum / (_M / 2);
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
