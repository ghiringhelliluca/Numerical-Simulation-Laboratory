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
#include <vector>
#include <iomanip>
#include "population.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Process ID (continent)
  MPI_Comm_size(MPI_COMM_WORLD, &size);  // Number of processes (continents)

  ofstream coutf("../OUTPUT/output.dat", ios::app);



  // Initialize population for this process
  Population pop;
  pop.initialize();
  
  if (rank == 0) {
    coutf << "Starting MPI Genetic Algorithm with " << size << " processes" << endl;
    coutf << "Migration every " << pop.get_N_migr() << " generations" << endl;
  }
  // Initial sorting
  pop.sort(0);
  
  if (rank == 0) {
    coutf << "Initial populations created and sorted" << endl;
  }

  // Main evolution loop
  for (int gen = 1; gen <= pop.get_generation(); gen++) {
    // Evolution step
    pop.evolve();
    pop.sort(gen);

    // Migration step
    if (gen % pop.get_N_migr() == 0 && size > 1) {

      // Get best individual from current population
      Individual best = pop.get_best();
      vector<int> data = best.encode();
      
      // Prepare buffer for receiving data
      vector<int> received(data.size());
      
      // Define communication partners (ring topology)
      int next = (rank + 1) % size;
      int prev = (rank - 1 + size) % size;
      
      // Exchange best individuals with neighboring processes
      MPI_Sendrecv(data.data(), data.size(), MPI_INT, next, 0,
                  received.data(), data.size(), MPI_INT, prev, 0,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      // Decode received individual and replace worst in population
      Individual immigrant;
      immigrant.set_base_individual(best.get_base_individual());
      immigrant.decode(received);
      pop.replace_worst(immigrant);
      
      // Re-sort population after migration
      pop.sort(gen, false);
        
      if (rank == 0) {
        coutf << "Migration completed at generation " << gen << endl;
      }
    }
    
    // Progress reporting (only rank 0, every 50 generations)
    if (rank == 0 && gen % 50 == 0) {
      double current_best = pop.get_best_fitness();
      coutf << "Generation " << gen << ", Best fitness: " 
          << fixed << setprecision(2) << current_best << " km" << endl;
    }
  }

  // Final global reduction to find the best solution across all processes
  double my_fitness = pop.get_best_fitness();
  double global_best;
  int best_rank;
  
  // Find minimum fitness and which rank has it
  struct {
    double fitness;
    int rank;
  } local, global;
  
  local.fitness = my_fitness;
  local.rank = rank;
  
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  
  // All processes now know the global best and which rank has it
  global_best = global.fitness;
  best_rank = global.rank;
  
  // Output results
  if (rank == 0) {
    coutf << "\n" << string(50, '=') << endl;
    coutf << "FINAL RESULTS:" << endl;
    coutf << "Global best fitness: " << fixed << setprecision(2) 
        << global_best << " km" << endl;
    coutf << "Best solution found by process: " << best_rank << endl;
    coutf << string(50, '=') << endl;
  }
  
  // Each process reports its final best fitness
  coutf << "Process " << rank << " final best: " 
    << fixed << setprecision(2) << my_fitness << " km" << endl;
  
  // Optionally, broadcast the best solution to all processes and save it
  if (rank == best_rank) {
    Individual global_best_individual = pop.get_best();
    
    // Save the globally best solution
    ofstream best_solution("../OUTPUT/global_best_solution.dat");
    best_solution << "Global best solution (distance: " << global_best << " km):" << endl;
    best_solution << "Province order: ";
    global_best_individual.save_path(best_solution);
    best_solution.close();

    coutf << "Global best solution saved to ../OUTPUT/global_best_solution.dat" << endl;
  }
  
  // Synchronize all processes before finalizing
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (rank == 0) {
    coutf << "All processes completed successfully." << endl;
  }

  // Finalize MPI
  MPI_Finalize();

  coutf.close();

  return 0;
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