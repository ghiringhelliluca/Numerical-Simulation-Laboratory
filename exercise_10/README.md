<!-- sudo apt update
sudo apt install mpich

Potrebbe non funzionare il prossimo codice,  Open MPI non riesce a rilevare il numero di cpu
COMPILARE
mpiexec -np 4 ./a.out

Vedo quante cpu sono disponibili
lscpu | grep '^CPU(s):'

Potrebbe non vedere le cpu disponibili, con questo comando lo forzo l'esecuzione su 4 processi:
ESEGUIRE
mpiexec --oversubscribe -np 4 ./a.out

mpirun --oversubscribe -np 11 ./main.exe -->





# Traveling Salesman Problem (TSP) solved with Genetic Algorithm (GA) and MPI

This project extends the Genetic Algorithm (GA) implementation for solving the Traveling Salesman Problem (TSP) by parallelizing the computation using MPI (Message Passing Interface). The objective is to efficiently find the shortest possible route that visits a set of cities exactly once and returns to the starting point.

The problem is solved for the case of the 110 Italian provincial capitals, using their longitude and latitude coordinates.

The project is organized into three main directories corresponding to different simulation setups:

## Directory Structure



.
├── MIGRATION/        # Parallel GA with inter-process migration of best individuals
│   ├── SOURCE/       # C++ source code (main.cpp, individual.cpp, population.cpp)
│   ├── INPUT/        # Simulation parameters and city positions
│   └── OUTPUT/       # Results of the simulation
├── NO\_MIGRATION/     # Parallel GA with independent populations (no migration)
│   ├── SOURCE/       # C++ source code
│   ├── INPUT/        # Simulation parameters
│   └── OUTPUT/       # Results of the simulation
├── INDEPENDENT/      # Single-process GA with larger population for comparison
│   ├── SOURCE/       # C++ source code
│   ├── INPUT/        # Simulation parameters
│   └── OUTPUT/       # Results of the simulation
└── Countries         # File to plot Italy borders in Python


### MIGRATION/

- **Parallel GA with migration:** 
  - Four processes connected in a ring topology.
  - Every `N_MIGR = 15` generations, each process exchanges its best individual with a neighbor.
  - Helps the processes share good solutions during evolution.

### NO_MIGRATION/

- **Parallel GA without migration:** 
  - Four completely independent GA populations.
  - No communication between processes during the search (`N_MIGR = 1001`).

### INDEPENDENT/

- **Single-process GA** with a population 4 times larger than in MIGRATION and NO_MIGRATION setups, to match the total computational effort.

## Compilation and Execution

Navigate to the `SOURCE/` directory of the desired setup and compile the code with:

```bash
make
````

Run the parallel simulation with MPI (example with 4 processes) and save execution time to a file:

```bash
/usr/bin/time -v mpirun --oversubscribe -np 4 ./main.exe 2> time.txt
```

Or for single-process execution:

```bash
/usr/bin/time -v mpirun --oversubscribe -np 1 ./main.exe 2> time.txt
```

## Configuration

Simulation parameters are defined in the `INPUT/` directory:

* `input.dat`: Contains parameters such as number of generations, population size, selection exponent `p`, crossover and mutation probabilities, and migration frequency `N_MIGR`.
* `provinces.dat`: Contains the dictionary of province IDs and their corresponding geographic coordinates (longitude, latitude).

## Output 

Each simulation generates the following output files inside `OUTPUT/`:

* `provinces.dat`: Dictionary of province IDs and corresponding coordinates.
* `best_path_rank{i}.dat`: Best path at each generation for process with rank `i`.
* `loss_rank{i}.dat`: Evolution of the best loss and the average loss over the best half of the population for process `i`.
* `global_best_solution.dat`: The globally best path found across all processes along with its final loss.


## Execution Time

Execution time for each simulation is measured and saved to `time.txt` in folder `SOURCE` using the `time` command as described above, allowing for direct comparison of computational efficiency between different setups.


