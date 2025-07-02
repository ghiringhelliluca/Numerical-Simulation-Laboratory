# Traveling Salesman Problem (TSP) solved with Genetic Algorithm (GA)

This project implements a Genetic Algorithm (GA) to solve the Traveling Salesman Problem (TSP). The objective is to find the shortest possible route that visits a set of cities exactly once and returns to the starting point.

The project is organized into two main directories corresponding to two different spatial distributions of the cities:

## Directory Structure



.
├── CIRCUMFERENCE/   # Cities uniformly distributed along a circumference
│   ├── SOURCE/       # C++ source code (main.cpp, individual.cpp, population.cpp)
│   ├── INPUT/        # Simulation parameters and city positions
│   └── OUTPUT/       # Results of the simulation
└── SQUARE/          # Cities uniformly distributed inside a square
│   ├── SOURCE/       # C++ source code (main.cpp, individual.cpp, population.cpp)
│   ├── INPUT/        # Simulation parameters and city positions
│   └── OUTPUT/       # Results of the simulation


### CIRCUMFERENCE/

- **Cities uniformly distributed on a circumference of radius 1.**
- Used to test the GA in a symmetric configuration of cities.

### SQUARE/

- **Cities uniformly distributed inside a square of side 1.**
- Used to test the GA with randomly positioned cities in 2D space.

## Compilation and Execution

Navigate to the `SOURCE/` directory of either **CIRCUMFERENCE/** or **SQUARE/** and compile the code with:

```bash
make
````

Execute the program:

```bash
./main.exe
```

## Configuration

Simulation parameters are defined in the `INPUT/` directory:

* `input.dat`: Contains parameters such as number of generations, population size, probabilities for crossover and mutation, and other settings.

## Output and Analysis

The simulation generates the following output files inside `OUTPUT/`:

* `best_path.dat`: Final optimized path as a sequence of city indices.
* `city.dat`: Dictionary matching city indices to their corresponding coordinates.
* `loss.dat`: Contains two columns: the best loss (shortest path length) and the average loss of the best half of the population for each generation. Useful for plotting convergence trends.

