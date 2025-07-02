
# Monte Carlo - Quantum Wave Function

This repository contains the simulation code and data for Exercises 5, which sample quantum wave function of a particle using Metropolis algorithm. 

---

## Directory Structure

The `es_05.1` exercise is organized into two main subdirectories:

- `EQ/`: used to study the equilibration behavior of the system

- `SIM/`: used for production sampling and statistical analysis

Each of these folders contains the following structure:

```
.
├── SOURCE/          # Source files (main.cpp, libraries, functions)
├── INPUT/           # Simulation configuration file (input.dat)
├── OUTPUT/          # Output results from the simulation
```

---

## Compilation Instructions

To compile the code, navigate to the `SOURCE` folder of either `EQ/` or `SIM/` and run:

```bash
make
```

This will generate an executable file named:

```bash
main.exe
```

---

## How to Run

Once compiled, run the simulation by executing:

```bash
./main.exe
```

All output files will be written to the `OUTPUT/` directory automatically.

---

## Input 

Each subdirectory (`EQ/` and `SIM/`) contains an `INPUT/input.dat` file which specifies the parameters for the Metropolis algorithm and parameters for the Blocking analysis. 


## Output

Simulation results are written to the `OUTPUT/` directory, and typically include sampled positions and/or observables and averages and statistical uncertainties.

