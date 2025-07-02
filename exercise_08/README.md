
# Variational Monte Carlo (VMC) & Simulated Annealing (SA) for Quantum Mechanics

This project is organized into two main directories, corresponding to Exercise 08.1 and Exercise 08.2.

The exercise is structured to separate equilibration studies, production simulations, and optimization routines.

## Directory Structure

.
├── es_08.1/
│   ├── EQ/        # Equilibration and block size analysis
│   │   ├── SOURCE/  # Source code for equilibration study
│   │   ├── INPUT/   # input.dat, properties.dat
│   │   └── OUTPUT/  # Results of equilibration analysis
│   └── SIM/       # Production simulation for variational energy
│       ├── SOURCE/  # Source code for simulation
│       ├── INPUT/   # input.dat, properties.dat
│       └── OUTPUT/  # Results of VMC simulation
├── es_08.2/       # Variational optimization with Simulated Annealing
│   ├── SOURCE/     # Source code for SA optimization
│   ├── INPUT/      # input.dat, parameters.dat, SA settings
│   └── OUTPUT/     # Results of SA optimization and analysis


### es_08.1

- **EQ/**: used to study the equilibration behavior of the Metropolis sampling and determine a suitable block size for data blocking.
- **SIM/**: performs the full Variational Monte Carlo simulation to estimate the energy using the right parameters found in EQ.

Both folders have the same internal structure with separate SOURCE, INPUT, and OUTPUT directories.

### es_08.2

- Contains the implementation of the Simulated Annealing (SA)** algorithm to optimize the variational parameters (μ, σ).
- Results include the trajectory of the optimization, final optimized energy, and probability density comparisons.

## Compilation and Execution

Navigate into the directory you want to work with (`EQ/`, `SIM/`, or `es_08.2/`) and compile in the `SOURCE/` folder using:

```bash
make
````

Run the compiled executable:

```bash
./main.exe   
```

## Configuration

All simulation and optimization parameters are configured in the following files inside each `INPUT/` directory:

* `input.dat`: Monte Carlo parameters (e.g., number of steps, Metropolis step size), SA algorithm settings (e.g., temperature, alpha, number of SA steps), variational parameters $(\mu, \sigma)$ for VMC

## Output and Analysis

* Outputs include instantaneous energy values, block averages with statistical errors, optimization trajectories, and probability density distributions.
* Data analysis scripts and plotting utilities are provided separately to generate the figures discussed in the report.



