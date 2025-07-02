# Molecular Dynamics (MD) & Monte Carlo (MC) Simulations


The project is organized into two main directories, corresponding to **Exercise 07.2** and **Exercise 07.4**.

Each directory is divided into **MD** and **MC** (`M(RT)^2/` folder name) simulations to allow separate study of the two methods.

```
.
├── es_07.2/
│ ├── MD/
│ │ ├── EQUILIBRATION/
│ │ │ ├── SOURCE/ # Source code
│ │ │ ├── INPUT/ # input.dat, properties.dat
│ │ │ └── OUTPUT/ # Results of equilibration study
│ │ ├── SIMULATION/
│ │ │ ├── SOURCE/ # Source code
│ │ │ ├── INPUT/ # input.dat, properties.dat
│ │ │ └── OUTPUT/ # Results of production simulation
│ └── M(RT)^2/ # Same structure as MD
```

- **EQUILIBRATION/**: used to study the equilibration behavior of the system.
- **SIMULATION/**: used for full production simulations after equilibration is verified.

In **07.2**, the outputs are **instantaneous values** of physical observables like energy, temperature, and pressure.

```
.
└── es_07.4/
│ ├── MD/
│ │ ├── SOURCE/
│ │ ├── INPUT/ # input.dat, properties.dat
│ │ └── OUTPUT/ # Simulation results
│ └── M(RT)^2/ # Same structure as MD
```

- Here **equilibration is fixed** and performed **within the simulation** itself.
- Output data is analyzed using **data blocking** to estimate statistical errors.




## Compilation and Execution

Navigate into the `MD/` or `M(RT)^2/` directory you want to work with and then compile in the `SOURCE/` folder with:

```bash
make
```

Run the compiled executable:

```bash
./equilibrator.exe  # For equilibration runs for EQUILIBRATION/ folder
./simulator.exe     # For production simulations
```


## Configuration

All simulation parameters and properties to be measured are set in the following files inside `INPUT/`:

* **`input.dat`** → Simulation parameters (e.g., temperature, number of steps, cutoff radius)
* **`properties.dat`** → List of physical quantities to measure during the simulation

