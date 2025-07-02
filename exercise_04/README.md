
# Molecular Dynamics Simulation Exercises – Lennard-Jones Model

This repository contains the simulation code and data for Exercises `es_04.1`, `es_04.2`, and `es_04.3`, which explore molecular dynamics using the Lennard-Jones potential. 

---

## Directory Structure

Each exercise directory (`es_04.1`, `es_04.2`, `es_04.3`) has the following structure:

```
.
├── SOURCE/          # Source files (NSL_SIMULATOR.cpp, libraries, functions)
├── INPUT/           # Input files
├── OUTPUT/          # Output files
```

---

## Compilation Instructions

To compile the code, navigate to the `SOURCE` folder of the desired exercise and run:

```bash
make
```

This will generate an executable file named:

```bash
simulator.exe
```

---

## How to Run

Once compiled, run the simulation by executing:

```bash
./simulator.exe
```

All output files will be written to the `OUTPUT/` directory automatically.

---

## Exercise Details

### `es_04.1`

- **Input**: Reads particle positions from `INPUT/CONFIG/config.xyz`.
- **Output**: 
  - `OUTPUT/CONFIG/config.xyz`: Final particle positions.
  - `OUTPUT/CONFIG/conf-1.xyz`: Positions at final time minus one timestep.

### `es_04.2`

- **Input**: Reads initial positions from `INPUT/CONFIG/config.fcc`.
- **Output**:
  - `OUTPUT/CONFIG/config.fcc`: Final particle positions.
  - `OUTPUT/CONFIG/conf-1.fcc`: Positions at final time minus one timestep.
  - `OUTPUT/CONFIG/config.vel`: Final particle velocities.

### `es_04.3`

- **Input**:
  - Uses the final positions and velocities from `es_04.2`, saved in `INPUT/CONFIG/`.
  - The only required manual input is setting the temperature (`TEMP`) in `INPUT/input.dat` as the equilibrium temperature (mean of last five blocks values) of `es_04.2` simulation.
  - The `restart` flag in `input.dat` must be set to `1`.
- **Output**:
  - `OUTPUT/CONFIG/`: Final positions and velocities after time reversal simulation.

---