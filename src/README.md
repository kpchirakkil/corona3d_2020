# Corona3D: 3D Monte Carlo Hot Atom Transport Model

A Monte Carlo simulation for modeling hot atom (O, H) transport and escape in planetary atmospheres (Mars, Venus).

## Quick Start

### Prerequisites
- g++ compiler
- Eigen3 library (`sudo apt install libeigen3-dev` on Ubuntu, `brew install eigen` on Mac)

### Build & Run
```bash
cd src/
make
./corona3d_2020
```

The simulation reads `corona3d_2020.cfg` and writes output to the specified output directory.

## What It Does

Corona3D tracks energetic ("hot") atoms through a planetary atmosphere:

1. **Spawns particles** from physical sources (e.g., O2+ dissociative recombination)
2. **Integrates trajectories** under planetary gravity (Verlet scheme)
3. **Simulates collisions** with background species (O, CO2, CO, N2) using energy-dependent cross-sections
4. **Tracks outcomes**: escape, thermalization, or boundary crossing

## Directory Structure

```
src/
├── corona3d_2020.cfg      # Main configuration
├── Hot_O.cfg / Hot_H.cfg  # Distribution configs
├── inputs/
│   ├── Mars/ Venus/       # Atmospheric profiles
│   └── collisions/        # Cross-section data (from Gacesa, Kharchenko, Kumar)
├── output/                # Default simulation output
└── scripts/               # Analysis + batch/HPC run tools
```

## Key Configuration (corona3d_2020.cfg)

| Parameter | Description |
|-----------|-------------|
| `num_testparts` | Number of test particles |
| `timesteps` / `dt` | Simulation duration and timestep (seconds) |
| `part_type` | H, O, N2, CO, or CO2 |
| `dist_type` | Hot_H, Hot_O, MB, or Import |
| `sim_upper_bound` / `sim_lower_bound` | Altitude boundaries (cm) |
| `bg_part*_config` | Background species config files |

## Output Files

| File | Contents |
|------|----------|
| `density1d_day.out`, `density1d_night.out` | Radial density profiles |
| `column_density_day.out` | Integrated column density |
| `EDF_day_*km.out`, `EDF_night_*km.out` | Energy distribution functions |
| `loss_rates.out` | Escape/loss rates summary |
