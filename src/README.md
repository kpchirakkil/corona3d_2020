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

## Inelastic Collision Notes

Background species config files can set `rot_population_model` for state-resolved inelastic sampling:

- `ji0` (default): only sample channels originating from `ji=0`
- `thermal`: use Boltzmann populations with nuclear-spin weights. Requires excited-initial-state integral data; currently available for CO2 only. Also requires `thermal_angular_model ji0_proxy` to explicitly accept that only ji=0 angular data exist. Each reversible pair shares one angular distribution, looked up at total energy K + E(ji) in both directions: the 0→j channel distribution for pairs containing j=0 (when available), the aggregate inelastic distribution otherwise. Optional `rot_population_max_j N` truncates the Boltzmann populations to j<=N and renormalizes (sensitivity option).

The Mars O–CO2, O–CO and O–N2 configs now use `rotational_cross_sections_file` with columns `total_energy_eV,ji,jf,sigma_cm2`, including elastic diagonals. Rates and event selection use the same open state-to-state channels. The energy argument is relative kinetic energy plus the initial rotational energy. Ground-state elastic rates retain the independent `total_sigma_file` reference. CO2 populations exclude odd j; missing thermal data cause an error instead of silently falling back to ji=0.

Use `inelastic_model state_resolved` (the default) with `rotational_cross_sections_file`. An enabled state model without its file is an error. Older averaged configurations must explicitly set `inelastic_model legacy_average` and supply all three keys: `total_sigma_file_total`, `elastic_fraction_file`, and `avg_energy_loss_file`. Combining that legacy choice with a rotational state file is an error. The removed `missing_deltaE_use_avg` option has no replacement; state energy transfers come from rotational levels.

Optional `energyN_inelastic_angles_file PATH` entries select transition-specific angular files on the configured energy grid. Paths are explicit and work for any species; there is no directory/filename inference. The supplied CO2 config lists its 41 raw files. Missing, malformed, duplicate, or out-of-range configured entries fail at startup. Files use blocks with a header `E_eV 0 jf`, followed by `theta_degrees DCS Te` rows; DCS and Te must use decimal/scientific notation to distinguish rows from integer state labels. Only ji=0 blocks are supported. Te is checked for finiteness but never used as rotational energy loss. Zero-weight blocks supply no CDF. For uncovered channels or energies, the configured `energyN_inelastic_file` aggregate distribution is used; those aggregate files remain required. Startup messages identify transition-specific coverage, aggregate fallback, and any accepted thermal proxy. The historical `inelastic_channels_iEngNN.csv` files are retained as source artifacts and are no longer loaded.

The signed rigid-rotor transfer is `B*[jf(jf+1)-ji(ji+1)]`. Scattering rotates the incoming relative velocity; closed transitions are excluded before selection. A null-collision sampler accounts for the full relative-speed weighting of Maxwellian targets and processes multiple collisions per timestep. Positions are held fixed during the collision substep, so spatial/timestep convergence still needs testing. Build without `-ffast-math`; finite-value checks are part of the physics validation.

**Data limits:** N2 integral channels come from Kumar et al. (2023), Table 4. CO channels come from Chhabra et al. (online 2022; MNRAS 519, 2023), Table 2, DOI `10.1093/mnras/stac3057`. That CO table contains only jf≤30 although the underlying calculation used j≤70: its inelastic rate is a **truncated estimate**, not a complete collision model. The CO inelastic angular distribution still uses the supplied 18O aggregate proxy, held at 0.3 eV below that angular-data limit. CO/N2 aggregate angles do not resolve correlations between final rotational state and angle. Vibrational, electronic and reactive channels are not included.

For a reversible pair, interpolate a shared reduced cross section `R(E) = (2*ji+1)*(E-Eji)/(E-Ejf)*sigma_up(E)` from the excitation data and obtain the reverse channel by detailed balance. Hold `R` outside the energy grid; excitation vanishes at its rotational opening, while the reverse cross section stays finite. Elastic cross sections are interpolated directly and held beyond their grid. These are explicit threshold/extrapolation approximations, not new quantum-scattering data. A thermal basis with available elastic diagonals must cover at least 99.9% of the population; the supplied CO2 data lack the j=100 diagonal, so that state's population is excluded subject to this tolerance. For elastic DCSs that integrate below the independent integral cross section, unresolved cross section is represented by a zero-angle component; the resolved angular contribution is retained. Its omitted small-angle transport requires angular-resolution sensitivity checks.

`total_cross_section`, `elastic_fraction` and `avg_energy_loss` CSVs in the inelastic directories are now diagnostic ji=0 summaries. The supplied configs sample the state data directly. Regenerate summaries with `python tools/rebuild_collision_summaries.py`; use `--rebuild-co2` to rebuild CO2 state data from the raw files. The older local `scripts/compute_inelastic_tables_*` scripts are legacy generators and must not be used to replace these summaries.

Run the physics regression suite from the repository root with `bash src/tests/run_collision_tests.sh`. It checks conservation, Galilean invariance, zero-angle and zero-speed limits, excitation thresholds, thermal symmetry, detailed balance, total-energy lookup, independent integral rates, thermal-gas/Poisson collision statistics, production-sampler channel frequencies, and explicit/malformed configuration handling.

Run a reproducible finite-duration timestep and CO2 angular-correlation study with `python3 src/tools/check_collision_sensitivity.py --output /tmp/corona-sensitivity-new` from the repository root after building. It retains isolated configs, seeds, logs, and outputs; defaults compare three timesteps and aggregate versus transition-specific CO2 angles over 10 seconds with three seeds. Survival and collision counts at that horizon do not establish converged escape probabilities or test every angular/extrapolation approximation.

The current `Atmosphere.cpp` stops particles below local escape speed and labels them “thermalized”. This is an escape-probability approximation, not a definition of thermal equilibrium; the collision fixes do not change that stopping rule.

## Output Files

| File | Contents |
|------|----------|
| `density1d_day.out`, `density1d_night.out` | Radial density profiles |
| `column_density_day.out` | Integrated column density |
| `EDF_day_*km.out`, `EDF_night_*km.out` | Energy distribution functions |
| `loss_rates.out` | Escape/loss rates summary |
