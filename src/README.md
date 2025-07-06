# Corona3D: 3D Monte Carlo Hot Atom Transport Model

## Overview

Corona3D is a sophisticated 3D Monte Carlo simulation framework for modeling hot atom transport in planetary atmospheres. The model tracks the trajectories of energetic atoms (hot oxygen and hot hydrogen) as they undergo collisions with atmospheric constituents, providing insights into atmospheric escape processes and exospheric dynamics at Mars and Venus.

## Table of Contents

1. [Physical Model](#physical-model)
2. [Mathematical Framework](#mathematical-framework)
3. [Code Architecture](#code-architecture)
4. [Core Classes](#core-classes)
5. [Collision Physics](#collision-physics)
6. [Input Configuration](#input-configuration)
7. [Output and Analysis](#output-and-analysis)
8. [Usage Guide](#usage-guide)

## Physical Model

### Hot Atom Transport

The model simulates the transport of "hot" atoms—atoms with kinetic energies significantly above thermal equilibrium—through planetary atmospheres at Mars and Venus. These hot atoms are typically produced by:

**Hot Oxygen Production Mechanisms:**
- **Dissociative recombination**: O₂⁺ + e⁻ → O + O*
- **Charge exchange**: O⁺ + CO₂ → O* + CO₂⁺
- **Photodissociation**: CO₂ + hν → CO + O*

**Hot Hydrogen Production Mechanisms:**
- **H₂ Photodissociation**: H₂ + hν → H + H*
- **Charge exchange**: H⁺ + atmospheric neutrals → H* + ions
- **HCO⁺ Dissociative recombination**: HCO⁺ + e⁻ → H* + CO

### Key Physical Processes

1. **Gravitational Forces**: Hot atoms experience planetary gravity according to Newton's law of universal gravitation, with acceleration: **a = -GM/r²**
2. **Elastic Collisions**: Kinetic energy and momentum are conserved during collisions with background atmospheric species
3. **Inelastic Collisions**: Energy is transferred to internal degrees of freedom (rotation/vibration) of collision partners
4. **Atmospheric Escape**: Particles escape when their total energy (kinetic + gravitational potential) exceeds zero
5. **Thermalization**: High-energy atoms lose energy through collisions until they reach thermal equilibrium with the atmosphere

### Physical Environment

#### Atmospheric Structure
- **Density Profiles**: Exponentially decreasing number density with altitude: n(h) = n₀ exp(-h/H)
- **Scale Heights**: Characteristic height over which density drops by factor e: H = kT/(mg)
- **Temperature Profiles**: Altitude-dependent atmospheric temperature affecting collision cross sections

#### Collision Partners

The model supports comprehensive collision physics between hot atoms and background atmospheric species for both Mars and Venus:

**Hot Oxygen Collision Systems:**
- **O-O**: Oxygen-oxygen elastic scattering (Kharchenko et al. 2000)
- **O-CO₂**: Oxygen-carbon dioxide interactions (elastic + inelastic channels, Gacesa et al. 2020)
- **O-CO**: Oxygen-carbon monoxide interactions (elastic + inelastic channels, Kumar et al. 2022)
- **O-N₂**: Oxygen-nitrogen interactions (elastic + inelastic channels, Kumar et al. 2023)

**Hot Hydrogen Collision Systems:**
- **H-O**: Hydrogen-oxygen elastic scattering
- **H-CO₂**: Hydrogen-carbon dioxide interactions
- **H-CO**: Hydrogen-carbon monoxide interactions
- **H-N₂**: Hydrogen-nitrogen interactions

## Mathematical Framework

### Classical Mechanics and Trajectory Integration

Hot atoms follow classical trajectories under planetary gravity using Newton's second law:

```
d²r/dt² = -GM_planet/r³ * r + F_collision/m
```

Where:
- **r** is the 3D position vector [x, y, z] in cm
- **G** is the gravitational constant (6.67430×10⁻⁸ cm³/g/s²)
- **M_planet** is the planetary mass (kg converted to g)
- **F_collision** represents impulsive collision forces
- **m** is the particle mass (g)

#### Numerical Integration Scheme
The code uses a modified Verlet integration scheme with fixed time steps:

```cpp
// Calculate acceleration at current position
double inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
Array<double, 3, 1> a = k_g*position.array()*inv_r_cube;

// Position update: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt²
position.array() = position.array() + (velocity.array()*dt) + (0.5*a*dt*dt);
radius = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
inverse_radius = 1.0 / radius;

// Calculate acceleration at new position
inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
a = a + k_g*position.array()*inv_r_cube;

// Velocity update: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
velocity.array() = velocity.array() + 0.5*a*dt;
```

### Collision Probability and Cross Sections

#### Collision Probability Calculation

The probability of collision during time step dt follows Poisson statistics:

```
P_collision = 1 - exp(-λ*dt)
```

Where the collision rate λ is:
```
λ = Σᵢ nᵢ(r) * σᵢ(E) * v_rel
```

- **nᵢ(r)** = number density of background species i at position r [cm⁻³]
- **σᵢ(E)** = energy-dependent total collision cross section [cm²]
- **v_rel** = relative velocity between particles [cm/s]

#### Energy-Dependent Cross Sections

Total cross sections σ_total(E) are provided as lookup tables and interpolated:
```cpp
double collision_energy_eV = 0.5 * μ * v_rel² / constants::ergev;
double total_sigma = interpolate_cross_section(collision_energy_eV);
```

Where μ is the reduced mass: **μ = m₁m₂/(m₁ + m₂)**

#### Differential Cross Sections

Angular-dependent differential cross sections determine post-collision scattering angles through cumulative distribution functions (CDFs):

```
dσ/dΩ (E,θ) = f(E, θ) [cm²/steradian]
```

The scattering angle θ is sampled from the normalized CDF:
```cpp
double u = uniform_random_0_to_1();
double theta = sample_from_CDF(u, collision_energy);
```

### Energy and Momentum Conservation

#### Center-of-Mass Transformation

For collision calculations, the code transforms to the center-of-mass (CM) frame:

```cpp
// Center-of-mass velocity
Matrix<double, 3, 1> v_cm = (m₁*v₁ + m₂*v₂) / (m₁ + m₂);

// Relative velocity in CM frame  
Matrix<double, 3, 1> v_rel = v₁ - v₂;
double v_rel_magnitude = sqrt(v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]);
```

#### Elastic Collision Dynamics

For elastic collisions, kinetic energy is conserved in the CM frame:

```cpp
// Post-collision relative velocity (same magnitude, new direction)
Matrix<double, 3, 1> v_rel_new = v_rel_magnitude * Matrix<double, 3, 1>(cos(θ), sin(θ)*cos(φ), sin(θ)*sin(φ));

// Transform back to lab frame
v₁_new = v_cm + (m₂/(m₁ + m₂)) * v_rel_new;
v₂_new = v_cm - (m₁/(m₁ + m₂)) * v_rel_new;
```

Where:
- **θ** = scattering angle sampled from differential cross section
- **φ** = azimuthal angle (uniformly random from 0 to 2π)

#### Inelastic Collision Energy Loss

For inelastic collisions, energy is transferred to internal degrees of freedom of the target molecule (rotation, vibration, electronic excitation):

```
E_kinetic_final = E_kinetic_initial - ΔE_internal
```

The energy loss ΔE_internal represents the specific quantum energy difference between initial and final rovibrational states:

```
ΔE_internal = E_internal(v',J') - E_internal(v,J)
```

Where (v,J) and (v',J') are the initial and final vibrational and rotational quantum numbers. This energy transfer directly affects the post-collision kinematics and must be properly accounted for in momentum conservation.

### Escape Condition and Energy Analysis

#### Gravitational Binding Energy

A particle escapes when its total mechanical energy exceeds zero:

```
E_total = E_kinetic + E_gravitational > 0
```

Where:
```
E_kinetic = ½mv²
E_gravitational = -GMm/r
```

#### Escape Velocity

The minimum velocity for escape at radius r is:
```
v_escape(r) = √(2GM/r)
```

#### Energy Distribution Functions (EDFs)

The code calculates energy distribution functions at specific altitudes:

```cpp
double energy_eV = 0.5 * mass * velocity² / constants::ergev;
int energy_bin = (int)(20.0 * energy_eV);  // 0.05 eV bins

double cos_theta = v_radial / v_total;  // Direction cosine
int angle_bin = (int)(100.0 * |cos_theta|);  // 0.01 cosine bins
```

### Statistical Sampling and Monte Carlo Method

#### Initial Condition Sampling

Hot atom initial conditions are sampled from physical production mechanisms:

**Dissociative Recombination**: O₂⁺ + e⁻ → O + O*
```cpp
double excess_energy_eV = 1.5;  // Typical excess energy
double velocity = sqrt(2.0 * excess_energy_eV * constants::ergev / mass);

// Isotropic velocity direction
double phi = 2π * uniform_random();
double cos_theta = 2.0 * uniform_random() - 1.0;
Vector3d v_initial = velocity * [sin(acos(cos_theta))*cos(phi), 
                                sin(acos(cos_theta))*sin(phi), 
                                cos_theta];
```

#### Background Density Profiles

Atmospheric density varies with altitude according to hydrostatic equilibrium:

```cpp
// Exponential atmosphere model
double scale_height = k_B * temperature / (mass_background * g);
double density_at_altitude = reference_density * exp(-(altitude - reference_altitude) / scale_height);
```

#### Collision Target Selection

When multiple background species are present, the collision target is selected probabilistically:

```cpp
double total_collision_rate = Σᵢ nᵢ * σᵢ * v_rel;
double selection_probability = (nᵢ * σᵢ * v_rel) / total_collision_rate;
```

## Code Architecture

### Directory Structure

```
src/
├── README.md                    # This file
├── main.cpp                     # Main simulation driver
├── corona3d_2020.cfg           # Main configuration file
├── makefile                     # Build configuration
│
├── Core Physics Classes:
├── Particle.hpp/.cpp           # Base particle class
├── Particle_H.hpp/.cpp         # Hydrogen particle implementation
├── Particle_O.hpp/.cpp         # Oxygen particle implementation
├── Particle_CO.hpp/.cpp        # Carbon monoxide particle
├── Particle_CO2.hpp/.cpp       # Carbon dioxide particle
├── Particle_N2.hpp/.cpp        # Nitrogen particle
│
├── Environment Classes:
├── Planet.hpp/.cpp              # Planetary properties and gravity
├── Atmosphere.hpp/.cpp          # Main simulation controller
├── Background_Species.hpp/.cpp  # Atmospheric background species
│
├── Distribution Classes:
├── Distribution.hpp/.cpp        # Base distribution class
├── Distribution_Hot_H.hpp/.cpp  # Hot hydrogen distributions
├── Distribution_Hot_O.hpp/.cpp  # Hot oxygen distributions
├── Distribution_MB.hpp/.cpp     # Maxwell-Boltzmann distributions
├── Distribution_Import.hpp/.cpp # Imported velocity distributions
│
├── Utility Classes:
├── Common_Functions.hpp/.cpp    # Shared utility functions
├── Interpolator.hpp/.cpp        # Data interpolation routines
├── vtally.hpp/.cpp             # Velocity tallying and statistics
│
└── Input Data:
    ├── inputs/                  # Configuration and data files
    │   ├── *.cfg               # Species-specific configuration
    │   ├── collisions/         # Collision cross section data
    │   └── */                  # Atmospheric profile data
    └── model_output_data/      # Reference output data
```

### Build System

The project uses a standard makefile with automatic dependency tracking:

```bash
make          # Build the executable
make clean    # Clean build artifacts
make depend   # Update dependencies
```

## Core Classes

### 1. Particle Class Hierarchy

#### `Particle` (Base Class)
- **Purpose**: Abstract base class for all particle types
- **Key Methods**:
  - `do_collision()`: Handles collision physics
  - `do_timestep()`: Advances particle position/velocity
  - `get_energy_in_eV()`: Returns kinetic energy
  - `is_active()`: Checks if particle is still being tracked

#### Derived Classes
- **`Particle_O`**: Oxygen atoms (hot oxygen transport)
- **`Particle_H`**: Hydrogen atoms (hot hydrogen transport)
- **`Particle_CO`**: Carbon monoxide molecules (background species)
- **`Particle_CO2`**: Carbon dioxide molecules (background species)
- **`Particle_N2`**: Nitrogen molecules (background species)

Each derived class implements:
```cpp
virtual double get_mass() const = 0;
virtual string get_name() const = 0;
```

### 2. Atmosphere Class

The `Atmosphere` class serves as the main simulation controller:

#### Key Responsibilities:
- **Particle Management**: Creates and tracks test particles
- **Time Evolution**: Advances simulation through timesteps
- **Collision Detection**: Determines when collisions occur
- **Output Generation**: Records trajectories and statistics

### Core Algorithm:
```cpp
for (timestep = 0; timestep < max_timesteps; timestep++) {
    for (int j = 0; j < active_parts; j++) {
        int particle_index = active_indices[j];
        
        // 1. Update particle statistics and velocity tallying
        update_stats(dt, particle_index);
        my_vtally.update_vtally(my_parts[particle_index]);
        
        // 2. Integrate equations of motion under gravity (leapfrog scheme)
        my_parts[particle_index]->do_timestep(dt, k_g);
        
        // 3. Check for collisions with background species
        if (bg_species.check_collision(my_parts[particle_index], dt)) {
            // Execute collision dynamics with selected target and scattering angle
            my_parts[particle_index]->do_collision(
                bg_species.get_collision_target(), 
                bg_species.get_collision_theta(), 
                i*dt, 
                my_planet.get_radius()
            );
        }
        
        // 4. Check deactivation conditions
        double v_esc_current = sqrt(2.0 * constants::G * my_planet.get_mass() / my_parts[particle_index]->get_radius());
        double v_thermal = v_esc_current;  // Thermalization threshold
        
        if (my_parts[particle_index]->get_total_v() < v_thermal) {
            my_parts[particle_index]->deactivate("Particle was thermalized");
            active_parts--;
            active_indices.erase(active_indices.begin() + j);
            j--;
        } else if (my_parts[particle_index]->get_radius() >= upper_bound && 
                   my_parts[particle_index]->get_total_v() >= v_esc_upper) {
            // Escape condition with hemispheric tracking
            if (my_parts[particle_index]->get_x() > 0.0) {
                my_parts[particle_index]->deactivate("Reached upper bound on day side with escape velocity");
                day_escape_count++;
            } else {
                my_parts[particle_index]->deactivate("Reached upper bound on night side with escape velocity");
                night_escape_count++;
            }
            active_parts--;
            active_indices.erase(active_indices.begin() + j);
            j--;
        } else if (my_parts[particle_index]->get_radius() <= lower_bound) {
            my_parts[particle_index]->deactivate("Dropped below lower bound");
            active_parts--;
            active_indices.erase(active_indices.begin() + j);
            j--;
        }
    }
    
    // 5. Output diagnostics and trace data
    if (print_status_freq > 0 && (i+1) % print_status_freq == 0) {
        output_simulation_status(i, dt, active_parts, day_escape_count, night_escape_count);
    }
    if (num_traced > 0) {
        output_trace_data();
    }
}
```

### 3. Background_Species Class

Manages atmospheric background constituents and collision physics:

#### Key Responsibilities:
- **Density Profile Management**: Altitude-dependent number densities from observational data
- **Temperature Profiles**: Atmospheric temperature structure affecting collision dynamics  
- **Cross Section Integration**: Energy-dependent collision cross sections with interpolation
- **Collision Detection**: Probabilistic collision occurrence using Monte Carlo sampling
- **Scattering Dynamics**: Post-collision velocity updates using conservation laws

#### Collision Detection Algorithm:
```cpp
bool Background_Species::check_collision(shared_ptr<Particle> p, double dt) {
    // 1. Calculate background densities at particle location
    double r = p->get_radius();
    double alt = r - my_planet.get_radius();
    vector<double> dens(num_species);
    
    if (use_dens_profile) {
        // Get density from imported profile
        for (int i = 0; i < num_species; i++) {
            dens[i] = get_density(alt, i);
        }
    } else {
        // Calculate density using exponential scale height
        double r_moved = my_planet.get_radius() + ref_height - r;
        for (int i = 0; i < num_species; i++) {
            dens[i] = calc_new_density(bg_densities[i][0], bg_scaleheights[i][0], r_moved);
        }
    }
    
    // 2. Calculate collision energies and cross sections
    vector<double> energy(num_species);
    vector<double> total_sig(num_species);
    double my_total_v = p->get_total_v();
    
    for (int i = 0; i < num_species; i++) {
        if (bg_sigma_defaults[i] == 0.0) {
            // Use energy-dependent cross section table
            double avg_v = (use_temp_profile) ? avg_v_interp[i]->loglinterp(alt) : bg_avg_v[i][0];
            my_dist->init_vonly(bg_parts[i], avg_v);
            energy[i] = calc_collision_e(p, bg_parts[i]);
            total_sig[i] = sigma_interp[i]->linterp(energy[i]);
        } else {
            // Use default cross section
            total_sig[i] = bg_sigma_defaults[i];
        }
    }
    
    // 3. Calculate total collision probability using Beer-Lambert law
    double tau = 0.0;
    for (int i = 0; i < num_species; i++) {
        tau += my_total_v * dt * total_sig[i] * dens[i];
    }
    double collision_probability = 1.0 - exp(-tau);
    
    // 4. Determine if collision occurs
    double u = common::get_rand();
    if (u > collision_probability) {
        collision_target = -1;
        return false;  // No collision
    }
    
    // 5. Select collision partner probabilistically
    num_collisions++;
    u = common::get_rand();
    double total_dens = 0.0;
    for (int i = 0; i < num_species; i++) {
        total_dens += dens[i];
    }
    
    double frac = 0.0;
    collision_target = 0;
    do {
        frac += dens[collision_target] / total_dens;
        collision_target++;
    } while (u >= frac && collision_target < num_species);
    collision_target--;
    
    // 6. Initialize collision target and sample scattering angle
    if (bg_sigma_defaults[collision_target] != 0.0) {
        double avg_v = (use_temp_profile) ? avg_v_interp[collision_target]->loglinterp(alt) : bg_avg_v[collision_target][0];
        my_dist->init_vonly(bg_parts[collision_target], avg_v);
        energy[collision_target] = calc_collision_e(p, bg_parts[collision_target]);
    }
    collision_theta = find_new_theta(collision_target, energy[collision_target]);
    return true;
}
```

#### Cross Section Data Integration:
- **Total Cross Sections**: Energy-dependent σ_total(E) from quantum scattering calculations
- **Differential Cross Sections**: Angular distributions dσ/dΩ(E,θ) for realistic scattering
- **Inelastic Channels**: State-to-state cross sections for rotational/vibrational excitation
- **Statistical Averaging**: Multiple potential energy surfaces combined for O-CO₂ collisions

#### Atmospheric Density Models:
```cpp
double Background_Species::get_density_at_altitude(double altitude, int species_index) {
    if (use_imported_profiles) {
        // Use spacecraft observations (MAVEN, Venus Express, etc.)
        return log_interpolate_density_profile(altitude, species_index);
    } else {
        // Use exponential atmosphere model
        double scale_height = k_B * temperature / (background_mass[species_index] * gravity);
        return reference_density * exp(-(altitude - reference_altitude) / scale_height);
    }
}
```

### 4. Distribution Classes

Generate initial conditions for hot atoms:

#### `Distribution_Hot_O` (Oxygen)
Production mechanisms for hot oxygen atoms at Mars and Venus:

**1. O₂⁺ Dissociative Recombination**:
```
O₂⁺ + e⁻ → O + O* (with ~1.5 eV excess energy)
```
Implementation:
```cpp
void Distribution_Hot_O::init_O2plus_DR_particle(shared_ptr<Particle> p) {
    // 1. Sample spatial location from O₂⁺ density profile
    double radius = sample_from_ion_density_CDF();
    
    // 2. Hemispherical constraint (dayside production)
    double phi = 2π * uniform_random();
    double cos_theta = uniform_random();  // 0 to 1 for hemisphere
    double x = radius * sqrt(1 - cos_theta²) * cos(phi);
    double y = radius * sqrt(1 - cos_theta²) * sin(phi);  
    double z = radius * cos_theta;
    if (x > 0) x = -x;  // Force to nightside hemisphere
    
    // 3. Calculate excess kinetic energy from quantum chemistry
    double excess_energy_eV = 1.5;  // Representative value from branching ratios
    double kinetic_energy_ergs = excess_energy_eV * constants::ergev;
    double velocity_magnitude = sqrt(2.0 * kinetic_energy_ergs / oxygen_mass);
    
    // 4. Isotropic velocity distribution
    double velocity_phi = 2π * uniform_random();
    double velocity_cos_theta = 2.0 * uniform_random() - 1.0;
    double vx = velocity_magnitude * sqrt(1 - velocity_cos_theta²) * cos(velocity_phi);
    double vy = velocity_magnitude * sqrt(1 - velocity_cos_theta²) * sin(velocity_phi);
    double vz = velocity_magnitude * velocity_cos_theta;
    
    p->init_particle(x, y, z, vx, vy, vz);
}
```

**2. Charge Exchange Reactions**:
```
O⁺ + CO₂ → O* + CO₂⁺ (with ~2-3 eV excess energy)
```

**3. Photodissociation**:
```
CO₂ + hν → CO + O* (UV photolysis with wavelength-dependent energy)
```

#### `Distribution_Hot_H` (Hydrogen)
Production mechanisms for hot hydrogen atoms at Mars and Venus:

**1. H₂ Photodissociation**:
```cpp
void Distribution_Hot_H::init_H2_photodiss_particle(shared_ptr<Particle> p) {
    // Sample from H₂ density profile and apply Lyman-α photolysis cross section
    double photodissociation_energy_eV = calculate_wavelength_dependent_energy();
    // ... similar spatial and velocity sampling
}
```

**2. HCO⁺ Dissociative Recombination**:
```cpp
void Distribution_Hot_H::init_HCOplus_DR_particle(shared_ptr<Particle> p) {
    // Sample from HCO⁺ density profile and apply DR cross section
    double excess_energy_eV = 1.5;  // Representative value
    // ... spatial and velocity sampling
}
```

**3. Charge Exchange**: H⁺ + atmospheric neutrals → H* + ions

**4. Ion-Neutral Reactions**: Complex multi-step processes with intermediate energy transfer

#### `Distribution_MB` (Maxwell-Boltzmann)
Thermal equilibrium distributions for validation and background particle initialization:

```cpp
void Distribution_MB::init(shared_ptr<Particle> p) {
    // 1. Sample spatial location on reference sphere
    double phi = 2π * uniform_random();
    double cos_theta = 2.0 * uniform_random() - 1.0;
    Vector3d position = reference_radius * Vector3d(
        sqrt(1 - cos_theta²) * cos(phi),
        sqrt(1 - cos_theta²) * sin(phi), 
        cos_theta
    );
    
    // 2. Sample velocity from Maxwell-Boltzmann distribution
    double thermal_velocity = sqrt(k_B * reference_temperature / particle_mass);
    Vector3d velocity = sample_maxwell_boltzmann_3D(thermal_velocity);
    
    p->init_particle(position.x(), position.y(), position.z(), 
                    velocity.x(), velocity.y(), velocity.z());
}

Vector3d Distribution_MB::sample_maxwell_boltzmann_3D(double v_thermal) {
    // Box-Muller transformation for Gaussian velocity components
    static bool spare_available = false;
    static double spare_value;
    
    auto gaussian_random = [&]() -> double {
        if (spare_available) {
            spare_available = false;
            return spare_value;
        }
        spare_available = true;
        double u1 = uniform_random();
        double u2 = uniform_random();
        double magnitude = sqrt(-2.0 * log(u1));
        spare_value = magnitude * cos(2π * u2);
        return magnitude * sin(2π * u2);
    };
    
    return v_thermal * Vector3d(gaussian_random(), gaussian_random(), gaussian_random());
}
```

#### `Distribution_Import` (External Data)
For importing observational or theoretical velocity distributions:

```cpp
class Distribution_Import : public Distribution {
    vector<Vector3d> imported_positions;
    vector<Vector3d> imported_velocities;
    
    void init(shared_ptr<Particle> p) override {
        int random_index = uniform_random_int(0, imported_positions.size() - 1);
        Vector3d pos = imported_positions[random_index];
        Vector3d vel = imported_velocities[random_index];
        p->init_particle(pos.x(), pos.y(), pos.z(), vel.x(), vel.y(), vel.z());
    }
};
```

## Collision Physics and Cross Section Implementation

### Current Code Status vs. Planned Implementation

**CRITICAL ACCURACY NOTE**: The extensive inelastic collision framework described in the sections below represents **PLANNED IMPLEMENTATION** based on literature best practices. Currently, the Corona3D codebase implements **ONLY elastic collision physics**. 

**Current Implementation Status**:
- ✅ **IMPLEMENTED**: Elastic collisions with energy-dependent total cross sections σ(E)
- ✅ **IMPLEMENTED**: Angular differential cross sections dσ/dΩ(E,θ) for realistic scattering
- ✅ **IMPLEMENTED**: Center-of-mass collision dynamics with conservation laws
- ✅ **IMPLEMENTED**: O-O, O-CO₂, O-CO, O-N₂ collision systems (elastic only)
- ❌ **NOT IMPLEMENTED**: InelasticCollisionHandler class
- ❌ **NOT IMPLEMENTED**: StateResolvedCollisionManager class  
- ❌ **NOT IMPLEMENTED**: State-resolved inelastic collision channels
- ❌ **NOT IMPLEMENTED**: Energy transfer to internal molecular modes

The detailed inelastic collision algorithms and classes described below are part of the development roadmap outlined in `AGENTS.md` but do not exist in the current codebase.

**Current Capabilities**:
- ✅ Elastic collisions with energy-dependent total cross sections σ(E)
- ✅ Angular differential cross sections dσ/dΩ(E,θ) for realistic scattering
- ✅ Center-of-mass collision dynamics with exact conservation laws
- ✅ O-O, O-CO₂, O-CO, O-N₂ collision systems (elastic only)
- ✅ Maxwell-Boltzmann thermal velocity distributions for background species
- ✅ Exponential atmospheric density profiles with scale heights
- ✅ Gravitational trajectory integration with leapfrog scheme

**Planned Implementation** (described in detail below):
- 🔄 State-resolved inelastic collision channels
- 🔄 Collision branching algorithm (elastic vs. inelastic selection)
- 🔄 Energy transfer to internal molecular modes (rotation, vibration)
- 🔄 Temperature-dependent thermal state populations
- 🔄 Detailed balance and conservation law validation

**Scientific Motivation for Inelastic Implementation**: Studies indicate that elastic-only treatments may underestimate energy loss, potentially leading to overestimated atmospheric escape rates. The planned inelastic implementation will provide more realistic physics.

### Cross Section Data Structure

The model implements a comprehensive collision physics framework using quantum mechanically calculated cross sections:

#### 1. Total Cross Sections σ_total(E)
Energy-dependent total cross sections determine collision probability:
```cpp
class Background_Species {
    vector<Interpolator*> sigma_interpolators;  // One per species
    
    double get_total_cross_section(int species, double energy_eV) {
        return sigma_interpolators[species]->linear_interpolate(energy_eV);
    }
};
```

**Data Format**: Two-column CSV files
```csv
# Column 1: Collision energy (eV)
# Column 2: Total cross section (cm²)
energy_eV,cross_section_cm2
0.030000,3.571670e-14
0.153740,2.189012e-14
...
```

#### 2. Differential Cross Sections dσ/dΩ(E,θ)
Angular-dependent cross sections for realistic scattering distributions:

```cpp
class Background_Species {
    // 4D array: [species][energy_index][0=CDF_values, 1=angles][angle_index]
    vector<vector<vector<vector<double>>>> differential_cross_section_CDFs;
    
    double sample_scattering_angle(int species, double energy_eV) {
        int energy_index = find_nearest_energy_index(species, energy_eV);
        double random_value = uniform_random_0_to_1();
        
        // Binary search through cumulative distribution function
        return interpolate_angle_from_CDF(species, energy_index, random_value);
    }
};
```

**Physical Significance**: Angular distributions encode the quantum mechanical scattering physics:
- **Forward scattering** (θ ≈ 0°): Grazing collisions, small energy transfer
- **Backward scattering** (θ ≈ 180°): Head-on collisions, large energy transfer  
- **Intermediate angles**: Mixed elastic/inelastic character

#### 3. Inelastic Cross Sections: Literature-Based Implementation Plan

**Current Status**: The codebase currently implements **elastic-only** collision physics. Inelastic collision channels represent a critical physics enhancement that will significantly improve model accuracy.

**Scientific Motivation**: Studies indicate that elastic-only treatments may underestimate energy loss compared to full elastic + inelastic models. Inelastic collisions transfer kinetic energy to internal molecular modes (rotation, vibration), leading to enhanced thermalization and more realistic atmospheric escape predictions.

**Implementation Plan**: State-Resolved Inelastic Collision Physics

The implementation follows established Monte Carlo practices from the quantum scattering literature:
<!-- (Balakrishnan & Dalgarno 2001; Cecchi-Pestellini et al. 2009) -->

**Step 1: Collision Channel Branching Algorithm**

```cpp
// Monte Carlo collision processing framework
enum CollisionType { ELASTIC, INELASTIC_ROTATIONAL, INELASTIC_VIBRATIONAL };

struct CollisionOutcome {
    CollisionType type;
    double energy_transfer_eV;     // ΔE = E_kinetic_final - E_kinetic_initial
    double scattering_angle_rad;   // Lab frame scattering angle
    int final_v, final_J;          // Target molecule final quantum state
};

class InelasticCollisionHandler {
public:
    CollisionOutcome process_collision(double collision_energy_eV, 
                                     double background_temperature_K,
                                     const string& target_species) {
        
        // 1. Calculate thermal state populations (Boltzmann distribution)
        auto state_populations = calculate_thermal_populations(background_temperature_K, target_species);
        
        // 2. Determine energetically accessible collision channels
        vector<CollisionChannel> open_channels;
        vector<double> weighted_cross_sections;
        
        // Loop over all possible initial states (thermally populated)
        for (auto& initial_state : state_populations) {
            int v_i = initial_state.v;
            int J_i = initial_state.J;
            double population = initial_state.boltzmann_factor;
            
            // Loop over all possible final states
            for (auto& final_state : get_accessible_states(collision_energy_eV, v_i, J_i)) {
                int v_f = final_state.v;
                int J_f = final_state.J;
                
                double energy_threshold = get_internal_energy_difference(v_i, J_i, v_f, J_f);
                
                if (collision_energy_eV >= energy_threshold) {
                    double sigma_partial = get_state_to_state_cross_section(v_i, J_i, v_f, J_f, collision_energy_eV);
                    double thermal_weight = population * sigma_partial;
                    
                    CollisionChannel channel{v_i, J_i, v_f, J_f, energy_threshold, sigma_partial};
                    open_channels.push_back(channel);
                    weighted_cross_sections.push_back(thermal_weight);
                }
            }
        }
        
        // 3. Add elastic channel (v_i,J_i → v_i,J_i for all thermal states)
        double sigma_elastic_total = calculate_thermal_elastic_cross_section(collision_energy_eV, state_populations);
        double sigma_total = sigma_elastic_total + std::accumulate(weighted_cross_sections.begin(), 
                                                                 weighted_cross_sections.end(), 0.0);
        
        // 4. Stochastic channel selection (branching)
        double random_selector = uniform_random() * sigma_total;
        
        if (random_selector <= sigma_elastic_total) {
            // Elastic collision: energy and internal states unchanged
            double theta = sample_elastic_scattering_angle(collision_energy_eV, target_species);
            return CollisionOutcome{ELASTIC, 0.0, theta, -1, -1};
        } else {
            // Inelastic collision: select specific state-to-state transition
            double cumulative_sigma = sigma_elastic_total;
            for (size_t i = 0; i < open_channels.size(); ++i) {
                cumulative_sigma += weighted_cross_sections[i];
                if (random_selector <= cumulative_sigma) {
                    CollisionChannel& selected = open_channels[i];
                    
                    double energy_transfer = get_internal_energy_difference(selected.v_initial, selected.J_initial,
                                                                          selected.v_final, selected.J_final);
                    double theta = sample_inelastic_scattering_angle(selected, collision_energy_eV);
                    
                    CollisionType type = (selected.v_initial != selected.v_final) ? 
                                       INELASTIC_VIBRATIONAL : INELASTIC_ROTATIONAL;
                    
                    return CollisionOutcome{type, energy_transfer, theta, selected.v_final, selected.J_final};
                }
            }
        }
    }
};
```

**Step 2: Data Requirements and Sources**

1. **State-to-State Cross Sections**: σ(v,J → v',J') in eV
   - **O-CO₂**: Gacesa et al. (2020) quantum scattering calculations
     - Repository: https://github.com/mgacesa66/O-CO2_cross-sections
     - Includes full rovibrational state resolution for ground electronic state
   - **O-CO**: Kumar et al. 2022 multi-surface calculations on 3A', 3A'', and 2-3A'' potential energy surfaces
     - Repository: https://github.com/mgacesa66/Cross-sections-O-CO
   - **O-N₂**: Kumar et al. (2023) high-level quantum chemistry with rotational coupling
     - Repository: https://github.com/snchtchhbr/n2_o_cross_section
   - **O-O**: Kharchenko et al. (2000)
        <!-- - elastic + inelastic channels -->

2. **Molecular Energy Level Data**:
   - **CO₂**: Rovibrational energy levels E(v₁,v₂,v₃,J) in cm⁻¹ or eV
   - **CO**: Vibrational levels E(v) and rotational constants B(v)
   - **N₂**: Dunham coefficients for accurate energy level calculation

3. **Angular Differential Cross Sections**: dσ/dΩ(E,θ) for each inelastic channel
   - Forward/backward scattering preferences depend on collision mechanism
   - Large energy transfers typically correlate with large scattering angles

**Step 3: Conservation Laws and Physics Validation**

```cpp
// Energy conservation check for inelastic collisions
void validate_inelastic_collision(const CollisionOutcome& outcome, 
                                double E_kinetic_initial, double E_kinetic_final) {
    double energy_balance = (E_kinetic_initial - E_kinetic_final) - outcome.energy_transfer_eV;
    
    // Energy transferred to internal modes must equal kinetic energy loss
    assert(abs(energy_balance) < 1e-12 * E_kinetic_initial); 
    
    // Kinetic energy cannot become negative
    assert(E_kinetic_final >= 0.0);
    
    // Energy transfer cannot exceed available kinetic energy
    assert(outcome.energy_transfer_eV <= E_kinetic_initial);
}

// Detailed balance validation (microscopic reversibility)
void validate_detailed_balance(double temperature_K) {
    // Forward rate: σ(v,J → v',J') * f(v,J) 
    // Reverse rate: σ(v',J' → v,J) * f(v',J')
    // Must satisfy: σ_forward * exp(-E_initial/kT) = σ_reverse * exp(-E_final/kT)
    
    for (auto& transition : all_transitions) {
        double sigma_forward = transition.cross_section_forward;
        double sigma_reverse = transition.cross_section_reverse;
        double delta_E = transition.energy_difference_eV;
        
        double balance_ratio = (sigma_forward / sigma_reverse) * exp(delta_E / (k_boltzmann * temperature_K));
        assert(abs(balance_ratio - 1.0) < 0.01); // 1% tolerance for numerical accuracy
    }
}
```

**Step 4: Implementation Priorities**

1. **Phase 1**: O-CO₂ inelastic channels (highest scientific priority)
   - Implement vibrational excitation: (0,0,0) → (0,0,1), (1,0,0), (0,1,0)
   - Add rotational transitions within vibrational levels
   - Validate against Gacesa et al. (2020) rate coefficients

2. **Phase 2**: O-CO and O-N₂ channels  
   - Simpler diatomic molecules: fewer states, cleaner implementation
   - Test detailed balance and energy conservation

3. **Phase 3**: Multi-quantum transitions and combination bands
   - (0,0,0) → (0,0,2), (1,0,1), etc. for CO₂
   - Verify against experimental rate measurements where available

**Key Physics Considerations**:

1. **Detailed Balance**: Satisfies microscopic reversibility (forward/reverse rates consistent with equilibrium)
2. **Energy Conservation**: Total energy (kinetic + internal) rigorously conserved in each collision
3. **Angular-Energy Correlation**: Inelastic processes exhibit different angular distributions than elastic scattering
4. **Temperature Coupling**: Background thermal state populations directly affect collision dynamics
5. **Threshold Effects**: Inelastic channels open only above specific collision energies

**Expected Scientific Impact**: Incorporating inelastic collision physics will:
- Reduce predicted escape rates due to enhanced energy loss
- Enable study of atmospheric temperature effects on escape efficiency
- Provide more realistic coupling between hot and thermal atmospheric populations

#### Implementation Summary and Next Steps

**To implement the inelastic collision framework described above**, the following code changes are required:

1. **Extend Background_Species class** with inelastic collision handling methods
2. **Add molecular energy level databases** for CO₂, CO, and N₂ internal states  
3. **Import state-to-state cross section data** from Gacesa et al. and Kumar et al. repositories
4. **Implement thermal state population calculations** using Boltzmann distributions
5. **Add collision branching logic** to select between elastic and inelastic channels
6. **Update collision dynamics** to handle energy transfer to internal modes
7. **Implement conservation law verification** for energy and momentum
8. **Add configuration options** to enable/disable inelastic physics

**See `AGENTS.md` in the project root** for a detailed implementation plan.

### Collision Dynamics Implementation

#### Center-of-Mass Frame Calculations

The collision algorithm uses center-of-mass transformations to ensure exact conservation of energy and momentum. The actual implementation in `Particle.cpp` uses a rotation matrix approach:

```cpp
void Particle::do_collision(shared_ptr<Particle> target, double theta, double time, double planet_radius) {
    // 1. Store pre-collision velocities for diagnostics
    double v_before = get_total_v();
    
    // 2. Calculate masses and center-of-mass velocity
    double my_mass = get_mass();
    double targ_mass = target->get_mass();
    Matrix<double, 3, 1> targ_v = {target->get_vx(), target->get_vy(), target->get_vz()};
    
    // 3. Calculate center-of-mass velocity
    Matrix<double, 3, 1> vcm = (my_mass*velocity.array() + targ_mass*targ_v.array()) / (my_mass + targ_mass);
    
    // 4. Transform to center-of-mass frame
    Matrix<double, 3, 1> v1v = velocity.array() - vcm.array();  // particle 1 c-o-m velocity
    double v1 = sqrt(v1v[0]*v1v[0] + v1v[1]*v1v[1] + v1v[2]*v1v[2]);  // particle 1 c-o-m scalar velocity
    
    // 5. Apply scattering transformation using rotation matrix
    // Unit vector parallel to particle 1 velocity
    Matrix<double, 3, 1> r = velocity.array() / sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
    
    double alpha = atan2(velocity[1], velocity[0]);
    double phi = atan2(velocity[2], sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1]));
    double gamma = constants::twopi*common::get_rand();  // random azimuthal angle
    
    // 6. Construct post-collision velocity in CM frame
    Matrix<double, 3, 1> vp;
    vp[0] = v1*cos(alpha)*cos(phi-theta);
    vp[1] = v1*sin(alpha)*cos(phi-theta);
    vp[2] = v1*sin(phi-theta);
    
    // 7. Apply rotation matrix transformation
    Matrix<double, 3, 3> Rrg;
    double Cg = cos(gamma);
    double Sg = sin(gamma);
    double Vg = 1.0-Cg;
    
    Rrg(0, 0) = r[0]*r[0]*Vg+Cg;
    Rrg(0, 1) = r[0]*r[1]*Vg+r[2]*Sg;
    Rrg(0, 2) = r[0]*r[2]*Vg-r[1]*Sg;
    Rrg(1, 0) = r[0]*r[1]*Vg-r[2]*Sg;
    Rrg(1, 1) = r[1]*r[1]*Vg+Cg;
    Rrg(1, 2) = r[1]*r[2]*Vg+r[0]*Sg;
    Rrg(2, 0) = r[0]*r[2]*Vg+r[1]*Sg;
    Rrg(2, 1) = r[1]*r[2]*Vg-r[0]*Sg;
    Rrg(2, 2) = r[2]*r[2]*Vg+Cg;
    
    Matrix<double, 3, 1> vrel1 = Rrg * vp;
    
    // 8. Transform back to laboratory frame and update velocity
    velocity = vcm.array() + vrel1.array();
    
    // 9. Log collision event for traced particles
    if (traced) {
        v_after = get_total_v()*1e-5;
        double alt_in_km = 1e-5*(radius - planet_radius);
        collision_log.push_back(to_string(time) + "\t\t" + to_string(alt_in_km) + "\t" + 
                               target->get_name() + "\t" + to_string(theta * (180.0/constants::pi)) + 
                               "\t" + to_string(v_before*1e-5) + "\t" + to_string(v_after));
    }
}
```

#### Energy and Momentum Verification

The code includes built-in conservation checks:
```cpp
// Verify energy conservation (for debugging)
double E_initial = 0.5 * m1 * v1_initial.squaredNorm() + 0.5 * m2 * v2_initial.squaredNorm();
double E_final = 0.5 * m1 * v1_final.squaredNorm() + 0.5 * m2 * v2_final.squaredNorm();
assert(abs(E_final - E_initial) / E_initial < 1e-12);

// Verify momentum conservation
Vector3d p_initial = m1 * v1_initial + m2 * v2_initial;
Vector3d p_final = m1 * v1_final + m2 * v2_final;
assert((p_final - p_initial).norm() / p_initial.norm() < 1e-12);
```

### Cross Section Data Sources and Validation

#### Quantum Scattering Calculations
- **O-CO₂**: Gacesa et al. (2020) - Ab initio potential energy surfaces with statistical averaging
- **O-CO**: Kumar et al. (2022) - Multi-surface scattering calculations (3A', 3A'', 2-3A'')  
- **O-N₂**: Kumar et al. (2023) - High-level quantum chemistry with rotational coupling
- **O-O**: Kharchenko et al. (2000) - Validated against experimental measurements

#### Data Processing Pipeline
1. **Unit Conversion**: Original data (often in cm⁻¹, radians) → eV, degrees
2. **Statistical Averaging**: Multiple potential energy surfaces → Single averaged cross section
3. **CDF Generation**: Differential cross sections → Cumulative distribution functions for sampling
4. **Interpolation Setup**: Energy grids → Continuous functions via linear/log interpolation

#### Validation Against Experimental Data
The cross sections are benchmarked against:
- **Laboratory measurements**: Molecular beam scattering experiments
- **Atmospheric observations**: MAVEN in-situ measurements  
- **Analytical models**: Hard sphere and classical trajectory approximations

### Cross Section Data Format

All cross section files follow a standardized 2-column CSV format:

#### Total Cross Sections:
```csv
# Collision Energy (eV), Cross Section (cm²)
energy_eV,cross_section_cm2
0.030000,3.571670e-14
0.153740,2.189012e-14
...
```

#### Differential Cross Sections:
```csv
# Scattering Angle (degrees), Cross Section (cm²)
angle_degrees,cross_section_cm2
0.000000,4.257504e-11
0.100000,4.162266e-11
...
```

## Input Configuration

### Main Configuration (`corona3d_2020.cfg`)

```ini
# Simulation Parameters
num_testparts     10000        # Number of test particles
part_type         O            # Particle type (H, O, N2, CO, CO2)
dist_type         Hot_O        # Initial distribution type (Hot_O, Hot_H, MB, Import)
timesteps         1000000      # Number of simulation timesteps
dt                1.0          # Time step (seconds)

# Planet Properties (Mars or Venus)
planet_mass       6.39e23      # Mars mass (kg) - uncomment Venus values for Venus
planet_radius     3.39e6       # Mars radius (m)
#planet_mass       4.87e24      # Venus mass (kg)
#planet_radius     6.05e6       # Venus radius (m)
ref_height        200e5        # Reference height (cm)
ref_temp          200.0        # Reference temperature (K)

# Simulation Domain
sim_lower_bound   100e5        # Lower boundary (cm)
sim_upper_bound   1000e5       # Upper boundary (cm)

# Background Species (Mars or Venus configurations)
num_bgparts       4            # Number of background species
./inputs/O_Mars_HotH.cfg       # O configuration (Mars)
./inputs/N2_Mars_HotH.cfg      # N₂ configuration (Mars)
./inputs/CO_Mars_HotH.cfg      # CO configuration (Mars)
./inputs/CO2_Mars_HotH.cfg     # CO₂ configuration (Mars)
# Use Venus configurations for Venus simulations

# Output Configuration
output_dir        ./output/    # Output directory
print_status_freq 10000        # Status print frequency
```

### Species Configuration Files

Each background species has its own configuration file for both Mars and Venus:

```ini
# CO2_Mars_HotO.cfg (Mars configuration)
type                CO2
ref_dens           6.68e7                    # Reference density (cm⁻³)
total_sigma_file   ./inputs/collisions/...   # Total cross section file
num_diff_energies  41                        # Number of differential energies

# Energy grid (eV)
energy1    0.029756
energy2    0.153740
...

# Differential cross section files
energy1_file    ./inputs/collisions/.../DCS_01.csv
energy2_file    ./inputs/collisions/.../DCS_02.csv
...

# CO2_Venus_HotH.cfg (Venus configuration)
type                CO2
ref_dens           1.2e8                     # Reference density (cm⁻³) - Venus values
total_sigma_file   ./inputs/collisions/...   # Total cross section file
num_diff_energies  41                        # Number of differential energies
...
```

## Output and Analysis

### Primary Simulation Outputs

#### 1. Escape Statistics and Energy Distributions
```cpp
class Atmosphere {
    // Escape counters by hemisphere
    int day_escape_count;   // Particles escaping from sunward side  
    int night_escape_count; // Particles escaping from anti-sunward side
    
    // Energy Distribution Functions (EDFs) at specified altitudes
    // Dimensions: [hemisphere][altitude_index][energy_bin][cos_theta_bin]
    vector<vector<vector<vector<double>>>> stats_EDFs;
    
    void calculate_escape_statistics() {
        double day_escape_fraction = (double)day_escape_count / (double)total_particles;
        double night_escape_fraction = (double)night_escape_count / (double)total_particles;
        double total_escape_probability = day_escape_fraction + night_escape_fraction;
        
        // Convert to physical escape flux using production rate
        double global_production_rate = my_distribution->get_global_rate(); // particles/s
        double hemispherical_rate = global_production_rate / 2.0;
        double escape_flux = total_escape_probability * hemispherical_rate; // particles/s
    }
};
```

**Energy Distribution Function Analysis**:
- **Energy Resolution**: 0.05 eV bins (20 bins per eV)
- **Angular Resolution**: 0.01 cosine bins (cos θ from -1 to +1)
- **Altitude Coverage**: User-specified atmospheric levels
- **Physical Interpretation**: 
  - cos θ > 0: Upward-moving particles (potential escape candidates)
  - cos θ < 0: Downward-moving particles (bound trajectories)

#### 2. Trajectory and Collision Data

**Particle Trajectory Logs**:
```cpp
struct TrajectoryPoint {
    double time_seconds;
    Vector3d position_cm;
    Vector3d velocity_cm_per_s;
    double altitude_km;
    double energy_eV;
    string status;  // "ACTIVE", "ESCAPED", "THERMALIZED", "ABSORBED"
};

vector<TrajectoryPoint> particle_trajectory_log;
```

**Collision Event Records**:
```cpp
struct CollisionEvent {
    double time_seconds;
    double altitude_km;
    string collision_partner;      // "CO2", "CO", "N2", "O", etc.
    double scattering_angle_deg;
    double velocity_before_km_s;
    double velocity_after_km_s;
    double energy_transfer_eV;
};

vector<CollisionEvent> collision_log;
```

#### 3. Density and Column Integration

**Radial Density Profiles**:
```cpp
void Atmosphere::update_density_statistics(int particle_index) {
    double radius = my_particles[particle_index]->get_radius();
    double altitude_km = (radius - planet_radius) * 1e-5;
    
    // Bin particles by altitude (1 km resolution)
    int altitude_bin = (int)altitude_km;
    if (altitude_bin >= 0 && altitude_bin < max_altitude_bins) {
        stats_density_counts[altitude_bin]++;
    }
    
    // Calculate column density contribution
    double column_density_factor = 1.0 / (4π * radius * radius);  // Geometric factor
    stats_column_density += column_density_factor;
}
```

**Limb Observation Geometry**:
```cpp
void Atmosphere::calculate_limb_viewing_geometry() {
    // For spacecraft limb observations, calculate line-of-sight column densities
    for (int altitude_index = 0; altitude_index < num_altitudes; altitude_index++) {
        double tangent_altitude = EDF_altitudes[altitude_index] * 1e5; // km to cm
        double impact_parameter = planet_radius + tangent_altitude;
        
        // Abel transform for limb geometry
        for (particle : active_particles) {
            double x = particle->get_x();  // Distance from planet center
            double y = particle->get_y();
            double z = particle->get_z();
            double projected_distance = sqrt(y*y + z*z);
            
            if (abs(x - impact_parameter) < slab_width/2.0 && 
                projected_distance <= limb_width/2.0) {
                stats_limb_density[altitude_index]++;
            }
        }
    }
}
```

### Data Analysis and Visualization Tools

#### Python Analysis Scripts (in `plotting_scripts/`)

**1. Escape Probability Analysis** (`plot_escape_prob.py`):
```python
def calculate_escape_probabilities(simulation_output_directory):
    """Calculate energy-dependent escape probabilities."""
    # Load EDF data from simulation outputs
    edf_data = load_edf_files(simulation_output_directory)
    
    escape_probabilities = {}
    for altitude in edf_data.keys():
        for energy_bin in range(len(edf_data[altitude])):
            # Sum upward-moving particles (cos_theta > 0)
            upward_flux = np.sum(edf_data[altitude][energy_bin][100:])  # cos_theta > 0
            total_flux = np.sum(edf_data[altitude][energy_bin][:])
            
            if total_flux > 0:
                escape_prob = upward_flux / total_flux
                escape_probabilities[(altitude, energy_bin)] = escape_prob
    
    return escape_probabilities

def plot_energy_dependent_escape(escape_probabilities):
    """Generate escape probability vs energy plots."""
    energies = np.arange(0, 10.0, 0.05)  # 0.05 eV bins
    
    for altitude in unique_altitudes:
        plt.figure(figsize=(10, 6))
        escape_probs = [escape_probabilities.get((altitude, i), 0) 
                       for i in range(len(energies))]
        plt.plot(energies, escape_probs, label=f'{altitude} km')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Escape Probability')
        plt.title(f'Energy-Dependent Escape Probability at {altitude} km')
        plt.grid(True)
        plt.savefig(f'escape_prob_{altitude}km.png', dpi=300)
```

**2. Velocity Distribution Analysis** (`plot_limb_vel_dist.py`):
```python
def analyze_limb_velocity_distributions(edf_data):
    """Analyze velocity distributions for limb observation comparison."""
    # Convert energy and angle bins to velocity space
    for altitude, energy_angle_matrix in edf_data.items():
        velocity_distribution = np.zeros((len(velocity_bins), len(angle_bins)))
        
        for e_bin, energy_eV in enumerate(energy_bins):
            for a_bin, cos_theta in enumerate(angle_bins):
                # Convert energy to velocity: v = sqrt(2E/m)
                velocity_cm_s = sqrt(2.0 * energy_eV * erg_per_eV / oxygen_mass_g);
                velocity_km_s = velocity_cm_s * 1e-5;
                
                # Map to velocity-angle grid
                v_bin = find_velocity_bin(velocity_km_s);
                velocity_distribution[v_bin, a_bin] = energy_angle_matrix[e_bin, a_bin];
        
        # Generate limb brightness profiles
        plot_limb_brightness_profile(altitude, velocity_distribution);

def compare_with_observations(simulation_results, observational_data):
    """Compare simulation outputs with spacecraft observations."""
    # Load MAVEN, Venus Express, or other observational datasets
    obs_altitudes, obs_densities, obs_uncertainties = load_observational_data(observational_data);
    
    # Interpolate simulation results to observational altitudes
    sim_densities = interpolate_simulation_to_observations(simulation_results, obs_altitudes);
    
    # Statistical comparison
    chi_squared = calculate_chi_squared(sim_densities, obs_densities, obs_uncertainties);
    correlation_coefficient = calculate_correlation(sim_densities, obs_densities);
    
    return {'chi_squared': chi_squared, 'correlation': correlation_coefficient};
```

**3. Weighted Escape Flux Calculations** (`calculate_weighted_escape_estimate.py`):
```python
def calculate_global_escape_flux(escape_statistics, production_mechanisms):
    """Calculate planet-wide escape flux from simulation statistics."""
    
    # Combine escape probabilities with production rates
    total_escape_flux = 0.0;
    
    for mechanism, rate_s in production_mechanisms.items():
        mechanism_escape_prob = escape_statistics[mechanism]['total_escape_probability'];
        mechanism_flux = rate_s * mechanism_escape_prob;
        total_escape_flux += mechanism_flux;
        
        print(f'{mechanism}: {rate_s:.2e} s⁻¹ × {mechanism_escape_prob:.3f} = {mechanism_flux:.2e} s⁻¹');
    
    # Convert to standard units
    escape_flux_per_cm2_per_s = total_escape_flux / (4 * np.pi * planet_radius_cm**2);
    escape_rate_kg_per_s = total_escape_flux * oxygen_atomic_mass_kg;
    
    return {
        'total_flux_particles_per_s': total_escape_flux,
        'flux_per_cm2_per_s': escape_flux_per_cm2_per_s,
        'mass_loss_rate_kg_per_s': escape_rate_kg_per_s,
        'mass_loss_rate_kg_per_year': escape_rate_kg_per_s * seconds_per_year
    };
```

### Comparison with Observational Data

#### MAVEN Deep Dip Campaign Analysis (Mars)
```python
def analyze_maven_deep_dip_comparison(simulation_directory, maven_data_file):
    """Compare simulation results with MAVEN deep dip observations."""
    
    # Load MAVEN in-situ measurements
    maven_data = load_maven_observations(maven_data_file);
    
    # Extract simulation results at MAVEN observation altitudes
    simulation_densities = extract_simulation_at_altitudes(
        simulation_directory, maven_data['altitudes_km']);
    
    # Statistical analysis
    residuals = simulation_densities - maven_data['measured_densities'];
    relative_errors = residuals / maven_data['measured_densities'];
    
    # Generate comparison plots
    plt.figure(figsize=(12, 8));
    plt.errorbar(maven_data['altitudes_km'], maven_data['measured_densities'],
                yerr=maven_data['uncertainties'], label='MAVEN Observations',
                fmt='o', capsize=5);
    plt.plot(maven_data['altitudes_km'], simulation_densities, 
            'r-', linewidth=2, label='Corona3D Simulation');
    plt.xlabel('Altitude (km)');
    plt.ylabel('Hot O Density (cm⁻³)');
    plt.yscale('log');
    plt.legend();
    plt.title('Simulation vs MAVEN Deep Dip Observations');
    plt.grid(True, alpha=0.3);
    plt.savefig('maven_comparison.png', dpi=300, bbox_inches='tight');
    
    return {'residuals': residuals, 'relative_errors': relative_errors};
```

<!-- #### Venus Express Analysis (Venus)
```python
def analyze_venus_express_comparison(simulation_directory, venus_data_file):
    """Compare simulation results with Venus Express observations."""
    
    # Load Venus Express in-situ measurements
    venus_data = load_venus_express_observations(venus_data_file);
    
    # Extract simulation results at Venus Express observation altitudes
    simulation_densities = extract_simulation_at_altitudes(
        simulation_directory, venus_data['altitudes_km']);
    
    # Statistical analysis
    residuals = simulation_densities - venus_data['measured_densities'];
    relative_errors = residuals / venus_data['measured_densities'];
    
    # Generate comparison plots
    plt.figure(figsize=(12, 8));
    plt.errorbar(venus_data['altitudes_km'], venus_data['measured_densities'],
                yerr=venus_data['uncertainties'], label='Venus Express Observations',
                fmt='s', capsize=5);
    plt.plot(venus_data['altitudes_km'], simulation_densities, 
            'b-', linewidth=2, label='Corona3D Simulation');
    plt.xlabel('Altitude (km)');
    plt.ylabel('Hot H Density (cm⁻³)');
    plt.yscale('log');
    plt.legend();
    plt.title('Simulation vs Venus Express Observations');
    plt.grid(True, alpha=0.3);
    plt.savefig('venus_express_comparison.png', dpi=300, bbox_inches='tight');
    
    return {'residuals': residuals, 'relative_errors': relative_errors};
``` -->

#### Model Validation Metrics
- **Chi-squared goodness of fit**: Quantitative comparison with observational data
- **Correlation analysis**: Statistical relationship between model and measurements  
- **Energy spectrum validation**: Comparison of predicted vs observed velocity distributions
- **Seasonal/temporal variations**: Model predictions vs time-dependent observations
- **Planetary comparison**: Validation against Mars (MAVEN) observations

## Goals

1. **Recalculation of Hot O (and Hot H) Escape Rates**: Utilize new doubly differential elastic cross-sections (O-CO2, O-CO, O-N2) and revised MAVEN data to improve accuracy in escape rate calculations for Mars (and Venus in the future).
2. **Automated Escape Probability Calculations**: Develop tools to automate the calculation of escape probabilities for each MAVEN orbit (Mars) during inbound periapsis in deep-dip campaigns.
3. **Inclusion of Inelastic Collision Physics**: Integrate state-resolved inelastic collision physics and cross-sections to enhance the physical realism of atmospheric escape modeling for hot oxygen (and hot hydrogen in the future).

## Usage Guide

### Basic Simulation

1. **Compile the Code**:
   ```bash
   cd src/
   make
   ```

2. **Configure Simulation**:
   - Edit `corona3d_2020.cfg` for main parameters
   - Modify species `.cfg` files for collision data

3. **Run Simulation**:
   ```bash
   ./corona3d_2020
   ```

4. **Analyze Results**:
   ```bash
   cd plotting_scripts/
   python plot_escape_prob.py
   ```

### Advanced Configuration and Future Development

#### Custom Initial Distribution Implementation
```cpp
// Example: Adding a new production mechanism
class Distribution_Custom_Photochemistry : public Distribution {
private:
    vector<double> altitude_profile;
    vector<double> production_rates;
    double characteristic_energy_eV;
    
public:
    Distribution_Custom_Photochemistry(Planet p, double ref_h, double ref_T, 
                                     string profile_file, double energy) 
        : Distribution(p, ref_h, ref_T), characteristic_energy_eV(energy) {
        load_production_profile(profile_file);
    }
    
    void init(shared_ptr<Particle> p) override {
        // 1. Sample altitude from production rate profile
        double altitude = sample_from_production_CDF();
        double radius = my_planet.get_radius() + altitude;
        
        // 2. Sample spherical coordinates
        double phi = 2π * uniform_random();
        double cos_theta = 2.0 * uniform_random() - 1.0;
        Vector3d position = radius * Vector3d(
            sqrt(1 - cos_theta²) * cos(phi),
            sqrt(1 - cos_theta²) * sin(phi),
            cos_theta
        );
        
        // 3. Energy distribution with characteristic energy
        double energy_eV = sample_energy_distribution();
        double velocity = sqrt(2.0 * energy_eV * constants::ergev / p->get_mass());
        
        // 4. Isotropic velocity direction
        Vector3d velocity_direction = sample_isotropic_direction();
        Vector3d velocity_vector = velocity * velocity_direction;
        
        p->init_particle(position.x(), position.y(), position.z(),
                        velocity_vector.x(), velocity_vector.y(), velocity_vector.z());
    }
    
    double get_global_rate() override {
        // Integrate production rate over all altitudes
        return integrate_production_profile();
    }
};
```

#### New Collision Partner Integration
```cpp
// Example: Adding a new atmospheric species
class Particle_H2O : public Particle {
private:
    static constexpr double mass_amu = 18.015;  // Water molecular mass
    static constexpr double mass_grams = mass_amu * constants::amu;
    
public:
    double get_mass() const override { return mass_grams; }
    string get_name() const override { return "H2O"; }
    
    // Additional methods for H2O-specific physics
    double get_rotational_temperature() const { return rotational_temp; }
    double get_vibrational_temperature() const { return vibrational_temp; }
};

// Integration into Background_Species configuration
void Background_Species::add_new_species(string species_config_file) {
    // Parse new species configuration
    string species_type = parse_species_type(species_config_file);
    
    // Create appropriate particle instance
    shared_ptr<Particle> new_species = set_particle_type(species_type);
    bg_particles.push_back(new_species);
    
    // Load cross section data
    load_cross_section_data(species_config_file, bg_particles.size() - 1);
    
    // Update collision detection arrays
    resize_collision_arrays();
}
```

#### Inelastic Collision Implementation: Literature-Based Approach

**Theoretical Foundation**: Implementation follows the rigorous quantum mechanical framework established by Gacesa et al. (2020), Kumar et al. (2022, 2023), Kharchenko et al. (2000), and related computational studies for atmosphere-relevant collision systems.

**Core Algorithm**: The collision branching algorithm implements a three-stage process consistent with Monte Carlo collision theory:

```cpp
// Literature-based inelastic collision implementation
class StateResolvedCollisionManager {
private:
    struct CollisionChannel {
        string channel_identifier;           // e.g., "O+CO2(v=0,j=0)->O+CO2(v=1,j=2)"
        int v_initial, J_initial;           // Initial rovibrational quantum numbers
        int v_final, J_final;               // Final rovibrational quantum numbers
        double energy_threshold_eV;         // Channel opening energy
        double internal_energy_change_eV;   // Exact quantum energy difference
        Interpolator* partial_cross_section; // σ_channel(E) from ab initio calculations
        Interpolator* angular_distribution;  // dσ/dΩ(E,θ) for this specific channel
    };
    
    vector<vector<CollisionChannel>> species_channels;  // [species_index][channel_index]
    vector<Interpolator*> elastic_cross_sections;       // Reference elastic σ(E)
    
public:
    /**
     * Processes collision using rigorous branching algorithm
     * @param species_index Background species involved in collision
     * @param collision_energy_eV Center-of-mass collision energy
     * @param target_temperature_K Background temperature for thermal state weighting
     * @return CollisionResult containing channel type, energy transfer, and scattering angle
     */
    CollisionResult process_collision(int species_index, double collision_energy_eV,
                                    double target_temperature_K) {
        // Step 1: Calculate thermal state distribution of target molecule
        vector<double> boltzmann_weights = calculate_thermal_populations(species_index, target_temperature_K);
        
        // Step 2: Identify all energetically accessible inelastic channels
        vector<int> open_channel_indices;
        vector<double> weighted_partial_cross_sections;
        
        for (size_t ch = 0; ch < species_channels[species_index].size(); ++ch) {
            const auto& channel = species_channels[species_index][ch];
            
            if (collision_energy_eV >= channel.energy_threshold_eV) {
                double sigma_partial = channel.partial_cross_section->interpolate(collision_energy_eV);
                double thermal_factor = boltzmann_weights[get_state_index(channel.v_initial, channel.J_initial)];
                
                open_channel_indices.push_back(ch);
                weighted_partial_cross_sections.push_back(sigma_partial * thermal_factor);
            }
        }
        
        // Step 3: Calculate total cross section including elastic contribution
        double sigma_elastic = elastic_cross_sections[species_index]->interpolate(collision_energy_eV);
        double sigma_total_inelastic = std::accumulate(weighted_partial_cross_sections.begin(), 
                                                     weighted_partial_cross_sections.end(), 0.0);
        double sigma_total = sigma_elastic + sigma_total_inelastic;
        
        // Step 4: Stochastic channel selection weighted by cross sections
        double random_selector = uniform_random_0_to_1() * sigma_total;
        
        // Check for elastic collision
        if (random_selector <= sigma_elastic) {
            double theta_elastic = sample_elastic_scattering_angle(species_index, collision_energy_eV);
            return CollisionResult{ELASTIC, 0.0, theta_elastic};
        }
        
        // Select inelastic channel
        double cumulative_cross_section = sigma_elastic;
        for (size_t i = 0; i < open_channel_indices.size(); ++i) {
            cumulative_cross_section += weighted_partial_cross_sections[i];
            if (random_selector <= cumulative_cross_section) {
                int selected_channel = open_channel_indices[i];
                const auto& channel = species_channels[species_index][selected_channel];
                
                double theta_inelastic = sample_inelastic_scattering_angle(selected, collision_energy_eV);
                return CollisionResult{INELASTIC, channel.internal_energy_change_eV, theta_inelastic};
            }
        }
        
        // Fallback to elastic (should not occur with proper normalization)
        logger->warn("Channel selection fallback to elastic collision");
        return CollisionResult{ELASTIC, 0.0, sample_elastic_scattering_angle(species_index, collision_energy_eV)};
    }
    
    /**
     * Calculates Boltzmann thermal population distribution for rovibrational states
     */
    vector<double> calculate_thermal_populations(int species_index, double temperature_K) {
        vector<double> populations(get_max_states(species_index), 0.0);
        double kT_eV = constants::boltzmann_eV * temperature_K;
        double partition_function = 0.0;
        
        // Calculate unnormalized populations and partition function
        for (int v = 0; v <= v_max[species_index]; ++v) {
            for (int J = 0; J <= J_max[species_index][v]; ++J) {
                double energy_eV = get_rovibrational_energy(species_index, v, J);
                double statistical_weight = (2*J + 1);  // Rotational degeneracy
                double boltzmann_factor = statistical_weight * exp(-energy_eV / kT_eV);
                
                int state_index = get_state_index(v, J);
                populations[state_index] = boltzmann_factor;
                partition_function += boltzmann_factor;
            }
        }
        
        // Normalize populations
        for (auto& pop : populations) {
            pop /= partition_function;
        }
        
        return populations;
    }
};

// Enhanced post-collision dynamics with proper energy-momentum conservation
void Particle::execute_inelastic_collision(shared_ptr<Particle> target, 
                                          double scattering_angle_rad,
                                          double internal_energy_transfer_eV,
                                          double simulation_time) {
    
    // Convert energy transfer to ergs
    double internal_energy_transfer_erg = internal_energy_transfer_eV * constants::ergev;
    
    // Store initial state for conservation checks
    Vector3d initial_momentum = get_mass() * velocity + target->get_mass() * target->velocity;
    double initial_kinetic_energy = 0.5 * get_mass() * velocity.squaredNorm() +
                                   0.5 * target->get_mass() * target->velocity.squaredNorm();
    
    // Apply elastic scattering kinematics first
    perform_elastic_collision_kinematics(target, scattering_angle_rad);
    
    // Calculate post-elastic kinetic energy
    double post_elastic_kinetic_energy = 0.5 * get_mass() * velocity.squaredNorm() +
                                        0.5 * target->get_mass() * target->velocity.squaredNorm();
    
    // Apply inelastic energy loss while conserving momentum
    double final_kinetic_energy = post_elastic_kinetic_energy - internal_energy_transfer_erg;
    
    if (final_kinetic_energy < 0) {
        logger->warn("Inelastic collision: insufficient kinetic energy for internal excitation");
        final_kinetic_energy = 0.1 * post_elastic_kinetic_energy;  // Retain small fraction
    }
    
    // Scale velocities to achieve correct final kinetic energy
    double energy_scaling_factor = sqrt(final_kinetic_energy / post_elastic_kinetic_energy);
    velocity *= energy_scaling_factor;
    target->velocity *= energy_scaling_factor;
    
    // Apply momentum conservation correction
    Vector3d final_momentum = get_mass() * velocity + target->get_mass() * target->velocity;
    Vector3d momentum_error = initial_momentum - final_momentum;
    double total_mass = get_mass() + target->get_mass();
    
    Vector3d momentum_correction = momentum_error / total_mass;
    velocity += momentum_correction * (target->get_mass() / total_mass);
    target->velocity += momentum_correction * (get_mass() / total_mass);
    
    // Update collision statistics
    increment_inelastic_collision_count();
    record_energy_transfer(internal_energy_transfer_eV);
}
```

**Implementation Validation Strategy**:

1. **Energy-Momentum Conservation Tests**: Verify that total energy and momentum are conserved to machine precision
2. **Detailed Balance Verification**: Confirm that forward/reverse collision rates satisfy microscopic reversibility  
3. **Cross Section Integration**: Validate that partial cross sections sum to total cross sections from literature
4. **Temperature Dependence**: Test Boltzmann distribution calculations against analytical results
5. **Benchmark Comparisons**: Compare escape rates with and without inelastic channels against published results

**Literature Data Integration**:

- **O-CO₂ system**: State-to-state cross sections from Gacesa et al. (2020) quantum calculations
- **O-N₂ system**: Rotationally-resolved data from Kumar et al. (2023)
- **O-CO system**: Vibrational excitation cross sections from Kumar et al. (2022)
- **Energy level data**: NIST/HITRAN rovibrational constants for accurate internal energy calculations

**Performance Considerations**:

- Pre-compute thermal state distributions for standard atmospheric temperature profiles
- Cache frequently accessed cross section interpolations
- Optimize channel selection algorithm for systems with many open channels
- Implement parallel collision processing with thread-safe random number generation

### Performance Optimization and Parallel Processing

#### OpenMP Parallelization
```cpp
// Parallel particle loop in Atmosphere::run_simulation()
void Atmosphere::run_simulation_parallel(double dt, int num_steps, ...) {
    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < active_parts; j++) {
        int particle_index = active_indices[j];
        
        // Thread-local random number generators
        thread_local static mt19937 local_rng(omp_get_thread_num() + base_seed);
        
        // Update particle dynamics
        my_parts[particle_index]->do_timestep(dt, gravitational_constant);
        
        // Check collisions (thread-safe)
        bool collision_occurred = false;
        #pragma omp critical(collision_check)
        {
            collision_occurred = bg_species.check_collision(my_parts[particle_index], dt);
        }
        
        if (collision_occurred) {
            #pragma omp critical(collision_execution)
            {
                my_parts[particle_index]->do_collision(
                    bg_species.get_collision_target(),
                    bg_species.get_collision_theta(),
                    current_time,
                    my_planet.get_radius()
                );
            }
        }
        
        // Thread-local deactivation tracking
        #pragma omp critical(deactivation_update)
        {
            check_deactivation_conditions(particle_index);
        }
    }
}
```

#### Memory-Efficient Trajectory Storage
```cpp
// Compressed trajectory storage for large simulations
class CompressedTrajectoryStorage {
private:
    struct TrajectoryPoint {
        float time_s;           // 4 bytes vs 8 for double
        float position_km[3];   // Position in km for better precision
        float velocity_km_s[3]; // Velocity in km/s
        uint16_t status_code;   // Encoded status
    };
    
    vector<vector<TrajectoryPoint>> particle_trajectories;
    int compression_factor = 10;  // Store every 10th timestep
    
public:
    void add_trajectory_point(int particle_id, double time, Vector3d pos, Vector3d vel, string status) {
        if (time_step_counter % compression_factor == 0) {
            TrajectoryPoint point{
                .time_s = (float)time,
                .position_km = {(float)(pos.x()*1e-5), (float)(pos.y()*1e-5), (float)(pos.z()*1e-5)},
                .velocity_km_s = {(float)(vel.x()*1e-5), (float)(vel.y()*1e-5), (float)(vel.z()*1e-5)},
                .status_code = encode_status(status)
            };
            particle_trajectories[particle_id].push_back(point);
        }
    }
};
```

#### Adaptive Time Stepping
```cpp
// Adaptive time step based on collision frequency and dynamics
double Atmosphere::calculate_adaptive_timestep(shared_ptr<Particle> p) {
    // Base time step from configuration
    double base_dt = configuration_timestep;
    
    // Adjust based on collision probability
    double collision_rate = bg_species.calculate_total_collision_rate(p);
    double collision_timescale = 1.0 / collision_rate;
    
    // Adjust based on gravitational dynamics
    double orbital_period = 2π * sqrt(pow(p->get_radius(), 3) / (constants::G * my_planet.get_mass()));
    double dynamical_timescale = orbital_period / 100.0;  // 1% of orbital period
    
    // Take minimum of all relevant timescales
    double adaptive_dt = std::min({base_dt, collision_timescale/10.0, dynamical_timescale});
    
    // Ensure numerical stability
    return std::max(adaptive_dt, minimum_allowed_timestep);
}
```

## Physics Validation and Benchmarking

### Validation Against Experimental Data

#### Laboratory Cross Section Measurements
The model's collision cross sections are validated against experimental measurements:

**Molecular Beam Scattering Experiments**:
- **O-O collisions**: Kharchenko et al. (2000) cross sections validated against beam experiments
- **O-CO₂ collisions**: Gacesa et al. (2020) ab initio calculations benchmarked vs experimental total cross sections
- **Energy range validation**: Cross sections tested from thermal (0.03 eV) to superthermal (5+ eV) energies

**Validation Metrics**:
```cpp
// Cross section validation against experimental data
double validate_cross_sections(string experimental_data_file, string model_data_file) {
    vector<double> exp_energies, exp_cross_sections, exp_uncertainties;
    vector<double> model_energies, model_cross_sections;
    
    load_experimental_data(experimental_data_file, exp_energies, exp_cross_sections, exp_uncertainties);
    load_model_data(model_data_file, model_energies, model_cross_sections);
    
    // Interpolate model to experimental energy grid
    vector<double> interpolated_model = interpolate_to_grid(model_energies, model_cross_sections, exp_energies);
    
    // Calculate chi-squared goodness of fit
    double chi_squared = 0.0;
    for (int i = 0; i < exp_energies.size(); i++) {
        double residual = interpolated_model[i] - exp_cross_sections[i];
        chi_squared += (residual * residual) / (exp_uncertainties[i] * exp_uncertainties[i]);
    }
    
    return chi_squared / (exp_energies.size() - 1);  // Reduced chi-squared
}
```

#### Spacecraft Observation Comparisons

**MAVEN Mission Validation**:
- **Deep Dip Campaigns**: Direct comparison with in-situ hot O measurements at 130-200 km altitude
- **Limb Observations**: Validation of column density predictions vs MAVEN IUVS observations
- **Energy Spectra**: Comparison of predicted vs observed velocity distributions

```cpp
// MAVEN comparison analysis
struct MavenComparisonResult {
    double altitude_km;
    double observed_density_cm3;
    double observed_uncertainty;
    double model_density_cm3;
    double relative_error;
    double chi_squared_contribution;
};

vector<MavenComparisonResult> validate_against_maven(string maven_data_file, 
                                                   string simulation_output_dir) {
    // Load MAVEN Deep Dip observations
    auto maven_data = load_maven_observations(maven_data_file);
    
    // Extract simulation results at MAVEN altitudes
    auto simulation_data = load_simulation_densities(simulation_output_dir);
    
    vector<MavenComparisonResult> comparison_results;
    
    for (auto& obs : maven_data) {
        double model_density = interpolate_simulation_density(simulation_data, obs.altitude_km);
        double relative_error = (model_density - obs.density_cm3) / obs.density_cm3;
        double chi_sq_contrib = pow((model_density - obs.density_cm3) / obs.uncertainty, 2);
        
        comparison_results.push_back({
            .altitude_km = obs.altitude_km,
            .observed_density_cm3 = obs.density_cm3,
            .observed_uncertainty = obs.uncertainty,
            .model_density_cm3 = model_density,
            .relative_error = relative_error,
            .chi_squared_contribution = chi_sq_contrib
        });
    }
    
    return comparison_results;
}
```

#### Analytical Model Benchmarks

##### Hard Sphere Collision Model
```cpp
// Analytical hard sphere cross section for validation
double hard_sphere_cross_section(double radius1_cm, double radius2_cm) {
    double collision_radius = radius1_cm + radius2_cm;
    return π * collision_radius * collision_radius;  // cm²
}

// Compare with quantum mechanical cross sections
void validate_hard_sphere_limit() {
    double O_radius = 1.4e-8;   // cm, approximate atomic radius
    double CO2_radius = 2.0e-8; // cm, approximate molecular radius    
    double hard_sphere_sigma = hard_sphere_cross_section(O_radius, CO2_radius);
    double quantum_sigma_low_energy = interpolate_quantum_cross_section(0.03); // eV
    
    cout << "Hard sphere limit: " << hard_sphere_sigma << " cm²" << endl;
    cout << "Quantum mechanical (0.03 eV): " << quantum_sigma_low_energy << " cm²" << endl;
    cout << "Ratio: " << quantum_sigma_low_energy / hard_sphere_sigma << endl;
}
```

##### Energy Conservation Tests
```cpp
// Rigorous energy conservation validation
void test_energy_conservation(int num_test_collisions = 10000) {
    double total_energy_error = 0.0;
    double max_energy_error = 0.0;
    
    for (int test = 0; test < num_test_collisions; test++) {
        // Create test particles with random velocities
        auto particle1 = make_shared<Particle_O>();
        auto particle2 = make_shared<Particle_CO2>();
        
        // Random initial velocities
        Vector3d v1_initial = generate_random_velocity(5e5);  // 5 km/s
        Vector3d v2_initial = generate_random_velocity(1e5);  // 1 km/s thermal
        particle1->init_particle_vonly(v1_initial.x(), v1_initial.y(), v1_initial.z());
        particle2->init_particle_vonly(v2_initial.x(), v2_initial.y(), v2_initial.z());
        
        // Calculate initial energy
        double E_initial = 0.5 * particle1->get_mass() * v1_initial.squaredNorm() +
                          0.5 * particle2->get_mass() * v2_initial.squaredNorm();
        
        // Perform collision
        double random_theta = π * uniform_random();
        particle1->do_collision(particle2, random_theta, 0.0, 3.39e8);
        
        // Calculate final energy
        double E_final = 0.5 * particle1->get_mass() * particle1->get_velocity().squaredNorm() +
                        0.5 * particle2->get_mass() * particle2->get_velocity().squaredNorm();
        
        // Check conservation
        double energy_error = abs(E_final - E_initial) / E_initial;
        total_energy_error += energy_error;
        max_energy_error = max(max_energy_error, energy_error);
    }
    
    cout << "Average energy conservation error: " << total_energy_error / num_test_collisions << endl;
    cout << "Maximum energy conservation error: " << max_energy_error << endl;
    assert(max_energy_error < 1e-12);  // Machine precision limit
}
```

##### Escape Velocity Validation
```cpp
// Test escape condition implementation
void validate_escape_calculations() {
    Planet mars;
    mars.init(6.39e26, 3.39e8);  // Mars mass and radius
    
    // Test particle at various altitudes
    vector<double> test_altitudes = {150e5, 200e5, 300e5, 500e5, 1000e5}; // cm
    
    for (double altitude : test_altitudes) {
        double radius = mars.get_radius() + altitude;
        double escape_velocity = sqrt(2.0 * constants::G * mars.get_mass() / radius);
        
        // Create test particle with escape velocity
        auto test_particle = make_shared<Particle_O>();
        test_particle->init_particle(radius, 0, 0, escape_velocity, 0, 0);
        
        // Check total energy
        double kinetic_energy = 0.5 * test_particle->get_mass() * escape_velocity * escape_velocity;
        double potential_energy = -constants::G * mars.get_mass() * test_particle->get_mass() / radius;
        double total_energy = kinetic_energy + potential_energy;
        
        cout << "Altitude: " << altitude*1e-5 << " km, ";
        cout << "Escape velocity: " << escape_velocity*1e-5 << " km/s, ";
        cout << "Total energy: " << total_energy << " ergs" << endl;
        
        // Total energy should be approximately zero for escape velocity
        assert(abs(total_energy) < 1e-6 * kinetic_energy);
    }
}
```

### Monte Carlo Statistical Validation

#### Convergence Testing
```cpp
// Test Monte Carlo convergence with increasing particle numbers
void test_monte_carlo_convergence() {
    vector<int> particle_numbers = {1000, 5000, 10000, 50000, 100000};
    vector<double> escape_probabilities;
    vector<double> statistical_uncertainties;
    
    for (int num_particles : particle_numbers) {
        // Run multiple simulations for statistics
        vector<double> escape_prob_samples;
        int num_runs = 10;
        
        for (int run = 0; run < num_runs; run++) {
            auto simulation_result = run_simulation_with_n_particles(num_particles);
            escape_prob_samples.push_back(simulation_result.escape_probability);
        }
        
        // Calculate mean and standard deviation
        double mean_escape_prob = calculate_mean(escape_prob_samples);
        double std_dev = calculate_standard_deviation(escape_prob_samples);
        double statistical_uncertainty = std_dev / sqrt(num_runs);
        
        escape_probabilities.push_back(mean_escape_prob);
        statistical_uncertainties.push_back(statistical_uncertainty);
        
        cout << "N = " << num_particles << ": ";
        cout << "Escape prob = " << mean_escape_prob << " ± " << statistical_uncertainty << endl;
    }
    
    // Check convergence (uncertainty should decrease as 1/sqrt(N))
    validate_statistical_convergence(particle_numbers, statistical_uncertainties);
}
```

### Benchmarking Against Lillis et al. (2017) MAVEN Analysis

```cpp
// Comprehensive comparison with published MAVEN analysis
struct LillisComparisonMetrics {
    double hot_O_scale_height_km;
    double exobase_density_cm3;
    double escape_flux_particles_cm2_s;
    double average_escape_energy_eV;
    vector<double> altitude_profile_km;
    vector<double> density_profile_cm3;
};

LillisComparisonMetrics compare_with_lillis_2017(string simulation_output_dir) {
    // Load Lillis et al. (2017) reference data
    auto lillis_data = load_lillis_reference_data();
    
    // Extract comparable quantities from simulation
    auto sim_data = analyze_simulation_output(simulation_output_dir);
    
    LillisComparisonMetrics comparison{
        .hot_O_scale_height_km = sim_data.calculate_scale_height(),
        .exobase_density_cm3 = sim_data.get_exobase_density(),
        .escape_flux_particles_cm2_s = sim_data.calculate_escape_flux(),
        .average_escape_energy_eV = sim_data.calculate_average_escape_energy(),
        .altitude_profile_km = sim_data.altitude_grid,
        .density_profile_cm3 = sim_data.density_profile
    };
    
    // Statistical comparison
    double scale_height_agreement = abs(comparison.hot_O_scale_height_km - lillis_data.scale_height_km) 
                                  / lillis_data.scale_height_uncertainty_km;
    
    cout << "Scale height agreement: " << scale_height_agreement << " sigma" << endl;
    
    return comparison;
}
```

## References and Technical Documentation

### Inelastic Collision Implementation Literature

**References for State-Resolved Collision Physics:**

1. **Gacesa, M., Lewkow, N., & Kharchenko, V.** (2020). "O(³P)+CO₂ collisions at hyperthermal energies: dynamics simulations and semiclassical rate coefficients for planetary aeronomy." *Monthly Notices of the Royal Astronomical Society*, 491(4), 5650-5664. DOI: 10.1093/mnras/stz3366
   - State-to-state inelastic cross sections for O-CO₂ rovibrational excitation
   - Quantum mechanical calculations on ab initio potential energy surfaces
   - **Data Repository**: https://github.com/mgacesa66/O-CO2_cross-sections

2. Kumar et al. (2022): Sanchit Kumar, Marko Gacesa, Malathe S Khalil, Amal Al Ghaferi, Nayla El-Kork, A quantum-mechanical investigation of O(3P) + CO scattering cross sections at superthermal collision energies, Monthly Notices of the Royal Astronomical Society, Volume 519, Issue 1, February 2023, Pages 1253–1260, https://doi.org/10.1093/mnras/stac3057

**O-CO Cross Sections**: https://github.com/mgacesa66/Cross-sections-O-CO
   - Multi-surface quantum scattering calculations
   - Energy range: 1.0×10⁻⁴ to 4.75 eV
   - 298 collision energies with angular distributions

3. Kumar et al. (2023): Sanchit Kumar, Sumit Kumar, Marko Gacesa, Nayla El-Kork, Sharma S R K C Yamijala, Quantum scattering cross-sections for O(3P) + N2 collisions for planetary aeronomy, Monthly Notices of the Royal Astronomical Society, Volume 526, Issue 4, December 2023, Pages 5675–5681, https://doi.org/10.1093/mnras/stad3149

 **O-N₂ Cross Sections**: https://github.com/snchtchhbr/n2_o_cross_section  
   - Rotationally resolved cross sections
   - Elastic and inelastic channels
   - Temperature-dependent collision dynamics

#### Data Format Specifications:
```
Total Cross Sections:
- Column 1: Collision Energy (eV)
- Column 2: Cross Section (cm²)
- Energy resolution: ~0.001-0.1 eV depending on collision system

Differential Cross Sections:
- Column 1: Scattering Angle (degrees)
- Column 2: Differential Cross Section (cm²/steradian)
- Angular resolution: 0.1-1.0 degrees
```

### Atmospheric Data Sources

#### Mars Atmospheric Profiles:
- **MAVEN Mission**: [https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/MAVEN/maven.html](https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/MAVEN/maven.html)
  - Neutral density profiles from NGIMS
  - Ion density profiles from STATIC and SWIA
  - Temperature profiles from LPW and NGIMS

- **Mars Climate Database**: [http://www-mars.lmd.jussieu.fr/mcd_python/](http://www-mars.lmd.jussieu.fr/mcd_python/)
  - Global circulation model outputs
  - Seasonal and diurnal variations
  - Dust storm impacts on atmospheric structure

#### Venus Atmospheric Profiles:
- **Venus Express Mission**: Neutral and ion density measurements
- **Pioneer Venus**: Historic atmospheric composition data
- **Akatsuki Mission**: Modern atmospheric dynamics observations

### Theoretical Background

#### Quantum Scattering Theory:
- **Born-Oppenheimer approximation**: Separation of electronic and nuclear motion
- **Close-coupling methods**: Inclusion of rotational and vibrational coupling
- **Statistical averaging**: Treatment of multiple potential energy surfaces

#### Classical Trajectory Methods:
- **Hamilton's equations**: Classical mechanics foundation
- **Symplectic integrators**: Conservation of phase space volume
- **Collision dynamics**: Center-of-mass transformations

#### Monte Carlo Statistical Mechanics:
- **Ergodic hypothesis**: Time averages equal ensemble averages
- **Detailed balance**: Microscopic reversibility in collision processes
- **Central limit theorem**: Convergence of statistical estimates

### Software Dependencies and Tools

#### Required Libraries:
```cpp
// Linear algebra and numerical methods
#include <eigen3/Eigen/Core>        // Matrix operations and vector algebra

// Standard C++ libraries
#include <vector>                   // Dynamic arrays
#include <memory>                   // Smart pointers for memory management
#include <random>                   // Mersenne Twister random number generation
#include <fstream>                  // File I/O operations
#include <algorithm>                // STL algorithms
```

**Build System**: 
- **Make**: Standard Unix build system (no automatic dependency tracking as claimed)
- **GCC**: Basic C++ compilation without OpenMP or advanced optimizations by default
- **Eigen3**: Linear algebra library with conditional includes for different platforms

#### Python Analysis Tools:
```python
# Required packages for data analysis
import numpy as np                  # Numerical arrays and operations
import matplotlib.pyplot as plt     # Plotting and visualization
import scipy.interpolate           # Interpolation and fitting
import scipy.optimize              # Optimization and curve fitting
import pandas as pd                # Data manipulation and analysis
```

### Code Documentation Standards

#### Header File Documentation:
```cpp
/**
 * @file Particle.hpp
 * @brief Base class for all particle types in corona3d simulation
 * @version 2.0
 * 
 * This header defines the abstract base class for all particle types.
 * Derived classes must implement get_mass() and get_name() methods.
 * The class handles collision dynamics, gravitational motion, and
 * particle lifecycle management.
 */
```

#### Function Documentation:
```cpp
/**
 * @brief Performs elastic collision between two particles
 * @param target Shared pointer to collision partner particle
 * @param theta Scattering angle in radians (0 = forward scattering)
 * @param time Current simulation time in seconds
 * @param planet_r Planet radius in cm for altitude calculations
 * 
 * @details Implements center-of-mass collision dynamics with exact
 * conservation of energy and momentum. The collision algorithm:
 * 1. Transforms to center-of-mass frame
 * 2. Applies scattering rotation by angle theta
 * 3. Transforms back to laboratory frame
 * 4. Updates both particle velocities
 * 
 * @note Only elastic collisions are currently implemented. Inelastic
 * channels are planned for future development.
 * 
 * @see Background_Species::check_collision() for collision detection
 * @see calc_collision_e() for collision energy calculation
 */
void do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r);
```

## Project Goals and Development Roadmap

### Primary Objectives

1. **Re-calculation of Hot O Escape Rates**: Re-do the Lillis et al. (2017) analysis using new doubly differential elastic cross-sections for O-CO₂, O-CO, and O-N₂ interactions, combined with revised MAVEN in-situ data to provide improved escape rate calculations

2. **Automated MAVEN Data Processing**: Write comprehensive output to file and automate the calculation of escape probabilities for each MAVEN orbit in-situ data during inbound periapsis passes and deep-dip campaigns

3. **Inelastic Collision Physics Integration**: Include state-resolved inelastic collision physics and cross-sections to provide more realistic energy transfer modeling and enhanced atmospheric escape predictions

### Scientific Impact

This enhanced model will provide:
- Updated atmospheric escape rate calculations using improved collision physics and latest MAVEN (Mars) dataset
- Systematic automated analysis of MAVEN observational data with new cross-section databases
- Quantitative assessment of the role of inelastic processes in hot atom thermalization for hot oxygen (and hot hydrogen in the future)
- Enhanced understanding of Mars (and Venus in the future) atmospheric evolution and current escape processes through improved physics implementation

---

*This documentation reflects the current state of the Corona3D model. The model continues to evolve with new cross section data, improved physics implementations, and inclusion of spacecraft observations.*
