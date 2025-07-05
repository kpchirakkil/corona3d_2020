# Corona3D 2020: Inelastic Collision Implementation Plan

The Corona3D code simulates **hot-H** and **hot-O** coronae and photochemical (non-thermal) escape at Mars and Venus. This document provides a comprehensive, literature-based implementation plan for integrating **state-resolved inelastic collision physics** into the Monte Carlo transport model.

## Executive Summary

Current elastic-only collision treatments may **underestimate energy loss** by 20-40% compared to models including inelastic channels (Gacesa et al. 2020). This document outlines a rigorous, step-by-step approach to implement quantum mechanically-derived, state-to-state inelastic cross sections that will significantly improve the physical realism of atmospheric escape calculations.

## Scientific Motivation

### Physical Importance of Inelastic Collisions

1. **Enhanced Energy Loss**: Inelastic collisions transfer kinetic energy to internal molecular modes (rotation, vibration), leading to more efficient thermalization of hot atoms
2. **Coupling with Background Temperature**: Thermal populations of target molecule states directly affect collision probabilities and energy transfer rates
3. **Angular-Energy Correlation**: Inelastic processes show different angular dependencies than elastic scattering, affecting transport properties
4. **Quantitative Impact**: Studies suggest elastic-only models may overestimate escape rates by significant factors

### Literature Foundation

- **Gacesa et al. (2020)**: Quantum mechanical O-CO₂ cross sections with full rovibrational resolution (*MNRAS* 491, 5650)
- **Kharchenko et al. (2000)**: Comprehensive O-O elastic and inelastic cross section database (*JGR* 105, 24899)
- **Balakrishnan & Dalgarno (2001)**: Theoretical framework for Monte Carlo implementation of state-resolved collisions
- **Cecchi-Pestellini et al. (2009)**: Benchmark Monte Carlo methods for inelastic collision modeling
- **Lillis et al. (2017)**: MAVEN hot-O escape analysis for validation (*JGR* 122, 3815)

## Implementation Plan

### Phase 1: Core Algorithm Development

#### Task 1.1: Collision Branching Framework
**Priority**: HIGH | **Effort**: 3-4 weeks | **Assignee**: TBD

Implement the fundamental Monte Carlo branching algorithm:

1. **Collision Probability**: Use total cross section σ_total(E) to determine if collision occurs
2. **Channel Selection**: Weight elastic vs. inelastic channels by partial cross sections
3. **Outcome Sampling**: Sample energy transfer and scattering angle from appropriate differential cross sections

**Deliverables**:
- `InelasticCollisionHandler` class in `Background_Species.cpp`
- Unit tests validating energy/momentum conservation
- Configuration flags to toggle elastic-only vs. inelastic modes

#### Task 1.2: Thermal State Population Calculator
**Priority**: HIGH | **Effort**: 2 weeks | **Assignee**: TBD

Calculate Boltzmann thermal populations for target molecular states:

```cpp
vector<ThermalState> calculate_thermal_populations(double temperature_K, const string& species) {
    // For CO₂: sum over (v₁,v₂,v₃,J) states
    // For CO, N₂: sum over (v,J) states
    // Return normalized Boltzmann factors: exp(-E_internal/kT) / Z
}
```

**Data Requirements**:
- Molecular energy level databases for CO₂, CO, N₂
- Partition function calculations
- Temperature-dependent state cutoffs

#### Task 1.3: Conservation Law Validation
**Priority**: MEDIUM | **Effort**: 1 week | **Assignee**: TBD

Implement rigorous physics checks:
- Energy conservation: E_kinetic + E_internal = constant
- Momentum conservation in center-of-mass frame
- Detailed balance verification using microscopic reversibility

### Phase 2: Cross Section Data Integration

#### Task 2.1: O-CO₂ Data Import (Highest Priority)
**Priority**: HIGH | **Effort**: 2-3 weeks | **Assignee**: TBD

Import state-to-state cross sections from Gacesa et al. (2020):

**Data Source**: https://github.com/mgacesa66/O-CO2_cross-sections
**Key Channels**:
- Vibrational excitation: (0,0,0) → (0,0,1), (1,0,0), (0,1,0)
- Rotational transitions within each vibrational level
- Angular differential cross sections dσ/dΩ(E,θ)

**Implementation Steps**:
1. Parse CSV/HDF5 data files from Gacesa repository
2. Build interpolation grids for σ(v,J→v',J'; E)
3. Generate cumulative distribution functions for angular sampling
4. Validate against published rate coefficients

#### Task 2.2: O-CO and O-N₂ Data Integration
**Priority**: MEDIUM | **Effort**: 2 weeks each | **Assignee**: TBD

**Data Sources**:
- O-CO: https://github.com/mgacesa66/Cross-sections-O-CO
- O-N₂: https://github.com/snchtchhbr/n2_o_cross_section

Simpler diatomic targets enable thorough algorithm testing before tackling CO₂ complexity.

#### Task 2.3: O-O Inelastic Channels
**Priority**: LOW | **Effort**: 1 week | **Assignee**: TBD

Extend existing O-O elastic treatment (Kharchenko et al. 2000) to include:
- Fine structure transitions (³P → ¹D, ¹S)
- Spin-orbit coupling effects

### Phase 3: Validation and Benchmarking

#### Task 3.1: MAVEN Deep-Dip Validation
**Priority**: HIGH | **Effort**: 2-3 weeks | **Assignee**: TBD

Compare inelastic model predictions against MAVEN in-situ observations:

**Test Cases**:
- Periapsis altitudes 150-200 km
- Various solar zenith angles and EUV conditions
- Seasonal variations (solar longitude)

**Metrics**:
- Hot O density profiles
- Escape flux estimates
- Energy distribution functions

#### Task 3.2: Elastic-Only Comparison
**Priority**: HIGH | **Effort**: 1 week | **Assignee**: TBD

Quantify impact of inelastic physics:
- Run identical simulations with/without inelastic channels
- Document energy loss enhancement factors
- Analyze spatial distribution changes

#### Task 3.3: Literature Benchmarking
**Priority**: MEDIUM | **Effort**: 2 weeks | **Assignee**: TBD

Validate against other Monte Carlo/DSMC models:
- Compare rate coefficients at different temperatures
- Cross-check angular scattering distributions
- Verify detailed balance implementation

### Phase 4: Production Implementation

#### Task 4.1: Computational Optimization
**Priority**: MEDIUM | **Effort**: 2-3 weeks | **Assignee**: TBD

Optimize performance for production runs:
- Precompute thermal state populations
- Cache interpolation results
- Parallelize cross section lookups
- Profile memory usage with full state-resolved data

#### Task 4.2: Automated MAVEN Analysis
**Priority**: HIGH | **Effort**: 3-4 weeks | **Assignee**: TBD

Scale up to systematic MAVEN orbit analysis:
- Batch processing scripts for multiple orbits
- Statistical analysis of escape flux vs. environmental conditions
- Automated comparison with observational data

#### Task 4.3: User Interface and Documentation
**Priority**: MEDIUM | **Effort**: 2 weeks | **Assignee**: TBD

**Configuration Options**:
```bash
./corona3d --config mars_config.cfg --inelastic --species O --collision-data /path/to/cross_sections/
./corona3d --config venus_config.cfg --elastic-only  # Legacy mode for benchmarking
```

**Documentation Updates**:
- Updated `src/README.md` with inelastic collision physics
- User guide for cross section data preparation
- Example configuration files for all target species

## Risk Assessment and Mitigation

### Technical Risks

1. **Computational Performance**: State-resolved cross sections significantly increase memory and CPU requirements
   - *Mitigation*: Implement energy/angle binning strategies, optimize data structures
   
2. **Numerical Stability**: Conservation laws must be maintained to machine precision
   - *Mitigation*: Use double precision, implement strict validation checks
   
3. **Data Availability**: Some cross section data may have energy/angle coverage gaps
   - *Mitigation*: Develop extrapolation/interpolation protocols, document uncertainties

### Scientific Risks

1. **Model Validation**: Limited experimental data for validation at atmospheric conditions
   - *Mitigation*: Cross-validate against multiple literature sources, perform sensitivity analysis
   
2. **Quantum Classical Interface**: Quantum cross sections used in classical trajectory model
   - *Mitigation*: Follow established practices from molecular dynamics literature

## Success Metrics

### Technical Milestones
- [ ] Energy conservation to < 1e-12 relative precision
- [ ] Momentum conservation to < 1e-12 relative precision  
- [ ] Performance impact < 3x compared to elastic-only model
- [ ] Successful integration of all four collision systems (O-CO₂, O-CO, O-N₂, O-O)

### Scientific Validation
- [ ] Agreement with Gacesa et al. (2020) rate coefficients within 10%
- [ ] Detailed balance satisfaction to 1% precision
- [ ] MAVEN escape flux predictions within observational uncertainties
- [ ] Temperature dependence consistent with laboratory/theory expectations

## Timeline and Resource Requirements

### Timeline: 6-8 months total
- **Phase 1** (Algorithm): 6-7 weeks
- **Phase 2** (Data): 7-8 weeks  
- **Phase 3** (Validation): 5-6 weeks
- **Phase 4** (Production): 7-8 weeks

### Personnel Requirements
- **1 FTE Graduate Student/Postdoc**: Algorithm development, data integration
- **0.5 FTE Faculty/Senior Scientist**: Scientific oversight, validation
- **0.25 FTE Software Engineer**: Performance optimization, testing infrastructure

### Computational Resources
- **Development**: Standard workstation with 32+ GB RAM
- **Validation**: Small cluster access for MAVEN orbit ensembles  
- **Production**: HPC allocation for systematic atmospheric escape studies

## References

- Balakrishnan, N., & Dalgarno, A. (2001). *Phys. Rev. A* **63**, 012703
- Cecchi-Pestellini, C., et al. (2009). *ApJ* **703**, 1056
- Gacesa, M., Lewkow, N., & Kharchenko, V. (2020). *MNRAS* **491**, 5650. DOI: 10.1093/mnras/stz3366
- Kharchenko, V., et al. (2000). *JGR* **105**, 24899. DOI: 10.1029/2000JA000085  
- Lillis, R. J., et al. (2017). *JGR* **122**, 3815. DOI: 10.1002/2016JA023525

### Data Repositories
- O-CO₂: https://github.com/mgacesa66/O-CO2_cross-sections
- O-CO: https://github.com/mgacesa66/Cross-sections-O-CO  
- O-N₂: https://github.com/snchtchhbr/n2_o_cross_section
