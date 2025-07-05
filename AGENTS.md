# Corona3D 2020: Inelastic Collision Implementation Plan

The Corona3D code simulates **hot-H** and **hot-O** coronae and photochemical (non-thermal) escape at Mars and Venus. This document provides a comprehensive, literature-based implementation plan for integrating **state-resolved inelastic collision physics** into the Monte Carlo transport model.

## Executive Summary

Current elastic-only collision treatments may **underestimate energy loss** by 20-40% compared to models including inelastic channels (Gacesa et al. 2020). This document outlines a rigorous, step-by-step approach to implement quantum mechanically-derived, state-to-state inelastic cross sections that will significantly improve the physical realism of atmospheric escape calculations. The current Corona3D implementation includes sophisticated elastic collision physics with energy-dependent cross sections and differential scattering, but lacks the inelastic energy transfer mechanisms described in this plan.

## Project Goals

### Primary Objectives

1. **Re-calculation of Hot O Escape Rates**: Re-do the Lillis et al. (2017) analysis using new doubly differential elastic cross-sections for O-CO₂, O-CO, and O-N₂ interactions, combined with revised MAVEN in-situ data to provide improved escape rate calculations

2. **Automated MAVEN Data Processing**: Write comprehensive output to file and automate the calculation of escape probabilities for each MAVEN orbit in-situ data, including inbound periapsis passes and deep-dip campaigns

3. **Inelastic Collision Physics Integration**: Include state-resolved inelastic collision physics and cross-sections to provide more realistic energy transfer modeling and enhanced atmospheric escape predictions

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

## Implementation Approach

### Phase 1: Enhanced Cross-Section Physics

#### Collision Branching Framework
Implement the fundamental Monte Carlo branching algorithm:

1. **Collision Probability**: Use total cross section σ_total(E) to determine if collision occurs
2. **Channel Selection**: Weight elastic vs. inelastic channels by partial cross sections
3. **Outcome Sampling**: Sample energy transfer and scattering angle from appropriate differential cross sections

**Key Components**:
- `InelasticCollisionHandler` class in `Background_Species.cpp`
- Energy/momentum conservation validation
- Configuration flags for elastic-only vs. inelastic modes

#### Thermal State Population Calculator
Calculate Boltzmann thermal populations for target molecular states across different atmospheric conditions and temperatures.

### Phase 2: Data Integration and Automation

#### Cross-Section Data Integration
Integrate quantum mechanically-derived cross-section data for O-CO₂, O-CO, and O-N₂ interactions from literature sources:

- **O-CO₂**: https://github.com/mgacesa66/O-CO2_cross-sections
- **O-CO**: https://github.com/mgacesa66/Cross-sections-O-CO
- **O-N₂**: https://github.com/snchtchhbr/n2_o_cross_section

#### MAVEN Data Processing Automation
Develop comprehensive automation for MAVEN orbit analysis:
- Automated processing of inbound periapsis passes
- Deep-dip campaign data analysis
- Systematic file output generation for escape probability calculations
- Batch processing capabilities for multiple orbits

### Phase 3: Scientific Analysis and Validation

#### Hot O Escape Rate Recalculation
Recalculate escape rates using updated cross-sections and revised MAVEN data to improve upon the Lillis et al. (2017) analysis with enhanced physics.

#### Model Validation and Comparison
- Compare results with previous studies to validate improvements
- Assess the quantitative impact of inelastic collision physics
- Validate against MAVEN observational data

## Scientific Impact

This enhanced model will provide:
- Updated atmospheric escape rate calculations using improved collision physics and latest MAVEN datasets
- Systematic automated analysis of MAVEN observational data with new cross-section databases  
- Quantitative assessment of the role of inelastic processes in hot atom thermalization
- Enhanced understanding of Mars atmospheric evolution and current escape processes through improved physics implementation
