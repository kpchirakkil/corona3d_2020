# Collision Cross Section File Format Summary

This document summarizes the standardized file formats for collision cross section data in the corona3d_2020 project.

## Overview

All collision cross section files have been standardized to use 2-column CSV format:
- **Total cross sections**: Column 1 = collision energy (eV), Column 2 = cross section (cm²)
- **Differential cross sections**: Column 1 = angle (degrees), Column 2 = cross section (cm²)

## O-CO Collision Data

### Total Cross Sections
- **Location**: `src/inputs/collisions/O-CO_full/elastic_cross_sections/`
- **Main file**: `total_elastic_cross_section_averaged.csv`
  - This is a **statistically averaged** total cross section file
  - Combines data from three potential energy surfaces: 2-3Adp, 3Adp, and 3Ap
  - Energy range: ~1.0e-4 to 4.75 eV
  - 300 data points with linear interpolation between surfaces
- **Individual surface files** (for reference):
  - `elastic_16_O_2-3Adp.csv` (converted from .dat format)
  - Original .dat files still available for comparison

### Differential Cross Sections
- **Location**: `src/inputs/collisions/O-CO_full/16O-CO_DCS-el_3pes/csv_format/`
- **Format**: `DCS_elastic_{energy:.6f}_eV.csv`
- **Energy range**: 0.000100 to 4.75 eV (298 energies)
- **Example files**:
  - `DCS_elastic_0.000100_eV.csv`
  - `DCS_elastic_0.000110_eV.csv`
  - `DCS_elastic_4.750000_eV.csv`

## O-N2 Collision Data

### Total Cross Sections
- **Location**: `src/inputs/collisions/O-N2_full/n2_o_data/`
- **File**: `total_elastic_cross_section_isotope16.csv`
- **Source**: Converted from original isotope16.dat file
- **Format**: 2-column CSV (energy in eV, cross section in cm²)

### Differential Cross Sections
- **Location**: `src/inputs/collisions/O-N2_full/n2_o_data/DCS/isotope16/elastic_DCS_isotope_16/csv_format/`
- **Format**: `DCS_elastic_{energy:.6f}_eV.csv`
- **Available energies** (10 files):
  - `DCS_elastic_0.309960_eV.csv` (from original elastic_2500.dat)
  - `DCS_elastic_0.619920_eV.csv` (from original elastic_5000.dat)
  - `DCS_elastic_1.249883_eV.csv` (from original elastic_10081.dat)
  - `DCS_elastic_1.499958_eV.csv` (from original elastic_12098.dat)
  - `DCS_elastic_2.249938_eV.csv` (from original elastic_18147.dat)
  - `DCS_elastic_2.749965_eV.csv` (from original elastic_22180.dat)
  - `DCS_elastic_3.249993_eV.csv` (from original elastic_26213.dat)
  - `DCS_elastic_3.500068_eV.csv` (from original elastic_28230.dat)
  - `DCS_elastic_3.749896_eV.csv` (from original elastic_30245.dat)
  - `DCS_elastic_3.999972_eV.csv` (from original elastic_32262.dat)

## Configuration Files Updated

### CO_Mars_HotO.cfg
- **Total sigma file**: Now references `total_elastic_cross_section_averaged.csv`
- **Differential files**: All 298 energy files use correct energy-based CSV filenames
- **Energy range**: 0.000100 to 4.75 eV

### N2_Mars_HotO.cfg
- **Total sigma file**: References `total_elastic_cross_section_isotope16.csv`
- **Differential files**: All 10 energy files use correct energy-based CSV filenames
- **Energy values match**: Configuration energies correspond to actual CSV file energies

## Key Benefits of Standardization

1. **Consistency**: All files use the same 2-column CSV format
2. **Statistical Averaging**: O-CO total cross sections now represent the average across all surfaces used for differential data
3. **Energy-based Naming**: Filenames clearly indicate the collision energy
4. **Proper Units**: All energies in eV, all cross sections in cm²
5. **Header Documentation**: Each file includes descriptive headers

## Conversion Scripts

- `convert_O_CO_files.py`: Converts O-CO cross section files
- `convert_O_N2_files.py`: Converts O-N2 cross section files
- `fix_co_config.py`: Creates averaged total sigma file and updates CO configuration
- `fix_n2_config.py`: Updates N2 configuration for correct file references

All original files are preserved for reference and comparison.
