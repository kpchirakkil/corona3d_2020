# Archive of Data Processing Scripts

This directory contains Python scripts used for converting and standardizing collision cross section data files for the corona3d_2020 project.

## Script Status and Purpose

### Active/Maintenance Scripts
- **`fix_co_config.py`** - Creates statistically averaged total sigma files and updates CO configuration
  - **Keep**: May be needed if collision data is updated in the future
  - **Purpose**: Generates averaged cross sections from multiple potential energy surfaces

### Conversion Templates  
- **`convert_O_CO_files.py`** - Converts O-CO collision data from .dat to CSV format
  - **Keep**: Useful template for future O-CO data conversions
  - **Purpose**: Standardizes file format and units

- **`convert_O_N2_files.py`** - Converts O-N2 collision data from .dat to CSV format  
  - **Keep**: Useful template for future O-N2 data conversions
  - **Purpose**: Standardizes file format and units

### Completed/Obsolete Scripts
- **`fix_n2_config.py`** - Updates N2 configuration file paths
  - **Status**: Task complete, can be archived
  - **Purpose**: One-time fix for filename mismatches

- **`update_co_config.py`** - Basic CO configuration updates
  - **Status**: Superseded by fix_co_config.py, can be archived  
  - **Purpose**: Simple regex-based file path updates

## Usage Notes

1. **For new collision data**: Adapt the `convert_*.py` scripts as templates
2. **For data updates**: Use `fix_co_config.py` to regenerate averaged files
3. **For configuration issues**: Reference the archived scripts for patterns

## File Locations After Processing

- **O-CO data**: `src/inputs/collisions/O-CO_full/`
  - Total: `elastic_cross_sections/total_elastic_cross_section_averaged.csv`
  - Differential: `16O-CO_DCS-el_3pes/csv_format/DCS_elastic_{energy}_eV.csv`

- **O-N2 data**: `src/inputs/collisions/O-N2_full/`
  - Total: `n2_o_data/total_elastic_cross_section_isotope16.csv`  
  - Differential: `n2_o_data/DCS/isotope16/elastic_DCS_isotope_16/csv_format/DCS_elastic_{energy}_eV.csv`
