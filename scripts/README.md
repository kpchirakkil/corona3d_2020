# Collision Cross Section Conversion Scripts

This directory contains Python scripts for converting and managing collision cross section data files in the corona3d_2020 project.

## Scripts Overview

### Data Conversion Scripts

- **`convert_O_CO_files.py`** - Converts O-CO collision cross section files from original .dat format to standardized 2-column CSV format
- **`convert_O_N2_files.py`** - Converts O-N2 collision cross section files from original .dat format to standardized 2-column CSV format  
- **`verify_and_fix_co2.py`** - Verifies and fixes O-CO2 collision cross section files to ensure proper 2-column CSV format

### Configuration Update Scripts

- **`fix_co_config.py`** - Creates statistically averaged total sigma file for O-CO collisions and updates CO_Mars_HotO.cfg configuration
- **`fix_n2_config.py`** - Updates N2_Mars_HotO.cfg configuration to use correct energy-based CSV filenames
- **`update_co_config.py`** - Earlier version of CO configuration update script (kept for reference)

## Required File Format

All scripts ensure collision cross section files follow this standardized format:

### Total Cross Section Files
```csv
# Header comments
# Column 1: Collision energy (eV)
# Column 2: Cross section (cm²)
energy_eV,cross_section_cm2
1.000000e-04,1.234567e-14
...
```

### Differential Cross Section Files
```csv
# Header comments  
# Column 1: Angle (degrees)
# Column 2: Cross section (cm²)
angle_degrees,cross_section_cm2
0.000000,1.234567e-14
...
```

## Usage Instructions

### To convert new collision data:

1. Place original .dat files in the appropriate `*_full` directories
2. Run the conversion script for the specific collision type:
   ```bash
   python scripts/convert_O_CO_files.py
   python scripts/convert_O_N2_files.py
   python scripts/verify_and_fix_co2.py
   ```

### To update configuration files:

```bash
python scripts/fix_co_config.py
python scripts/fix_n2_config.py
python scripts/update_co2_config.py
```

### To verify file formats:

```bash
python scripts/verify_and_fix_co2.py
```

## Script Categories

### Active Conversion Tools
- `convert_O_CO_files.py` - Template for future O-CO data conversions
- `convert_O_N2_files.py` - Template for future O-N2 data conversions  
- `verify_and_fix_co2.py` - Active tool for O-CO2 data verification

### Configuration Management
- `fix_co_config.py` - Creates averaged total sigma files (reusable for data updates)
- `fix_n2_config.py` - Completed N2 configuration fix (reference for similar issues)

### Legacy/Reference
- `update_co_config.py` - Earlier CO configuration script (superseded by fix_co_config.py)

## Important Notes

1. **Backup Original Data**: All scripts preserve original .dat files - converted CSV files are created alongside or in `csv_format/` subdirectories

2. **Unit Conversions**: Scripts automatically handle unit conversions:
   - Energy: cm⁻¹ → eV (using factor 1.24e-4)
   - Angles: radians → degrees  
   - Cross sections: Various formats → cm²

3. **Statistical Averaging**: For O-CO data, the `fix_co_config.py` script creates a statistically averaged total cross section file combining data from multiple potential energy surfaces

4. **Energy-based Naming**: Differential cross section files use energy values in filenames for clarity

## Future Maintenance

Keep these scripts for:
- Converting new collision data as it becomes available (conversion templates)
- Verifying file formats after any manual edits (`verify_and_fix_co2.py`)
- Updating configuration files when energy grids change (`fix_co_config.py`)
- Regenerating averaged cross sections (`fix_co_config.py`)
- Reference for similar configuration issues (`fix_n2_config.py`, `update_co_config.py`)
- Documenting the exact conversion methodology used

## Organization Notes

All scripts are consolidated in this single `scripts/` directory for simplicity:
- **Active tools**: Use regularly for data management and verification
- **Templates**: Adapt for new collision types or data formats  
- **Reference/Legacy**: Older versions kept for documentation and troubleshooting
- **Configuration tools**: Scripts for managing .cfg files and generating averaged data

This unified structure eliminates confusion and maintains all tools in one accessible location.

## Dependencies

All scripts use only Python standard library modules:
- `os`, `csv`, `math` - No external dependencies required
