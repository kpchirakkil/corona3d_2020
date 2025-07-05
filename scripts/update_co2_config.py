#!/usr/bin/env python3
"""
Script to update CO2_Mars_HotO.cfg to use the new CSV file paths in O-CO2_full.
"""

import re

def update_co2_config():
    """Update CO2_Mars_HotO.cfg to use new file paths in O-CO2_full."""
    
    config_file = "src/inputs/CO2_Mars_HotO.cfg"
    
    # Read the config file
    with open(config_file, 'r') as f:
        content = f.read()
    
    # Update differential cross section file paths
    # Replace: ./inputs/collisions/O-CO2/O-CO2_DCS-3pes_iEngXX.csv
    # With:    ./inputs/collisions/O-CO2_full/csv_format/differential_cross_sections/O-CO2_DCS-3pes_iEngXX.csv
    
    pattern = r'(\./inputs/collisions/)O-CO2(/O-CO2_DCS-3pes_iEng\d+\.csv)'
    replacement = r'\1O-CO2_full/csv_format/differential_cross_sections\2'
    
    updated_content = re.sub(pattern, replacement, content)
    
    # Write the updated content back to the file
    with open(config_file, 'w') as f:
        f.write(updated_content)
    
    print("Updated CO2_Mars_HotO.cfg to use new file paths in O-CO2_full")
    print("- Total cross section: O-CO2_full/csv_format/total_cross_sections/")
    print("- Differential cross sections: O-CO2_full/csv_format/differential_cross_sections/")

if __name__ == "__main__":
    update_co2_config()
