#!/usr/bin/env python3
"""
Script to update CO_Mars_HotO.cfg to use the new CSV format differential cross section files.
"""

import re

# Read the current configuration file
with open('src/inputs/CO_Mars_HotO.cfg', 'r') as f:
    content = f.read()

# Replace all energy file references to use the new CSV format
# Pattern matches: energyN_file      ./inputs/collisions/O-CO_full/16O-CO_DCS-el_3pes/O-CO_DCS-el_3pes_iE_N.dat
# Replace with: energyN_file      ./inputs/collisions/O-CO_full/16O-CO_DCS-el_3pes/csv_format/DCS_elastic_N_eV.csv

pattern = r'(energy\d+_file\s+)\.\/inputs\/collisions\/O-CO_full\/16O-CO_DCS-el_3pes\/O-CO_DCS-el_3pes_iE_(\d+)\.dat'
replacement = r'\1./inputs/collisions/O-CO_full/16O-CO_DCS-el_3pes/csv_format/DCS_elastic_\2_eV.csv'

updated_content = re.sub(pattern, replacement, content)

# Write the updated content back to the file
with open('src/inputs/CO_Mars_HotO.cfg', 'w') as f:
    f.write(updated_content)

print("Updated CO_Mars_HotO.cfg to use CSV format differential cross section files")
