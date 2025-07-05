#!/usr/bin/env python3
"""
MAINTENANCE SCRIPT - KEEP FOR FUTURE USE

Script to create a statistically averaged total sigma file for O-CO collisions
and update the CO_Mars_HotO.cfg configuration file.

This script should be retained because:
1. It may be needed if new O-CO collision data is provided
2. It generates the averaged total cross section from multiple surfaces
3. It ensures consistency between total and differential cross section data

Usage: Run this script whenever O-CO collision data is updated to regenerate
the statistically averaged total cross section file.

This script:
1. Reads the three O-CO total sigma files (2-3Adp, 3Adp, 3Ap surfaces)
2. Creates a statistically averaged file based on the three surfaces
3. Updates the CO configuration file to use the new averaged file

Author: Auto-generated script for corona3d_2020 project
"""

import os
import csv

def read_dat_file(filepath):
    """Read a .dat file and return energy and cross section arrays."""
    energies = []
    cross_sections = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # Energy in 1/cm, convert to eV (multiply by 1.24e-4)
                    energy_eV = float(parts[0]) * 1.24e-4
                    # Cross section in 1E-16 cm^2, convert to cm^2
                    cross_section_cm2 = float(parts[1]) * 1e-16
                    energies.append(energy_eV)
                    cross_sections.append(cross_section_cm2)
                except ValueError:
                    continue
    
    return energies, cross_sections

def interpolate(x_new, x_old, y_old):
    """Simple linear interpolation."""
    y_new = []
    for x in x_new:
        if x <= x_old[0]:
            y_new.append(y_old[0])
        elif x >= x_old[-1]:
            y_new.append(y_old[-1])
        else:
            # Find the two points to interpolate between
            for i in range(len(x_old) - 1):
                if x_old[i] <= x <= x_old[i + 1]:
                    # Linear interpolation
                    t = (x - x_old[i]) / (x_old[i + 1] - x_old[i])
                    y = y_old[i] + t * (y_old[i + 1] - y_old[i])
                    y_new.append(y)
                    break
    return y_new

def create_averaged_total_sigma():
    """Create statistically averaged total sigma file for O-CO."""
    
    base_dir = "src/inputs/collisions/O-CO_full/elastic_cross_sections"
    
    # File paths for the three surfaces
    files = {
        '2-3Adp': os.path.join(base_dir, 'elastic_16_O_2-3Adp.dat'),
        '3Adp': os.path.join(base_dir, 'elastic_16_O_3Adp.dat'), 
        '3Ap': os.path.join(base_dir, 'elastic_16_O_3Ap.dat')
    }
    
    # Read all three files
    data = {}
    for surface, filepath in files.items():
        energies, cross_sections = read_dat_file(filepath)
        # Sort by energy
        sorted_pairs = sorted(zip(energies, cross_sections))
        energies, cross_sections = zip(*sorted_pairs)
        data[surface] = {
            'energies': list(energies),
            'cross_sections': list(cross_sections)
        }
        print(f"Read {surface}: {len(energies)} data points")
    
    # Find common energy range
    all_energies = []
    for surface_data in data.values():
        all_energies.extend(surface_data['energies'])
    
    energy_min = min(all_energies)
    energy_max = max(all_energies)
    
    # Create a common energy grid with 300 points
    common_energies = []
    for i in range(300):
        energy = energy_min + i * (energy_max - energy_min) / 299
        common_energies.append(energy)
    
    # Interpolate each surface to the common grid
    interpolated_cross_sections = []
    
    for surface, surface_data in data.items():
        interp_cs = interpolate(common_energies, 
                               surface_data['energies'],
                               surface_data['cross_sections'])
        interpolated_cross_sections.append(interp_cs)
        print(f"Interpolated {surface} to common grid")
    
    # Calculate statistical average (mean of the three surfaces)
    avg_cross_sections = []
    for i in range(len(common_energies)):
        avg = sum(surface_cs[i] for surface_cs in interpolated_cross_sections) / 3
        avg_cross_sections.append(avg)
    
    # Create output directory if it doesn't exist
    output_file = os.path.join(base_dir, 'total_elastic_cross_section_averaged.csv')
    
    # Write the averaged data
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['# Collision Energy (eV)', 'Cross Section (cm^2)'])
        writer.writerow(['# Statistically averaged elastic cross section for O-CO collisions'])
        writer.writerow(['# Averaged over three potential energy surfaces: 2-3Adp, 3Adp, 3Ap'])
        
        for energy, cross_section in zip(common_energies, avg_cross_sections):
            writer.writerow([f'{energy:.6e}', f'{cross_section:.6e}'])
    
    print(f"Created averaged total sigma file: {output_file}")
    print(f"Energy range: {energy_min:.6e} to {energy_max:.6e} eV")
    print(f"Number of data points: {len(common_energies)}")
    
    return output_file

def update_co_config():
    """Update CO_Mars_HotO.cfg to use the new averaged total sigma file."""
    
    config_file = "src/inputs/CO_Mars_HotO.cfg"
    
    # Read the config file
    with open(config_file, 'r') as f:
        lines = f.readlines()
    
    # Update the total_sigma_file line
    for i, line in enumerate(lines):
        if line.strip().startswith('total_sigma_file') and not line.strip().startswith('#'):
            # Replace with the new averaged file
            lines[i] = 'total_sigma_file    ./inputs/collisions/O-CO_full/elastic_cross_sections/total_elastic_cross_section_averaged.csv\n'
            break
    
    # Write the updated config file
    with open(config_file, 'w') as f:
        f.writelines(lines)
    
    print(f"Updated {config_file} to use averaged total sigma file")

if __name__ == "__main__":
    print("Creating statistically averaged total sigma file for O-CO...")
    create_averaged_total_sigma()
    
    print("\nUpdating CO_Mars_HotO.cfg configuration...")
    update_co_config()
    
    print("\nDone! The configuration has been updated to use the statistically averaged total sigma file.")
