#!/usr/bin/env python3
"""
Script to verify and correct O-CO2 collision cross section files format.

This script:
1. Compares original O-CO2_full .dat files with converted O-CO2 .csv files
2. Identifies format issues and corrects them
3. Ensures proper 2-column CSV format with correct units

Expected format:
- Total cross sections: energy (eV), cross section (cm²)
- Differential cross sections: angle (degrees), cross section (cm²)

Author: Auto-generated script for corona3d_2020 project
"""

import os
import csv
import math

def convert_radians_to_degrees(rad):
    """Convert radians to degrees."""
    return rad * 180.0 / math.pi

def check_total_cross_section():
    """Check and verify the total cross section file format."""
    
    original_file = "src/inputs/collisions/O-CO2_full/CSs_3pes/O-CO2_elastic-cross-sections_3pes.dat"
    csv_file = "src/inputs/collisions/O-CO2/total_cross_section_O(3P)_CO2_elastic.csv"
    
    print("Checking total cross section file...")
    
    # Read original file
    original_data = []
    with open(original_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                energy_eV = float(parts[0])
                # Original units: 1e+16 cm^2, convert to cm^2
                cs_cm2 = float(parts[1]) * 1e-16
                original_data.append((energy_eV, cs_cm2))
    
    # Read CSV file
    csv_data = []
    with open(csv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split(',')
            if len(parts) >= 2:
                try:
                    energy_eV = float(parts[0])
                    cs_cm2 = float(parts[1])
                    csv_data.append((energy_eV, cs_cm2))
                except ValueError:
                    continue
    
    print(f"Original file: {len(original_data)} data points")
    print(f"CSV file: {len(csv_data)} data points")
    
    # Compare first few points
    print("\nFirst 3 data points comparison:")
    print("Original (eV, cm²) vs CSV (eV, cm²):")
    for i in range(min(3, len(original_data), len(csv_data))):
        orig = original_data[i]
        csv_point = csv_data[i]
        print(f"  {orig[0]:.6f}, {orig[1]:.6e} vs {csv_point[0]:.6f}, {csv_point[1]:.6e}")
        
        # Check if values match (within tolerance)
        energy_match = abs(orig[0] - csv_point[0]) < 1e-6
        cs_match = abs(orig[1] - csv_point[1]) / orig[1] < 0.01  # 1% tolerance
        
        if energy_match and cs_match:
            print(f"    ✓ Point {i+1} matches")
        else:
            print(f"    ✗ Point {i+1} does not match")
    
    return len(original_data) == len(csv_data)

def check_differential_cross_sections():
    """Check and verify differential cross section files."""
    
    print("\nChecking differential cross section files...")
    
    # Read energy table
    energy_table = {}
    with open("src/inputs/collisions/O-CO2_full/DCSs_3pes/energy_table.dat", 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                index = int(parts[0])
                energy = float(parts[1])
                energy_table[index] = energy
    
    print(f"Found {len(energy_table)} energies in energy table")
    
    issues_found = []
    
    # Check a few files
    for file_index in [1, 2, 5, 10]:
        if file_index not in energy_table:
            continue
            
        original_file = f"src/inputs/collisions/O-CO2_full/DCSs_3pes/O-CO2_DCS-3pes_iEng{file_index:02d}.dat"
        csv_file = f"src/inputs/collisions/O-CO2/O-CO2_DCS-3pes_iEng{file_index:02d}.csv"
        
        if not os.path.exists(original_file) or not os.path.exists(csv_file):
            continue
        
        print(f"\nChecking file {file_index} (Energy: {energy_table[file_index]:.6f} eV)...")
        
        # Read original data
        original_data = []
        with open(original_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    angle_rad = float(parts[0])
                    angle_deg = convert_radians_to_degrees(angle_rad)
                    # Original units: 1E-16 cm^2/Sr, convert to cm²/Sr
                    dcs_cm2 = float(parts[1]) * 1e-16
                    original_data.append((angle_deg, dcs_cm2))
        
        # Read CSV data
        csv_data = []
        with open(csv_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split(',')
                if len(parts) >= 2:
                    try:
                        angle_deg = float(parts[0])
                        dcs_cm2 = float(parts[1])
                        csv_data.append((angle_deg, dcs_cm2))
                    except ValueError:
                        continue
        
        print(f"  Original: {len(original_data)} points, CSV: {len(csv_data)} points")
        
        # Compare first few points
        if len(original_data) > 0 and len(csv_data) > 0:
            orig = original_data[0]
            csv_point = csv_data[0]
            
            angle_match = abs(orig[0] - csv_point[0]) < 0.1  # 0.1 degree tolerance
            cs_match = abs(orig[1] - csv_point[1]) / max(orig[1], 1e-20) < 0.01  # 1% tolerance
            
            print(f"  First point: {orig[0]:.2f}°, {orig[1]:.6e} vs {csv_point[0]:.2f}°, {csv_point[1]:.6e}")
            if angle_match and cs_match:
                print(f"    ✓ File {file_index} format appears correct")
            else:
                print(f"    ✗ File {file_index} has issues")
                issues_found.append(file_index)
        else:
            print(f"    ✗ File {file_index} has missing data")
            issues_found.append(file_index)
    
    return issues_found

def fix_differential_files(problematic_files):
    """Fix any problematic differential cross section files."""
    
    if not problematic_files:
        print("\nNo differential files need fixing.")
        return
    
    print(f"\nFixing {len(problematic_files)} differential files...")
    
    # Read energy table
    energy_table = {}
    with open("src/inputs/collisions/O-CO2_full/DCSs_3pes/energy_table.dat", 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                index = int(parts[0])
                energy = float(parts[1])
                energy_table[index] = energy
    
    for file_index in problematic_files:
        if file_index not in energy_table:
            continue
            
        original_file = f"src/inputs/collisions/O-CO2_full/DCSs_3pes/O-CO2_DCS-3pes_iEng{file_index:02d}.dat"
        csv_file = f"src/inputs/collisions/O-CO2/O-CO2_DCS-3pes_iEng{file_index:02d}.csv"
        
        if not os.path.exists(original_file):
            continue
        
        print(f"  Fixing file {file_index}...")
        
        # Read and convert original data
        data_points = []
        with open(original_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    angle_rad = float(parts[0])
                    angle_deg = convert_radians_to_degrees(angle_rad)
                    # Original units: 1E-16 cm^2/Sr, convert to cm²/Sr
                    dcs_cm2 = float(parts[1]) * 1e-16
                    data_points.append((angle_deg, dcs_cm2))
        
        # Write corrected CSV file
        with open(csv_file, 'w') as f:
            f.write("# O(3P)+CO2(000; j=0) -> O(3P)+CO2(000; j') scattering, statistically averaged (3A2+3B2+3B1) potential energy surfaces\n")
            f.write("# Differential cross sections; packaged by M. Gacesa (marko.gacesa@nasa.gov), Apr 30 2019\n")
            f.write(f"# E = {energy_table[file_index]:.6f} eV\n")
            f.write("# theta(deg),elastic DCS (cm²)\n")
            
            for angle_deg, dcs_cm2 in data_points:
                f.write(f"{angle_deg:.6f},{dcs_cm2:.12e}\n")
        
        print(f"    ✓ Fixed file {file_index} with {len(data_points)} data points")

def fix_all_differential_files():
    """Fix all differential cross section files to ensure proper format."""
    
    print("\nConverting all differential cross section files...")
    
    # Read energy table
    energy_table = {}
    with open("src/inputs/collisions/O-CO2_full/DCSs_3pes/energy_table.dat", 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                index = int(parts[0])
                energy = float(parts[1])
                energy_table[index] = energy
    
    converted_count = 0
    
    for file_index in range(1, 42):  # Files 1-41
        if file_index not in energy_table:
            continue
            
        original_file = f"src/inputs/collisions/O-CO2_full/DCSs_3pes/O-CO2_DCS-3pes_iEng{file_index:02d}.dat"
        csv_file = f"src/inputs/collisions/O-CO2/O-CO2_DCS-3pes_iEng{file_index:02d}.csv"
        
        if not os.path.exists(original_file):
            continue
        
        # Read and convert original data
        data_points = []
        with open(original_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    angle_rad = float(parts[0])
                    angle_deg = convert_radians_to_degrees(angle_rad)
                    # Original units: 1E-16 cm^2/Sr, convert to cm²/Sr
                    dcs_cm2 = float(parts[1]) * 1e-16
                    data_points.append((angle_deg, dcs_cm2))
        
        # Write corrected CSV file
        with open(csv_file, 'w') as f:
            f.write("# O(3P)+CO2(000; j=0) -> O(3P)+CO2(000; j') scattering, statistically averaged (3A2+3B2+3B1) potential energy surfaces\n")
            f.write("# Differential cross sections; packaged by M. Gacesa (marko.gacesa@nasa.gov), Apr 30 2019\n")
            f.write(f"# E = {energy_table[file_index]:.6f} eV\n")
            f.write("# Column 1: Angle (degrees)\n")
            f.write("# Column 2: Cross section (cm²)\n")
            f.write("angle_degrees,cross_section_cm2\n")
            
            for angle_deg, dcs_cm2 in data_points:
                f.write(f"{angle_deg:.6f},{dcs_cm2:.12e}\n")
        
        converted_count += 1
        if converted_count % 10 == 0:
            print(f"  Converted {converted_count} files...")
    
    print(f"Successfully converted {converted_count} differential cross section files")

if __name__ == "__main__":
    print("O-CO2 Collision Cross Section File Format Verification")
    print("=" * 60)
    
    # Check total cross sections
    total_ok = check_total_cross_section()
    
    # Check differential cross sections
    issues = check_differential_cross_sections()
    
    if not total_ok:
        print("\n⚠️  Total cross section file needs attention")
    else:
        print("\n✓ Total cross section file format is correct")
    
    if issues:
        print(f"\n⚠️  Found issues with {len(issues)} differential files")
        fix_all_differential_files()
    else:
        print("\n✓ Differential cross section files format appears correct")
        # Still run the fix to ensure all files are properly formatted
        fix_all_differential_files()
    
    print("\nDone! All O-CO2 files should now be in proper 2-column CSV format.")
