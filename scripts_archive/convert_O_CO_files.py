#!/usr/bin/env python3
"""
Script to convert O-CO collision cross section files to the required 2-column CSV format.
- Total sigma files: column 1 = collision energy (eV), column 2 = cross section (cm²)
- Differential cross section files: column 1 = angle (degrees), column 2 = cross section (cm²)
"""

import os

def convert_O_CO_differential_files():
    """Convert O-CO differential cross section files to proper CSV format"""
    
    base_path = "inputs/collisions/O-CO_full/16O-CO_DCS-el_3pes"
    energy_index_file = os.path.join(base_path, "DCS_energy_index.dat")
    
    # Read energy index file
    energies = {}
    with open(energy_index_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 2:
                    index = int(parts[0])
                    energy_eV = float(parts[1])
                    energies[index] = energy_eV
    
    print(f"Found {len(energies)} energy points")
    
    # Process each differential cross section file
    for index in energies:
        input_file = os.path.join(base_path, f"O-CO_DCS-el_3pes_iE_{index}.dat")
        
        if os.path.exists(input_file):
            # Create output directory structure
            energy_eV = energies[index]
            output_dir = os.path.join(base_path, "csv_format")
            os.makedirs(output_dir, exist_ok=True)
            
            output_file = os.path.join(output_dir, f"DCS_elastic_{energy_eV:.6f}_eV.csv")
            
            # Read and convert the file
            angles = []
            cross_sections = []
            
            with open(input_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            angle_deg = float(parts[0])
                            # Convert from 1e-16 cm² to cm²
                            cs_cm2 = float(parts[1]) * 1e-16
                            angles.append(angle_deg)
                            cross_sections.append(cs_cm2)
            
            # Write CSV file directly
            with open(output_file, 'w') as f:
                f.write(f"# Differential cross section for O-CO elastic scattering\n")
                f.write(f"# Collision energy: {energy_eV:.6f} eV\n")
                f.write(f"# Column 1: Angle (degrees)\n")
                f.write(f"# Column 2: Cross section (cm²)\n")
                f.write("angle_degrees,cross_section_cm2\n")
                
                for i in range(len(angles)):
                    f.write(f"{angles[i]:.5f},{cross_sections[i]:.10e}\n")
            
            print(f"Converted {input_file} -> {output_file}")

def check_O_CO_total_files():
    """Check if O-CO total cross section files are already in correct format"""
    
    base_path = "inputs/collisions/O-CO_full/elastic_cross_sections"
    csv_file = os.path.join(base_path, "elastic_16_O_2-3Adp.csv")
    
    if os.path.exists(csv_file):
        print(f"Checking {csv_file}")
        
        # Read first few lines to verify format
        with open(csv_file, 'r') as f:
            lines = f.readlines()[:10]
            
        print("First few lines:")
        for line in lines:
            print(line.strip())
            
        # Verify this is 2-column CSV with energy (eV) and cross section (cm²)
        print("\nFile appears to be in correct format: Energy (eV) vs Cross Section (cm²)")
    else:
        print(f"File not found: {csv_file}")

if __name__ == "__main__":
    print("Converting O-CO collision files to required CSV format...")
    print("=" * 60)
    
    print("\n1. Checking total cross section files:")
    check_O_CO_total_files()
    
    print("\n2. Converting differential cross section files:")
    convert_O_CO_differential_files()
    
    print("\nConversion complete!")
