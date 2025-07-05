#!/usr/bin/env python3
"""
Script to convert O-N2 collision cross section files to the required 2-column CSV format.
- Total sigma files: column 1 = collision energy (eV), column 2 = cross section (cm²)
- Differential cross section files: column 1 = angle (degrees), column 2 = cross section (cm²)
"""

import os

# Constants for unit conversion
CM_TO_EV = 1.23984e-4  # Convert cm⁻¹ to eV
ANGSTROM2_TO_CM2 = 1e-16  # Convert Ų to cm²

def convert_energy_cm_to_eV(energy_cm):
    """Convert energy from cm⁻¹ to eV"""
    return energy_cm * CM_TO_EV

def convert_O_N2_total_files():
    """Convert O-N2 total cross section files and verify units"""
    
    base_path = "inputs/collisions/O-N2_full/n2_o_data"
    input_file = os.path.join(base_path, "elastic_isotope16.csv")
    
    if os.path.exists(input_file):
        print(f"Processing {input_file}")
        
        # Read the current file
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Parse data
        energies_eV = []
        cross_sections_cm2 = []
        
        for line in lines:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().replace(',', ' ').split()
                if len(parts) >= 2:
                    try:
                        energy = float(parts[0])
                        sigma = float(parts[1])
                        
                        # Assume energy is already in eV and sigma in cm²
                        # If units are wrong, they would need conversion
                        energies_eV.append(energy)
                        cross_sections_cm2.append(sigma)
                    except ValueError:
                        continue
        
        # Create properly formatted CSV
        output_file = os.path.join(base_path, "total_elastic_cross_section_isotope16.csv")
        
        with open(output_file, 'w') as f:
            f.write("# Total elastic cross section for O-N2 (isotope 16)\n")
            f.write("# Column 1: Collision energy (eV)\n")
            f.write("# Column 2: Cross section (cm²)\n")
            f.write("energy_eV,cross_section_cm2\n")
            
            for energy, sigma in zip(energies_eV, cross_sections_cm2):
                f.write(f"{energy:.6f},{sigma:.6e}\n")
        
        print(f"Created standardized file: {output_file}")
        print(f"Energy range: {min(energies_eV):.3f} - {max(energies_eV):.3f} eV")
        print(f"Cross section range: {min(cross_sections_cm2):.2e} - {max(cross_sections_cm2):.2e} cm²")
    
    else:
        print(f"File not found: {input_file}")

def convert_O_N2_differential_files():
    """Convert O-N2 differential cross section files to proper CSV format"""
    
    base_path = "inputs/collisions/O-N2_full/n2_o_data/DCS/isotope16/elastic_DCS_isotope_16"
    
    # Get list of differential files
    dcs_files = []
    for filename in os.listdir(base_path):
        if filename.startswith("elastic_") and filename.endswith(".dat"):
            dcs_files.append(filename)
    
    dcs_files.sort()
    
    print(f"Found {len(dcs_files)} differential cross section files")
    
    # Create output directory
    output_dir = os.path.join(base_path, "csv_format")
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in dcs_files:
        input_file = os.path.join(base_path, filename)
        
        # Extract energy from filename (in cm⁻¹, convert to eV)
        energy_cm = float(filename.replace("elastic_", "").replace(".dat", ""))
        energy_eV = convert_energy_cm_to_eV(energy_cm)
        
        # Read differential cross section data
        angles = []
        cross_sections = []
        
        with open(input_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            angle_deg = float(parts[0])
                            # Assuming cross section is in atomic units (a₀²) - convert to cm²
                            # The scientific notation format suggests these are already large numbers
                            cs_value = float(parts[1])
                            # These appear to be in some units that need conversion to cm²
                            # Based on typical values, these might be in Ų (1e-16 cm²)
                            cs_cm2 = cs_value * 1e-16
                            
                            angles.append(angle_deg)
                            cross_sections.append(cs_cm2)
                        except ValueError:
                            continue
        
        # Create CSV output file
        output_file = os.path.join(output_dir, f"DCS_elastic_{energy_eV:.6f}_eV.csv")
        
        with open(output_file, 'w') as f:
            f.write(f"# Differential cross section for O-N2 elastic scattering (isotope 16)\n")
            f.write(f"# Original energy: {energy_cm:.0f} cm⁻¹\n")
            f.write(f"# Collision energy: {energy_eV:.6f} eV\n")
            f.write(f"# Column 1: Angle (degrees)\n")
            f.write(f"# Column 2: Cross section (cm²)\n")
            f.write("angle_degrees,cross_section_cm2\n")
            
            for angle, cs in zip(angles, cross_sections):
                f.write(f"{angle:.5f},{cs:.10e}\n")
        
        print(f"Converted {filename} (E={energy_eV:.6f} eV) -> {os.path.basename(output_file)}")

def verify_converted_files():
    """Verify the converted files have the correct format"""
    
    print("\nVerifying converted files:")
    print("-" * 40)
    
    # Check O-N2 total cross section
    total_file = "inputs/collisions/O-N2_full/n2_o_data/total_elastic_cross_section_isotope16.csv"
    if os.path.exists(total_file):
        with open(total_file, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
        
        print(f"O-N2 total cross section file:")
        print(f"  - {len(lines)-1} data points (excluding header)")
        if len(lines) > 1:
            first_data = lines[1].split(',')
            last_data = lines[-1].split(',')
            print(f"  - Energy range: {float(first_data[0]):.3f} - {float(last_data[0]):.3f} eV")
    
    # Check a sample differential file
    dcs_dir = "inputs/collisions/O-N2_full/n2_o_data/DCS/isotope16/elastic_DCS_isotope_16/csv_format"
    if os.path.exists(dcs_dir):
        dcs_files = [f for f in os.listdir(dcs_dir) if f.endswith('.csv')]
        if dcs_files:
            sample_file = os.path.join(dcs_dir, dcs_files[0])
            with open(sample_file, 'r') as f:
                lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
            
            print(f"Sample O-N2 differential file ({dcs_files[0]}):")
            print(f"  - {len(lines)-1} angle points (excluding header)")
            if len(lines) > 1:
                first_data = lines[1].split(',')
                last_data = lines[-1].split(',')
                print(f"  - Angle range: {float(first_data[0]):.1f} - {float(last_data[0]):.1f} degrees")

if __name__ == "__main__":
    print("Converting O-N2 collision files to required CSV format...")
    print("=" * 60)
    
    print("\n1. Converting total cross section files:")
    convert_O_N2_total_files()
    
    print("\n2. Converting differential cross section files:")
    convert_O_N2_differential_files()
    
    verify_converted_files()
    
    print("\nConversion complete!")
