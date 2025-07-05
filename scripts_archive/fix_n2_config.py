#!/usr/bin/env python3
"""
Script to update N2_Mars_HotO.cfg to use the correct CSV format differential cross section files.
Maps the file numbers to the actual energy-based CSV filenames.
"""

import os

# Conversion factor from cm⁻¹ to eV
CM_TO_EV = 1.24e-4

def map_filename_to_energy():
    """Create mapping from original filenames to energy-based CSV filenames."""
    
    base_path = "src/inputs/collisions/O-N2_full/n2_o_data/DCS/isotope16/elastic_DCS_isotope_16"
    csv_dir = os.path.join(base_path, "csv_format")
    
    # Map from original file numbers to CSV filenames
    original_to_csv = {}
    
    # Get all original .dat files
    for filename in os.listdir(base_path):
        if filename.startswith("elastic_") and filename.endswith(".dat"):
            # Extract energy in cm⁻¹
            energy_cm = int(filename.replace("elastic_", "").replace(".dat", ""))
            
            # Convert to eV
            energy_eV = energy_cm * CM_TO_EV
            
            # Find corresponding CSV file
            csv_filename = f"DCS_elastic_{energy_eV:.6f}_eV.csv"
            
            original_to_csv[filename] = csv_filename
            
    return original_to_csv

def update_n2_config():
    """Update N2_Mars_HotO.cfg to use correct CSV format files."""
    
    config_file = "src/inputs/N2_Mars_HotO.cfg"
    filename_mapping = map_filename_to_energy()
    
    # Read the config file
    with open(config_file, 'r') as f:
        content = f.read()
    
    # Update each energy file reference
    for original_filename, csv_filename in filename_mapping.items():
        # Replace references to old filenames with new CSV format
        old_pattern = f"csv_format/{original_filename}"
        new_pattern = f"csv_format/{csv_filename}"
        
        content = content.replace(old_pattern, new_pattern)
    
    # Write the updated content back to the file
    with open(config_file, 'w') as f:
        f.write(content)
    
    print(f"Updated {config_file} to use correct energy-based CSV format differential cross section files")
    
    # Print the mappings for verification
    print("\nFilename mappings:")
    for original, csv in sorted(filename_mapping.items()):
        print(f"  {original} -> {csv}")

if __name__ == "__main__":
    update_n2_config()
