# Functions for plotting scripts for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

###########################
# Read input files to return profiles of background species' densities
def find_bg_dens(SA,planet):
    import pandas as pd
    if planet == 'Mars':
        if SA == 'LSA':
            # LSA data files #
            bg_z = pd.read_csv('../inputs/Mars/bg_densities_LSA_Fox2015.csv')['#alt(cm)']
            CO2_density = pd.read_csv('../inputs/Mars/bg_densities_LSA_Fox2015.csv')['CO2(cm-3)']
            O_density = pd.read_csv('../inputs/Mars/bg_densities_LSA_Fox2015.csv')['O(cm-3)']
            CO_density = pd.read_csv('../inputs/Mars/bg_densities_LSA_Fox2015.csv')['CO(cm-3)']
            N2_density = pd.read_csv('../inputs/Mars/bg_densities_LSA_Fox2015.csv')['N2(cm-3)']
        elif SA == 'HSA':
            # HSA data files #
            bg_z = pd.read_csv('../inputs/Mars/bg_densities_HSA_Fox2015.csv')['#alt(cm)']
            CO2_density = pd.read_csv('../inputs/Mars/bg_densities_HSA_Fox2015.csv')['CO2(cm-3)']
            O_density = pd.read_csv('../inputs/Mars/bg_densities_HSA_Fox2015.csv')['O(cm-3)']
            CO_density = pd.read_csv('../inputs/Mars/bg_densities_HSA_Fox2015.csv')['CO(cm-3)']
            N2_density = pd.read_csv('../inputs/Mars/bg_densities_HSA_Fox2015.csv')['N2(cm-3)']

    elif planet == 'Venus':
        if SA == 'LSA':
            # LSA data files #
            bg_z = pd.read_csv('../inputs/Venus/bg_densities_LSA_FoxSung01_edit.csv')['#alt(cm)']
            CO2_density = pd.read_csv('../inputs/Venus/bg_densities_LSA_FoxSung01_edit.csv')['CO2(cm-3)']
            O_density = pd.read_csv('../inputs/Venus/bg_densities_LSA_FoxSung01_edit.csv')['O(cm-3)']
            CO_density = pd.read_csv('../inputs/Venus/bg_densities_LSA_FoxSung01_edit.csv')['CO(cm-3)']
            N2_density = pd.read_csv('../inputs/Venus/bg_densities_LSA_FoxSung01_edit.csv')['N2(cm-3)']
        elif SA == 'HSA':
            # HSA data files #
            bg_z = pd.read_csv('../inputs/Venus/bg_densities_HSA_FoxSung01_edit.csv')['#alt(cm)']
            CO2_density = pd.read_csv('../inputs/Venus/bg_densities_HSA_FoxSung01_edit.csv')['CO2(cm-3)']
            O_density = pd.read_csv('../inputs/Venus/bg_densities_HSA_FoxSung01_edit.csv')['O(cm-3)']
            CO_density = pd.read_csv('../inputs/Venus/bg_densities_HSA_FoxSung01_edit.csv')['CO(cm-3)']
            N2_density = pd.read_csv('../inputs/Venus/bg_densities_HSA_FoxSung01_edit.csv')['N2(cm-3)']

    return ([bg_z,CO2_density,O_density,CO_density,N2_density])

##########################
# Interpolate background densities and calculate total column density above each z
def find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,SA,planet): 
    import bisect
    import numpy as np
    import pdb
    CO2_dens_new = []
    O_dens_new = []
    CO_dens_new = []
    N2_dens_new = []
    total_dens = []

    for z,alt in enumerate(z_new):
        if alt <= bg_z[0]:
            CO2_dens_new.append(CO2_density[0])
            O_dens_new.append(O_density[0])
            CO_dens_new.append(CO_density[0])
            N2_dens_new.append(N2_density[0])
        elif alt >= bg_z[len(bg_z)-1]: 
            CO2_dens_new.append(CO2_density[len(bg_z)-1])
            O_dens_new.append(O_density[len(bg_z)-1])
            CO_dens_new.append(CO_density[len(bg_z)-1])
            N2_dens_new.append(N2_density[len(bg_z)-1])
        elif alt < bg_z[len(bg_z)-1]:               
            i = bisect.bisect_right(bg_z,alt)
            dndz = (CO2_density[i] - CO2_density[i-1])/(bg_z[i] - bg_z[i-1])
            CO2_dens_new.append(CO2_density[i-1] + dndz*(alt-bg_z[i-1]))
            dndz = (O_density[i] - O_density[i-1])/(bg_z[i] - bg_z[i-1])
            O_dens_new.append(O_density[i-1] + dndz*(alt-bg_z[i-1]))
            dndz = (CO_density[i] - CO_density[i-1])/(bg_z[i] - bg_z[i-1])
            CO_dens_new.append(CO_density[i-1] + dndz*(alt-bg_z[i-1]))
            dndz = (N2_density[i] - N2_density[i-1])/(bg_z[i] - bg_z[i-1])
            N2_dens_new.append(N2_density[i-1] + dndz*(alt-bg_z[i-1]))
       
        total_dens.append(CO2_dens_new[z] + O_dens_new[z] + CO_dens_new[z] + N2_dens_new[z])
        

    # estimate total density above 400 km (=N_400 x scale height_400)
    if planet == 'Mars':
        scale_heights_LSA = [3.47023e6,1.98196e6,1.9822e6,1.26158e6] # scale heights for O, N2, CO, CO2 calculated in model for LSA at 400 km
        scale_heights_HSA = [5.21514e6,2.97854e6,2.97889e6,1.89593e6] # scale heights for O, N2, CO, CO2 calculated in model for HSA at 400 km
    elif planet == 'Venus':
        scale_heights_LSA = [1.7474E+06,9.9802E+05,9.9814E+05,6.3527E+05]
        scale_heights_HSA = [2.1034E+06,1.2013E+06,1.2015E+06,7.6468E+05] # scale heights for O, N2, CO and CO2 calculated by me using 'Venus NH calculation.xlsx'
    densities_at_400 = [O_dens_new[len(z_new)-1],N2_dens_new[len(z_new)-1],CO_dens_new[len(z_new)-1],CO2_dens_new[len(z_new)-1]] # O, N2, CO, CO2 densities at 400 km
    total_dens_above_400 = 0.0
    for i in [0,1,2,3]:
        if SA == 'LSA':
            total_dens_above_400 += (densities_at_400[i]*scale_heights_LSA[i])
        elif SA == 'HSA':
            total_dens_above_400 += (densities_at_400[i]*scale_heights_HSA[i])
        
    # calculate column density above z
    col_dens_above = []
    for z1,alt1 in enumerate(z_new):
        col_dens_z = total_dens_above_400
        for z2,alt2 in enumerate(z_new):
            if z2 > z1:
                col_dens_z += total_dens[z2]*bin_size
        col_dens_above.append(col_dens_z)


    return (col_dens_above)
