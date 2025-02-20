# for energy = 5 eV
def func1(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*3.89e-15))

# for energy = 0.2 eV
def func2(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*5.52e-15))

# for energy = 1 eV
def func4(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*4.63E-15))

# for energy = 10 eV
def func5(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*3.68E-15))

# for energy = 0.5 eV
def func6(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*5.03E-15))

def fit_escape_prob():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import bisect
    from scipy.optimize import curve_fit
    import matplotlib as mpl
    import os
    from Common_Plotting_Functions import find_bg_dens, find_col_dens_above

    planet = 'Mars'

 
    energies = {0.2: 5.52e-15,
                0.5:5.03E-15,
                1:4.63E-15,
                5:3.89e-15,
                10:3.68E-15}

    z_new = range(8000000,40025000,25000)
    bin_size = 25000

    escape_rates_to_return = []
    escape_z_to_return = []
    calc_esc_dens_to_return = []


    esc_prob_prob_0_2eV_LSA = []
    esc_prob_prob_0_2eV_HSA = []
    esc_prob_prob_5eV_HSA = []
    esc_prob_prob_1eV_LSA = []
    esc_prob_prob_10eV_LSA = []
    esc_prob_prob_0_5eV_LSA = []
    esc_prob_prob_5eV_comb = [] # for combined data
    esc_prob_prob_1eV_comb = [] # for combined data in future
    esc_prob_prob_10eV_comb = [] # for combined data in future
    esc_prob_prob_0_5eV_comb = [] # for combined data in future
    esc_prob_prob_0_2eV_comb = [] # for combined data
    esc_prob_col_dens_5eV = [] # for combined data
    esc_prob_col_dens_1eV = [] # for combined data in future
    esc_prob_col_dens_10eV = [] # for combined data in future
    esc_prob_col_dens_0_5eV = [] # for combined data in future
    esc_prob_col_dens_0_2eV = [] # for combined data

    
    for SA in ['LSA','HSA']:
        esc_prob_z_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 5 eV %s' % (SA)]
        esc_prob_prob_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 5 eV %s' % (SA)]
        esc_prob_prob_1_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 1' % (SA)]
        esc_prob_prob_2_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 2' % (SA)]
        esc_prob_prob_3_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 3' % (SA)]
        esc_prob_prob_4_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 4' % (SA)]
        esc_prob_z_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 0.2 eV %s' % (SA)]
        esc_prob_prob_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 0.2 eV %s' % (SA)]
        esc_prob_prob_1_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 1' % (SA)]
        esc_prob_prob_2_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 2' % (SA)]
        esc_prob_prob_3_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 3' % (SA)]
        esc_prob_prob_4_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 4' % (SA)]
                                           
        bg_z, CO2_density, O_density, CO_density, N2_density = find_bg_dens(SA,planet) # call function to find bg density data from file

        if SA == 'LSA': # as these energies don't have HSA data yet
            esc_prob_z_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 1 eV %s' % (SA)]
            esc_prob_prob_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 1 eV %s' % (SA)]
            esc_prob_prob_1_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['1 eV %s 1' % (SA)]
            esc_prob_prob_2_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['1 eV %s 2' % (SA)]
            esc_prob_prob_3_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['1 eV %s 3' % (SA)]
            esc_prob_prob_4_1eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['1 eV %s 4' % (SA)]

            esc_prob_z_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 10 eV %s' % (SA)]
            esc_prob_prob_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 10 eV %s' % (SA)]
            esc_prob_prob_1_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['10 eV %s 1' % (SA)]
            esc_prob_prob_2_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['10 eV %s 2' % (SA)]
            esc_prob_prob_3_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['10 eV %s 3' % (SA)]
            esc_prob_prob_4_10eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['10 eV %s 4' % (SA)]

            esc_prob_z_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 0.5 eV %s' % (SA)]
            esc_prob_prob_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 0.5 eV %s' % (SA)]
            esc_prob_prob_1_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.5 eV %s 1' % (SA)]
            esc_prob_prob_2_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.5 eV %s 2' % (SA)]
            esc_prob_prob_3_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.5 eV %s 3' % (SA)]
            esc_prob_prob_4_0_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.5 eV %s 4' % (SA)]

        # Interpolate bg densities to match new z grid
        if SA == 'LSA':
            col_dens_above_LSA = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,'LSA',planet)
            col_dens_above = col_dens_above_LSA
        elif SA == 'HSA':
            col_dens_above_HSA = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,'HSA',planet)
            col_dens_above = col_dens_above_HSA


        for energy in energies:
            calculated_escape_z = []
            calculated_escape_dens = []
            for z,alt in enumerate(z_new):
                calculated_escape_z.append(alt)

            if SA == 'LSA':
                if energy == 5:
                    for i,item in enumerate(esc_prob_prob_5eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_5eV[i]*1e5:
                                esc_prob_col_dens_5eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_5eV_comb.append(item)
                    
                elif energy == 1:
                    for i,item in enumerate(esc_prob_prob_1eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_1eV[i]*1e5:
                                esc_prob_col_dens_1eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_1eV_comb.append(item)

                elif energy == 10:
                    for i,item in enumerate(esc_prob_prob_10eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_10eV[i]*1e5:
                                esc_prob_col_dens_10eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_10eV_comb.append(item)


                elif energy == 0.5:
                    for i,item in enumerate(esc_prob_prob_0_5eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_0_5eV[i]*1e5:
                                esc_prob_col_dens_0_5eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_0_5eV_comb.append(item)
                    
                elif energy == 0.2:
                    for i,item in enumerate(esc_prob_prob_0_2eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_0_2eV[i]*1e5:
                                esc_prob_col_dens_0_2eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_0_2eV_comb.append(item)
                    
            elif SA == 'HSA':
                if energy == 5:
                    for i,item in enumerate(esc_prob_prob_5eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_5eV[i]*1e5:
                                esc_prob_col_dens_5eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_5eV_comb.append(item)

                elif energy == 0.2:
                    for i,item in enumerate(esc_prob_prob_0_2eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_0_2eV[i]*1e5:
                                esc_prob_col_dens_0_2eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_0_2eV_comb.append(item) # for combined data

    # Combined data 5 eV
    N_5eV = esc_prob_col_dens_5eV #uses 0s
    N_5eV = np.array(N_5eV, dtype=float)
    popt_M, pcov_M = curve_fit(func1,N_5eV,esc_prob_prob_5eV_comb) # uses 0s
    print('combined 5 eV',popt_M,pcov_M)

    calculated_escape_5eV_LSA = []
    calculated_escape_5eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_5eV_LSA.append(popt_M[0]*np.exp(popt_M[1]*-1*col_dens_above_LSA[z]*energies[5]))
        calculated_escape_5eV_HSA.append(popt_M[0]*np.exp(popt_M[1]*-1*col_dens_above_HSA[z]*energies[5]))

  
    # Combined data 0.2 eV
    N_0_2eV =esc_prob_col_dens_0_2eV #uses 0s
    N_0_2eV = np.array(N_0_2eV, dtype=float)
    popt_M_LE, pcov_M_LE = curve_fit(func2,N_0_2eV,esc_prob_prob_0_2eV_comb) # uses 0s
    print('combined 0.2 eV',popt_M_LE,pcov_M_LE)

    calculated_escape_0_2eV_LSA = []
    calculated_escape_0_2eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_0_2eV_LSA.append(popt_M_LE[0]*np.exp(popt_M_LE[1]*-1*col_dens_above_LSA[z]*energies[0.2]))
        calculated_escape_0_2eV_HSA.append(popt_M_LE[0]*np.exp(popt_M_LE[1]*-1*col_dens_above_HSA[z]*energies[0.2]))
                
    # Data 1 eV
    N_1eV = esc_prob_col_dens_1eV #uses 0s
    N_1eV = np.array(N_1eV, dtype=float)
    popt_M_1eV, pcov_M_1eV = curve_fit(func4,N_1eV,esc_prob_prob_1eV_comb) # uses 0s
    print('combined 1 eV',popt_M_1eV,pcov_M_1eV)

    calculated_escape_1eV_LSA = []
    calculated_escape_1eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_1eV_LSA.append(popt_M_1eV[0]*np.exp(popt_M_1eV[1]*-1*col_dens_above_LSA[z]*energies[1]))
        calculated_escape_1eV_HSA.append(popt_M_1eV[0]*np.exp(popt_M_1eV[1]*-1*col_dens_above_HSA[z]*energies[1]))
   
    # Data 0.5 eV
    N_0_5eV = esc_prob_col_dens_0_5eV #uses 0s
    N_0_5eV = np.array(N_0_5eV, dtype=float)
    popt_M_0_5eV, pcov_M_0_5eV = curve_fit(func6,N_0_5eV,esc_prob_prob_0_5eV_comb) # uses 0s
    print('combined 0.5 eV',popt_M_0_5eV,pcov_M_0_5eV)

    calculated_escape_0_5eV_LSA = []
    calculated_escape_0_5eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_0_5eV_LSA.append(popt_M_0_5eV[0]*np.exp(popt_M_0_5eV[1]*-1*col_dens_above_LSA[z]*energies[0.5]))
        calculated_escape_0_5eV_HSA.append(popt_M_0_5eV[0]*np.exp(popt_M_0_5eV[1]*-1*col_dens_above_HSA[z]*energies[0.5]))

    # Data 10 eV
    N_10eV = esc_prob_col_dens_10eV #uses 0s
    N_10eV = np.array(N_10eV, dtype=float)
    popt_M_10eV, pcov_M_10eV = curve_fit(func5,N_10eV,esc_prob_prob_10eV_comb) # uses 0s
    print('combined 10 eV',popt_M_10eV,pcov_M_10eV)

    calculated_escape_10eV_LSA = []
    calculated_escape_10eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_10eV_LSA.append(popt_M_10eV[0]*np.exp(popt_M_10eV[1]*-1*col_dens_above_LSA[z]*energies[10]))
        calculated_escape_10eV_HSA.append(popt_M_10eV[0]*np.exp(popt_M_10eV[1]*-1*col_dens_above_HSA[z]*energies[10]))
     
    for escape_rates in [calculated_escape_5eV_LSA,calculated_escape_5eV_HSA,calculated_escape_0_2eV_LSA,calculated_escape_0_2eV_HSA,calculated_escape_1eV_LSA,calculated_escape_1eV_HSA,calculated_escape_10eV_LSA,calculated_escape_10eV_HSA,calculated_escape_0_5eV_LSA,calculated_escape_0_5eV_HSA]:
        escape_rates_to_return.append(escape_rates)
    for escape_z in [np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5),np.divide(calculated_escape_z,1e5)]:
        escape_z_to_return.append(escape_z)

    return (escape_rates_to_return,escape_z_to_return)
    
fit_escape_prob()
            
            
