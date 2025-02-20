# Script to calculate escaping H production rates and plot figures for Gregory et al. (2023b)
# Bethan Gregory
# 17/06/2023

def calculate_weighted_escape_estimate(SA_list):
        
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import pandas as pd
    import bisect
    import math
    import pdb
    import sys
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                   AutoMinorLocator, LogLocator, LinearLocator)
    import matplotlib.ticker as ticker
    from scipy.signal import butter,lfilter,freqz

    from Common_Plotting_Functions import find_bg_dens,find_col_dens_above
    from fit_escape_prob import fit_escape_prob

    # Make directories for script outputs and figures
    if os.path.exists('../model_output_data/figures') == False:
        os.mkdir('../model_output_data/figures')
    if os.path.exists('../model_output_data/escape_estimates') == False:
        os.mkdir('../model_output_data/escape_estimates')
    if os.path.exists('../model_output_data/escape_estimates/process_figs') == False:
        os.mkdir('../model_output_data/escape_estimates/process_figs')

    planet = 'Mars' 
    
    fig_choice = 14 # do you want to...
    #(1) calculate escape rates from and plot all the mechanisms; Note --- this will break if ecsape_estimates/escape_estimates_LSA and escape_estimates/escape_estimates_HSA already exist.
    #(8) Plot escape estimate method (Gregory et al. (2023) Figure 1);
    #12) Plot reactant densities, production rates, escape probabilities, and escaping H production rates for top 5 mechanisms (Figure 3 in paper);
    #13) Plot accumulated production and escaping production plots, with pie charts (Figure 4 in paper);
    #14) 8 but with an additional panel for thermalisation fraction profile and only charge exchange reactions shown (Figure S1);

    temp_greying_out = 'yes'


    mechanism_labels = ['HCO+ + e- \u2192 CO + H',
                        'OCOH+ DR \u2192 CO2 + H',
                        'OCOH+ DR \u2192 O + CO + H',
                        'OH+ DR \u2192 O(1D) + H',
                        'H3+ DR \u2192 H2 + H',
                        'H3+ DR \u2192 H + H + H',## this one x3
                        'HNO+ DR \u2192 NO + H',
                        'N2H+ DR \u2192 N2 + H',
                        'H2+ DR \u2192 H + H', ## this one x2
                        'CH+ DR \u2192 C(1D) + H',
                        'CH+ DR \u2192 C(1S) + H',
                        'HO2+ DR \u2192 O2 + H',
                        'HO2+ DR \u2192 O + O + H',
                        'ArH+ DR \u2192 Ar + H', 
                        
                        'H + H+ \u2192 H+ + H',
                        'H+ + O \u2192 O+ + H',
                        'CO2 + H+ \u2192 CO2+ + H', # Note that this reaction is endothermic at room temperature so has been neglected

                        'N2+ + H2 \u2192 N2H+ + H',
                        'O+(2D) + H2 \u2192 OH+ + H', # from Fox et al. (1996)
                        'O+(2P) + H2 \u2192 OH+ + H',  # from Fox et al. (1996)
                        'O+(2P) + H2 \u2192 O + H+ + H', # from Fox et al. (1996)
                        'O+(4S) + H2 \u2192 OH+ + H',  # from Fox (2015) -  O(2D)+ - O(2P)+ from Fox et al. (1996)
                        'CO2+ + H2 \u2192 OCOH+ + H',
                        'H2+ + H2 \u2192 H3+ + H',
                        'H2+ + CO2 \u2192 OCOH+ + H',
                        'H2+ + O \u2192 OH+ + H', 
                        'CO+ + H2 \u2192 HCO+ + H', 
                        
                        'CO+ + H2 \u2192 HOC+ + H',
                        'N+ + H2 \u2192 NH+ + H',
                        'OH+ + H2 \u2192 H2O+ + H',
                        'OH+ + N \u2192 NO+ + H', 
                        'OH+ + O \u2192 O2+ + H',
                        'HNO+ + O \u2192 NO2+ + H',
                        'H2+ + N2 \u2192 N2H+ + H',
                        'H2+ + CO \u2192 HCO+ + H',
                        'H2+ + CO \u2192 HOC+ + H', 
                        'H2+ + O2 \u2192 HO2+ + H', 
 
                        'H2+ + N \u2192 NH+ + H', 
                        'H3+ + O \u2192 H2O+ + H',
                        'CH+ + H2 \u2192 CH2+ + H',
                        'CH+ + N \u2192 CN+ + H', 
                        'CH+ + O \u2192 CO+ + H',

                        'H2+ + C \u2192 CH+ + H', 
                        'CH+ + C \u2192 C2+ + H',
                        'Ar+ + H2 \u2192 ArH+ + H', 
                        'H2+ + Ar \u2192 ArH+ + H',

                        'HO2+ + N \u2192 NO2+ + H',
                        'He+ + H2 \u2192 H+ + He + H',
                        'H3O+ + e \u2192 OH + H + H']  # last from Krasnopolsky (2019)

                        
                        
    # Reactant 1 for each mechanism
    mechanism_R1s = ['HCO+','OCOH+','OCOH+','OH+','H3+','H3+','HNO+','N2H+','H2+','CH+','CH+','HO2+','HO2+','ArH+',
                     'H','H+','CO2',
                     'N2+','O+2D','O+2P','O+2P','O+','CO2+','H2+','H2+','H2+','CO+',
                     'CO+','N+','OH+','OH+','OH+','HNO+','H2+','H2+','H2+','H2+',
                     'H2+','H3+','CH+','CH+','CH+','H2+','CH+','Ar+','H2+',
                     'HO2+','He+','H3O+']
    # Reactant 2 for each mechanism
    mechanism_R2s = ['e','e','e','e','e','e','e','e','e','e','e','e','e','e',
                     'H+','O','H+',
                     'H2','H2','H2','H2','H2','H2','H2','CO2','O','H2',
                     'H2','H2','H2','N','O','O','N2','CO','CO','O2',
                     'N','O','H2','N','O','C','C','H2','Ar',
                     'N','H2','e']
    # Specify temperature-dependence (Te = electron temperature, Ti = ion T, Tn = neutral T, T_none = no T dependence)
    mechanism_Ts =  ['Te','Te','Te','Te','Te','Te','Te','Te','Te','Te','Te','Te','Te','Te',
                     'Ti_and_Tn','Ti','T_none',
                     'T_none','T_none','T_none','Ti','T_none','Ti','Ti','T_none','T_none','T_none',
                     'T_none','Ti','T_none','T_none','T_none','T_none','T_none','T_none','T_none','T_none',
                     'T_none','Ti','T_none','T_none','T_none','T_none','T_none','T_none','T_none',
                     'T_none','T_none','Te'] # options are 'Tn' or 'Te' or 'Ti' (or 'T_none' if not required for rate)

    total_esc_prod_order = [1,14,6,26,39,35,36,19,30,32,34,40,41,42,
                        9,7,4,
                            8,3,100,100,100,2,21,13,10,11,
                        12,18,16,22,5,43,15,23,24,29,
                        28,27,31,37,20,33,38,17,25,
                            44]

    # Specify excess energies (eV)
    excess_energies = [7.31,7.95,2.42,6.63,9.23,4.76,7.97,8.47,10.95,5.92,4.50,9.25,4.13,9.68,
                       'see text','see text','see text',
                       2.61,3.86,5.56,0.52,0.539,1.35,1.71,3.00,2.35,2.04,
                       0.61,3.09,1.02,5.84,1.65,1.33,2.43,3.47,2.05,1.67,
                       3.98,1.66,0.18,0.80,4.27,3.78,1.27,1.60,1.27,
                       4.00,6.51,1.31]
                       

    if fig_choice == 12 or fig_choice == 13: # most important reactions, for plot options 12 and 13 - nb. this is v. hardcoded
        LSA_5eV_top = ['HCO+ + e- \u2192 CO + H',
                       'CO2+ + H2 \u2192 OCOH+ + H',
                       'O+(4S) + H2 \u2192 OH+ + H', 
                       'OH+ + O \u2192 O2+ + H',
                       'OCOH+ DR \u2192 O + CO + H'] 
        LSA_0_2eV_top = LSA_5eV_top
        HSA_5eV_top = ['HCO+ + e- \u2192 CO + H',
                       'CO2+ + H2 \u2192 OCOH+ + H',
                       'O+(4S) + H2 \u2192 OH+ + H',
                       'OH+ + O \u2192 O2+ + H',
                       'N2+ + H2 \u2192 N2H+ + H']
        HSA_0_2eV_top = HSA_5eV_top

        LSA_colors = {'HCO+':'k',
                      'H':'k',
                  #    'H+':'k',
                      'e':'k',
                      'O':'gray',
                      'CO2+':'silver',
                      'H2':'k', # silver when H+ + H included
                      'OH+':'gray',
                      'O+':'gray',
                      'OCOH+':'k',
                      'N2+':'silver',}

        HSA_colors = {'HCO+':'k',
                      'H':'k',
                   #   'H+':'k',
                      'e':'k',
                      'O':'gray',
                      'OCOH+':'k',
                      'N2+':'silver',
                      'H2':'k', # silver when H+ + H included
                      'OH+':'gray',
                      'O+':'gray',
                      'CO2+':'silver'}
        colors = ['#c51160','#f7c1bb','#84b082','#353a47','#885a5a','#4e8069'] # brighter pink = #dc136c
        
    elif fig_choice == 8: # only plot 2 mechanisms
        LSA_5eV_top = ['HCO+ + e- \u2192 CO + H',
                       'O+(4S) + H2 \u2192 OH+ + H']
        LSA_0_2eV_top = LSA_5eV_top
        HSA_5eV_top = ['HCO+ + e- \u2192 CO + H',
                       'O+(4S) + H2 \u2192 OH+ + H']
        HSA_0_2eV_top = HSA_5eV_top
        colors = ['#c51160','#84b082'] # dark pink and light green
        gradient_colors_HCOpl = ['#E288AF','#D0407F','#c51160','#B10F56','#890B43'] #'#E79FBF',
        gradient_colors_Opl = ['#b3ceb2','#9cbf9a','#84b082','#5a8b58','#4b7349']
        
    elif fig_choice == 14: # only plot charge exchange mechanisms
        LSA_5eV_top = ['H + H+ \u2192 H+ + H',
                        'H+ + O \u2192 O+ + H']
   #                     'CO2 + H+ \u2192 CO2+ + H']
        LSA_0_2eV_top = LSA_5eV_top
        HSA_5eV_top = ['H + H+ \u2192 H+ + H',
                        'H+ + O \u2192 O+ + H']
                   #     'CO2 + H+ \u2192 CO2+ + H']
        HSA_0_2eV_top = HSA_5eV_top
        colors = ['rebeccapurple','mediumseagreen','darkorange'] # dark pink and light green
        gradient_colors_HCOpl = ['#E288AF','#D0407F','#c51160','#B10F56','#890B43'] #'#E79FBF',
        gradient_colors_Opl = ['#b3ceb2','#9cbf9a','#84b082','#5a8b58','#4b7349']
        gradient_colors = ['rebeccapurple','mediumseagreen','darkorange'] ##E288AF','#D0407F','#c51160','#B10F56','#890B43']
        
        
    else:
        colors = ['#c51160','#f7c1bb','#84b082','#353a47','#885a5a','#4e8069'] # brighter pink = #dc136c

    gradient_colors_grey = ['0.75','0.55','0.4','0.25','0']

    accumulated_escape_LSA = []
    accumulated_escape_HSA = []
    accumulated_prod_LSA = []
    accumulated_prod_HSA = []
    accumulated_colors_LSA = []
    accumulated_colors_HSA = []


    if fig_choice == 12:
        fig2 = plt.figure(figsize=(18,8))
        axb4 = fig2.add_subplot(2,4,1) # LSA - reactant dens. vs alt.
        axb5 = fig2.add_subplot(2,4,2) # LSA - production rate vs. alt.
        axb6 = fig2.add_subplot(2,4,3) # LSA - esc. prob. vs. alt. (0.2 eV, 0.5 eV, 1 eV, 5 eV, 10 eV)
        axb7 = fig2.add_subplot(2,4,4) # LSA - escaping H production rate
        axb8 = fig2.add_subplot(2,4,5) # HSA - reactant dens. vs. alt.
        axb9 = fig2.add_subplot(2,4,6) # HSA - prod. rate vs. alt.
        axb10 = fig2.add_subplot(2,4,7) # HSA - esc. prob. vs. alt (0.2 eV, 0.5 eV, 1 eV, 5 eV, 10 eV)
        axb11 = fig2.add_subplot(2,4,8) # HSA - escaping H production rate
    elif fig_choice == 8:
        fig8 = plt.figure(figsize=(14,5))
        axb5 = fig8.add_subplot(1,3,1) # LSA - production rate vs. alt.
        axb6 = fig8.add_subplot(1,3,2) # LSA - esc. prob. vs. alt. (0.2 eV, 0.5 eV, 1 eV, 5 eV, 10 eV)
        axb7 = fig8.add_subplot(1,3,3) # LSA - escaping H production rate
    elif fig_choice == 13:
        fig13 = plt.figure(figsize=(12,9))
        axb7 = fig13.add_subplot(2,2,1) # LSA - cumulative production rate vs. alt.
        axb11 = fig13.add_subplot(2,2,3) # LSA - cumulative production rate vs. alt.
        axb12 = fig13.add_subplot(2,2,2) # HSA - cumulative escaping H production rate vs. alt.
        axb13 = fig13.add_subplot(2,2,4) # HSA - cumulative escaping H production rate vs. alt.
    elif fig_choice == 14:
        fig14 = plt.figure(figsize=(16,5))
        axb5 = fig14.add_subplot(1,4,1) # LSA - reactant dens. vs alt.
        axb6 = fig14.add_subplot(1,4,2) # LSA - production rate vs. alt.
        axb7 = fig14.add_subplot(1,4,4) # LSA - escaping H production rate
        axb12 = fig14.add_subplot(1,4,3) # LSA - fraction of particles produced above escape energy vs. alt.

    color = 'tab:blue'
    color2 = 'rebeccapurple'
    markers = ['x','o','^','d']
    linestyles = ['-','--']

    reactants_LSA = []
    reactants_HSA = []
    mech_count_LSA = 0
    mech_count_HSA = 0

    # Iterate over low and high solar activity
    for s,SA in enumerate(SA_list):
    
        ### Escape probability data ### No longer utilised
        esc_prob_file = '../model_output_data/escape_probabilities_%s.csv' % (SA)
        esc_prob_z_5eV = pd.read_csv(esc_prob_file)['Altitude (km) 5 eV %s' % (SA)]
        esc_prob_prob_5eV = pd.read_csv(esc_prob_file)['Mean probability 5 eV %s' % (SA)]
        esc_prob_1_5eV = pd.read_csv(esc_prob_file)['5 eV %s 1' % (SA)]
        esc_prob_2_5eV = pd.read_csv(esc_prob_file)['5 eV %s 2' % (SA)]
        esc_prob_3_5eV = pd.read_csv(esc_prob_file)['5 eV %s 3' % (SA)]
        esc_prob_4_5eV = pd.read_csv(esc_prob_file)['5 eV %s 4' % (SA)]

        ### Create new z grid ###                                                                                       
        z_new = range(8000000,40025000,25000) # new z grid 80-400 km in cm, with 0.25 km spacing
        bin_size = 25000 # 0.25 km

        ### Import T profile data ###
        T_file = '../inputs/Mars/MarsTemp%s_Fox2015.csv' % (SA)
        z_orig_T = pd.read_csv(T_file)['#alt(cm)']
        Tn_orig_T = pd.read_csv(T_file)['Tn(K)']
        Te_orig_T = pd.read_csv(T_file)['Te(K)']
        Ti_orig_T = pd.read_csv(T_file)['Ti(K)']

        remnant_altitude = pd.read_csv('Charge_exchange_reactions_energy_equation.csv')['Altitude %s' % (SA)]
        CO2_remnant_fraction = pd.read_csv('Charge_exchange_reactions_energy_equation.csv')['CO2 %s: Fraction above 4991 m/s' % (SA)]
        O_remnant_fraction = pd.read_csv('Charge_exchange_reactions_energy_equation.csv')['O %s using E: Fraction above 0.13 eV' % (SA)]
        H_remnant_fraction = pd.read_csv('Charge_exchange_reactions_energy_equation.csv')['H %s: Fraction above 4991 m/s' % (SA)]
    
        # Check that the input profiles increase in altitude (as this script will only work that way)
        if z_orig_T[0] > z_orig_T[len(z_orig_T)-1]:
            print('Please use T profiles with increasing altitude')
            stop

        ### Interpolate T for new grid - Tn, Ti and Te ###
        Tn_new = []
        Ti_new = []
        Te_new = []
        # - if the new altitude is below the lowest altitude in the dataset, set T at the new altitude to the lowest T in the dataset                                                                                                                  
        # - if the new altitude is above the highest altitude in the dataset, set T at the new altitude to the highest T in the dataset                                                                                                                
        # - otherwise, do a linear interpolation using the gradient between the two existing points either side of the new point and the difference in z between the new z and the lowest adjacent existing z
        for z,alt in enumerate(z_new):
            if alt <= z_orig_T[0]: # I don't think this significantly happens
                Tn_new.append(Tn_orig_T[0]) 
                Te_new.append(Te_orig_T[0]) 
                Ti_new.append(Ti_orig_T[0]) 
            elif alt >= z_orig_T[len(z_orig_T)-1]: # I don't think this significantly happens
                Tn_new.append(Tn_orig_T[len(z_orig_T)-1])
                Te_new.append(Te_orig_T[len(z_orig_T)-1])
                Ti_new.append(Ti_orig_T[len(z_orig_T)-1])
            else:
                i = bisect.bisect_right(z_orig_T,alt)
                dndz = (Tn_orig_T[i] - Tn_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Tn_new.append(Tn_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))
                dndz = (Te_orig_T[i] - Te_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Te_new.append(Te_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))
                dndz = (Ti_orig_T[i] - Ti_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Ti_new.append(Ti_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))

        # Interpolate non-thermalised remnant fractions for charge exchange reactions
        CO2_remnant_fraction_new = []
        O_remnant_fraction_new = []
        H_remnant_fraction_new = []
        for z,alt in enumerate(z_new):
            if alt <= remnant_altitude[0]: # I don't think this significantly happens
                CO2_remnant_fraction_new.append(CO2_remnant_fraction[0])
                O_remnant_fraction_new.append(O_remnant_fraction[0])
                H_remnant_fraction_new.append(H_remnant_fraction[0])
            elif alt >= remnant_altitude[len(remnant_altitude)-1]:
                CO2_remnant_fraction_new.append(CO2_remnant_fraction[len(remnant_altitude)-1])
                O_remnant_fraction_new.append(O_remnant_fraction[len(remnant_altitude)-1])
                H_remnant_fraction_new.append(H_remnant_fraction[len(remnant_altitude)-1])
            else:
                i = bisect.bisect_right(remnant_altitude,alt)
                dfdz = (CO2_remnant_fraction[i] - CO2_remnant_fraction[i-1])/(remnant_altitude[i] - remnant_altitude[i-1])
                CO2_remnant_fraction_new.append(CO2_remnant_fraction[i-1] + dfdz*(alt-remnant_altitude[i-1]))
                dfdz = (O_remnant_fraction[i] - O_remnant_fraction[i-1])/(remnant_altitude[i] - remnant_altitude[i-1])
                O_remnant_fraction_new.append(O_remnant_fraction[i-1] + dfdz*(alt-remnant_altitude[i-1]))
                dfdz = (H_remnant_fraction[i] - H_remnant_fraction[i-1])/(remnant_altitude[i] - remnant_altitude[i-1])
                H_remnant_fraction_new.append(H_remnant_fraction[i-1] + dfdz*(alt-remnant_altitude[i-1]))
        if fig_choice ==14:
            linewidths12=[2,0.75]
            axb12.plot(O_remnant_fraction_new,np.divide(z_new,1e5),color=gradient_colors[1],linewidth=linewidths12[s])
            axb12.plot(H_remnant_fraction_new,np.divide(z_new,1e5),color=gradient_colors[0],linewidth=linewidths12[s])
 
        returned_escape_rates = fit_escape_prob()[0]
        returned_escape_z = fit_escape_prob()[1]
        returned_escape_rates_LSA_5eV = returned_escape_rates[0]
        returned_escape_rates_LSA_0_2eV = returned_escape_rates[2]
        returned_escape_rates_HSA_5eV = returned_escape_rates[1]
        returned_escape_rates_HSA_0_2eV = returned_escape_rates[3]
        returned_escape_rates_LSA_1eV = returned_escape_rates[4]
        returned_escape_rates_LSA_10eV = returned_escape_rates[6]
        returned_escape_rates_HSA_1eV = returned_escape_rates[5]
        returned_escape_rates_HSA_10eV = returned_escape_rates[7]
        returned_escape_rates_LSA_0_5eV = returned_escape_rates[8]
        returned_escape_rates_HSA_0_5eV = returned_escape_rates[9]
        returned_escape_z_LSA_5eV = returned_escape_z[0]
        returned_escape_z_LSA_0_2eV = returned_escape_z[2]
        returned_escape_z_HSA_5eV = returned_escape_z[1]
        returned_escape_z_HSA_0_2eV = returned_escape_z[3]
        returned_escape_z_LSA_1eV = returned_escape_z[4]
        returned_escape_z_LSA_10eV = returned_escape_z[6]
        returned_escape_z_LSA_0_5eV = returned_escape_z[8]
        returned_escape_z_HSA_1eV = returned_escape_z[5]
        returned_escape_z_HSA_10eV = returned_escape_z[7]
        returned_escape_z_HSA_0_5eV = returned_escape_z[9]    

        ### Interpolate escape probability for new grid, if 'smooth_interp' chosen
        esc_prob_new = []

        if fig_choice == 1:
            # Create output csv file for weighted escape estimates - NB. will throw an error if the file already exists
            with open('../model_output_data/escape_estimates/escape_estimates_%s.csv' % (SA),'x') as f1: # 'x' option means that it will fail if the file already exists
                if SA == 'LSA':
                    f1.write('R1,R2,Mechanism,Rough escape estimate 5 eV,Rough escape estimate 0.2 eV,Esc flux 2dp 5eV, Esc flux 2dp 0.2eV,Peak prod alt,Total prod,90 percent 5 eV,90 percent 0.2 eV,Esc flux 1 eV,Esc flux 10 eV,Esc flux 0.5 eV,Excess energy (eV)\n')
                elif SA == 'HSA':
                    f1.write('R1,R2,Mechanism,Rough escape estimate 5 eV,Rough escape estimate 0.2 eV,esc flux 2dp 5eV, Esc flux 2dp 0.2eV,Peak prod alt,Total prod,90 percent 5 eV,90 percent 0.2 eV,Esc flux 1 eV,Esc flux 10 eV,Esc flux 0.5 eV,Excess energy (eV)\n')
    
        ### Important reactant profiles. ###

        for m, mech in enumerate(mechanism_labels):
      #  for m,mech in enumerate(['HCO+ + e- \u2192 CO + H','OCOH+ DR \u2192 CO2 + H']): # only show 2 plots if debugging

            if fig_choice == 1:
                ### Set up plot ###
                fig1 = plt.figure(figsize=(30,5))
                ax1 = fig1.add_subplot(1,6,1) # R1 density with z
                ax2 = fig1.add_subplot(1,6,2) # R2 density with z
                ax3 = fig1.add_subplot(1,6,3) # T (chosen) with z
                ax4 = fig1.add_subplot(1,6,4) # production rate with z
                ax5 = fig1.add_subplot(1,6,5) # escape probability with z
                ax6 = fig1.add_subplot(1,6,6) # weighted escape with z

            ### PLot escape probability data ###
              # Interpolate escape probability data for 'smooth' style
                if SA == 'LSA':
                    ax5.plot(returned_escape_rates_LSA_5eV,returned_escape_z_LSA_5eV,label='fitted curve 5 eV')
                    ax5.plot(returned_escape_rates_LSA_0_2eV,returned_escape_z_LSA_0_2eV,label='fitted curve 0.2 eV')
                elif SA == 'HSA':
                    ax5.plot(returned_escape_rates_HSA_5eV,returned_escape_z_HSA_5eV,label='fitted curve 5 eV')
            
                ax5.set_xlabel('Escape probability')
                ax5.set_xlim([-0.10,1.1])
    
            elif fig_choice == 8  or fig_choice == 12 or fig_choice == 14:
                if SA == 'LSA' and m == 0:
                    if fig_choice == 8 or fig_choice == 12: # Plot 0.5 eV, 1 eV and 10 eV too
                        axb6.plot(returned_escape_rates_LSA_10eV,returned_escape_z_LSA_10eV,label='fitted curve 10 eV',linewidth=2.75,color=gradient_colors_grey[4])
                        axb6.plot(returned_escape_rates_LSA_5eV,returned_escape_z_LSA_5eV,label='fitted curve 5 eV',linewidth=2.25,color=gradient_colors_grey[3])
                        axb6.plot(returned_escape_rates_LSA_1eV,returned_escape_z_LSA_1eV,label='fitted curve 1 eV',linewidth=1.75,color=gradient_colors_grey[2])
                        axb6.plot(returned_escape_rates_LSA_0_5eV,returned_escape_z_LSA_0_5eV,label='fitted curve 0.5 eV',linewidth=1.25,color=gradient_colors_grey[1])
                        axb6.plot(returned_escape_rates_LSA_0_2eV,returned_escape_z_LSA_0_2eV,label='fitted curve 0.2 eV',linewidth=0.75,color=gradient_colors_grey[0])
                    elif fig_choice == 14:
                        axb6.plot(returned_escape_rates_LSA_0_2eV,returned_escape_z_LSA_0_2eV,label='fitted curve 0.2 eV',linewidth=2,color='k')
                        
                          
                    else:
                        axb6.plot(returned_escape_rates_LSA_5eV,returned_escape_z_LSA_5eV,label='fitted curve 5 eV',linewidth=2,color='k')
                        axb6.plot(returned_escape_rates_LSA_0_2eV,returned_escape_z_LSA_0_2eV,label='fitted curve 0.2 eV',linewidth=0.5,color='k')
                else:
                    if fig_choice == 14:
                         axb6.plot(returned_escape_rates_HSA_0_2eV,returned_escape_z_HSA_0_2eV,label='fitted curve 0.2 eV',linewidth=0.75,color='k')
                    elif  fig_choice == 12:# Plot 0.5 eV, 1 eV and 10 eV too
                        axb10.plot(returned_escape_rates_HSA_10eV,returned_escape_z_HSA_10eV,label='fitted curve 10 eV',linewidth=2.75,color=gradient_colors_grey[4])
                        axb10.plot(returned_escape_rates_HSA_5eV,returned_escape_z_HSA_5eV,label='fitted curve 5 eV',linewidth=2.25,color=gradient_colors_grey[3])
                        axb10.plot(returned_escape_rates_HSA_1eV,returned_escape_z_HSA_1eV,label='fitted curve 1 eV',linewidth=1.75,color=gradient_colors_grey[2])
                        axb10.plot(returned_escape_rates_HSA_0_5eV,returned_escape_z_HSA_0_5eV,label='fitted curve 0.5 eV',linewidth=1.25,color=gradient_colors_grey[1])
                        axb10.plot(returned_escape_rates_HSA_0_2eV,returned_escape_z_HSA_0_2eV,label='fitted curve 0.2 eV',linewidth=0.75,color=gradient_colors_grey[0])
                        
                        
                        
                    elif fig_choice != 8 and fig_choice != 14:
                        axb10.plot(returned_escape_rates_HSA_5eV,returned_escape_z_HSA_5eV,label='fitted curve 5 eV',linewidth=1.5,color='k')
                        axb10.plot(returned_escape_rates_HSA_0_2eV,returned_escape_z_HSA_0_2eV,label='fitted curve 0.2 eV',linewidth=0.5,color='k')
                

            T_type = mechanism_Ts[m]
            if fig_choice == 1:
                ### Plot T data (choose Tn, Te, Ti depending on mechanism rate dependence)
                if T_type == 'Tn':
                    ax3.plot(Tn_orig_T,np.divide(z_orig_T,1e5),color=color,linewidth=1)
                    ax3.set_xlabel('T$_{n}$ (K)')
                elif T_type == 'Te':
                    ax3.plot(Te_orig_T,np.divide(z_orig_T,1e5),color=color,linewidth=1)
                    ax3.set_xlabel('T$_{e}$ (K)')
                elif T_type == 'Ti':
                    ax3.plot(Ti_orig_T,np.divide(z_orig_T,1e5),color=color,linewidth=1)
                    ax3.set_xlabel('T$_{i}$ (K)')
                elif T_type == 'Ti_and_Tn':
                    ax3.plot(Tn_orig_T,np.divide(z_orig_T,1e5),color=color,linewidth=1)
                    ax3.plot(Ti_orig_T,np.divide(z_orig_T,1e5),color='tomato',linewidth=1)
                    ax3.set_xlabel('T$_{n}$ and T$_{i}$ (K)')
                    ax3.legend()


            if mechanism_R1s[m] == 'CO2' or mechanism_R1s[m] == 'CO' or mechanism_R1s[m] == 'O' or mechanism_R1s[m] == 'N2':
                r1_data_file = '../inputs/Mars/bg_densities_%s_Fox2015.csv' % (SA)
            elif mechanism_R1s[m] == 'H' or mechanism_R1s[m] == 'H2' or mechanism_R1s[m] == 'N' or mechanism_R1s[m] == 'O2' or mechanism_R1s[m] == 'C' or mechanism_R1s[m] == 'Ar':
                r1_data_file = '../inputs/Mars/%s_density_profile_%s_Fox2015.csv' % (mechanism_R1s[m],SA)
            elif mechanism_R1s[m] == 'HCO+':  # CHANGE THIS BACK!!
                r1_data_file = '../inputs/Mars/%s_density_profile_%s_eroded_Fox2015.csv' % (mechanism_R1s[m],SA)
            elif mechanism_R1s[m] == 'O+2D' or mechanism_R1s[m] == 'O+2P':
                r1_data_file = '../inputs/Mars/%s_density_profile_%s_FoxEtAl1996.csv' % (mechanism_R1s[m],SA) # Uses O+(2D) and O+(2P) profiles from Fox et al. (1996)
            elif mechanism_R1s[m] == 'O+':
                r1_data_file = '../inputs/Mars/O+_density_profile_%s_eroded_Fox2015.csv' % (SA)
            elif mechanism_R1s[m] == 'H3O+': # from Krasnopolsky (2019)
                r1_data_file = '../inputs/Mars/H3O+_density_profile_MSA_Krasnopolsky2019.csv'
            else:
                r1_data_file = '../inputs/Mars/%s_density_profile_%s_eroded_Fox2015.csv' % (mechanism_R1s[m],SA)
            
            if mechanism_R2s[m] == 'e':
                r2_data_file = '../inputs/Mars/electron_density_profile_%s_eroded_Fox2015.csv' % (SA)
            elif mechanism_R2s[m] == 'H' or mechanism_R2s[m] == 'H2' or mechanism_R2s[m] == 'N' or mechanism_R2s[m] == 'O2' or mechanism_R2s[m] == 'C' or mechanism_R2s[m] == 'Ar':
                r2_data_file = '../inputs/Mars/%s_density_profile_%s_Fox2015.csv' % (mechanism_R2s[m],SA)
            elif mechanism_R2s[m] == 'CO2' or mechanism_R2s[m] == 'CO' or mechanism_R2s[m] == 'O' or mechanism_R2s[m] == 'N2':
                r2_data_file = '../inputs/Mars/bg_densities_%s_Fox2015.csv' % (SA)
            else:
                r2_data_file = '../inputs/Mars/%s_density_profile_%s_eroded_Fox2015.csv' % (mechanism_R2s[m],SA)
            
            n_orig_r1 = pd.read_csv(r1_data_file)['%s(cm-3)' % (mechanism_R1s[m])]
            n_orig_r2 = pd.read_csv(r2_data_file)['%s(cm-3)' % (mechanism_R2s[m])]

            z_orig_r1 = pd.read_csv(r1_data_file)['#alt(cm)']
            z_orig_r2 = pd.read_csv(r2_data_file)['#alt(cm)']
                                                                                             
            # Check that the input profiles increase in altitude (as this script will only work that way)                
            for prof in [z_orig_r1,z_orig_r2]:
                if prof[0] > prof[len(prof)-1]:
                    print('Please use profiles with increasing altitude')
                    stop

    
            ### Interpolate density inputs and calculate production rate to have value for each altitude where there is an escape probability profile ###                                                                                         

            n_new_r1 = []
            n_new_r2 = []
            production_rate = []
        #    esc_estimate = []
            esc_estimate_5eV = []
            esc_estimate_0_2eV = []
            esc_estimate_1eV = []
            esc_estimate_10eV = []
            esc_estimate_0_5eV = []
            weighted_escape = 0.0
            weighted_escape_5eV = 0.0
            weighted_escape_0_2eV = 0.0
            weighted_escape_1eV = 0.0
            weighted_escape_10eV = 0.0
            weighted_escape_0_5eV = 0.0
            total_prod_LSA = 0.0
            total_prod_HSA = 0.0

            # Interpolate densities for new grid                                                                                       
            # - if the new altitude is below the lowest altitude in the dataset, set the density at the new altitude to the lowest density in the dataset                                                                                                                  
            # - if the new altitude is above the highest altitude in the dataset, set the density at the new altitude to the highest density in the dataset                                                                                                                
            # - otherwise, do a linear interpolation using the gradient between the two existing points either side of the new point and the difference in z between the new z and the lowest adjacent existing z
            # Interpolate R1 densities
            for z,alt in enumerate(z_new):
                if alt <= z_orig_r1[0]:
                    extrap_number_dens = extrapolation(mechanism_R1s[m],alt,n_orig_r1,z_orig_r1,'lower')
                    n_new_r1.append(extrap_number_dens)
                elif alt >= z_orig_r1[len(z_orig_r1)-1]:
                    extrap_number_dens = extrapolation(mechanism_R1s[m],alt,n_orig_r1,z_orig_r1,'upper')
                    n_new_r1.append(extrap_number_dens)
                else:
                    # nb. this is linear interpolation, not log interpolation. Consider which is best, especially for HCO. Might need to consider what is best for model too.
                    i = bisect.bisect_right(z_orig_r1,alt)
                    dndz = (n_orig_r1[i] - n_orig_r1[i-1])/(z_orig_r1[i] - z_orig_r1[i-1])
                    n_new_r1.append(n_orig_r1[i-1] + dndz*(alt-z_orig_r1[i-1]))
            if mech == 'O+(2P) + H2 \u2192 OH+ + H':
                n_new_r1_Opl2P = n_new_r1
            elif mech == 'O+(2D) + H2 \u2192 OH+ + H':
                n_new_r1_Opl2D = n_new_r1
            elif mech == 'O+(4S) + H2 \u2192 OH+ + H':
                n_new_r1 = np.subtract(n_new_r1,n_new_r1_Opl2P)
                n_new_r1 = np.subtract(n_new_r1,n_new_r1_Opl2D) # O+4S reaction must be after others for this to work
                    
            if fig_choice == 1:
                ax1.plot(n_new_r1,np.divide(z_new,1e5),color=color,linewidth=1)
                ax1.plot(n_orig_r1,np.divide(z_orig_r1,1e5),color=color2,linewidth=1)
                if mechanism_R1s[m] == 'ArH+' or mechanism_R1s[m] == 'CH+' or mechanism_R1s[m] == 'CO+' or mechanism_R1s[m] == 'H2+' or mechanism_R1s[m] == 'H3+' or mechanism_R1s[m] == 'He+' or mechanism_R1s[m] == 'HNO+' or mechanism_R1s[m] == 'HO2+' or mechanism_R1s[m] == 'N+' or mechanism_R1s[m] == 'N2+' or mechanism_R1s[m] == 'N2H+' or mechanism_R1s[m] == 'OCOH+' or mechanism_R1s[m] == 'OH+' or mechanism_R1s[m] == 'Ar+' or mechanism_R1s[m] == 'C+' or mechanism_R1s[m] == 'C' or mechanism_R1s[m] == 'N' or mechanism_R1s[m] == 'O+' or mechanism_R1s[m] == 'O+2D' or mechanism_R1s[m] == 'O+2P':
                    ax1.scatter(n_orig_r1[2],np.divide(z_orig_r1[2],1e5),color='red',marker='x')
                    ax1.scatter(n_orig_r1[0],np.divide(z_orig_r1[0],1e5),color='red',marker='x')
                if mechanism_R1s[m] == 'C' or mechanism_R1s[m] == 'N' or mechanism_R1s[m] == 'Ar' or mechanism_R1s[m] == 'O2' or mechanism_R2s[m] == 'O+2D' or mechanism_R2s[m] == 'O+2P':
                    ax1.scatter(n_orig_r1[len(n_orig_r1)-1],np.divide(z_orig_r1[len(z_orig)-1],1e5),color='red',marker='x')
                    ax1.scatter(n_orig_r1[len(n_orig_r1)-3],np.divide(z_orig_r1[len(n_orig)-3],1e5),color='red',marker='x')
                ax1.set_xlabel('%s number density (cm$^{-3}$)' % (mechanism_R1s[m]))
                ax1.set_ylabel('Altitude (km)')
                ax1.set_title('%s %s' % (mechanism_labels[m],SA))



            # Interpolate R2 densities for new grid                                                        
            for z,alt in enumerate(z_new):
                if alt <= z_orig_r2[0]:
                    extrap_number_dens = extrapolation(mechanism_R2s[m],alt,n_orig_r2,z_orig_r2,'lower')
                    n_new_r2.append(extrap_number_dens)
                elif alt >= z_orig_r2[len(z_orig_r2)-1]:
                    extrap_number_dens = extrapolation(mechanism_R2s[m],alt,n_orig_r2,z_orig_r2,'upper')
                    n_new_r2.append(extrap_number_dens)
                else:
                    i = bisect.bisect_right(z_orig_r2,alt)
                    dndz = (n_orig_r2[i] - n_orig_r2[i-1])/(z_orig_r2[i] - z_orig_r2[i-1])
                    n_new_r2.append(n_orig_r2[i-1] + dndz*(alt-z_orig_r2[i-1]))
                
                    
            if fig_choice == 1:
                ax2.plot(n_new_r2,np.divide(z_new,1e5),color=color,linewidth=1) 
                ax2.plot(n_orig_r2,np.divide(z_orig_r2,1e5),color=color2,linewidth=1)
                if mechanism_R2s[m] == 'ArH+' or mechanism_R2s[m] == 'CH+' or mechanism_R2s[m] == 'CO+' or mechanism_R2s[m] == 'H2+' or mechanism_R2s[m] == 'H3+' or mechanism_R2s[m] == 'He+' or mechanism_R2s[m] == 'HNO+' or mechanism_R2s[m] == 'HO2+' or mechanism_R2s[m] == 'N+' or mechanism_R2s[m] == 'N2+' or mechanism_R2s[m] == 'N2H+' or mechanism_R2s[m] == 'OCOH+' or mechanism_R2s[m] == 'OH+' or mechanism_R2s[m] == 'Ar+' or mechanism_R2s[m] == 'C+' or mechanism_R2s[m] == 'C' or mechanism_R2s[m] == 'N' or mechanism_R2s[m] == 'O+':
                    ax2.scatter(n_orig_r2[2],np.divide(z_orig_r2[2],1e5),color='red',marker='x')
                    ax2.scatter(n_orig_r2[0],np.divide(z_orig_r2[0],1e5),color='red',marker='x')
                if mechanism_R2s[m] == 'C' or mechanism_R2s[m] == 'N' or mechanism_R2s[m] == 'Ar' or mechanism_R2s[m] == 'O2':
                    ax2.scatter(n_orig_r2[len(n_orig_r2)-1],np.divide(z_orig_r2[len(z_orig_r2)-1],1e5),color='red',marker='x')
                    ax2.scatter(n_orig_r2[len(n_orig_r2)-3],np.divide(z_orig_r2[len(n_orig_r2)-3],1e5),color='red',marker='x')
                ax2.set_xlabel('%s number density (cm$^{-3}$)' % (mechanism_R2s[m]))


            elif fig_choice == 12:
                if SA == 'LSA' and (mech in LSA_5eV_top or mech in LSA_0_2eV_top):
                    if mechanism_R1s[m] not in reactants_LSA: # and mechanism_R1s[m] not in ['O','CO','CO2','N2']: # uncomment if including bg densities elsewhere in the paper
                        axb4.plot(n_new_r1,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R1s[m],color=LSA_colors[mechanism_R1s[m]])
                        reactants_LSA.append(mechanism_R1s[m])
                    if mechanism_R2s[m] not in reactants_LSA: # and mechanism_R2s[m] not in ['O','CO','CO2','N2']: # uncomment if including bg densities elsewhere in the paper
                        if mechanism_R2s[m] == 'e':
                            axb4.plot(n_new_r2,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R2s[m],color=LSA_colors[mechanism_R2s[m]],linestyle='--')
                        else:
                            axb4.plot(n_new_r2,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R2s[m],color=LSA_colors[mechanism_R2s[m]])          
                        reactants_LSA.append(mechanism_R2s[m])
                if SA in 'HSA' and (mech in HSA_5eV_top or mech in HSA_0_2eV_top):
                    if mechanism_R1s[m] not in reactants_HSA: # and mechanism_R1s[m] not in ['O','CO','CO2','N2']: # uncomment if including bg densities elsewhere in the paper                
                        axb8.plot(n_new_r1,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R1s[m],color=HSA_colors[mechanism_R1s[m]])
                        reactants_HSA.append(mechanism_R1s[m])
                    if mechanism_R2s[m] not in reactants_HSA: # and mechanism_R2s[m] not in ['O','CO','CO2','N2']:  # uncomment if including bg densities elsewhere in the paper
                        if mechanism_R2s[m] == 'e':
                            axb8.plot(n_new_r2,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R2s[m],color=HSA_colors[mechanism_R2s[m]],linestyle='--')
                        else:
                            axb8.plot(n_new_r2,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R2s[m],color=HSA_colors[mechanism_R2s[m]])
                            if fig_choice == 9:
                                axb4.plot(n_new_r2,np.divide(z_new,1e5),linewidth=0.5,label=mechanism_R2s[m],color=HSA_colors[mechanism_R2s[m]])
                        reactants_HSA.append(mechanism_R2s[m])
                # Also plot OCOH+ DR on HSA panel and N2+ in LSA panel:
                if SA == 'HSA' and  mechanism_R1s[m] == 'OCOH+' and mechanism_R1s[m] not in reactants_HSA:
                    axb8.plot(n_new_r1,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R1s[m],color=LSA_colors[mechanism_R1s[m]])
                    reactants_HSA.append('OCOH+')
                if SA == 'LSA' and mechanism_R1s[m] == 'N2+' and mechanism_R1s[m] not in reactants_LSA:
                    axb4.plot(n_new_r1,np.divide(z_new,1e5),linewidth=1.5,label=mechanism_R1s[m],color=HSA_colors[mechanism_R1s[m]])
                    reactants_LSA.append('N2+')

                    

            ### Calculate production rate for new grid ###
            for z, alt in enumerate(z_new):
                rate_coefficient = find_rate_coefficient(mechanism_labels[m],T_type,Tn_new,Te_new,Ti_new,z)
                if mech == 'H3+ DR \u2192 H + H + H':
                    production_rate.append(3*rate_coefficient*n_new_r1[z]*n_new_r2[z]) # because 3 H atoms produced (NB. each has the chosen energy - it isn't divided)
                elif mech == 'H2+ DR \u2192 H + H':
                    production_rate.append(2*rate_coefficient*n_new_r1[z]*n_new_r2[z]) # because 2 H atoms produced (NB. each has the chosen energy = it isn't divided)
                else:
                    production_rate.append(rate_coefficient*n_new_r1[z]*n_new_r2[z])

                if SA == 'LSA':
                    if mech == 'H+ + O \u2192 O+ + H':
                        esc_estimate_5eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z])
                        esc_estimate_0_2eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z])
                        esc_estimate_1eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z])
                        esc_estimate_10eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z])
                        esc_estimate_0_5eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z])
                        weighted_escape_5eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z]*bin_size)
                        total_prod_LSA += O_remnant_fraction_new[z]*production_rate[z]*bin_size
                    elif mech == 'CO2 + H+ \u2192 CO2+ + H':
                        esc_estimate_5eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z])
                        esc_estimate_0_2eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z])
                        esc_estimate_1eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z])
                        esc_estimate_10eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z])
                        esc_estimate_0_5eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z])
                        weighted_escape_5eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z]*bin_size)
                        total_prod_LSA += CO2_remnant_fraction_new[z]*production_rate[z]*bin_size
                    elif mech == 'H + H+ \u2192 H+ + H':
                        esc_estimate_5eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z])
                        esc_estimate_0_2eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z])
                        esc_estimate_1eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z])
                        esc_estimate_10eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z])
                        esc_estimate_0_5eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z])
                        weighted_escape_5eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_LSA_0_5eV[z]*bin_size)
                        total_prod_LSA += H_remnant_fraction_new[z]*production_rate[z]*bin_size
                        
                    else:
                        esc_estimate_5eV.append(production_rate[z]*returned_escape_rates_LSA_5eV[z])
                        esc_estimate_0_2eV.append(production_rate[z]*returned_escape_rates_LSA_0_2eV[z])
                        esc_estimate_1eV.append(production_rate[z]*returned_escape_rates_LSA_1eV[z])
                        esc_estimate_10eV.append(production_rate[z]*returned_escape_rates_LSA_10eV[z])
                        esc_estimate_0_5eV.append(production_rate[z]*returned_escape_rates_LSA_0_5eV[z])
                        weighted_escape_5eV += (production_rate[z]*returned_escape_rates_LSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (production_rate[z]*returned_escape_rates_LSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (production_rate[z]*returned_escape_rates_LSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (production_rate[z]*returned_escape_rates_LSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (production_rate[z]*returned_escape_rates_LSA_0_5eV[z]*bin_size)
                        total_prod_LSA += production_rate[z]*bin_size
                           
                elif SA == 'HSA':
                    if mech == 'H+ + O \u2192 O+ + H':
                        esc_estimate_5eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z])
                        esc_estimate_0_2eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z])
                        esc_estimate_1eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z])
                        esc_estimate_10eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z])
                        esc_estimate_0_5eV.append(O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z])
                        weighted_escape_5eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (O_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z]*bin_size)
                        total_prod_HSA += O_remnant_fraction_new[z]*production_rate[z]*bin_size
                    elif mech == 'CO2 + H+ \u2192 CO2+ + H':
                        esc_estimate_5eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z])
                        esc_estimate_0_2eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z])
                        esc_estimate_1eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z])
                        esc_estimate_10eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z])
                        esc_estimate_0_5eV.append(CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z])
                        weighted_escape_5eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (CO2_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z]*bin_size)
                        total_prod_HSA += CO2_remnant_fraction_new[z]*production_rate[z]*bin_size
                    elif mech == 'H + H+ \u2192 H+ + H':
                        esc_estimate_5eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z])
                        esc_estimate_0_2eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z])
                        esc_estimate_1eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z])
                        esc_estimate_10eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z])
                        esc_estimate_0_5eV.append(H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z])
                        weighted_escape_5eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (H_remnant_fraction_new[z]*production_rate[z]*returned_escape_rates_HSA_0_5eV[z]*bin_size)
                        total_prod_HSA += H_remnant_fraction_new[z]*production_rate[z]*bin_size
                    else:
                        esc_estimate_5eV.append(production_rate[z]*returned_escape_rates_HSA_5eV[z])
                        esc_estimate_0_2eV.append(production_rate[z]*returned_escape_rates_HSA_0_2eV[z])
                        esc_estimate_1eV.append(production_rate[z]*returned_escape_rates_HSA_1eV[z])
                        esc_estimate_10eV.append(production_rate[z]*returned_escape_rates_HSA_10eV[z])
                        esc_estimate_0_5eV.append(production_rate[z]*returned_escape_rates_HSA_0_5eV[z])
                        weighted_escape_5eV += (production_rate[z]*returned_escape_rates_HSA_5eV[z]*bin_size)
                        weighted_escape_0_2eV += (production_rate[z]*returned_escape_rates_HSA_0_2eV[z]*bin_size)
                        weighted_escape_1eV += (production_rate[z]*returned_escape_rates_HSA_1eV[z]*bin_size)
                        weighted_escape_10eV += (production_rate[z]*returned_escape_rates_HSA_10eV[z]*bin_size)
                        weighted_escape_0_5eV += (production_rate[z]*returned_escape_rates_HSA_0_5eV[z]*bin_size)
                        total_prod_HSA += production_rate[z]*bin_size
                        

            if fig_choice == 1:
                ax4.plot(production_rate,np.divide(z_new,1e5),color=color,linewidth=1)
                peak_prod_alt = peak_prod_altitude(production_rate,z_new)[0]
                peak_prod = peak_prod_altitude(production_rate,z_new)[1]
                if mech == 'CH+ + H2 \u2192 CH2+ + H':
                    peak_prod_alt = 22150000 # hardcoding because there are several altitudes that have the (same) max value. This is in the middle. The default, now overwritten, was the smallest.
                ax4.scatter(peak_prod,np.divide(peak_prod_alt,1e5),marker='x',linewidth=0.5)
                ax4.set_xlabel('Production rate (cm$^{-3}$ s$^{-1}$)')

                # Find altitude below which 90% of escaping H is produced
                if SA == 'LSA':
                    percent_escaping_alt_5eV = percent_escaping_alt(weighted_escape_5eV,returned_escape_rates_LSA_5eV,production_rate,0.9,z_new,bin_size)
                    percent_escaping_alt_0_2eV = percent_escaping_alt(weighted_escape_0_2eV,returned_escape_rates_LSA_0_2eV,production_rate,0.9,z_new,bin_size)
                elif SA == 'HSA':
                    percent_escaping_alt_5eV = percent_escaping_alt(weighted_escape_5eV,returned_escape_rates_HSA_5eV,production_rate,0.9,z_new,bin_size)
                    percent_escaping_alt_0_2eV = percent_escaping_alt(weighted_escape_0_2eV,returned_escape_rates_HSA_0_2eV,production_rate,0.9,z_new,bin_size)
    
                # Print to output csv file
                with open('../model_output_data/escape_estimates/escape_estimates_%s.csv' % (SA),'a') as f1:
                    if SA == 'HSA':
                        f1.write('%s,%s,%s,%s,%s,  %.2E,%.2E,%s,%s,%s,%s,%s,%s,%s,%s\n' % (mechanism_R1s[m],mechanism_R2s[m],mechanism_labels[m],weighted_escape_5eV,weighted_escape_0_2eV,  weighted_escape_5eV, weighted_escape_0_2eV,peak_prod_alt,total_prod_HSA,percent_escaping_alt_5eV,percent_escaping_alt_0_2eV,weighted_escape_1eV, weighted_escape_10eV,weighted_escape_0_5eV,excess_energies[m]))
                    if SA == 'LSA':
                        f1.write('%s,%s,%s,%s,%s,  %.2E,%.2E,%s,%s,%s,%s,%s,%s,%s,%s\n' % (mechanism_R1s[m],mechanism_R2s[m],mechanism_labels[m],weighted_escape_5eV,weighted_escape_0_2eV,  weighted_escape_5eV, weighted_escape_0_2eV,peak_prod_alt,total_prod_LSA,percent_escaping_alt_5eV,percent_escaping_alt_0_2eV,weighted_escape_1eV, weighted_escape_10eV,weighted_escape_0_5eV,excess_energies[m]))


            elif fig_choice == 8  or fig_choice == 12 or fig_choice == 13 or fig_choice == 14:
                if SA == 'LSA': # Use this if want to plot total hot H production from 47 mechs using fig_choice == 13
                    accumulated_prod_LSA.append(production_rate)
                elif SA == 'HSA':
                    accumulated_prod_HSA.append(production_rate)
                if SA == 'LSA' and (mech in LSA_5eV_top or mech in LSA_0_2eV_top):

                    if fig_choice != 13:
                        axb5.plot(production_rate,np.divide(z_new,1e5),linewidth=2,label=m,color=colors[mech_count_LSA]) #mech)
                    if fig_choice == 8:
                        axb5.axhspan(146.5,173.75,color='lightgrey')

                elif SA == 'HSA' and (mech in HSA_5eV_top or mech in HSA_0_2eV_top):
                    if mech == 'N2+ + H2 \u2192 N2H+ + H':
                        if fig_choice != 13:
                            axb9.plot(production_rate,np.divide(z_new,1e5),linewidth=2,label=m,color=colors[5]) #mech)
                        if fig_choice == 9:
                            axb5.plot(production_rate,np.divide(z_new,1e5),linewidth=0.5,label=m,color=colors[5]) #mech)
                    else:
                        if fig_choice != 13 and fig_choice != 8 and fig_choice != 14:
                            axb9.plot(production_rate,np.divide(z_new,1e5),linewidth=2,label=m,color=colors[mech_count_HSA]) #mech)
                        if fig_choice == 9:
                            axb5.plot(production_rate,np.divide(z_new,1e5),linewidth=0.5,label=m,color=colors[mech_count_HSA]) #mech)
                        if fig_choice == 14:
                            axb5.plot(production_rate,np.divide(z_new,1e5),linewidth=0.75,label=m,color=colors[mech_count_HSA]) #mech)

            if fig_choice == 1:
                ax6.plot(esc_estimate_5eV,np.divide(z_new,1e5),color='tab:orange',linewidth=1,label='5eV')
                ax6.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),color='tab:green',linewidth=1,label='0.2 eV')    
                ax6.set_xlabel('Weighted escape estimate (cm$^{-3}$ s$^{-1}$)')
                ax6.set_xlim([1e-9,1e1])
                
            elif  fig_choice == 8  or fig_choice == 12 or fig_choice == 13 or fig_choice == 14:
                if SA == 'LSA' and (mech in LSA_5eV_top or mech in LSA_0_2eV_top):
                    if LSA_5eV_top == LSA_0_2eV_top:
                        if fig_choice == 8:
                            if mech == 'HCO+ + e- \u2192 CO + H':
                                gradient_colors = gradient_colors_HCOpl
                            else:
                                gradient_colors = gradient_colors_Opl
                            axb7.plot(esc_estimate_10eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2.75,color=gradient_colors[4]) #color=colors[mech_count_LSA])
                            axb7.plot(esc_estimate_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2.25,color=gradient_colors[3]) #colors[mech_count_LSA])
                            axb7.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=1.75,color=gradient_colors[2]) #color=colors[mech_count_LSA])
                            axb7.plot(esc_estimate_0_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=1.25,color=gradient_colors[1]) #color=colors[mech_count_LSA])
                            axb7.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=0.75,color=gradient_colors[0]) #color=colors[mech_count_LSA])
                            
                        elif fig_choice == 14:
                            axb7.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=2,color=gradient_colors[mech_count_LSA]) #color=colors[mech_count_LSA])
                            
                        elif fig_choice == 12 or fig_choice == 13:
                            if mech == 'HCO+ + e- \u2192 CO + H':
                                if fig_choice == 12:
                                    axb7.plot(esc_estimate_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_LSA]) # 5eV version
                                accumulated_escape_LSA.append(esc_estimate_5eV)
                                accumulated_colors_LSA.append(colors[mech_count_LSA])
                            elif mech == 'CO2+ + H2 \u2192 OCOH+ + H':
                                if fig_choice == 12:
                                    axb7.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_LSA]) # 1 eV
                                accumulated_escape_LSA.append(esc_estimate_1eV)
                                accumulated_colors_LSA.append(colors[mech_count_LSA])
                            elif mech == 'O+(4S) + H2 \u2192 OH+ + H':
                                if fig_choice == 12:
                                    axb7.plot(esc_estimate_0_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_LSA]) # 0.5 eV
                                accumulated_escape_LSA.append(esc_estimate_0_5eV)
                                accumulated_colors_LSA.append(colors[mech_count_LSA])
                            elif mech == 'OH+ + O \u2192 O2+ + H':
                                if fig_choice == 12:
                                    axb7.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_LSA]) # 1 eV
                                accumulated_escape_LSA.append(esc_estimate_1eV)
                                accumulated_colors_LSA.append(colors[mech_count_LSA])
                            elif mech == 'OCOH+ DR \u2192 O + CO + H':
                                if fig_choice == 12:
                                    axb7.plot(esc_estimate_1eV,np.divide(z_new,1e5),linewidth=2,color=colors[mech_count_LSA]) # 1 eV
                                accumulated_escape_LSA.append(esc_estimate_1eV)
                                accumulated_colors_LSA.append(colors[mech_count_LSA])
                        else:
                            axb7.plot(esc_estimate_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_LSA])
                            if fig_choice != 11:
                                axb7.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=1,color=colors[mech_count_LSA])
                        

                    else:
                        axb7.plot(esc_estimate_5eV,np.divide(z_new,1e5),linewidth=1,label='%s 5eV' % (mech),linestyle='-',color=colors[mech_count_LSA])
                        axb7.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=1,label='%s 0.2 eV' % (mech),linestyle='--',color=colors[mech_count_LSA])
                    mech_count_LSA += 1
                elif SA == 'HSA' and (mech in HSA_5eV_top or mech in HSA_0_2eV_top):
                    if fig_choice == 14:
                        axb7.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=0.75,color=gradient_colors[mech_count_HSA]) #color=colors[mech_count_LSA])
                    if mech=='N2+ + H2 \u2192 N2H+ + H':
                        if fig_choice == 12 or fig_choice == 13:
                            if fig_choice == 12:
                                axb11.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[5]) # 1 eV
                            accumulated_escape_HSA.append(esc_estimate_1eV)
                            accumulated_colors_HSA.append(colors[5])
                        else:
                            axb11.plot(esc_estimate_5eV,np.divide(z_new,1e5),linewidth=2,label=mech,color=colors[5])
                            if fig_choice != 11:
                                axb11.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=1,label=mech,color=colors[5])
                            
                    else:
                                                       
                        if fig_choice == 12 or fig_choice == 13:
                            if mech == 'HCO+ + e- \u2192 CO + H':
                                if fig_choice == 12:
                                    axb11.plot(esc_estimate_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_HSA]) # 5eV version
                                accumulated_escape_HSA.append(esc_estimate_5eV)
                                accumulated_colors_HSA.append(colors[mech_count_HSA])
                            elif mech == 'CO2+ + H2 \u2192 OCOH+ + H':
                                if fig_choice == 12:
                                    axb11.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_HSA]) # 1 eV
                                accumulated_escape_HSA.append(esc_estimate_1eV)
                                accumulated_colors_HSA.append(colors[mech_count_HSA])
                            elif mech == 'O+(4S) + H2 \u2192 OH+ + H':
                                if fig_choice == 12:
                                    axb11.plot(esc_estimate_0_5eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_HSA]) # 0.5 eV
                                accumulated_escape_HSA.append(esc_estimate_0_5eV)
                                accumulated_colors_HSA.append(colors[mech_count_HSA])
                            elif mech == 'OH+ + O \u2192 O2+ + H':
                                if fig_choice == 12:
                                    axb11.plot(esc_estimate_1eV,np.divide(z_new,1e5),label='%s' % (mech),linewidth=2,color=colors[mech_count_HSA]) # 1 eV
                                accumulated_escape_HSA.append(esc_estimate_1eV)
                                accumulated_colors_HSA.append(colors[mech_count_HSA])
                                
                        elif fig_choice != 8 and fig_choice != 14:
                            axb11.plot(esc_estimate_5eV,np.divide(z_new,1e5),linewidth=2,label=mech,color=colors[mech_count_HSA])
                            if fig_choice != 11:
                                axb11.plot(esc_estimate_0_2eV,np.divide(z_new,1e5),linewidth=1,label=mech,color=colors[mech_count_HSA])
                            
                        
                    mech_count_HSA += 1

            ### Fig 1 plot aesthetics ###
            if fig_choice == 1:
                for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
                    ax.set_ylim([80,400])
                    ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
                    ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
                for ax in [ax1,ax2,ax3,ax4,ax6]:
                    ax.set_xscale('log')
                for ax in [ax2,ax3,ax4,ax5,ax6]:
                    for tick in ax.yaxis.get_majorticklabels():
                        tick.set_color("none")
                for ax in [ax5,ax6]:
                    ax.legend(loc='upper left')
                if mech == 'H + H+ \u2192 H+ + H':
                    ax3.legend()

                plt.gcf().subplots_adjust(hspace=0.18,wspace=0.02,top=0.95,bottom=0.12,left=0.05,right=0.99)
                plt.savefig('../model_output_data/escape_estimates/process_figs/%s_%s.png' % (mech,SA))
                
        ### Other Fig 12 bits
        if  fig_choice == 12:
            [bg_z,CO2_density,O_density,CO_density,N2_density] = find_bg_dens(SA,planet)
            col_dens_above = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,SA,planet)
        

        
        ### Fig 8, 12, 14 plot aesthetics
        if fig_choice == 8 or fig_choice == 12 or fig_choice == 14:       
            if fig_choice == 8 or fig_choice == 14:
                for ax in [axb5,axb6,axb7]:
                    ax.set_ylim([110,400])
            else:
                for ax in [axb4,axb5,axb6,axb7,axb8,axb9,axb10,axb11]:
                    ax.set_ylim([110,400])


            if fig_choice != 8 and fig_choice != 14:
                for ax in [axb4,axb8]:
                    ax.set_ylabel('Altitude (km)',fontsize=14)
                
                    ax.set_xscale('log')
                for ax in [axb8]:
                    ax.set_xlabel('Number density (cm$^{-3}$)',fontsize=14)

                for ax in [axb4,axb8]:
                    ax.set_xlim([1e-2,1e10]) # was [1e-7,1e9] before; u
                    ax.xaxis.set_minor_locator(ticker.FixedLocator([1e0,1e2,1e4,1e6,1e8,1e10]))
                    for tick in ax.xaxis.get_minorticklabels():
                        tick.set_color('none')
                        
                axb9.set_xscale('log')
                axb9.set_xlim([3e-3,20]) #1e2]) #1e-4 when H+ + H included

                axb10.set_xlim([-0.01,1.01])

            axb5.set_xscale('log')
            axb5.set_xlim([3e-3,20]) #1e2]) #1e-4 when H+ + H included

            axb6.set_xlim([-0.01,1.01])
                
            if fig_choice != 8 and fig_choice != 14:
                for ax in [axb7,axb11]:
                    ax.set_xscale('log')
                    ax.set_xlim([4e-4,6]) # was 3e-3 before

            if fig_choice !=8 and fig_choice != 12 and fig_choice != 14:
                    ax.set_xlabel('Escaping H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)

            if fig_choice == 12:
                axb11.set_xlabel('Escaping H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                for ax in [axb4,axb5,axb6,axb7,axb8,axb9,axb10,axb11]:
                    ax.axhspan(80,400, facecolor='w', alpha=1,zorder=-20)
                    

            if fig_choice == 8:
                axb5.set_ylabel('Altitude (km)',fontsize=14)
                axb7.set_xscale('log')
                axb7.set_xlim([3e-3,7])
                axb5.set_xlabel('H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                axb6.set_xlabel('Escape probability',fontsize=14)
                axb7.set_xlabel('Escaping H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                
            elif fig_choice == 14:
                axb5.set_ylabel('Altitude (km)',fontsize=14)
                axb7.set_xscale('log')
                axb7.set_xlim([3e-8,4e-2])
                axb5.set_xlim([3e-5,0.3])
                axb12.set_xlim([-0.01,0.88])
                axb12.set_ylim([110,400])
                axb6.set_xlim([-0.01,0.9])
                axb5.set_xlabel('H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                axb6.set_xlabel('Escape probability',fontsize=14)
                axb7.set_xlabel('Escaping H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                axb12.set_xlabel('Fraction produced above E$_{\mathrm{esc}}$',fontsize=14)
                axb5.text(0.05,0.92,'(a)',fontsize=16,fontweight='bold',transform=axb5.transAxes)
                axb6.text(0.05,0.92,'(b)',fontsize=16,fontweight='bold',transform=axb6.transAxes)
                axb7.text(0.05,0.92,'(d)',fontsize=16,fontweight='bold',transform=axb7.transAxes)
                axb12.text(0.05,0.92,'(c)',fontsize=16,fontweight='bold',transform=axb12.transAxes)

            else:
                axb9.set_xlabel('H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=14)
                axb10.set_xlabel('Escape probability',fontsize=14)
                
                axb5.set_title('Production rates\n(top 5 mechanisms)',fontsize=14)
                axb6.set_title('Escape probability\nprofiles',fontsize=14)
                axb7.set_title('Escaping H production rates \n(top 5 mechanisms)',fontsize=14)

            if fig_choice != 8 and fig_choice != 14:
                axb4.set_title('Reactant number\ndensity profiles',fontsize=14)

        
            if fig_choice == 12:
                for ax in [axb4,axb5,axb6,axb8,axb9,axb10]:
                    ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
                    ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
            elif fig_choice == 8 or fig_choice == 14:
                for ax in [axb5,axb7]:
                    ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
                    ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
                axb6.tick_params('both',width=1,length=8,direction='in',top=True,right=True,pad=5)
                axb6.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True,pad=5)
                if fig_choice == 14:
                    axb12.tick_params('both',width=1,length=8,direction='in',top=True,right=True,pad=5)
                    axb12.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True,pad=5)
    
            if fig_choice != 8 and fig_choice != 14:
                for ax in [axb7,axb11]:
          #      ax.legend(fontsize=10,bbox_to_anchor=(1,1.05))
                    ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
                    ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)

            if  fig_choice == 9 or fig_choice == 12:
                for ax in [axb5,axb6,axb7,axb9,axb10,axb11]:
                    for tick in ax.yaxis.get_majorticklabels():
                        tick.set_color('none')
            elif fig_choice == 8 or fig_choice ==14:
                for tick in axb5.yaxis.get_majorticklabels():
                    tick.set_fontsize(14)
                for ax in [axb6,axb7]:
                    for tick in ax.yaxis.get_majorticklabels():
                        tick.set_color('none')
                if fig_choice == 14:
                    for tick in axb12.yaxis.get_majorticklabels():
                        tick.set_color('none')
                    for tick in axb12.xaxis.get_majorticklabels():
                        tick.set_fontsize(14)
                for ax in [axb5,axb6,axb7]:
                    for tick in ax.xaxis.get_majorticklabels():
                        tick.set_fontsize(14)

            if fig_choice != 8 and fig_choice != 14:
                for tick in axb4.xaxis.get_majorticklabels():
                    tick.set_color('none')
                for ax in [axb4,axb8]:
                    for tick in ax.yaxis.get_majorticklabels():
                        tick.set_fontsize(14)
                for tick in axb8.xaxis.get_majorticklabels():
                    tick.set_fontsize(14)
                    
            if  fig_choice != 8 and fig_choice != 14:
                for ax in [axb5,axb6,axb7]:
                    for tick in ax.xaxis.get_majorticklabels():
                        tick.set_color('none')

                for ax in [axb9,axb10,axb11]:
                    for tick in ax.xaxis.get_majorticklabels():
                        tick.set_fontsize(14)


            if fig_choice == 12:
                axb4.text(1e-6,250,'Low\nsolar\nactivity',fontweight='bold',fontsize=18,ha='center')
                axb8.text(1e-6,250,'High\nsolar\nactivity',fontweight='bold',fontsize=18,ha='center')

            if fig_choice == 12:
                # Figure letters
                axb4.text(8e9,375,'(a)',fontweight='bold',fontsize=16,ha='right')
                axb5.text(17,375,'(b)',fontweight='bold',fontsize=16,ha='right')
                axb6.text(0.01,375,'(c)',fontweight='bold',fontsize=16)
                axb7.text(5.5,375,'(d)',fontweight='bold',fontsize=16,ha='right')

                axb8.text(8e9,375,'(e)',fontweight='bold',fontsize=16,ha='right')
                axb9.text(17,375,'(f)',fontweight='bold',fontsize=16,ha='right')
                axb10.text(0.01,375,'(g)',fontweight='bold',fontsize=16)
                axb11.text(5.5,375,'(h)',fontweight='bold',fontsize=16,ha='right')

                ### Line labels
                # Reactant density line labels
                axb4.text(2e-2,180,'OH$^+$',fontsize=14,color=LSA_colors['OH+'])
                axb4.text(1e3,375,'e$^-$',fontsize=14,color=LSA_colors['e'])
                axb4.text(2e4,310,'H$_2$',fontsize=14,color=LSA_colors['H2'])
                axb4.text(7e6,250,'O',fontsize=14,color=LSA_colors['O'])
                axb4.text(3e-2,230,'OCOH$^+$',fontsize=14,color=LSA_colors['OCOH+'])
                axb4.text(2e-1,160,'HCO$^+$',fontsize=14,color=LSA_colors['HCO+'])
                axb4.plot([4.62,613],[164.3,164.3],color=LSA_colors['HCO+'],linewidth=0.3)
                axb4.text(2e-2,125,'N$_2^+$',fontsize=14,color=LSA_colors['N2+'])
                axb4.text(25,115,'O$^+$',fontsize=14,color=LSA_colors['O+'])
                axb4.text(6e2,140,'CO$_2^+$',fontsize=14,color=LSA_colors['CO2+'])
                
                axb8.text(2e-2,300,'OCOH$^+$',fontsize=14,color=LSA_colors['OCOH+'])
                axb8.text(1e-1,220,'OH$^+$',fontsize=14,color=LSA_colors['OH+'])
                axb8.text(8e3,310,'e$^-$',fontsize=14,color=LSA_colors['e'])
                axb8.text(2e5,350,'H$_2$',fontsize=14,color=LSA_colors['H2'])
                axb8.text(3e7,280,'O',fontsize=14,color=LSA_colors['O'])
                axb8.text(1.1e-1,145,'HCO$^+$',fontsize=14,color=LSA_colors['HCO+'])
                axb8.plot([2.67683,178.4],[145.3,122.0],color=LSA_colors['HCO+'],linewidth=0.3)
                axb8.text(1.5e1,185,'N$_2^+$',fontsize=14,color=LSA_colors['N2+'])
                axb8.text(32,145,'O$^+$',fontsize=14,color=LSA_colors['O+'])
                axb8.text(6e2,138,'CO$_2^+$',fontsize=14,color=LSA_colors['CO2+'])

                # Production rate line labels
                axb5.text(2.2e-2,350,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2])
                axb5.text(1e-1,300,'CO$_2^+$ + H$_2$',fontsize=14,color=colors[3])
                axb5.plot([9e-2,0.0078],[310,310],color=colors[3],linewidth=0.4)
                axb5.text(1e0,195,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0])
                axb5.text(1e-0,115,'OCOH$^+$\n+ e$^-$',fontsize=14,color=colors[1])
                axb5.text(4e-3,142,'OH$^+$ + O',fontsize=14,color=colors[4])

                axb9.text(1.5e-1,350,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2])
                axb9.text(6e-1,130,'CO$_2^+$ + H$_2$',fontsize=14,color=colors[3])
                axb9.text(1.5e0,200,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0])
                axb9.text(9e-3,230,'N$_2^+$ + H$_2$',fontsize=14,color=colors[5])
                axb9.text(3.1e-3,183,'OH$^+$ + O',fontsize=14,color=colors[4])


            if fig_choice == 8:

                # Figure letters
                axb5.text(17,375,'(a)',fontweight='bold',fontsize=16,ha='right')
                axb6.text(0.01,375,'(b)',fontweight='bold',fontsize=16)
                axb7.text(5.5,375,'(c)',fontweight='bold',fontsize=16,ha='right')

                ### Line labels
                # Production rate line labels
                axb5.text(2.2e-2,350,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[1])
                axb5.text(1e0,195,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0])

                # Multiplication and equals signs
                fig8.text(0.37,0.6,'×',fontsize=80,va='center',ha='center',color='k',fontweight=540)
                fig8.text(0.68,0.6,'=',fontsize=80,va='center',ha='center',color='k',weight='light')

                # Escape probability line labels
                axb6.text(0.18,200,'0.2 eV',fontsize=14,color=gradient_colors_grey[0])
                axb6.plot([0.249,0.295],[194.9,163.9],color=gradient_colors_grey[0],linewidth=0.4)
                axb6.text(0.4,195,'0.5 eV',fontsize=14,color=gradient_colors_grey[1])
                axb6.plot([0.487,0.529],[191.1,163],color=gradient_colors_grey[1],linewidth=0.4)
                axb6.text(0.75,140,'1 eV',fontsize=14,color=gradient_colors_grey[2])
                axb6.plot([0.811,0.779],[158.6,174.4],color=gradient_colors_grey[2],linewidth=0.4)
                axb6.text(0.45,124,'5 eV',fontsize=14,color=gradient_colors_grey[3]) #,fontweight='bold')
                axb6.plot([0.504,0.504],[137.5,152],color=gradient_colors_grey[3],linewidth=0.4)
                axb6.text(0.2,120,'10 eV',fontsize=14,color=gradient_colors_grey[4])

                axb7.text(7e-2,350,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[1],ha='center',fontweight='bold')
                axb7.text(1e-2,200,'0.2 eV',fontsize=14,color=gradient_colors_Opl[0],ha='center')
                axb7.plot([0.0178,0.0534],[196,173.9],color=gradient_colors_Opl[0],linewidth=0.4)
                axb7.text(1e-2,170,'10 eV',fontsize=14,color=gradient_colors_Opl[4],ha='center',fontweight='bold')
                axb7.plot([0.0178,0.0510],[167.4,148.3],color=gradient_colors_Opl[4],linewidth=0.4)

                axb7.text(0.38,170,'0.2 eV',fontsize=14,color=gradient_colors_HCOpl[0],ha='center')
                axb7.plot([0.4011,0.4011],[165.6,151.1],color=gradient_colors_HCOpl[0],linewidth=0.4)
                axb7.text(3e-0,115,'10 eV',fontsize=14,color=gradient_colors_HCOpl[4],ha='center',fontweight='bold')
                axb7.plot([3.342,3.342],[127.4,134],color=gradient_colors_HCOpl[4],linewidth=0.4)
                
                axb7.text(1e0,200,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0],ha='center',fontweight='bold')

            if fig_choice == 14:
                axb6.text(0.15,140,'LSA',fontsize=14,color='k',fontweight='bold')
                axb6.text(0.15,195,'HSA',fontsize=14,color='k')

                axb5.text(4e-4,210,'H + H$^+$',fontsize=14,color='rebeccapurple')
                axb5.plot([0.00076,0.000422],[226,242],linewidth=0.5,color='k')
                axb5.plot([0.00076,0.00131],[206,195],linewidth=0.5,color='k')
                axb5.text(0.04,300,'O + H$^+$',fontsize=14,color='mediumseagreen')
                axb5.plot([0.0846,0.0547],[294,285],linewidth=0.5,color='k')
                axb5.plot([0.0846,0.0547],[294,261],linewidth=0.5,color='k')

                axb12.text(0.3,300,'O + H$^+$',fontsize=14,color='mediumseagreen')
                axb12.text(0.3,210,'H + H$^+$',fontsize=14,color='rebeccapurple')

                axb7.text(1e-6,240,'H + H$^+$',fontsize=14,color='rebeccapurple')
                axb7.plot([3.039e-6,3.039e-6],[234.3,216.3],linewidth=0.5,color='k')
                axb7.plot([3.039e-6,9.7e-7],[234.3,184],linewidth=0.5,color='k')
                axb7.text(1e-3,180,'O + H$^+$',fontsize=14,color='mediumseagreen')
                axb7.plot([0.00314,0.00314],[190,211],linewidth=0.5,color='k')
                axb7.plot([0.00314,0.0104],[190,238],linewidth=0.5,color='k')

                
            if fig_choice == 12:
                # Escape probability line labels
                axb6.text(0.15,170,'0.2 eV',fontsize=14,color=gradient_colors_grey[0])
                axb6.text(0.4,195,'0.5 eV',fontsize=14,color=gradient_colors_grey[1])
                axb6.plot([0.487,0.529],[191.1,163],color=gradient_colors_grey[1],linewidth=0.4)
                axb6.text(0.75,140,'1 eV',fontsize=14,color=gradient_colors_grey[2])
                axb6.plot([0.811,0.779],[158.6,174.4],color=gradient_colors_grey[2],linewidth=0.4)
                axb6.text(0.45,124,'5 eV',fontsize=14,color=gradient_colors_grey[3]) #,fontweight='bold')
                axb6.plot([0.504,0.504],[137.5,152],color=gradient_colors_grey[3],linewidth=0.4)
                axb6.text(0.2,120,'10 eV',fontsize=14,color=gradient_colors_grey[4])

                axb10.text(0.15,205,'0.2 eV',fontsize=14,color=gradient_colors_grey[0])
                axb10.text(0.4,240,'0.5 eV',fontsize=14,color=gradient_colors_grey[1],ha='center')
                axb10.plot([0.473,0.622],[236.4,201.7],color=gradient_colors_grey[1],linewidth=0.4)
                axb10.text(0.7,155,'1 eV',fontsize=14,color=gradient_colors_grey[2])
                axb10.plot([0.757,0.697],[173.6,200.9],color=gradient_colors_grey[2],linewidth=0.4)
                axb10.text(0.5,135,'5 eV',fontsize=14,color=gradient_colors_grey[3])
                axb10.plot([0.560,0.560],[152.1,172],color=gradient_colors_grey[3],linewidth=0.4)
                axb10.text(0.2,125,'10 eV',fontsize=14,color=gradient_colors_grey[4]) #,fontweight='bold')

                # Escaping H production rate line labels
                axb7.text(2e-2,350,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2])
                axb7.text(4.1e-4,375,'CO$_2^+$+ H$_2$',fontsize=14,color=colors[3])
                axb7.text(3e-1,205,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0])
                axb7.text(6e-4,193,'OCOH$^+$ + e$^-$',fontsize=14,color=colors[1])
                axb7.text(1e-1,280,'OH$^+$ + O',fontsize=14,color=colors[4])
                axb7.plot([0.1033,0.0468],[278,257.8],color=colors[4],linewidth=0.4)

                axb11.text(1.3e-1,370,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2])
                axb11.text(2.01e-1,170,'CO$_2^+$ + H$_2$',fontsize=14,color=colors[3],ha='center')
                axb11.text(2e-1,120,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0])
                axb11.text(6e-3,230,'N$_2^+$ + H$_2$',fontsize=14,color=colors[5])
                axb11.text(3e-1,300,'OH$^+$ + O',fontsize=14,color=colors[4])
                axb11.plot([0.2666,0.116677],[307.5,307.5],color=colors[4],linewidth=0.4)



    if fig_choice == 13:
        from weighted_esc_pie import weighted_esc_pie
        fractions_LSA = weighted_esc_pie()[0]
        fractions_HSA = weighted_esc_pie()[1]

        print(fractions_LSA,fractions_HSA)
        z_new_km = np.divide(z_new,1e5)
        accumulated_escape_LSA_old = np.zeros(len(z_new_km))
        accumulated_escape_LSA_new = np.zeros(len(z_new_km))
        for i in [0,3,2,4,1]:
            accumulated_escape_LSA_new += accumulated_escape_LSA[i]
            axb7.fill_betweenx(z_new_km,accumulated_escape_LSA_old,accumulated_escape_LSA_new, facecolor=accumulated_colors_LSA[i], alpha=1)
            accumulated_escape_LSA_old += accumulated_escape_LSA[i]
        axb7.set_xlim([0,7])
        axb7.set_ylim([110,400])

        accumulated_escape_HSA_old = np.zeros(len(z_new_km))
        accumulated_escape_HSA_new = np.zeros(len(z_new_km))
        for i in [0,3,2,4,1]:
            accumulated_escape_HSA_new += accumulated_escape_HSA[i]
            axb11.fill_betweenx(z_new_km,accumulated_escape_HSA_old,accumulated_escape_HSA_new, facecolor=accumulated_colors_HSA[i], alpha=1)
            accumulated_escape_HSA_old += accumulated_escape_HSA[i]
        axb11.set_xlim([0,7])
        axb11.set_ylim([110,400])

        # Add cumulative production rate profiles
        accumulated_prod_LSA_new = np.zeros(len(z_new_km)) # These four lines are for total production. Vector production needs to be outside top 5 mechanism selection if loop above
        for i in range(0,len(mechanism_labels)):
           accumulated_prod_LSA_new += accumulated_prod_LSA[i]
        axb7.plot(accumulated_prod_LSA_new, z_new_km, color='k', alpha=1)

        accumulated_prod_HSA_new = np.zeros(len(z_new_km))
        for i in range(0,len(mechanism_labels)):
            accumulated_prod_HSA_new += accumulated_prod_HSA[i]
        axb11.plot(accumulated_prod_HSA_new, z_new_km, color='k', alpha=1)

        # Pie charts
        colours_LSA = ['#c51160','#353a47','oldlace','#84b082','#885a5a','#f7c1bb','lightgrey']
        colours_HSA = ['#c51160','oldlace','#353a47','#84b082','#885a5a','#4e8069','lightgrey']
        axb12.pie(
            fractions_LSA,
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
            startangle=90,
            radius=1000,
            colors=colours_LSA,
            normalize=True,
            frame=False)

        axb13.pie(
            fractions_HSA,
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
            startangle=90,
            colors=colours_HSA,
            radius=938.36,
            normalize=True,
            frame=False)
       

        axb12.set_ylim([-1275,925])
        axb12.set_xlim([-1100,1100])
        axb13.set_ylim([-1000.38,1199.62])
        axb13.set_xlim([-1024.62,1175.38])

        axb7.text(0.9,260,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2],ha='left',va='center',fontweight='bold')
        axb7.plot([0.8,0.06],[260,260],color='k',linewidth=0.3)
        axb7.text(5,125,'CO$_2^+$ + H$_2$',fontsize=14,color=colors[3],ha='center',fontweight='bold')
        axb7.plot([5,5],[140,155],color='k',linewidth=0.3)
        axb7.text(0.5,120,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0],ha='left',fontweight='bold')
        axb7.text(1.9,210,'OCOH$^+$ + e$^-$',fontsize=14,color=colors[1],ha='left',va='center',fontweight='bold')
        axb7.plot([1.8,0.57],[210,210],color='k',linewidth=0.3)
        axb7.text(1.45,235,'OH$^+$ + O',fontsize=14,color=colors[4],ha='left',va='center',fontweight='bold')
        axb7.plot([1.35,0.27],[235,235],color='k',linewidth=0.3)
        axb7.text(5.7,185,'escaping & bound\ntotal production\n(R1–47)',color='k',fontsize=14,ha='center') #,fontweight='bold')
        axb7.plot([5.7,5.7],[183,170],color='k',linewidth=0.3)

        axb11.text(1.3,310,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2],ha='left',va='center',fontweight='bold')
        axb11.plot([1.2,0.1],[310,310],color='k',linewidth=0.3)
        axb11.text(2.5,165,'CO$_2^+$ + H$_2$',fontsize=14,color=colors[3],ha='left',fontweight='bold')
        axb11.plot([2.5,2],[172,172],color='k',linewidth=0.3)
        axb11.text(0.2,120,'HCO$^+$ + e$^-$',fontsize=14,color=colors[0],ha='left',fontweight='bold')
        axb11.plot([0.5,0.5],[135,155],color='k',linewidth=0.3)
        axb11.text(2.6,250,'N$_2^+$ + H$_2$',fontsize=14,color=colors[5],ha='left',va='center',fontweight='bold')
        axb11.plot([2.5,1],[250,250],color='k',linewidth=0.3)
        axb11.text(2.1,280,'OH$^+$ + O',fontsize=14,color=colors[4],ha='left',va='center',fontweight='bold')
        axb11.plot([2,0.5],[280,280],color='k',linewidth=0.3)
        axb11.text(5.7,185,'escaping & bound\ntotal production\n(R1–47)',color='k',fontsize=14,ha='center') #fontweight='bold')
        axb11.plot([5.7,5.7],[182,170],color='k',linewidth=0.3)

        axb12.text(664,785,'O$^+$($^4$S) + H$_2$',fontsize=14,color=colors[2],ha='left',fontweight='bold')
        axb12.text(242,-621,'CO$_2^+$\n + H$_2$',fontsize=14,color='w',ha='center',fontweight='bold')
        axb12.text(-838,36,'HCO$^+$ + e$^-$',fontsize=14,color='w',ha='left',fontweight='bold')
        
        axb12.text(260,1060,'OCOH$^+$ + e$^-$',fontsize=14,color=colors[1],ha='left',fontweight='bold')
        axb12.text(400,929,'OH$^+$ + O',fontsize=14,color=colors[4],ha='left',fontweight='bold')
        axb12.text(-130,1025,'Other',fontsize=14,color='grey',ha='left',fontweight='bold')
        axb12.text(594,160,r'O$_{\rm hot}$ + H',fontsize=14,color='#353a47',ha='center',fontweight='bold')
        axb12.text(594,-50,'Shematovich\n(2013)',fontsize=10,color='#353a47',ha='center',fontweight='bold')

        axb12.text(0,-1300,'Total escape (R1–47 & O$_{\mathrm{hot}}$ + H):\n3.42 × 10$^7$ cm$^{-2}$ s$^{-1}$',ha='center',va='bottom',fontsize=14,fontweight='bold')  #EDIT
        axb13.text(75.38,-1250.62,'Total escape (R1–47 & O$_{\mathrm{hot}}$ + H):\n3.01 × 10$^7$ cm$^{-2}$ s$^{-1}$',ha='center',va='bottom'
                   ,fontsize=14,fontweight='bold') # EDIT
        

        axb13.text(174,-65,'O$^+$($^4$S) + H$_2$',fontsize=14,color='w',fontweight='bold')
        axb13.text(250,-572,'CO$_2^+$\n+ H$_2$',fontsize=14,color='w',fontweight='bold',ha='center')
        axb13.text(-790,161,'HCO$^+$ + e$^-$',fontsize=14,color='w',fontweight='bold')
        axb13.text(500,849,'N$_2^+$ + H$_2$',fontsize=14,color=colors[5],fontweight='bold')
        axb13.text(749,669,'OH$^+$ + O',fontsize=14,color=colors[4],fontweight='bold')
        axb13.text(-50,1000,'Other',fontsize=14,color='grey',ha='left',fontweight='bold')
        axb13.text(-420,-360,r'O$_{\rm hot}$ + H',fontsize=14,color='#353a47',ha='center',fontweight='bold')
        axb13.text(-420,-560,'Shematovich\n(2013)',fontsize=10,color='#353a47',ha='center',fontweight='bold')

        for ax in [axb7,axb11]:
            ax.set_ylabel('Altitude km)',fontsize=16)

        axb11.set_xlabel('Escaping H production rate (cm$^{-3}$ s$^{-1}$)',fontsize=16)

        for ax in [axb7,axb11]:
            ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
            ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
            for tick in ax.yaxis.get_majorticklabels():
                tick.set_fontsize(14)
        for ax in [axb7]:
            for tick in ax.xaxis.get_majorticklabels():
                tick.set_color('none')
        for ax in [axb11]:
            for tick in ax.xaxis.get_majorticklabels():
                tick.set_fontsize(14)

        axb7.text(-1.8,250,'Low\nsolar\nactivity',fontweight='bold',fontsize=18,ha='center')
        axb11.text(-1.8,250,'High\nsolar\nactivity',fontweight='bold',fontsize=18,ha='center')

        # Figure letters
        axb7.text(6.3,375,'(a)',fontweight='bold',fontsize=16)
        axb11.text(6.3,375,'(c)',fontweight='bold',fontsize=16)
        axb7.text(7.2,375,'(b)',fontweight='bold',fontsize=16)
        axb11.text(7.2,375,'(d)',fontweight='bold',fontsize=16)
        
        fig13.subplots_adjust(hspace=0.02,wspace=0.02,top=0.95,bottom=0.1,left=0.11,right=0.95)
        l, b, w, h = axb7.get_position().bounds
        axb7.set_position([l+0.13*w,b,w,h])
        l, b, w, h = axb11.get_position().bounds
        axb11.set_position([l+0.13*w,b,w,h])
        
        fig13.savefig('../model_output_data/figures/4_Accumulation.pdf',dpi=1200)
        plt.show()

    print('fig_choice',fig_choice)
    if fig_choice == 1:
        plt.show()

    
    if fig_choice == 8:
        fig8.subplots_adjust(hspace=0.01,wspace=0.01,top=0.97,bottom=0.11,left=0.06,right=0.99)
        fig8.savefig('../model_output_data/figures/1_Method_figure.pdf',dpi=1200)
        plt.show()
    elif fig_choice == 12:
        fig2.subplots_adjust(hspace=0.01,wspace=0.01,top=0.93,bottom=0.1,left=0.08,right=0.99)
        l, b, w, h = axb6.get_position().bounds
        axb6.set_position([l+0.15*w,b,w*0.7,h])
        l, b, w, h = axb10.get_position().bounds
        axb10.set_position([l+0.15*w,b,w*0.7,h])
        for ax in [axb4,axb5,axb8,axb9]:
            l, b, w, h = ax.get_position().bounds
            ax.set_position([l+0.15*w,b,w,h])
        for ax in [axb7,axb11]:
            l, b, w, h = ax.get_position().bounds
            ax.set_position([l-0.15*w,b,w,h])

        fig2.savefig('../model_output_data/figures/3_Results_profiles_2.pdf',dpi=1200,transparent=True) # consider using pdf instead??
        plt.show()
        
    elif fig_choice == 14:
        fig14.subplots_adjust(hspace=0.02,wspace=0.01,top=0.96,bottom=0.12,left=0.06,right=0.99)
        l, b, w, h = axb6.get_position().bounds
        axb6.set_position([l,b,w*0.8,h])
        l, b, w, h = axb12.get_position().bounds
        axb12.set_position([l-w*0.2,b,w*0.8,h])
        l, b, w, h = axb7.get_position().bounds
        axb7.set_position([l-w*0.4,b,w,h])
        
        fig14.savefig('../model_output_data/figures/S4_CX.pdf',dpi=1200)
        plt.show()

def find_rate_coefficient(mechanism,T_type,Tn_new,Te_new,Ti_new,z):
    import math
    from numpy import sqrt

    if T_type == 'Tn':
        T_new = Tn_new
    elif T_type == 'Te':
        T_new = Te_new
    elif T_type == 'Ti':
        T_new = Ti_new
    
    if mechanism == 'HCO+ + e- \u2192 CO + H':
        if T_new[z] <= 300:
            rate_coefficient = 2e-7*(T_new[z]/300)**-1.25
        else:
            rate_coefficient = 2e-7*(300/T_new[z])

    elif mechanism == 'OCOH+ DR \u2192 CO2 + H':
        rate_coefficient = 1.75e-8*(T_new[z]/300)**-0.5

    elif mechanism == 'OCOH+ DR \u2192 O + CO + H':
        rate_coefficient = 2.38e-7*(T_new[z]/300)**-0.5    

    elif mechanism == 'OH+ DR \u2192 O(1D) + H':
        rate_coefficient = 3.94e-8*(T_new[z]/300)**-1.28

    elif mechanism == 'H3+ DR \u2192 H2 + H':
        rate_coefficient = 1.7e-8*(T_new[z]/300)**-0.52

    elif mechanism == 'H3+ DR \u2192 H + H + H':
        rate_coefficient = 5.1e-8*(T_new[z]/300)**-0.52

    elif mechanism == 'HNO+ DR \u2192 NO + H':
        rate_coefficient = 3e-7*(T_new[z]/300)**-0.5

    elif mechanism == 'N2H+ DR \u2192 N2 + H':
        rate_coefficient = 2.325e-7*(T_new[z]/300)**-0.84

    elif mechanism == 'H2+ DR \u2192 H + H':
        rate_coefficient = 1.75e-8*(T_new[z]/300)**-0.4

    elif mechanism == 'CH+ DR \u2192 C(1D) + H':
        rate_coefficient = 6.16e-8*(T_new[z]/300)**-0.4

    elif mechanism == 'CH+ DR \u2192 C(1S) + H':
        rate_coefficient = 1.64e-8*(T_new[z]/300)**-0.4

    elif mechanism == 'HO2+ DR \u2192 O2 + H':
        rate_coefficient = 6.0e-8*(T_new[z]/300)**-0.5

    elif mechanism == 'HO2+ DR \u2192 O + O + H':
        rate_coefficient = 6.0e-8*(T_new[z]/300)**-0.5

    elif mechanism == 'ArH+ DR \u2192 Ar + H':
        rate_coefficient = 1.0e-9*(T_new[z]/300)**-0.5


        
    elif mechanism == 'H + H+ \u2192 H+ + H':
        rate_coefficient = 8.696e-11*sqrt(Tn_new[z] + Ti_new[z]) # from Rodriguez+ 1984
        
    elif mechanism == 'H+ + O \u2192 O+ + H':
        rate_coefficient = 6.86e-10 * ((Ti_new[z]/300)**0.26) * math.exp(-224.3/Ti_new[z])  # from UMIST database apparently  

    elif mechanism == 'CO2 + H+ \u2192 CO2+ + H':
        rate_coefficient = 3e-9


        
    elif mechanism == 'O+(4S) + H2 \u2192 OH+ + H':
        rate_coefficient = 1.65e-9

    elif mechanism == 'O+(2D) + H2 \u2192 OH+ + H':
        rate_coefficient = 1.5e-9

    elif mechanism == 'O+(2P) + H2 \u2192 OH+ + H':
        rate_coefficient = 1.0e-9

    elif mechanism == 'O+(2P) + H2 \u2192 O + H+ + H':
        rate_coefficient = 2.16e-10*(T_new[z]/300)**-0.97*math.exp(-292/T_new[z])

    elif mechanism == 'CO2+ + H2 \u2192 OCOH+ + H':
        rate_coefficient = 9.5e-10*(T_new[z]/300)**-0.15

    elif mechanism == 'He+ + H2 \u2192 H+ + He + H':
        rate_coefficient = 8.3e-14

    elif mechanism == 'H2+ + H2 \u2192 H3+ + H':
        rate_coefficient = 2.24e-9*(T_new[z]/300)**-0.042*math.exp(T_new[z]/-46600)

    elif mechanism == 'H2+ + CO2 \u2192 OCOH+ + H':
        rate_coefficient = 2.35e-9

    elif mechanism == 'CO+ + H2 \u2192 HCO+ + H' or mechanism == 'CO+ + H2 \u2192 HOC+ + H':
        rate_coefficient = 7.5e-10

    elif mechanism == 'N+ + H2 \u2192 NH+ + H':
        rate_coefficient = 8.23e-10*math.exp(-209/T_new[z])

    elif mechanism == 'N2+ + H2 \u2192 N2H+ + H':
        rate_coefficient = 1.52e-9

    elif mechanism == 'OH+ + H2 \u2192 H2O+ + H':
        rate_coefficient = 9.7e-10

    elif mechanism == 'OH+ + N \u2192 NO+ + H':
        rate_coefficient = 8.9e-10

    elif mechanism == 'OH+ + O \u2192 O2+ + H':
        rate_coefficient = 7.2e-10

    elif mechanism == 'HNO+ + O \u2192 NO2+ + H':
        rate_coefficient = 1.0e-12

    elif mechanism == 'H2+ + N2 \u2192 N2H+ + H':
        rate_coefficient = 2.3e-9

    elif mechanism == 'H2+ + O \u2192 OH+ + H':
        rate_coefficient = 1.5e-9

    elif mechanism == 'H2+ + CO \u2192 HCO+ + H' or mechanism == 'H2+ + CO \u2192 HOC+ + H':
        rate_coefficient = 7.65e-10

    elif mechanism == 'H2+ + O2 \u2192 HO2+ + H':
        rate_coefficient = 1.53e-9

    elif mechanism == 'H2+ + N \u2192 NH+ + H':
        rate_coefficient = 1.9e-9

    elif mechanism == 'H3+ + O \u2192 H2O+ + H':
        rate_coefficient = 3.42e-10*(T_new[z]/300)**-0.156 * math.exp(-1.41/T_new[z])

    elif mechanism == 'HO2+ + N \u2192 NO2+ + H':
        rate_coefficient = 1e-12

    elif mechanism == 'CH+ + H2 \u2192 CH2+ + H':
        rate_coefficient = 1.2e-9

    elif mechanism == 'CH+ + N \u2192 CN+ + H':
        rate_coefficient = 1.9e-10

    elif mechanism == 'CH+ + O \u2192 CO+ + H':
        rate_coefficient = 3.5e-10

    elif mechanism == 'C+ + H2 \u2192 CH+ + H':
        rate_coefficient = 7.4e-10*math.exp(-4538/T_new[z])

    elif mechanism == 'H2+ + C \u2192 CH+ + H':
        rate_coefficient = 2.4e-9

    elif mechanism == 'CH+ + C \u2192 C2+ + H':
        rate_coefficient = 1.2e-9

    elif mechanism == 'Ar+ + H2 \u2192 ArH+ + H':
        rate_coefficient = 8.72e-10

    elif mechanism == 'H2+ + Ar \u2192 ArH+ + H':
        rate_coefficient = 2.3e-9

    elif mechanism == 'H3O+ + e \u2192 OH + H + H':
        rate_coefficient = 3.05e-7*(300/T_new[z])**0.5

    return rate_coefficient

def peak_prod_altitude(production_rate,z_new):
    print(z_new[production_rate.index(max(production_rate))])
    print(max(production_rate))
    return (z_new[production_rate.index(max(production_rate))], max(production_rate))

def percent_escaping_alt(weighted_escape,escape_probability,production_rate,n,z_new,bin_size):
    # Returns z_new index below which fraction n of escaping hot H is produced
    total = 0.0
    for i,item in enumerate(production_rate):
        total += item*escape_probability[i]*bin_size
        if total/weighted_escape > n:
            print('percent_escaping_alt 2',z_new[i])
            return(z_new[i])
            break


def extrapolation(species,alt,n_orig,z_orig,end):
    # Function to extrapolate density profiles down to 80 km, if the input profiles don't reach this low

    import numpy as np

    if end == 'lower':
        # For a lot of species, extrapolate using by assuming the z-log(n) gradient between the lowest and third lowest points in the dataset describes the z-log(n) gradient below the lowest available values.
        if species == 'ArH+' or species == 'CH+' or species == 'CO+' or species == 'H2+' or species == 'H3+' or species == 'He+' or species == 'HNO+' or species == 'HO2+' or species == 'N+' or species == 'N2+' or species == 'N2H+' or species == 'OCOH+' or species == 'OH+' or species == 'Ar+' or species == 'C+' or species == 'C' or species == 'N' or species == 'O+' or species == 'O+2D' or species == 'O+2P' or species == 'O+4S':

            dndz = (z_orig[2]-z_orig[0])/(np.log10(n_orig[2])-np.log10(n_orig[0]))
            number_density = 10** (np.log10(n_orig[0]) + (alt - z_orig[0])*(np.log10(n_orig[2]) - np.log10(n_orig[0]))/(z_orig[2] - z_orig[0]))
        else:
            number_density = n_orig[0]

    elif end == 'upper':
       if species == 'C' or species == 'N' or species == 'Ar' or species == 'O2' or species == 'O+2D' or species == 'O+2P':
            dndz = (z_orig[len(z_orig)-1]-z_orig[len(z_orig)-3])/(np.log10(n_orig[len(n_orig)-1])-np.log10(n_orig[len(n_orig)-3]))
            number_density = 10** (np.log10(n_orig[len(z_orig)-1]) + (alt - z_orig[len(z_orig)-1])*(np.log10(n_orig[len(n_orig)-3]) - np.log10(n_orig[len(n_orig)-1]))/(z_orig[len(z_orig)-3] - z_orig[len(z_orig)-1]))
       else:
            number_density = n_orig[len(n_orig)-1]
   
    return number_density

calculate_weighted_escape_estimate(['LSA','HSA'])


                       
                        
