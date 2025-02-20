# Script to plot Figure 1a-c for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

# This script plots density profiles of (a) background species and (b) reactants, as well as the production rate of H via HCO+ dissociative recombination (DR)
# The script also prints out escape fluxes of hot H via HCO+ DR when the particles initially have kinetic energies of 5 eV or 0.2 eV

def plot_Figure1abc(SA_list):

    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import pandas as pd
    import bisect
    import pdb
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                   AutoMinorLocator, LogLocator, LinearLocator)
    import matplotlib.ticker as ticker

    from Common_Plotting_Functions import find_bg_dens,find_col_dens_above
    from escape_prob_curves import escape_prob_curves

    planet = 'Mars'

    ### What you need to input: ###
    # 1) R1 data file
    # 2) R2 data file
    # 3) T data file, T type ('mechanism_Ts'; neutral, ion, electron)
    # 4) Reaction rate coefficient equation

    mechanism_labels = ['HCO+ + e- \u2192 CO + H']                       
    mechanism_R1s = ['HCO+']                   
    mechanism_R2s = ['e']
    mechanism_Ts =  ['Te'] # options are 'Tn' or 'Te' or 'Ti' (or 'T_none' if not required for rate)

    fig1 = plt.figure(figsize=(13,5))
    axa = fig1.add_subplot(1,3,1)
    axb = fig1.add_subplot(1,3,2)
    axc = fig1.add_subplot(1,3,3)

    # Import fitted escape probability values (escape probabilities and z) from escape_prob_curves.py
    returned_escape_rates,returned_escape_z = escape_prob_curves('Return')
    returned_escape_rates_LSA_5eV, returned_escape_rates_LSA_0_2eV, returned_escape_rates_HSA_5eV, returned_escape_rates_HSA_0_2eV  = returned_escape_rates
    returned_escape_z_LSA_5eV,returned_escape_z_LSA_0_2eV,returned_escape_z_HSA_5eV,returned_escape_z_HSA_0_2eV  = returned_escape_z

    # Create new z grid
    z_new = range(8000000,40025000,25000) # new z grid 80-400 km in cm, with 0.25 km spacing
    bin_size = 25000 # 0.25 km

    # Loop over low and high solar activity
    for s,SA in enumerate(SA_list):

        # Import T profile data
        T_file = '../inputs/Mars/MarsTemp%s_Fox2015.csv' % (SA)
        z_orig_T = pd.read_csv(T_file)['#alt(cm)']
        Tn_orig_T = pd.read_csv(T_file)['Tn(K)']
        Te_orig_T = pd.read_csv(T_file)['Te(K)']
        Ti_orig_T = pd.read_csv(T_file)['Ti(K)']
    
        # Check that the input profiles increase in altitude (as this script will only work that way)
        if z_orig_T[0] > z_orig_T[len(z_orig_T)-1]:
            print('Please use T profiles with increasing altitude')
            stop

        # Interpolate T for new grid - neutral (Tn), ion (Ti) and electron (Te)
        Tn_new = []
        Ti_new = []
        Te_new = []
        # - if the new altitude is below the lowest altitude in the dataset, set T at the new altitude to the lowest T in the dataset                                                                                   
        # - if the new altitude is above the highest altitude in the dataset, set T at the new altitude to the highest T in the dataset                                                                                 
        # - otherwise, do a linear interpolation using the gradient between the two existing points either side of the new point and the difference in z between the new z and the lowest adjacent existing z
        for z,alt in enumerate(z_new):
            if alt <= z_orig_T[0]: # Extrapolate below lowest data point using lowest-altitude data point
                Tn_new.append(Tn_orig_T[0]) 
                Te_new.append(Te_orig_T[0]) 
                Ti_new.append(Ti_orig_T[0]) 
            elif alt >= z_orig_T[len(z_orig_T)-1]: # Extrapolate above highest data point using highest-altitude data point (doesn't really happen)
                Tn_new.append(Tn_orig_T[len(z_orig_T)-1])
                Te_new.append(Te_orig_T[len(z_orig_T)-1])
                Ti_new.append(Ti_orig_T[len(z_orig_T)-1])
            else: # Linear interpolation
                i = bisect.bisect_right(z_orig_T,alt)
                dndz = (Tn_orig_T[i] - Tn_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Tn_new.append(Tn_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))
                dndz = (Te_orig_T[i] - Te_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Te_new.append(Te_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))
                dndz = (Ti_orig_T[i] - Ti_orig_T[i-1])/(z_orig_T[i] - z_orig_T[i-1])
                Ti_new.append(Ti_orig_T[i-1] + dndz*(alt-z_orig_T[i-1]))

    
        # Import reactant profiles
        for m, mech in enumerate(mechanism_labels):
            T_type = mechanism_Ts[m] # Te, Tn or Ti (on which rate coefficient is dependent)

            r1_data_file = '../inputs/Mars/%s_density_profile_%s_eroded_Fox2015.csv' % (mechanism_R1s[m],SA) # Fox (2015) eroded model
            r2_data_file = '../inputs/Mars/electron_density_profile_%s_eroded_Fox2015.csv' % (SA) # Fox (2015) eroded model
            
            n_orig_r1 = pd.read_csv(r1_data_file)['%s(cm-3)' % (mechanism_R1s[m])] # HCO+ densities
            n_orig_r2 = pd.read_csv(r2_data_file)['%s(cm-3)' % (mechanism_R2s[m])] # e- densities

            z_orig_r1 = pd.read_csv(r1_data_file)['#alt(cm)']
            z_orig_r2 = pd.read_csv(r2_data_file)['#alt(cm)']
                                                                                             
            # Check that the input profiles increase in altitude (as this script will only work that way)                
            for prof in [z_orig_r1,z_orig_r2]:
                if prof[0] > prof[len(prof)-1]:
                    print('Please use profiles with increasing altitude')
                    stop


            # Interpolate densities for new grid                                                                                       
            n_new_r1 = []
            n_new_r2 = []

            # Extrapolate/interpolate R1 (HCO+) densities
            for z,alt in enumerate(z_new):
                if alt <= z_orig_r1[0]:
                    extrap_number_dens = extrapolation(mechanism_R1s[m],alt,n_orig_r1,z_orig_r1,'lower') # extrapolate below lowest point using linear extrapolation in altitude-log(density) space, using gradient between lowest three points
                    n_new_r1.append(extrap_number_dens)
                elif alt >= z_orig_r1[len(z_orig_r1)-1]:
                    extrap_number_dens = extrapolation(mechanism_R1s[m],alt,n_orig_r1,z_orig_r1,'upper') # extrapolate above highest point using linear extrapolation in altitude-log(density) space, using gradient between highest three points
                    n_new_r1.append(extrap_number_dens)
                else:
                    # linear interpolation
                    i = bisect.bisect_right(z_orig_r1,alt)
                    dndz = (n_orig_r1[i] - n_orig_r1[i-1])/(z_orig_r1[i] - z_orig_r1[i-1])
                    n_new_r1.append(n_orig_r1[i-1] + dndz*(alt-z_orig_r1[i-1]))

            # Interpolate R2 (electron) densities for new grid                                                        
            for z,alt in enumerate(z_new):
                if alt <= z_orig_r2[0]:
                    extrap_number_dens = extrapolation(mechanism_R2s[m],alt,n_orig_r2,z_orig_r2,'lower') # extrapolate below lowest point using linear extrapolation in altitude-log(density) space, using gradient between lowest three points
                    n_new_r2.append(extrap_number_dens)
                elif alt >= z_orig_r2[len(z_orig_r2)-1]:
                    extrap_number_dens = extrapolation(mechanism_R2s[m],alt,n_orig_r2,z_orig_r2,'upper') # extrapolate above highest point using linear extrapolation in altitude-log(density) space, using gradient between highest three points
                    n_new_r2.append(extrap_number_dens)
                else:
                    # linear interpolation
                    i = bisect.bisect_right(z_orig_r2,alt)
                    dndz = (n_orig_r2[i] - n_orig_r2[i-1])/(z_orig_r2[i] - z_orig_r2[i-1])
                    n_new_r2.append(n_orig_r2[i-1] + dndz*(alt-z_orig_r2[i-1]))

            # Plot HCO+ and e- density profiles in panel (b)
            color_HCOpl = 'darkblue'
            color_e = '#42CAFD'
            if SA == 'LSA' and mechanism_R1s[m] == 'HCO+' and mechanism_R2s[m] == 'e':
                axb.plot(n_new_r1,np.divide(z_new,1e5),linewidth=3,label=mechanism_R1s[m],color=color_HCOpl) # HCO+
                axb.plot(n_new_r2,np.divide(z_new,1e5),linewidth=3,label=mechanism_R2s[m],color=color_e) # e-
            elif SA == 'HSA' and mechanism_R1s[m] == 'HCO+' and mechanism_R2s[m] == 'e':
                axb.plot(n_new_r1,np.divide(z_new,1e5),linewidth=0.5,label=mechanism_R1s[m],color=color_HCOpl) # HCO+
                axb.plot(n_new_r2,np.divide(z_new,1e5),linewidth=0.5,label=mechanism_R2s[m],color=color_e) # e-
                    
            production_rate = []
            escape_flux_5eV = 0.0
            escape_flux_0_2eV = 0.0
            # Calculate production rate and escape flux for new grid
            for z, alt in enumerate(z_new):
                rate_coefficient = find_rate_coefficient(mechanism_labels[m],T_type,Tn_new,Te_new,Ti_new,z)
                production_rate.append(rate_coefficient*n_new_r1[z]*n_new_r2[z])
                if SA == 'LSA':
                    escape_flux_5eV += (production_rate[z]*returned_escape_rates_LSA_5eV[z]*bin_size)
                    escape_flux_0_2eV += (production_rate[z]*returned_escape_rates_LSA_0_2eV[z]*bin_size)
                           
                elif SA == 'HSA':
                    escape_flux_5eV += (production_rate[z]*returned_escape_rates_HSA_5eV[z]*bin_size)
                    escape_flux_0_2eV += (production_rate[z]*returned_escape_rates_HSA_0_2eV[z]*bin_size)

            # The line below prints out the escape fluxes for HCO+ DR when all particles initially have 5 eV or 0.2 eV (calculation results described in Section 5)
            print('HCO+ DR escape flux for %s is %s for 5 eV, %s for 0.2 eV' % (SA,escape_flux_5eV,escape_flux_0_2eV))
            # Plot production rate profiles in panel (c)
            if  mechanism_R1s[m] == 'HCO+' and mechanism_R2s[m] == 'e':
                if SA == 'LSA':
                    axc.plot(production_rate,np.divide(z_new,1e5),linewidth=3,label=m,color='k')
                elif SA == 'HSA':
                    axc.plot(production_rate,np.divide(z_new,1e5),linewidth=0.5,label=m,color='k')
                           
        # Plot background species densities in panel (a)
        [bg_z,CO2_density,O_density,CO_density,N2_density] = find_bg_dens(SA,planet)                                    # Find background density profiles
        col_dens_above = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,SA,planet) # Find profiles of column density above z with z

        colors_bg = ['tab:blue','tab:red','tab:orange','darkslateblue']
        if SA == 'LSA':
            lw_panel1 = 3
        elif SA == 'HSA':
            lw_panel1 = 0.5
        axa.plot(CO2_density,np.divide(bg_z,1e5),label='CO2',linewidth=lw_panel1,color=colors_bg[0])
        axa.plot(O_density,np.divide(bg_z,1e5),label='O',linewidth=lw_panel1,color=colors_bg[1])
        axa.plot(N2_density,np.divide(bg_z,1e5),label='N2',linewidth=lw_panel1,color=colors_bg[2])
        axa.plot(CO_density,np.divide(bg_z,1e5),label='CO',linewidth=lw_panel1,color=colors_bg[3])
        axa.plot(col_dens_above,np.divide(z_new,1e5),color='k',label='Col. dens. above',linewidth=lw_panel1)


    # Plot aesthetics
    fig1.subplots_adjust(hspace=0.02,wspace=0.03,top=0.93,bottom=0.12,left=0.06,right=0.98)
    l, b, w, h = axb.get_position().bounds
  #  axb.set_position([l,b,w*0.7,h])
    axb.set_position([l+w*0.2,b,w*0.8,h])
    l, b, w, h = axc.get_position().bounds
  #  axc.set_position([l-0.3*w,b,w*0.7,h])
    axc.set_position([l,b,w*0.8,h])
    l, b, w, h = axa.get_position().bounds
    axa.set_position([l,b,w*1.2,h])
    for ax in [axa, axb, axc]:
        ax.set_ylim([80,400])
        ax.set_xscale('log')
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_fontsize(14)
        ax.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
        ax.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
    axa.set_ylabel('Altitude (km)',fontsize=16)
    axa.set_xlabel('Number (cm$^{-3}$) or column (cm$^{-2}$) density',fontsize=16)
    axb.set_xlabel('Number density (cm$^{-3}$)',fontsize=16)
    axc.set_xlabel('Production rate (cm$^{-3}$ s$^{-1}$)',fontsize=16)
        
    axb.set_xlim([1e1,3e5])
    axc.set_xlim([1e-3,3e1])
    
    # Labels
    axa.text(2e18,145,'col. dens.\nabove',fontsize=16,color='k',ha='center') #,ha='right')
    axa.text(1e7,380,'O',fontsize=16,color='tab:red')
    axa.text(9e4,280,'N$_2$',fontsize=16,color='tab:orange')
    axa.text(1e6,150,'CO',fontsize=16,color='darkslateblue')
    axa.text(1e14,90,'CO$_2$',fontsize=16,color='tab:blue')
    axa.text(1e17,340,'LSA',fontsize=16,fontweight='bold',color='k')
    axa.text(1e17,320,'HSA',fontsize=16,color='k')
    axa.text(0.9,0.93,'(a)',fontsize=16,color='k',transform=axa.transAxes)
        
    axb.text(4e1,180,'HCO$^+$',fontsize=16,color='darkblue')
    axb.text(3e4,130,'e$^-$',fontsize=16,color='#42CAFD')
    axb.text(4e4,340,'LSA',fontsize=16,fontweight='bold',color='k')
    axb.text(4e4,320,'HSA',fontsize=16,color='k')
    axb.text(0.88,0.93,'(b)',fontsize=16,color='k',transform=axb.transAxes)   
        
    axc.text(2.2e-2,190,'LSA',fontsize=16,fontweight='bold',color='k')
    axc.text(2.2e-2,300,'HSA',fontsize=16,color='k')
    axc.text(0.88,0.93,'(c)',fontsize=16,color='k',transform=axc.transAxes)

    for ax in [axb,axc]:
        for tick in ax.yaxis.get_majorticklabels():
            tick.set_color('none')
    for tick in axa.yaxis.get_majorticklabels():
        tick.set_fontsize(14)
    fig1.savefig('Figure1.pdf',dpi=1200)




def find_rate_coefficient(mechanism,T_type,Tn_new,Te_new,Ti_new,z):
    # Find rate coefficient for selected mechanism at altitude z
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

    return (rate_coefficient)

def extrapolation(species,alt,n_orig,z_orig,end):
    # Function to extrapolate density profiles down to 80 km (or above 400 km), if the input profiles don't reach this low (high)

    import numpy as np

    if end == 'lower':
        # For a lot of species, extrapolate using by assuming the z-log(n) gradient between the lowest and third lowest points in the dataset describes the z-log(n) gradient below the lowest available values.
        if species == 'ArH+' or species == 'CH+' or species == 'CO+' or species == 'H2+' or species == 'H3+' or species == 'He+' or species == 'HNO+' or species == 'HO2+' or species == 'N+' or species == 'N2+' or species == 'N2H+' or species == 'OCOH+' or species == 'OH+' or species == 'Ar+' or species == 'C+' or species == 'C' or species == 'N' or species == 'O+':
            dndz = (z_orig[2]-z_orig[0])/(np.log10(n_orig[2])-np.log10(n_orig[0]))
            number_density = 10** (np.log10(n_orig[0]) + (alt - z_orig[0])*(np.log10(n_orig[2]) - np.log10(n_orig[0]))/(z_orig[2] - z_orig[0]))
        else:
            number_density = n_orig[0]

    elif end == 'upper':
       if species == 'C' or species == 'N' or species == 'Ar' or species == 'O2':
            dndz = (z_orig[len(z_orig)-1]-z_orig[len(z_orig)-3])/(np.log10(n_orig[len(n_orig)-1])-np.log10(n_orig[len(n_orig)-3]))
            number_density = 10** (np.log10(n_orig[len(z_orig)-1]) + (alt - z_orig[len(z_orig)-1])*(np.log10(n_orig[len(n_orig)-3]) - np.log10(n_orig[len(n_orig)-1]))/(z_orig[len(z_orig)-3] - z_orig[len(z_orig)-1]))
       else:
            number_density = n_orig[len(n_orig)-1]
   
    return number_density


plot_Figure1abc(['LSA','HSA'])


                       
                        
