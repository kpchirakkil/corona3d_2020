# Script to produce fitted escape probability curves for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

# Called as "escape_prob_curves('Return')" by "plot_Figure1abc.py" to pass escape probability data to calculate estimated escape flux for HCO+ DR (see Section 5)
# Called as "escape_prob_curves('Plot')" by "plot_Figure3.py" to plot escape probability figure

def escape_prob_curves(return_or_plot):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    import os
    from Common_Plotting_Functions import find_bg_dens, find_col_dens_above

    # Set up figure
    include_data = 'yes' # plot output from Monte Carlo model as points
    dpi = 300
    fig3 = plt.figure(figsize=(6,5))
    ax1 = fig3.add_subplot(1,1,1)
    planet = 'Mars'

    # energies and total cross sections (Zhang et al., 2009)
    energies = {0.2: 5.52e-15,
                5:3.89e-15}

    # new altitude profile (cm)
    z_new = range(8000000,40025000,25000) # 80-400 km with 0.25 km alt bins
    bin_size = 25000

    escape_rates_to_return = []
    escape_z_to_return = []

    esc_prob_prob_5eV_comb = [] # for combined LSA/HSA escape probability data
    esc_prob_prob_0_2eV_comb = [] # for combined LSA/HSA escape probability data
    esc_prob_col_dens_5eV = [] # for combined LSA/HSA column density above data
    esc_prob_col_dens_0_2eV = [] # for combined LSA/HSA column density above data

    for SA in ['LSA','HSA']: # loop over solar activity conditions
        # Read in escape probability data from Monte Carlo model output
        esc_prob_z_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 5 eV %s' % (SA)] # altitude (km)
        esc_prob_prob_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 5 eV %s' % (SA)] # Mean escape probability
        esc_prob_prob_1_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 1' % (SA)] # Escape probability from run 1...
        esc_prob_prob_2_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 2' % (SA)] # ... run 2
        esc_prob_prob_3_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 3' % (SA)] # ... run 3
        esc_prob_prob_4_5eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['5 eV %s 4' % (SA)] # ... run 4
        esc_prob_z_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Altitude (km) 0.2 eV %s' % (SA)]
        esc_prob_prob_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['Mean probability 0.2 eV %s' % (SA)]
        esc_prob_prob_1_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 1' % (SA)]
        esc_prob_prob_2_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 2' % (SA)]
        esc_prob_prob_3_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 3' % (SA)]
        esc_prob_prob_4_0_2eV = pd.read_csv('../model_output_data/escape_probabilities_%s.csv' % (SA))['0.2 eV %s 4' % (SA)]
                                           
        bg_z, CO2_density, O_density, CO_density, N2_density = find_bg_dens(SA,planet) # call function to find background density data from file

        # Interpolate background densities to match new z grid
        if SA == 'LSA':
            col_dens_above_LSA = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,'LSA',planet)
            col_dens_above = col_dens_above_LSA
        elif SA == 'HSA':
            col_dens_above_HSA = find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,'HSA',planet)
            col_dens_above = col_dens_above_HSA

        # Produce vectors of: (1) column density above z; (2) escape probability at altitude z
        for energy in energies:
            if SA == 'LSA':
                if energy == 5:
                    for i,item in enumerate(esc_prob_prob_5eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_5eV[i]*1e5:
                                esc_prob_col_dens_5eV.append(col_dens_above[z]) # for combined data (LSA and HSA together)
                                esc_prob_prob_5eV_comb.append(item) # for combined data
                 
                elif energy == 0.2:
                    for i,item in enumerate(esc_prob_prob_0_2eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_0_2eV[i]*1e5:
                                esc_prob_col_dens_0_2eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_0_2eV_comb.append(item) # for combined data
                
            elif SA == 'HSA':             
                if energy == 5:
                    for i,item in enumerate(esc_prob_prob_5eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_5eV[i]*1e5:
                                esc_prob_col_dens_5eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_5eV_comb.append(item) # for combined data
               
                elif energy == 0.2:
                    for i,item in enumerate(esc_prob_prob_0_2eV):
                        for z,alt in enumerate(z_new):
                            if alt == esc_prob_z_0_2eV[i]*1e5:
                                esc_prob_col_dens_0_2eV.append(col_dens_above[z]) # for combined data
                                esc_prob_prob_0_2eV_comb.append(item) # for combined data


        # Plot original escape/z data
        colors = ['k','silver']
        if include_data == 'yes':
            if SA == 'LSA':
                ax1.scatter(esc_prob_prob_5eV,esc_prob_z_5eV,linewidth=2,marker='x',zorder=2,color=colors[0])
                ax1.scatter(esc_prob_prob_1_5eV,esc_prob_z_5eV,linewidth=1,marker='x',zorder=2,color=colors[0],s=5)
                ax1.scatter(esc_prob_prob_2_5eV,esc_prob_z_5eV,linewidth=1,marker='x',zorder=2,color=colors[0],s=5)
                ax1.scatter(esc_prob_prob_3_5eV,esc_prob_z_5eV,linewidth=1,marker='x',zorder=2,color=colors[0],s=5)
                ax1.scatter(esc_prob_prob_4_5eV,esc_prob_z_5eV,linewidth=1,marker='x',zorder=2,color=colors[0],s=5)
            
                ax1.scatter(esc_prob_prob_0_2eV,esc_prob_z_0_2eV,linewidth=2,marker='x',zorder=2,color=colors[1])
                ax1.scatter(esc_prob_prob_1_0_2eV,esc_prob_z_0_2eV,linewidth=1,marker='x',zorder=2,color=colors[1],s=5)
                ax1.scatter(esc_prob_prob_2_0_2eV,esc_prob_z_0_2eV,linewidth=1,marker='x',zorder=2,color=colors[1],s=5)
                ax1.scatter(esc_prob_prob_3_0_2eV,esc_prob_z_0_2eV,linewidth=1,marker='x',zorder=2,color=colors[1],s=5)
                ax1.scatter(esc_prob_prob_4_0_2eV,esc_prob_z_0_2eV,linewidth=1,marker='x',zorder=2,color=colors[1],s=5)
           
            elif SA == 'HSA':
                ax1.scatter(esc_prob_prob_5eV,esc_prob_z_5eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[0],color='none')
                ax1.scatter(esc_prob_prob_1_5eV,esc_prob_z_5eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[0],color='none',s=5)
                ax1.scatter(esc_prob_prob_2_5eV,esc_prob_z_5eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[0],s=5,color='none')
                ax1.scatter(esc_prob_prob_3_5eV,esc_prob_z_5eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[0],s=5,color='none')
                ax1.scatter(esc_prob_prob_4_5eV,esc_prob_z_5eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[0],s=5,color='none')
            
                ax1.scatter(esc_prob_prob_0_2eV,esc_prob_z_0_2eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[1],color='none')
                ax1.scatter(esc_prob_prob_1_0_2eV,esc_prob_z_0_2eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[1],s=5,color='none')
                ax1.scatter(esc_prob_prob_2_0_2eV,esc_prob_z_0_2eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[1],color='none',s=5)
                ax1.scatter(esc_prob_prob_3_0_2eV,esc_prob_z_0_2eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[1],color='none',s=5)
                ax1.scatter(esc_prob_prob_4_0_2eV,esc_prob_z_0_2eV,linewidth=0.5,marker='o',zorder=2,edgecolor=colors[1],color='none',s=5)
                
    ## Combined data 5 eV ##
    N_5eV = np.array(esc_prob_col_dens_5eV, dtype=float)
    popt_5eV, pcov_5eV = curve_fit(func1,N_5eV,esc_prob_prob_5eV_comb) # find constants A and b for fitted escape probability curve (p = A exp(- b N_5eV(z) sigma))
    print('Constants A and b for 5 eV:',popt_5eV)

    calculated_escape_5eV_LSA = []
    calculated_escape_5eV_HSA = []
    for z,alt in enumerate(z_new):
        calculated_escape_5eV_LSA.append(popt_5eV[0]*np.exp(popt_5eV[1]*-1*col_dens_above_LSA[z]*energies[5]))
        calculated_escape_5eV_HSA.append(popt_5eV[0]*np.exp(popt_5eV[1]*-1*col_dens_above_HSA[z]*energies[5]))

    # Plot fitted escape probability curves
    ax1.plot(calculated_escape_5eV_LSA,np.divide(z_new,1e5),color=colors[0],linewidth=2,label = '5 eV') # 5eV LSA
    for v,value in enumerate(calculated_escape_5eV_LSA): # plot 1/e escape level
        if value > (1/np.exp(1)) and calculated_escape_5eV_LSA[v-1] < (1/np.exp(1)):
            ax1.plot([0,1],[np.divide(z_new[v],1e5),np.divide(z_new[v],1e5)],color='k',linestyle='--',linewidth=1,zorder=-2) 
            ax1.text(0.99,np.divide(z_new[v],1e5)-9,'LSA 5 eV',color='k',fontsize=9,ha='right')
            
    ax1.plot(calculated_escape_5eV_HSA,np.divide(z_new,1e5),color='k',linewidth=0.5) # 5 eV HSA
    for v,value in enumerate(calculated_escape_5eV_HSA): # plot 1/e escape level
        if value > (1/np.exp(1)) and calculated_escape_5eV_HSA[v-1] < (1/np.exp(1)):
            ax1.plot([0,1],[np.divide(z_new[v],1e5),np.divide(z_new[v],1e5)],color='k',linestyle='--',linewidth=0.5,zorder=-2) 
            ax1.text(0.99,np.divide(z_new[v],1e5)-9,'HSA 5 eV',color='k',fontsize=9,ha='right')

    ## Combined data 0.2 eV ##
    N_0_2eV = np.array(esc_prob_col_dens_0_2eV, dtype=float)
    popt_0_2eV, pcov_0_2eV = curve_fit(func2,N_0_2eV,esc_prob_prob_0_2eV_comb) # uses 0s
    print('Constants A and b for 0.2 eV:',popt_0_2eV)

    calculated_escape_0_2eV_LSA = []
    calculated_escape_0_2eV_HSA = []
           
    for z,alt in enumerate(z_new):
        calculated_escape_0_2eV_LSA.append(popt_0_2eV[0]*np.exp(popt_0_2eV[1]*-1*col_dens_above_LSA[z]*energies[0.2]))
        calculated_escape_0_2eV_HSA.append(popt_0_2eV[0]*np.exp(popt_0_2eV[1]*-1*col_dens_above_HSA[z]*energies[0.2]))

    # Plot fitted escape probability curves
    ax1.plot(calculated_escape_0_2eV_LSA,np.divide(z_new,1e5),color=colors[1],linewidth=2,label='0.2 eV')
    for v,value in enumerate(calculated_escape_0_2eV_LSA): # plot 1/e escape level
        if value > (1/np.exp(1)) and calculated_escape_0_2eV_LSA[v-1] < (1/np.exp(1)):
            ax1.plot([0,1],[np.divide(z_new[v],1e5),np.divide(z_new[v],1e5)],color='silver',linestyle='--',linewidth=1,zorder=-2) # Plot 1/e escape level
            ax1.text(0.99,np.divide(z_new[v],1e5)+2,'LSA 0.2 eV',color='darkgrey',fontsize=9,ha='right')
    ax1.plot(calculated_escape_0_2eV_HSA,np.divide(z_new,1e5),color='silver',linewidth=0.5) # 0.2 eV
    for v,value in enumerate(calculated_escape_0_2eV_HSA): # plot 1/e escape level
        if value > (1/np.exp(1)) and calculated_escape_0_2eV_HSA[v-1] < (1/np.exp(1)):
            ax1.plot([0,1],[np.divide(z_new[v],1e5),np.divide(z_new[v],1e5)],color='silver',linestyle='--',linewidth=0.5,zorder=-2) # Plot 1/e escape level
            ax1.text(0.99,np.divide(z_new[v],1e5)+2,'HSA\n 0.2 eV',color='darkgrey',fontsize=9,ha='right')

    # Lists of escape probabilities and altitudes (from fitted curves), to return to calling script, if 'Return' is input
    for escape_rates in [calculated_escape_5eV_LSA,calculated_escape_5eV_HSA,calculated_escape_0_2eV_LSA,calculated_escape_0_2eV_HSA]:
        escape_rates_to_return.append(escape_rates)
    for escape_z in [np.divide(z_new,1e5),np.divide(z_new,1e5),np.divide(z_new,1e5),np.divide(z_new,1e5)]:
        escape_z_to_return.append(escape_z)

    # Figure 3 aesthetics
    ax1.set_xlim([0,1])
    ax1.set_ylim([80,400])
    ax1.set_xlabel('Escape probability',fontsize=16)
    ax1.set_ylabel('Altitude (km)',fontsize=16)
    ax1.tick_params('both',width=1,length=8,direction='in',top=True,right=True)
    ax1.tick_params('both',which='minor',width=1,length=3,direction='in',top=True,right=True)
    
    fig3.subplots_adjust(hspace=0.02,wspace=0.01,top=0.93,bottom=0.08,left=0.25,right=0.96)

    # Labels
    ax1.text(0.7,110,'5 eV',fontsize=16)
    ax1.text(0.4,240,'0.2 eV',fontsize=16,color='silver')
    ax1.text(0.05,360,'Low solar activity',fontweight='bold',fontsize=16,color='dimgrey')
    ax1.text(0.05,340,'High solar activity',fontsize=16,color='dimgrey')

    
    if return_or_plot == 'Plot':
        plt.gcf().subplots_adjust(hspace=0.3,wspace=0.2,top=0.94,bottom=0.1,left=0.12,right=0.98)
        fig3.savefig('Figure3.pdf',dpi=1200)
    else:
        return (escape_rates_to_return,escape_z_to_return)

# Forms of the escape probbility curves: p(z,E) = A exp (-b N(z) sigma(E))
# for energy = 5 eV
def func1(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*3.89e-15)) # total cross section from Zhang et al. (2009)

# for energy = 0.2 eV
def func2(N,A,b):
    import numpy as np
    return(A*np.exp(-1*b*N*5.52e-15)) # total cross section from Zhang et al. (2009)
            
            
