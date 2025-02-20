# Function to plot pie chart with relative importance of HCO+ DR, other top 5 mechanisms, and other mechanisms for total nonthermal.

#  Note hard-coded numbers

def weighted_esc_pie():
    import matplotlib.pyplot as plt
    import numpy as np

    fractions_to_plot = []

    # Taking into account right energy for mechanism
    for SA in ['LSA','HSA']:
        if SA == 'LSA': # LSA proper energies
            mechanisms = ['HCO$^+$ + e$^-$','CO2$+$ + H2','Momentum exchange (Shematovich, 2013)','O$^+$ + H2','OH$^+$ + O','OCOH$^+$ + e$^-$','Other']
            fluxes = [17079445.53,6712321.574,6e6,1585643.556,900464.1873,601039.86] # LSA proper energy
            total_nonthermal = 28156487.11 + 6e6 # without extra O+ + H2
            fluxes.append(total_nonthermal - sum(fluxes))
            colours = ['#c51160','#353a47','oldlace','#84b082','#885a5a','#f7c1bb','lightgrey']

        elif SA == 'HSA': # HSA proper energies
            mechanisms = ['HCO$^+$ + e$^-$',
                          'Momentum exchange (Shematovich, 2013)',
                          'CO2$+$ + H2',
                          'O$^+$ + H2',
                          'OH$^+$ + O',
                          'N2+ + H2','Other']
            fluxes = [9045711.694,6e6,5138867.68,3787462.966,2964338.48,690963.3437] # HSA proper energies
            total_nonthermal = 24045585.48 + 6e6 # without extra O+ + H2
            fluxes.append(total_nonthermal - sum(fluxes))
            colours = ['#c51160','oldlace','#353a47','#84b082','#885a5a','#4e8069','lightgrey']
            
        fractions = []
        for flux in fluxes:
            fractions.append(np.divide(flux,total_nonthermal)) # Other
        print(SA,fractions)
        print(np.sum(fractions))

        patches = plt.pie(
            fractions,
        #labels=mechanisms,
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
    #    textprops={'size': 'x-large'},
            startangle=90,
            colors=colours,
            normalize=True)
     #   plt.savefig('../model_output_data/figures/weighted_esc_pie.pdf' % (SA),dpi=2000,transparent=True)
        fractions_to_plot.append(fractions)
    return(fractions_to_plot)
    
weighted_esc_pie()
