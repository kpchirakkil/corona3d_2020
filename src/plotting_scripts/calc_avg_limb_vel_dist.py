# Script to plot Monte Carlo model output hot H velocity distributions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def calc_avg_limb_vel_dist(ax,altitudes):
    
    output_directories = ['../output/150km/vtally_1','../output/150km/vtally_3','../output/150km/vtally_8','../output/150km/vtally_9','../output/150km/vtally_4','../output/150km/vtally_5', '../output/150km/vtally_6','../output/150km/vtally_7','../output/150km/vtally_11','../output/150km/vtally_12','../output/150km/vtally_13','../output/150km/vtally_14','../output/150km/vtally_15','../output/150km/vtally_16']

    ax_twiny = ax.twiny()
    ax_twinx = ax.twinx()   

    colors = ['#c51160','hotpink','darkred','crimson','#440154FF','#2D708EFF','#3CBB75FF']
    
    for a,alt in enumerate(altitudes):
        cumulative_col_dens = []
        
        for d, output_directory in enumerate(output_directories):
            if d != 0:
                min_velocity_in_bin_prev = min_velocity_in_bin
                vel_bin_width_prev = vel_bin_width
            
            min_velocity_in_bin = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['min_velocity_in_bin_[km/s]']       
            vel_bin_width = min_velocity_in_bin[1] - min_velocity_in_bin[0]

            if d != 0:
                if min_velocity_in_bin[0] != min_velocity_in_bin_prev[0]:
                    print('WARNING: Mismatch in min_velocity_in_bin values')
                    stop

                if vel_bin_width != vel_bin_width_prev:
                    print('WARNING: Mismatch in vel_bin_width values')
                    stop
        

            col_dens = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['column_density_[cm-2]']

            if d == 0:
                cumulative_col_dens = col_dens
            else:
                cumulative_col_dens = cumulative_col_dens + col_dens        

        for d, output_directory in enumerate(output_directories):
            min_velocity_in_bin = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['min_velocity_in_bin_[km/s]']       
            vel_bin_width = min_velocity_in_bin[1] - min_velocity_in_bin[0]
            col_dens = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['column_density_[cm-2]']

            if d == 0:
                col_dens_for_sd = (col_dens - cumulative_col_dens/len(output_directories))**2
            else:
                col_dens_for_sd = col_dens_for_sd + ((col_dens - cumulative_col_dens/len(output_directories)))**2

        ax.errorbar(min_velocity_in_bin+vel_bin_width/2, cumulative_col_dens/(1e9*len(output_directories)*vel_bin_width),linewidth=3,yerr=1e-9*np.sqrt(col_dens_for_sd/len(output_directories)),elinewidth=1,capsize=2,color=colors[a])
        
    # Plot thermal velocity distribution Gaussian
    T = 200 # K
    m = 1.674e-27 # kg (mass of H atom)
    kB = 1.38e-23 # J/K (Boltzmann constant)

    velocities = []
    f_v = []
    lots_of_vel = np.linspace(-40,40,8001)

    for v, vel in enumerate(lots_of_vel):
        f_v.append(np.sqrt(m/(2*np.pi*kB*T)) * np.exp((-1*m*(lots_of_vel[v]*1e3)**2)/(2*kB*T)))
    
    ax.plot(lots_of_vel, np.multiply(f_v, 1e3*251000), color='deepskyblue',lw=3)
    

    # plot lines to mark 1) escape velocity, 2) echelle resolution, and 3) velocity equivalent to 7.31 eV
    ax.plot([5.03,5.03],[0,1e5],linewidth=1.5,color='k', zorder=-1)
    ax.plot([-5.03,-5.03],[0,1e5],linewidth=1.5,color='k', zorder=-1)

    ax.plot([10,10],[0,1e5],linewidth=1.5, color='dimgrey', zorder=-1, ls='--')
    ax.plot([-10,-10],[0,1e5],linewidth=1.5,color='dimgrey', zorder=-1, ls='--')
    
    ax.plot([37.45,37.45],[0,1e5],linewidth=1.5,color='darkgrey', zorder=-1, ls=':')
    ax.plot([-37.45,-37.45],[0,1e5],linewidth=1.5,color='darkgrey', zorder=-1, ls=':')

    ax.text(7,7e8,'300 km',color=colors[0],fontsize=16)
    ax.text(9,5e8,'1000 km',color=colors[1],fontsize=16)
    ax.text(29,1e8,'4000 km',color=colors[2],fontsize=16)

    for tick in ax.xaxis.get_majorticklabels():
        tick.set_fontsize(16)
    for tick in ax_twiny.xaxis.get_majorticklabels():
        tick.set_fontsize(16)
    for tick in ax.yaxis.get_majorticklabels():
        tick.set_fontsize(16)
    for tick in ax_twinx.yaxis.get_majorticklabels():
        tick.set_fontsize(16)

    ax.set_xlabel('Velocity (km/s)',fontsize=16)
    ax.set_ylabel('N$_v$ (10$^9$ atoms cm$^{-2}$ (km/s)$^{-1}$)', fontsize=16)
    ax.set_yscale('log')
    ax_twinx.set_yscale('log')
  #  ax.set_ylim([0,2])
    y_limits = [1e-2,4e3]
    ax.set_ylim(y_limits)
    ax_twinx.set_ylim(y_limits)
    ax.set_xlim([-40,40])
    ax_twiny.set_xlim([-40*121.56/(1e-3*3e5),40*121.56/(1e-3*3e5)])
    ax_twiny.set_xlabel('$\Delta \lambda$ (10$^{-3}$ nm)',fontsize=16)
    ax_twinx.set_ylabel('Brightness (R)',fontsize=16)
    plt.savefig('../figures/vel_dist.pdf',dpi=300,transparent=True)
    plt.show()
                
                  
                  
                                          

    return (ax)
# List altitudes at which we want to plot velocity distributions
fig1 = plt.figure(figsize=(18,5))
ax2 = fig1.add_subplot(1,1,1)
calc_avg_limb_vel_dist(ax2,[150]) #,300,1000]) #,1000,4000]))
#calc_avg_limb_vel_dist(ax2,[300,1000,4000])
