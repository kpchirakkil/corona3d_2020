# Script to plot Monte Carlo model output hot H velocity distributions
import pandas as pd
import matplotlib.pyplot as plt

def plot_limb_vel_dist(ax,altitudes):
    
    output_directories = ['../output/150km/vtally_1','../output/150km/vtally_3','../output/150km/vtally_8','../output/150km/vtally_9']

    for d,output_directory in enumerate(output_directories):
        
        for a,alt in enumerate(altitudes):
            min_velocity_in_bin = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['min_velocity_in_bin_[km/s]']
            vel_bin_width = min_velocity_in_bin[1] - min_velocity_in_bin[0]
            print(vel_bin_width)
            average_col_dens = pd.read_csv('%s/vtally_%skm.csv' % (output_directory, str(alt)),dtype=float)['column_density_[cm-2]']
            ax.plot(min_velocity_in_bin+vel_bin_width/2,average_col_dens/vel_bin_width,marker='x')

    ax.set_xlabel('Velocity (km/s)',fontsize=16)
    ax.set_ylabel('N$_v$ (atoms cm$^{-2}$ (km/s)$^{-1}$', fontsize=16)
    ax.set_yscale('log')
    plt.show()
                
    return (ax)
# List altitudes at which we want to plot velocity distributions
fig1 = plt.figure(figsize=(10,5))
ax2 = fig1.add_subplot(1,1,1)
plot_limb_vel_dist(ax2,([150,300,1000]))
