# Script to plot Figure 2 for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import os

# Set up axes
fig = plt.figure()
fig.set_figheight(10)
fig.set_figwidth(14)

ax1 = plt.subplot2grid(shape=(4, 40), loc=(0, 0), colspan=16, rowspan=2)
ax2a = plt.subplot2grid(shape=(4,40), loc=(0, 17), colspan=16, rowspan=2)
ax2b = ax2a.twiny()
ax3 = plt.subplot2grid(shape=(4, 40), loc=(2, 0), colspan=7, rowspan=1)
ax4 = plt.subplot2grid(shape=(4, 40), loc=(2, 7), colspan=7, rowspan=1)
ax5 = plt.subplot2grid(shape=(4, 40), loc=(2, 14), colspan=7, rowspan=1)
ax6 = plt.subplot2grid(shape=(4, 40), loc=(2, 21), colspan=7, rowspan=1)
ax7 = plt.subplot2grid(shape=(4, 40), loc=(3, 0), colspan=7, rowspan=1)
ax8 = plt.subplot2grid(shape=(4, 40), loc=(3, 7), colspan=7, rowspan=1)
ax9 = plt.subplot2grid(shape=(4, 40), loc=(3, 14), colspan=7, rowspan=1)
ax10 = plt.subplot2grid(shape=(4, 40), loc=(3, 21), colspan=7, rowspan=1)

output_directory_LSA = '../model_output_data/LSA/1/output'
output_directory_HSA = '../model_output_data/HSA/1/output'


# Read in density data
dens_data_day_LSA = np.genfromtxt('%s/density1d_day.out' % (output_directory_LSA))
dens_data_day_HSA = np.genfromtxt('%s/density1d_day.out' % (output_directory_HSA))

# Plot panel (a) ---  densities
markers_on = [145, 300, 1000, 4000]
ax1.plot(dens_data_day_LSA[:,1], dens_data_day_LSA[:,0], markevery=markers_on, marker='x',linewidth=1.5,mew=1,ms=8)#, mfc='r', mec='r')
ax1.plot(dens_data_day_HSA[:,1], dens_data_day_HSA[:,0], markevery=markers_on, marker='o',linewidth=1.5,markeredgecolor='tab:orange',markerfacecolor='none',ms=8)

# Panel (a) aesthetics
ax1.xaxis.set_tick_params(top='on', direction='in', which='both', labelsize=12)
ax1.yaxis.set_tick_params(right='on', direction='in', which='both', labelsize=12)
ax1.xaxis.set_tick_params(top='on', direction='in', which='major', labelsize=12,length=5)
ax1.yaxis.set_tick_params(top='on', direction='in', which='major', labelsize=12,length=5)
ax1.annotate('low solar\nactivity\n(LSA)', xy=(15, 800), xytext=(15, 800), color='tab:blue', fontsize=14,ha='center') #,fontweight='bold')
ax1.annotate('high\nsolar\nactivity\n(HSA)', xy=(1.2, 1000), xytext=(1.2, 1000), color='tab:orange', fontsize=14)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(90, 5000)
ax1.set_xlim(1, 200)
ax1.set_ylabel('Altitude (km)',fontsize=12)
ax1.set_xlabel('Dayside hot H density ($\mathrm{cm^{-3}}$)',fontsize=12)

output_directory_LSA = '../model_output_data/LSA/example/'
output_directory_HSA = '../model_output_data/HSA/example/'

# Plot panel (b) --- brightnesses
# Find nadir brightnesses
from density_integration import nadir_density
nadir_profile_LSA = []
nadir_profile_HSA = []
dens_prof = []
altitudes = range(100,5001,100)
for z in altitudes:
    I_LSA = nadir_density(z,'%s/density1d_day.out' % (output_directory_LSA))[1]
    I_HSA = nadir_density(z,'%s/density1d_day.out' % (output_directory_HSA))[1]
    nadir_profile_LSA.append(I_LSA)
    nadir_profile_HSA.append(I_HSA)
# Plot nadir brightnesses
ax2a.plot(nadir_profile_LSA,altitudes,label='Nadir LSA',color='tab:blue',linewidth=1)
ax2a.plot(nadir_profile_HSA,altitudes,label='Nadir HSA',color='tab:orange',linewidth=1)

# Find and plot limb brightnesses
from density_integration import limb_density
alt,col_dens,I = limb_density('%s/angleavg_dens.out' % (output_directory_LSA))
ax2a.plot(I,alt,label='Limb LSA',color='tab:blue',linewidth=2)
alt,col_dens,I = limb_density('%s/angleavg_dens.out' % (output_directory_HSA))
ax2a.plot(I,alt,label='Limb HSA',color='tab:orange',linewidth=2)

# Panel (b) aesthetics
ax2a.xaxis.set_tick_params(top='on', direction='in', which='both', labelsize=12)
ax2a.yaxis.set_tick_params(right='on', direction='in', which='both', labelsize=12)
ax2a.xaxis.set_tick_params(top='on', direction='in', which='major', labelsize=12,length=5)
ax2a.yaxis.set_tick_params(right='on', direction='in', which='major', labelsize=12,length=5)

ax2a.set_yscale('log')
ax2a.set_ylim(90, 5000)
ax2a.set_xlim(0, 12)
ax2a.set_xlabel('Brightness (R)',fontsize=12)
ax2b.xaxis.set_tick_params(top='on', direction='in', which='both', labelsize=12)
ax2b.yaxis.set_tick_params(right='on', direction='in', which='both', labelsize=12)
ax2b.set_xlim(0, 12)
ax2b.set_xlabel('Column density (10$^{9}$ cm$^{-3}$)',fontsize=12)

ax2a.text(0.1,2.5e3,'HSA\nnadir',color='tab:orange',fontsize=11)
ax2a.text(1.3,750,'LSA\nnadir',color='tab:blue',fontsize=14)
ax2a.text(3.5,320,'HSA\nlimb',color='tab:orange',fontsize=14,fontweight='bold')
ax2a.text(7.3,700,'LSA\nlimb',color='tab:blue',fontsize=14,fontweight='bold')

output_directory_LSA = '../model_output_data/LSA/1/output'
output_directory_HSA = '../model_output_data/HSA/1/output'

# Plot panels (g)--(j) --- energy distribution functions
v_min = 0.01
v_max = 100
my_cmap = mpl.cm.get_cmap('inferno')
my_cmap.set_bad('k')
my_norm = LogNorm(vmin=v_min, vmax=v_max)

EDF_day1_LSA = np.genfromtxt('%s/EDF_day_145km.out' % (output_directory_LSA))
EDF_day1_HSA = np.genfromtxt('%s/EDF_day_145km.out' % (output_directory_HSA))
im1 = ax7.imshow(EDF_day1_LSA.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day2_LSA = np.genfromtxt('%s/EDF_day_300km.out' % (output_directory_LSA))
EDF_day2_HSA = np.genfromtxt('%s/EDF_day_300km.out' % (output_directory_HSA))
im2 = ax8.imshow(EDF_day2_LSA.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day3_LSA = np.genfromtxt('%s/EDF_day_1000km.out' % (output_directory_LSA))
EDF_day3_HSA = np.genfromtxt('%s/EDF_day_1000km.out' % (output_directory_HSA))
im3 = ax9.imshow(EDF_day3_LSA.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day4_LSA = np.genfromtxt('%s/EDF_day_4000km.out' % (output_directory_LSA))
EDF_day4_HSA = np.genfromtxt('%s/EDF_day_4000km.out' % (output_directory_HSA))
im4 = ax10.imshow(EDF_day4_LSA.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)

# Panels (f)--(j) aesthetics
for ax in [ax7,ax8,ax9,ax10]:
    ax.annotate('LSA', xy=(6, -0.5), xytext=(6, -0.5), color='white', fontsize=12)


# Plot panels (c)--(f) --- energy distribution function, integrated over mu
EDF_LSAday_summed1 = np.sum(EDF_day1_LSA, axis=1)
EDF_LSAday_summed2 = np.sum(EDF_day2_LSA, axis=1)
EDF_LSAday_summed3 = np.sum(EDF_day3_LSA, axis=1)
EDF_LSAday_summed4 = np.sum(EDF_day4_LSA, axis=1)
EDF_HSAday_summed1 = np.sum(EDF_day1_HSA, axis=1)
EDF_HSAday_summed2 = np.sum(EDF_day2_HSA, axis=1)
EDF_HSAday_summed3 = np.sum(EDF_day3_HSA, axis=1)
EDF_HSAday_summed4 = np.sum(EDF_day4_HSA, axis=1)

EDF_LSAday_summed1 = EDF_LSAday_summed1*(.01)
EDF_LSAday_summed2 = EDF_LSAday_summed2*(.01)
EDF_LSAday_summed3 = EDF_LSAday_summed3*(.01)
EDF_LSAday_summed4 = EDF_LSAday_summed4*(.01)
EDF_HSAday_summed1 = EDF_HSAday_summed1*(.01)
EDF_HSAday_summed2 = EDF_HSAday_summed2*(.01)
EDF_HSAday_summed3 = EDF_HSAday_summed3*(.01)
EDF_HSAday_summed4 = EDF_HSAday_summed4*(.01)

x = np.arange(0, 10.05, 0.05)
ax3.plot(x, EDF_LSAday_summed1,linewidth=1.5)
ax3.plot(x, EDF_HSAday_summed1,linewidth=1.5)
ax4.plot(x, EDF_LSAday_summed2,linewidth=1.5)
ax4.plot(x, EDF_HSAday_summed2,linewidth=1.5)
ax5.plot(x, EDF_LSAday_summed3,linewidth=1.5)
ax5.plot(x, EDF_HSAday_summed3,linewidth=1.5)
ax6.plot(x, EDF_LSAday_summed4,linewidth=1.5)
ax6.plot(x, EDF_HSAday_summed4,linewidth=1.5)

# Panels (c)--(f) aesthetics
for ax in [ax3,ax4,ax5,ax6]:
    ax.set_yscale('log')
    ax.set_xlim(0, 10)
    ax.set_ylim(0.03, 300)
    ax.xaxis.set_major_locator(MultipleLocator(2))

ax3.text(3,1e1,'LSA',color='tab:blue',fontsize=14)
ax3.text(3,4e-1,'HSA',color='tab:orange',fontsize=14)

# General aesthetics
ax3.set_title('145 km')
ax4.set_title('300 km')
ax5.set_title('1000 km')
ax6.set_title('4000 km')
ax3.set_ylabel('Dist. ($\mathrm{cm^{-3} eV^{-1}}$)',fontsize=12)
ax7.set_ylabel('$\mu$',fontsize=12)

for ax in [ax7,ax8,ax9,ax10]:
    ax.set_xlim(0,10)
    ax.set_ylim(-1,1)
    ax.set_xlabel('Energy (eV)',fontsize=12)
    ax.xaxis.set_major_locator(MultipleLocator(2))
    if ax != ax10:
        for tick in ax.xaxis.get_majorticklabels()[6:8:2]:
            tick.set_color('none')

plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), fontsize=12)
plt.setp(ax7.get_yticklabels(), fontsize=12)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax6.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)
plt.setp(ax10.get_yticklabels(), visible=False)
plt.setp(ax8.get_xticklabels(), fontsize=12)
plt.setp(ax9.get_xticklabels(), fontsize=12)
plt.setp(ax10.get_xticklabels(), fontsize=12)

for ax in [ax3,ax4,ax5,ax6]:
    ax.tick_params('x',colors='k',which='both',direction='in',pad=5,top=True,right=True)
    ax.tick_params('y',colors='k',which='both',direction='in',pad=5,top=True,right=True)
ax3.tick_params('y',colors='k',which='major',direction='inout',pad=5,top=True,right=True,length=8)
ax3.tick_params('y',colors='k',which='minor',direction='inout',pad=5,top=True,right=True,length=4)

ax1.text(1.3e2,4e3,'(a)',fontsize=12)
ax2a.text(11,4e3,'(b)',fontsize=12)
ax3.text(8,100,'(c)',fontsize=12)
ax4.text(8,100,'(d)',fontsize=12)
ax5.text(8,100,'(e)',fontsize=12)
ax6.text(8,100,'(f)',fontsize=12)
ax7.text(8,0.7,'(g)',fontsize=12,color='w')
ax8.text(8,0.7,'(h)',fontsize=12,color='w')
ax9.text(8,0.7,'(i)',fontsize=12,color='w')
ax10.text(8,0.7,'(j)',fontsize=12,color='w')

# Colorbar
cax=plt.axes([0.675,0.08,0.015,0.385])
cbar = plt.colorbar(im1, cax=cax,label='Distribution ($\mathrm{cm^{-3} eV^{-1} \mu^{-1}}$)',pad=0.1) #,anchor=(0.5,0.1))
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Distribution ($\mathrm{cm^{-3} eV^{-1} \mu^{-1}}$)',fontsize=12)

fig.subplots_adjust(wspace=0, hspace=0,bottom=0.08,top=0.85)
fig.subplots_adjust(hspace=0)

pos = ax1.get_position()
ax1.set_position([pos.x0, pos.y0+(0.22*pos.height), pos.width*0.95, pos.height])

pos = ax2a.get_position()
ax2a.set_position([pos.x0-(0.03*pos.width), pos.y0+(0.22*pos.height), pos.width*0.95, pos.height])

plt.savefig('Figure2.pdf',dpi=1200)
plt.show()


