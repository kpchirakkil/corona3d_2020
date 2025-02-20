# Script to plot Figure 1, lower right panel for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pdb

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
plan_rad = 3397.0  # Mars radius (km)
ax_lim = 8400
x = plan_rad * np.outer(np.cos(u), np.sin(v))
y = plan_rad * np.outer(np.sin(u), np.sin(v))
z = plan_rad * np.outer(np.ones(np.size(u)), np.cos(v))

output_directory = '../../../safe/model_output_data/LSA/traced_positions_LSA'

# Plot initial positions of 1e5 test particles
datFile = '%s/positions1000.out' % (output_directory)
data = np.genfromtxt(datFile)
fig = plt.figure(figsize=[8,8])
ax = plt.axes(projection="3d")
ax.plot_wireframe(x, y, z, rstride=5, cstride=5, color='k') # Mars surface
ax.scatter(data[:,0] / 100000, data[:,1] / 100000, data[:,2] / 100000, marker='.',s=0.1,color='#c51160')
ax.set_xlabel('x (km)',fontsize=16)
ax.set_ylabel('y (km)',fontsize=16,labelpad=15)
ax.set_zlabel('z (km)',fontsize=16,labelpad=15)
ax.set_xlim(-ax_lim, ax_lim)
ax.set_ylim(-ax_lim, ax_lim)
ax.set_zlim(-ax_lim, ax_lim)
for tick in ax.xaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.xaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
for tick in ax.yaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.yaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
for tick in ax.zaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.zaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
ax.tick_params(axis='z', which='major', pad=9)
plt.savefig('initial_positions.png')
plt.show()

# Plot end positions of 1e5 test particles
datFile = '%s/positions4711000.out'  % (output_directory)
data = np.genfromtxt(datFile)
fig = plt.figure(figsize=[8,8])
ax = plt.axes(projection="3d")
ax.plot_wireframe(x, y, z, rstride=5, cstride=5, color='k') # Mars surface
ax.scatter(data[:,0] / 100000, data[:,1] / 100000, data[:,2] / 100000, marker='.',s=0.1,color='#c51160')
ax.set_xlabel('x (km)',fontsize=16)
ax.set_ylabel('y (km)',fontsize=16,labelpad=15)
ax.set_zlabel('z (km)',fontsize=16,labelpad=15)
ax.set_xlim(-ax_lim, ax_lim)
ax.set_ylim(-ax_lim, ax_lim)
ax.set_zlim(-ax_lim, ax_lim)
for tick in ax.xaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.xaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
for tick in ax.yaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.yaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
for tick in ax.zaxis.get_majorticklabels():
    tick.set_fontsize(14)
for tick in ax.zaxis.get_majorticklabels()[2:16:2]:
    tick.set_color('none')
ax.tick_params(axis='z', which='major', pad=9)
plt.savefig('end_positions.png',dpi=300)
plt.show()
