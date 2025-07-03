Differential cross sections for 16^O+CO(v=0,j=0) -> 16^O+CO(v',j') scattering as a functions of angle and collision energy.
DCSs are constructed as a statistical average of 3 potential energy surfaces: 3A', 3A'', 2-3A''.
Total elastic cross sections can be reconstructed as:
   2*Pi*Integrate[eldcs[th]*Sin[th/180*Pi],{th,0,180}]*Pi/180   (in Mathematica notation)

Energy grid is given in the file DCS_energy_index.dat
Files O-CO_DCS-el_3pes* are for elastic channel only: O+CO(v=0,j=0) -> O+CO(v=0,j=0). 
Files O-CO_DCS-total_3pes* are for all channels: O+CO(v=0,j=0) -> O+CO(v',j'), where v'=0-20, j'=0-70

by Sanchit Chhabra, Marko Gacesa, Malathe S. Khalil, Amal Al Ghaferi, and Nayla El-Kork on September 15 2022 