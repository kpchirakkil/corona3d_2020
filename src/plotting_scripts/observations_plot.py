def observations_plot():
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os
    import math
    import pdb
    import pandas as pd
    from matplotlib.ticker import FormatStrFormatter

    plt.rcParams["ytick.major.size"] = 0
    
    # Set up figure and axes
    fig1 = plt.figure(figsize=(10,5))
    ax1 = fig1.add_subplot(1,2,1)
    ax2 = fig1.add_subplot(2,2,2)
    ax3 = fig1.add_subplot(2,2,4)
    ax4 = ax1.twiny()
    ax5 = ax2.twinx().twiny()
    ax6 = ax3.twinx().twiny()

#    # Plot Chaffin+ 2018 brightnesses for comparison:
#    min_I = np.multiply(pd.read_csv('Chaffin2018_min.csv')['I'],1000)
#    min_alt = pd.read_csv('Chaffin2018_min.csv')['alt']
#    max_I = np.multiply(pd.read_csv('Chaffin2018_max.csv')['I'],1000)
#    max_alt = pd.read_csv('Chaffin2018_max.csv')['alt']
#    ax1.plot(min_I,min_alt,color='darkorange')
#    ax1.plot(max_I,max_alt,color='darkorange')
  #  ax1.fill_betweenx(z_new_km,accumulated_prod_LSA_old,accumulated_prod_LSA_new, facecolor=accumulated_colors_LSA[i], alpha=1)plot(min_I
    

    # 1) Plot limb and nadir brightness profiles
    from density_integration import nadir_density
    nadir_profile = []
    nadir_profile2 = []
    dens_prof = []
    altitudes = [1000,2000,3000] #range(100,5001,100) #5001
    for z in altitudes:
        col_dens = nadir_density(z,'/Users/begr3234/corona3d2023/safe/model_output_data/LSA/1/output/density1d_day.out')[0]
        I = nadir_density(z,'/Users/begr3234/corona3d2023/safe/model_output_data/LSA/1/output/density1d_day.out')[1]
        col_dens2 = nadir_density(z,'/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test4/density1d_day.out')[0]
        I2 = nadir_density(z,'/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test4/density1d_day.out')[1]
       # nadir_profile.append(nadir_density(z)) # argument is spacecraft altitude in km
        nadir_profile.append(I)
        dens_prof.append(col_dens)
        nadir_profile2.append(I2)
    ax1.scatter(nadir_profile,altitudes,label='Nadir',color='tab:blue',marker='x',s=80)
 #   ax1.plot(nadir_profile2,altitudes,label='Nadir EDF slab test 4')
 #   ax4.plot(dens_prof,altitudes,color='k',linestyle='--')

    from calculate_N_from2d import calculate_N_from2d
    z2  = calculate_N_from2d()[0]
    limb_I = calculate_N_from2d()[2]
    limb_N = calculate_N_from2d()[1]
    z3  = calculate_N_from2d()[3]
    limb_N3 = calculate_N_from2d()[4]
    limb_I3 = calculate_N_from2d()[5]
    z4  = calculate_N_from2d()[6]
    limb_N4 = calculate_N_from2d()[7]
    limb_I4 = calculate_N_from2d()[8]
    z5  = calculate_N_from2d()[9] 
    limb_N5 = calculate_N_from2d()[10]
    limb_I5 = calculate_N_from2d()[11]
    z6  = calculate_N_from2d()[12] # Every 1 km - very noisy line
    limb_N6 = calculate_N_from2d()[13]
    limb_I6 = calculate_N_from2d()[14]
    z7  = calculate_N_from2d()[15]  # EDFslab_test9, 1e4 parts, xy plane, w=1e8??
    limb_N7 = calculate_N_from2d()[16]
    limb_I7 = calculate_N_from2d()[17]
    z8  = calculate_N_from2d()[18] # EDFslab_test10, 1e5 parts, xy plane
    limb_N8 = calculate_N_from2d()[19]
    limb_I8 = calculate_N_from2d()[20]
    z9  = calculate_N_from2d()[21] # EDFslab_test10, 1e5, xz plane
    limb_N9 = calculate_N_from2d()[22]
    limb_I9 = calculate_N_from2d()[23]
    z10  = calculate_N_from2d()[24] # EDFslab_test11, 1e5 parts, xz plane
    limb_N10 = calculate_N_from2d()[25]
    limb_I10 = calculate_N_from2d()[26]
    z11  = calculate_N_from2d()[27] # EDFslab_test11, 1e5 parts, xy plane 
    limb_N11 = calculate_N_from2d()[28]
    limb_I11 = calculate_N_from2d()[29]

    z12  = calculate_N_from2d()[30] # 
    limb_N12 = calculate_N_from2d()[31]
    limb_I12 = calculate_N_from2d()[32]
    z13  = calculate_N_from2d()[33] # 
    limb_N13 = calculate_N_from2d()[34]
    limb_I13 = calculate_N_from2d()[35]
    z14  = calculate_N_from2d()[36] # 
    limb_N14 = calculate_N_from2d()[37]
    limb_I14 = calculate_N_from2d()[38]
    z15  = calculate_N_from2d()[39] #
    limb_N15 = calculate_N_from2d()[40]
    limb_I15 = calculate_N_from2d()[41]

#    ax1.plot(limb_I, z2,label='Limb: 513 grid',marker='x',color='r') # Two options to remember to think about which is best...
#    ax1.plot(limb_I3, z3,label='Limb: 1025 grid',marker='x',color='tab:blue') #... but this one is probably best, might just need smoothing
#    ax1.plot(limb_I4, z4,label='Limb: 1025 grid EDF slab test 4',marker='x',color='grey') #... but this one is probably best, might just need smoothing
#    ax1.plot(limb_I5, z5,label='Limb: from 2d count; single z=0 line',marker='x',color='tab:orange') #Limb: 1025 grid EDF slab test 6#... but this one is probably best, might just need smoothing
  #  ax1.plot(limb_I6, z6,label='Limb: from 2d count; single z=0 line,1km2 bins',marker='x')
#    ax1.plot(limb_I7, z7,label='Limb: 2d; y=0; 1e4; test9; w=1e8?',marker='o',color='tab:green')
#    ax1.plot(limb_I8, z8,label='Limb: 2d; y=0; 1e5; test10',marker='o',color='fuchsia')
#    ax1.plot(limb_I9, z9,label='Limb: 2d; z=0; 1e5; test10',marker='x',color='fuchsia')
#    ax1.plot(limb_I10, z10,label='Limb: 2d; z=0; 1e5; test11',marker='x',color='darkred')
#    ax1.plot(limb_I11, z11,label='Limb: 2d; y=0; 1e5; test11',marker='o',color='darkred')
#    ax1.plot(limb_I12, z12,label='Limb: 2d; z=0; 1e5 seed; test11',marker='x',color='darkblue')
#    ax1.plot(limb_I13, z13,label='Limb: 2d; z=0; 1e5 seed; test6',marker='+',color='darkblue')
#    ax1.plot(limb_I14, z14,label='Limb: 2d; z=0; 1e5 seed; test11, -1',marker='*',color='darkblue')
#    ax1.plot(limb_I15, z15,label='Limb: 2d; z=0; 1e5 seed; test11, +1',marker='^',color='darkblue')
#    ax4.plot(limb_N5,z5)
 ##   ax4.plot(limb_N,z2,color='k',linestyle='--')
##    ax4.plot(limb_N3,z3,color='k',linestyle='--')

    # Limb brightness profiles from alternative clever method
    from density_integration import limb_density
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test5/EDF_linear_dens_Wed.out')
#    ax1.plot(I,alt,color='darkblue',label='new limb Wed method 100',marker='x')
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test6/EDF_linear_dens_Wed.out')
#    ax1.plot(I,alt,label='Limb',color='skyblue') #,marker='x') #,color='skyblue' # Note - NOT  multiplying by 2 for now to account for fact that this run missed a factor of 2
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test7/EDF_linear_dens_Wed.out')
#    ax1.plot(I,alt,label='Limb w=1e4',color='r')
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test8/EDF_linear_dens_Wed.out')
#    ax1.plot(I,alt,label='Limb w=1e8',color='b')
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test9/EDF_linear_dens_Wed.out')
#    ax1.plot(I,alt,label='Limb: AA; 1e4; test9; w=1e8?',color='tab:green')
    alt,col_dens,I = limb_density('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test11/angleavg_dens.out')
    ax1.plot(I,alt,label='Limb: AA; 1e5; test11',color='darkred')

    # 2) Plot limb velocity distributions
    from Nv_v import Nv_v
    from plot_limb_vel_dist import plot_limb_vel_dist
    plot_limb_vel_dist(ax2,[300,1000,4000]) # send to function to plot average LOS velocity in limb observation at chosen altitudes
    for z in [300,1000,4000]:
        v, Nv = Nv_v(z,'limb')
  #      ax2.plot(v,np.divide(Nv,2*(3.397e8+5000e5)),label='%s km' % (str(z)))
#        ax2.plot([4.81,4.81],[0,6e22],linestyle='--',color='tab:blue') # escape velocity at 300 km
#        ax2.plot([4.41,4.41],[0,6e22],linestyle='--',color='tab:orange') # escape velocity at 1000 km
#        ax2.plot([3.40,3.40],[0,6e22],linestyle='--',color='tab:green') # escape velocity at 4000 km
        v,Nv = Nv_v(z,'limb_weighted')
        Nv = np.divide(Nv,1e5) # divide by 1e5 cm because the plane at constant x is 1 km thick
        d_lamda = np.divide(v,(1e-3*3e5/121.56)) # delta lamda in pm
        I = np.multiply(Nv,1e-9) # brightness - multiply Nv by 1e-3 and 1e-6 to get I in R
   ##     ax5.plot(d_lamda,I,color='k',linestyle='--')
#        ax2.plot(v,Nv,label='%s km' % (str(z))) 
       # n,Nv = Nv_v(z,'limb_single')

    # 3) Plot nadir velocity distributions
    for z in [300,1000,4000]:
        v,Nv = Nv_v(z,'nadir')
        v2,Nv2 = Nv_v(z,'nadir2')
        Nv = np.divide(Nv,math.pi*1000e5*1000e5) # divide by cross-section of cylinder (here radius=1000 km, but beware hardcoded)
        Nv2 = np.divide(Nv2,math.pi*1000e5*1000e5) # divide by cross-section of cylinder (here radius=1000 km, but beware hardcoded)
        I = np.multiply(Nv2,1e-9) # brightness - multiply Nv by 1e-3 and 1e-6 to get I in R
        d_lamda = np.divide(v2,(1e-3*3e5/121.56)) # delta lamda in pm
   #     ax3.plot(v,Nv,label='%s km' % (str(z)))
        ax3.plot(v2,Nv2,label='%s km' % (str(z)))
        ax6.plot(d_lamda,I)

        integrated_N = 0
        for n in Nv:
            integrated_N += n*2
        print(z,integrated_N)
        

    
    # Plot echelle channel limits: 20 km/s
# CHANGE BACK    ax2.plot([10,10],[0,6e22],color='k')
# CHANGE BACK    ax2.plot([-10,-10],[0,6e22],color='k')
    ax3.plot([-10,-10],[0,2.1e8],color='k')
    ax3.plot([10,10],[0,2.1e8],color='k')
#    ax2.set_yscale('log')
        
    
    # Plot aesthetics
    ax1.set_xlabel('I (R)')
    ax4.set_xlabel('N (10$^{10}$ cm$^{-2}$)')
    ax1.set_xlim([0,12]) #14000]) # [0,12] #[0,22])
    ax4.set_xlim([0,12e9]) #14e12]) #12e9]) #[0,22e9])
  #      ax.set_yscale('log')
    ax3.set_ylabel('v (km/s)')
    ax1.set_ylabel('Altitude (km)')
    for ax in [ax3]: #CHANGE BACK[ax2,ax3]:
        ax.legend()
        ax.set_ylabel('N$_v$ (10$^8$ cm$^{-2}$/(km/s))')
   
# CHANGE BACK    ax5.set_xlabel('$\Delta $\u03BB (10$^{-3}$ nm)')
#    ax6.set_xlabel('$\Delta $\u03BB (nm)')

    ax7 = ax2.twinx()
    ax8 = ax3.twinx()
    ax9 = ax2.twiny()
    ax10 = ax3.twiny()

    for ax in [ax3]: # CHANGE BACK[ax2,ax3]:
        ax.set_xlim([-40,40])
        ax.ticklabel_format(axis='y',style='sci')
      #  ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e1'))
    for ax in [ax6,ax10]: # CHANGE BACK[ax5,ax6,ax9,ax10]:
        ax.set_xlim([-40*121.56/(1e-3*3e5),40*121.56/(1e-3*3e5)])
    #ax2.set_ylim([0,1.5e17])
# CHANGE BACK    ax2.set_ylim([0,1.5e8])
# CHANGE BACK    for ax in [ax5,ax7]:
        #ax.set_ylim([0,1.5e17*1e-9])
# CHANGE BACK        ax.set_ylim([0,1.5e8*1e-9])
    ax3.set_ylim([0,2.1e8])

    for ax in [ax6,ax8]:
        ax.set_ylim([0,2.1e8*1e-9])
    ax5.set_ylabel('I here')
   # ax7.set_ylim([0,1.5e17*1e-9])
   # ax8.set_ylim([0,2.1e8*1e-9])
    
    ax3.set_xlabel('v (km/s)')
    ax2.text(-36,1.35e17,'Limb',fontweight='bold')
    ax3.text(-36,1.9e8,'Nadir',fontweight='bold')
  #  ax1.text(0,110,'(a)',fontweight='bold')
    ax2.text(-39,0.05e17,'(b)',fontweight='bold')
    ax3.text(-39,0.05e8,'(c)',fontweight='bold')
    


    ax1.set_ylim([90,5e3])
    ax1.set_yscale('log')


 #   ax1.legend()
    
    fig1.subplots_adjust(hspace=0.04,left=0.08)

 # CHANGE BACK   ax7.set_ylabel('I (R)')
    ax8.set_ylabel('I (R)')

    ax1.tick_params('both',width=1,length=8,direction='in',right=True)
    ax1.tick_params('both',which='minor',width=1,length=3,direction='in',right=True)

    for ax in [ax2,ax5,ax6,ax10]:
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_color('none')
    for ax in [ax5,ax6,ax7,ax8]:
        for tick in ax.yaxis.get_majorticklabels():
            tick.set_color('none')
 #   for ax in [ax6,ax7]:
 #       ax.tick_params('both',width=1,length=0,direction='in')
    for ax in [ax2,ax3,ax4,ax7,ax8,ax9,ax10]:
        ax.tick_params('both',width=1,length=8,direction='in') 
        ax.tick_params('both',which='minor',width=1,length=3,direction='in')
    for ax in [ax5,ax6]:
        ax.tick_params('both',width=1,length=0,direction='in')
  #  for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]:
  #      ax.tick_params('both',which='both',width=1,length=0,direction='in') 
        

    plt.savefig('obs_test.pdf')
    plt.show()
    
observations_plot()
                    
