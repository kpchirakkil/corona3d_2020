def calculate_N_from2d():
    import numpy as np
    import matplotlib.pyplot as plt
    
    grid_width = 200e5 #200e5 # cm
    mid =  257
    surface = 274
    
    g = 1e-3

    extra_factor = (500e-6*6.13e25)/(1e5*grid_width*grid_width)
    with open ('/Users/begr3234/corona3d2023/safe/model_output_data/LSA/1/output/density2d.out','r') as f1:
        lines = f1.readlines()
        z_0_line = lines[mid].split()
        N1 = []
        alt1 = []
        I1 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt1.append((grid_width/1e5)*(d-surface))
                N1.append(float(dens)*extra_factor)
                I1.append(float(dens)*g*1e-6*extra_factor)
 #       plt.plot(I1,alt1)



    grid_width = 100e5 #200e5 # cm
    mid = 513 # 257
    surface = 547 #274
 #   extra_factor = (500e-6*6.13e25)/(1e5*grid_width*grid_width)
    
    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFlinear_r1500_1e5test/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N2 = []
        alt2 = []
        I2 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt2.append((grid_width/1e5)*(d-surface))
                N2.append(float(dens))
                I2.append(float(dens)*g*1e-6)
#        plt.plot(I2,alt2)
#        plt.ylim([0,6000])
#        plt.xlabel('I (R)')
#        plt.ylabel('Altitude (km)')
    #    plt.xscale('log')

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test4/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N3 = []
        alt3 = []
        I3 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt3.append((grid_width/1e5)*(d-surface))
                N3.append(float(dens))
                I3.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test6/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N4 = []
        alt4 = []
        I4 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt4.append((grid_width/1e5)*(d-surface))
                N4.append(float(dens))
                I4.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test9/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N7 = []
        alt7 = []
        I7 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt7.append((grid_width/1e5)*(d-surface))
                N7.append(float(dens))
                I7.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test10/density2d_xy.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N8 = []
        alt8 = []
        I8 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt8.append((grid_width/1e5)*(d-surface))
                N8.append(float(dens))
                I8.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test10/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N9 = []
        alt9 = []
        I9 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt9.append((grid_width/1e5)*(d-surface))
                N9.append(float(dens))
                I9.append(float(dens)*g*1e-6)


    grid_width = 100e5  # cm
    mid =  91
    surface = 125
    
    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test11/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N10 = []
        alt10 = []
        I10 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt10.append((grid_width/1e5)*(d-surface))
                N10.append(float(dens))
                I10.append(float(dens)*g*1e-6)


    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test11/density2d_xy.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N11 = []
        alt11 = []
        I11 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt11.append((grid_width/1e5)*(d-surface))
                N11.append(float(dens))
                I11.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src3/output/test11_reproduce/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N12 = []
        alt12 = []
        I12 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt12.append((grid_width/1e5)*(d-surface))
                N12.append(float(dens))
                I12.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src3/output/test11_reproduce/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid-1].split()
        N14 = []
        alt14 = []
        I14 = []
        for d, dens in enumerate(z_0_line):
            if d > surface -1:
                alt14.append((grid_width/1e5)*(d-surface-1))
                N14.append(float(dens))
                I14.append(float(dens)*g*1e-6)

    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src3/output/test11_reproduce/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid+1].split()
        N15 = []
        alt15 = []
        I15 = []
        for d, dens in enumerate(z_0_line):
            if d > surface +1:
                alt15.append((grid_width/1e5)*(d-surface+1))
                N15.append(float(dens))
                I15.append(float(dens)*g*1e-6)

   # grid_width = 100e5 #200e5 # cm
   # mid = 513 # 257
   # surface = 547 #274
    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src3/output/test6_reproduce/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N13 = []
        alt13 = []
        I13 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt13.append((grid_width/1e5)*(d-surface))
                N13.append(float(dens))
                I13.append(float(dens)*g*1e-6)


    grid_width = 1e5  # cm
    mid =  9000
    surface = 12397
    with open('/Users/begr3234/corona_github_bethang/corona3d_2020/src/output/EDFslab_test8/density2d.out','r') as f2:
        lines = f2.readlines()
        z_0_line = lines[mid].split()
        N5 = []
        alt5 = []
        I5 = []
        for d, dens in enumerate(z_0_line):
            if d > surface:
                alt5.append((grid_width/1e5)*(d-surface))
                N5.append(float(dens))
                I5.append(float(dens)*g*1e-6)
                
    return(alt1,N1,I1, alt2,N2,I2, alt3,N3,I3, alt4,N4,I4, alt5,N5,I5, alt7,N7,I7,alt8,N8,I8,alt9,N9,I9, alt10,N10,I10,alt11,N11,I11, alt12,N12,I12,alt13,N13,I13, alt14,N14,I14,alt15,N15,I15)
    
#    plt.show()


#calculate_N_from2d()
            
