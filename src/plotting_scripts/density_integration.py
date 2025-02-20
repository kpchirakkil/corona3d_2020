# Script for "HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars"
# Bethan Gregory, Rodney Elliott, Justin Deighan, Hannes Groeller, Michael Chaffin
# November 2022

# Functions to integrate densities of output hot H density (dayside) values, looking straight down (nadir), straight up (zenith), or a limb view from a chosen spacecraft altitude.

# Funtions print integrated density (N, in cm-2) and resulting brightness (I, in Rayleighs), assuming I = gN. g = 1e-3 photons particle-1 s-1 (Anderson & Hord, 1971).

def nadir_density(sc_alt,filename):
    g = 1e-3 # g at Mars
    f1 = open(filename,'r')
    col_dens = 0
    for l,line in enumerate(f1.readlines()):
        if l > 0 and l < (sc_alt + 2):
            col_dens += float(line.split()[1])*1e5 # Assuming output file gives density values every 1 km
    return(col_dens,col_dens*g*1e-6)
    
def zenith_density(sc_alt):
    g = 1e-3 # g at Mars
    f1 = open('LSA/1/output/density1d_day.out','r')
    col_dens = 0
    for l,line in enumerate(f1.readlines()):
        if l > sc_alt:
            col_dens += float(line.split()[1])*1e5 # Assuming output file gives density values every 1 km
    return(col_dens,col_dens*g*1e-6)

def limb_density(filename): # densities integrated along the line of sight
    g = 1e-3 # g at Mars
    f1 = open(filename,'r')
    alt = []
    col_dens = []
    I = []
    for l,line in enumerate(f1.readlines()):
        line = line.split()
        alt.append(int(line[0]))
        col_dens.append(float(line[1]))
        I.append(float(line[1])*g*1e-6)
    return(alt,col_dens,I)

