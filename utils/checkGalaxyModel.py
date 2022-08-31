#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import h5py

# CONFIGURATION SECTION
# ADJUST PARAMETERS BELOW TO MATCH THOSE IN compound_galaxy.info

# DISK 
M_d = 1.
h = 1.
z_0 = .2

# BULGE 
M_b = .33333
a = .2
c = .1

# HALO
M_h = 5.8
r_c = 10.
gamma = 1.
q = gamma/r_c
alpha = 1./(1.-np.sqrt(np.pi)*q*np.exp(q*q)*(1.-erf(q)))

# END CONFIGURATION

MAT_COLORS = ["#00ff00", "#0000ff", "#ff007f"]

# analytical functions

#cylindrical radius
def cylR(x):
    return np.sqrt(x[0]**2.+x[1]**2.)

#bulge
def m(x):
    return np.sqrt((x[0]**2.+x[1]**2.)/a**2. + x[2]**2./c**2.)

def rho_b(m):
    return M_b/(2.*np.pi*a**2.*c*m*(1.+m)**3.)
    #return M_b/(2.*np.pi*a*c**2.*m*(1.+m)**3.)


# halo
def rho_h(r):
    return M_h*alpha*np.exp(-r**2./(r_c**2.))/(2*np.pi**(3./2.)*r_c*(r**2.+gamma**2.))


# binning particles by sorted arrays
def binBy(radii, matIds, matId2bin, ppBin):
    rMin = 0.
    pCounter = 0
    rBins = []
    binsizes = []
    binRadii = []
    rMin = 0.

    print("Binning material ", matId2bin, " ...")
    
    for i, r in enumerate(radii):
        #print("Radius2bin ", r, ", index ", i, ", material ", matIds[i])
        if matIds[i] == matId2bin:
            #if pCounter == 0:
            #    rMin = r
            pCounter = pCounter+1
            binRadii.append(r)
            if pCounter == ppBin:
                rBins.append(np.mean(binRadii))
                binsizes.append(r - rMin)
                print("    ... bin ", rMin, " -- ", r, " full ...")
                pCounter = 0
                binRadii = []
                rMin = r

    print("... done.")
    
    binsizes = np.array(binsizes)
    rBins = np.array(rBins)
    return rBins, binsizes

# calculate mean velocity in bins
def binVelBy(radii, vels, matIds, matId2bin, ppBin):
    rMin = 0.
    pCounter = 0
    rBins = []
    vBins = []
    binRadii = []
    binVels = []
    rMin = 0.

    print("Binning material and calculating mean velocities ", matId2bin, " ...")
    
    for i, r in enumerate(radii):
        #print("Radius2bin ", r, ", index ", i, ", material ", matIds[i])
        #if pCounter == 0:
        #    rMin = r
        if matIds[i] == matId2bin:
            pCounter = pCounter+1
            binRadii.append(r)
            binVels.append(vels[i])
            if pCounter == ppBin:
                rBins.append(np.mean(binRadii))
                vBins.append(np.mean(binVels))
                print("    ... bin ", rMin, " -- ", r, " full with veocity ", vBins[-1], " ...")
                pCounter = 0
                binRadii = []
                rMin = r

    print("... done.")

    rBins = np.array(rBins)        
    vBins = np.array(vBins)
    return rBins, vBins


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Check density and velocity profiles of galaxy components.")
    parser.add_argument("--ic_file", "-i", metavar="str", help="initial conditions file to read velocities from", nargs="?", default="../output/ic.h5")
    parser.add_argument("--bins", "-b", type=int, metavar="int", help="number of bins to sample profiles", nargs="?", default=100)

    args = parser.parse_args()

    print("Reading '", args.ic_file, "' ...")
    data = h5py.File(args.ic_file, 'r')
    pos = data["x"][:]

    vel = data["v"][:]
    speeds = np.linalg.norm(vel, axis=1)
    matIds = data["materialId"][:]
    numParticles = len(data["m"])
    print("... done.")

    bins = args.bins

    unique, counts = np.unique(matIds, return_counts=True)
    N_mat = dict(zip(unique, counts))
    N_d = N_mat[0]
    N_b = N_mat[1]
    N_h = N_mat[2]
    
    #N_c = numParticles/3 # divisible by three when constructed by ParticleDistributor
    #print("Number of component particles N_c = ", N_c)
    
    if N_d % bins != 0: print("WARNING: ", N_d, " disk particles not binnable without remainder.")
    if N_b % bins != 0: print("WARNING: ", N_b, " bulge particles not binnable without remainder.")
    if N_h % bins != 0: print("WARNING: ", N_h, " halo particles not binnable without remainder.")

    # particles per bin
    ppBin_d = int(N_d/bins)
    ppBin_b = int(N_b/bins)
    ppBin_h = int(N_h/bins)
    
    # Bulge density profile

    # calculate elliptical radii m
    mRadii = np.apply_along_axis(m, 1, pos)

    mPerm = mRadii.argsort()
    mRadii = mRadii[mPerm]
    matIds_b = matIds[mPerm]

    mBins_b, binsizes_b = binBy(mRadii, matIds_b, 1, ppBin_b)

    # Halo density profile
    radii = np.linalg.norm(pos, axis=1)
    rPerm = radii.argsort()
    radii = radii[rPerm]
    matIds_h = matIds[rPerm]
    
    rBins_h, binsizes_h = binBy(radii, matIds_h, 2, ppBin_h)

    # Velocity profiles in cylindrical radius
    cylRadii = np.apply_along_axis(cylR, 1, pos)

    cylRPerm = cylRadii.argsort()
    cylRadii = cylRadii[cylRPerm]
    matIdsVel = matIds[cylRPerm]
    vels = speeds[cylRPerm]

    matLabels = ["Disk", "Bulge", "Halo"]
    ppBins = [ppBin_d, ppBin_b, ppBin_h]
    
    plt.figure(0, figsize=(8, 6), dpi=200)
    for matId in [0, 1, 2]:
        cylR, cylSpeed = binVelBy(cylRadii, vels, matIdsVel, matId, ppBins[matId])
        plt.plot(cylR, cylSpeed, 'x', color=MAT_COLORS[matId], label=matLabels[matId])

    plt.legend()
    plt.grid()
    plt.show()
    exit(0)

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{siunitx}')

    plt.figure(1, figsize=(8, 6), dpi=200)
    # BULGE
    # see Binney & Tremaine 2008, p. 86
    plt.loglog(mBins_b, M_b/N_b*ppBin_b/(binsizes_b*4.*np.pi*a**2.*c*mBins_b**2.), 'x', color=MAT_COLORS[1], label=r'sampled')

    mPlt = np.linspace(mBins_b[0], mBins_b[-1], 100*bins)
    plt.loglog(mPlt, rho_b(mPlt), 'k-', label=r"analytical")

    plt.title(r"Density profile of \textbf{Bulge} component")
    plt.xlabel(r"Elliptical radius $m$")
    plt.ylabel(r"Density $\varrho_b(m)$")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("output/bulgeDensity.png")
    
    plt.figure(2, figsize=(8, 6), dpi=200)
    # HALO
    # TODO: check why factor of 3./(2.*r_c**2.) fixes density profile
    plt.loglog(rBins_h, M_h/N_h*ppBin_h/(binsizes_h*4.*np.pi*rBins_h**2.), 'x', color=MAT_COLORS[2], label=r'sampled')

    rPlt = np.linspace(rBins_h[0], rBins_h[-1], 100*bins)
    plt.loglog(rPlt, rho_h(rPlt), 'k-', label=r"analytical")

    plt.title(r"Density profile of \textbf{Halo} component")
    plt.xlabel(r"Radius $r$")
    plt.ylabel(r"Density $\varrho_h(r)$")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("output/haloDensity.png")
    #plt.show() 

    
            
        
