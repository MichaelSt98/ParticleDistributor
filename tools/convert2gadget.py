#!/usr/bin/env python3

import argparse
import h5py
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert output files from ParticleDistributor to GADGET-4 format")
    parser.add_argument("--file", "-f", metavar="str", type=str, help="path to the h5 initial conditions file to convert", required=True)

    args = parser.parse_args()

    #reading input
    inpH5 = h5py.File(args.file, "r")

    posData = inpH5["x"]
    velData = inpH5["v"]
    massData = inpH5["m"]
    N = len(massData)
    
    outH5 = h5py.File("gadgetIC.h5", "w")

    # ADD HEADER
    header = outH5.create_group("Header")

    header.attrs.create("NumPart_ThisFile", data=[N], dtype="u4")
    header.attrs.create("NumPart_Total", data=[N], dtype="u8")
    print("WARNING: Only works with equal mass particles.")
    header.attrs.create("MassTable", data=[massData[0]], dtype="f8") # assuming equal mass particles
    header.attrs.create("NumFilesPerSnapshot", data=1, dtype="i4")
    header.attrs.create("Time", data=0., dtype="f8")
    header.attrs.create("Redshift", data=0., dtype="f8")
    header.attrs.create("BoxSize", data=0., dtype="f8")
    
    
    # ADD PARTICLE DATA
    particleType0 = outH5.create_group("PartType0")
    particleType0.create_dataset("Coordinates", data=posData, dtype="f8")
    particleType0.create_dataset("Velocities", data=velData, dtype="f8")
    particleType0.create_dataset("ParticleIDs", data=np.arange(len(massData)), dtype="u8")
    
