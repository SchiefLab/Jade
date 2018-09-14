#!/usr/bin/env python
from argparse import ArgumentParser
import os
import sys


if __name__ == "__main__":
    parser = ArgumentParser('Create Maps for use in Rosetta from CIF files and PDBs using Phenix default options')
    parser.add_argument("--cif_dir", help = "Directory with CIF files downloaded from rcsb", default = "cif")
    parser.add_argument("--pdb_dir", help = "Directory with PDB files (Symmetrized)", default = "pdbs")
    parser.add_argument("--pdb_roots", help = "pdb_roots_density file", default="pdb_roots_density_sym.txt")
    parser.add_argument("--pdb_suffix", help = "Suffix of the PDB", default = ".pdb")

    options = parser.parse_args()

    FILE = open(options.pdb_roots, 'r')
    pdbs = []
    for line in FILE:
        line = line.strip()
        if not line: continue
        if line.startswith('#'): continue
        lineSP = line.split()
        #print(lineSP)
        pdb = lineSP[0]
        pdbs.append(pdb)
    FILE.close()


    for pdb in pdbs:
        pdb_path = options.pdb_dir+"/"+pdb+options.pdb_suffix
        cif_path = options.cif_dir+"/"+pdb+"-sf.cif"
        cmd = 'phenix.maps '+pdb_path+' '+cif_path
        print cmd
        os.system(cmd)

print "done!"