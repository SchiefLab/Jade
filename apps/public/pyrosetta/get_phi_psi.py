#!/usr/bin/env python

#Simple Script to get Phi/Psi using PyRosetta.
from __future__ import print_function
from argparse import ArgumentParser

import rosetta
from pyrosetta import *
pyrosetta.init()







if __name__ == "__main__":
    parser = ArgumentParser("Get Phi/Psi of all residues in protein or a range of residues")


    parser.add_argument('-s', help="Input Structure", required=True)
    parser.add_argument('--start', help = "Starting resnum (pose/PDB - EX:24L)")
    parser.add_argument('--span', help = "Number of residues to print from start")


    options = parser.parse_args()


    p = pose_from_pdb(options.s)

    start = 1
    if options.start:
        start = rosetta.core.pose.parse_resnum(options.start, p)

    end = p.total_residue()
    if options.span:
        end = start+int(options.span)

    out = "res resPDB phi psi"
    print(out)
    for i in range(start, end+1):
        out = str(i)+" "+p.pdb_info().pose2pdb(i)+" "+str(180+p.phi(i))+" "+str(180+p.psi(i))
        print(out)
