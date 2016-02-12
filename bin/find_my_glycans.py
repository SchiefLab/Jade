#!/usr/bin/env python

from rosetta import *
rosetta.init("-include_sugars -read_pdb_link_records -write_pdb_link_records")
import sys
import os

#This code finds rosetta numbering of a glycan structure.
#Useful while we figure out this PDBInfo problem...

if len(sys.argv) == 1:
    sys.exit("Please specify a PDB on the command line as the only argument")

pose = pose_from_pdb(sys.argv[1])

for resnum in range(1, pose.total_residue() +1 ):
    if pose.residue( resnum ).is_carbohydrate():
        bp = pose.residue( resnum ).is_branch_point()
        print "Carbohydrate: "+repr( resnum ) +" BP: "+repr(bp)
    elif pose.residue( resnum ).is_branch_point():
        print "Branch Point: "+pose.residue( resnum ).name3()+" "+repr(resnum)

