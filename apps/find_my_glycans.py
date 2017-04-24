#!/usr/bin/env python

from rosetta import *
from pyrosetta import *

import sys
import os

#This code finds rosetta numbering of a glycan structure.
#Useful while we figure out this PDBInfo problem...

init("-include_sugars -write_pdb_link_records")
if len(sys.argv) == 1:
    sys.exit("Please specify a PDB on the command line as the only argument")

if sys.argv[1] == "--help":
    print "\n\nThis simple script aims to identify glycosylated positions in a PDB and thier associated rosetta Resnums."
    print "Please specifiy a PDB as the only argument to this script.\n"
    sys.exit()

pose = pose_from_pdb(sys.argv[1])

for resnum in range(1, pose.total_residue() +1 ):
    if pose.residue( resnum ).is_carbohydrate():
        bp = pose.residue( resnum ).is_branch_point()
        print "Carbohydrate: "+repr( resnum ) +" BP: "+repr(bp)
    elif pose.residue( resnum ).is_branch_point():
        print "Branch Point: "+pose.residue( resnum ).name3()+" "+repr(resnum)

