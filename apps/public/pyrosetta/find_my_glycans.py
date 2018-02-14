#!/usr/bin/env python

import pyrosetta
import rosetta

import sys
import os

from argparse import ArgumentParser

#This code finds rosetta numbering of a glycan structure.
#Useful while we figure out this PDBInfo problem...

def get_parser():
    parser = ArgumentParser(description="This app is the PyRosetta equivalent of GlycanInfo.  "
                            "Print carbohydrate info about the pose. Pass the pose in as an argument")

    return parser


if __name__ == "__main__":

    pyrosetta.init("-include_sugars -write_pdb_link_records")
    if len(sys.argv) == 1:
        sys.exit("Please specify a PDB on the command line as the only argument")

    if sys.argv[1] == "--help":
        print "\n\nThis simple script aims to identify glycosylated positions in a PDB and thier associated rosetta Resnums."
        print "Please specifiy a PDB as the only argument to this script.\n"
        sys.exit()

    pose = pyrosetta.pose_from_pdb(sys.argv[1])

    for resnum in range(1, pose.total_residue() +1 ):
        if pose.residue( resnum ).is_carbohydrate():
            bp = pose.residue( resnum ).is_branch_point()
            print "Carbohydrate: "+repr( resnum ) +" BP: "+repr(bp)
        elif pose.residue( resnum ).is_branch_point():
            print "Branch Point: "+pose.residue( resnum ).name3()+" "+repr(resnum)

