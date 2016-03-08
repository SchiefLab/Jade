#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser

from rosetta import *
rosetta.init("-include_sugars -read_pdb_link_records -ignore_unrecognized_res -ignore_zero_occupancy false")


parser = ArgumentParser("Simple app to scan a PDB file and print PDB info and Rosetta understood chains and resnums.")

parser.add_argument("pdb_file", help = "The PDB file to scan.", required = True)

options = parser.parse_args()


pose = pose_from_pdb(options.pdb_file)

print "#resnum chain_num chain pdb_num"
for i in range(1, pose.total_residue() +1 ):
    print repr(i)+" "+repr(pose.chain(i))+" "+pose.pdb_info().pose2pdb(i)

