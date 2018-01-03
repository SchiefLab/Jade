#! /usr/bin/env python
#Author: Jared Adolf-Bryfogle

import os,sys,re
from rosetta import *
from rosetta.core.conformation.carbohydrates import *
from rosetta.core.pose.carbohydrates import *
from pyrosetta import *


rosetta_options="""

-include_sugars 
-ignore_unrecognized_res

-ignore_zero_occupancy false

-cryst::crystal_refine

-alternate_3_letter_codes pdb_sugar
-auto_detect_glycan_connections
-min_bond_length 1.30
-max_bond_length 1.50


"""

init(rosetta_options.replace('\n', " "))

if len(sys.argv) == 1: sys.exit("Pass a structure.  "
                                "This simple code checks to make sure all dihedrals of "
                                "carbohydrates can be set in PyRosetta.  It is a debugging tool.")

pose = pose_from_pdb(sys.argv[1])

for i in range(1, pose.total_residue() + 1):
    if (pose.residue_type(i).is_carbohydrate()):
        #print "Residue: "+repr(i)+" "+pose.pdb_info().pose2pdb(i)
        n_torsions = get_n_glycosidic_torsions_in_res(pose, i)
        dih_angle = 5.0
        for torsion in range(1, n_torsions + 1):
            try:
                #print "  Torsion "+repr(torsion)
                dih_angle = get_glycosidic_torsion(torsion, pose, i)
            except Exception:
                print "Could not get torsion for torsion "+repr(torsion)+" res: "+repr(i)+" "+pose.pdb_info().pose2pdb(i)

            try:
                set_glycosidic_torsion(torsion, pose, i, dih_angle + 10.0)
            except Exception:
                print "Could not set torsion for torsion " + repr(torsion) + " res: " + repr(
                    i) + " " + pose.pdb_info().pose2pdb(i)

print "Done"