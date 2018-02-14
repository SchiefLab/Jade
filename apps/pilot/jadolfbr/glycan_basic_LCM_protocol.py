#!/usr/bin/env python

import os,sys
from argparse import ArgumentParser

import rosetta
import pyrosetta

pyrosetta.init("-include_sugars -write_pdb_link_records")


def get_parser():
    parser = ArgumentParser("This is a sample demonstration of the LinkageConformerMover in PyRosetta.  "
                            "For more, use GlycanRelax, which uses the LCM. ")

    parser.add_argument("--infile", "-s",
                        help = "Input PDB.",
                        required = True)

    parser.add_argument("--glycosylation_position", "-g",
                        help = "Glycosylation site. Rosetta resnum or resnumChain, ex: 463G",
                        required = True)

    parser.add_argument("--glycosylation_name", "-n",
                        help = "Glycosylation name",
                        default = "man5")

    parser.add_argument("--nstruct",
                        default = 1,
                        help = "Number of output structures")

    parser.add_argument("--cycles", "-c",
                        default = 75,
                        help = "Total number of cycles to attempt using the LCM")

    return parser

if __name__ == "__main__":


    parser = get_parser()
    options = parser.parse_args()




    p = pyrosetta.pose_from_pdb(options.infile)
    scorefxn = pyrosetta.get_score_function()
    KT = 1.0




    LCM = rosetta.protocols.carbohydrates.LinkageConformerMover()
    LCM.set_x_standard_deviations(2)
    SGM = rosetta.protocols.carbohydrates.SimpleGlycosylateMover()

    resnum = rosetta.core.pose.parse_resnum(options.glycosylation_position, p)
    SGM.set_position(resnum)
    SGM.set_glycosylation(options.glycosylation_name)
    SGM.apply(p)

    ### This will be unnessessary very soon!
    mm = pyrosetta.MoveMap()
    for i in range(1, p.total_residue() + 1):
        if (p.residue_type(i).is_carbohydrate()):
            mm.set_bb(i, True)

    LCM.set_movemap(mm)


    for i in range(1, int(options.nstruct) + 1 ):
        p_copy = pyrosetta.Pose(p)

        MC = pyrosetta.MonteCarlo(p_copy, scorefxn, KT)
        outname = os.path.basename(options.infile).replace(".pdb", "")+"_out_"+repr(i)+".pdb"
        for x in range(1, int(options.cycles)+1):
            print "LCM Round "+repr(x)

            LCM.apply(p_copy)
            MC.boltzmann(p_copy)

            print "Result: "+repr(scorefxn(p_copy))

        p_copy = MC.lowest_score_pose()

        print "End Score: "+repr(scorefxn(p_copy))
        print "Writing "+outname

        p_copy.dump_pdb(outname)
    print "Complete"





