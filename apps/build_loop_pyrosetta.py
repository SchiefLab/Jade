#!/usr/bin/env python

from __future__ import print_function
import os,re,sys
from argparse import ArgumentParser

from jade.rosetta_jade.flag_util import get_common_flags_string_for_init

from rosetta import *
from pyrosetta import *

init( get_common_flags_string_for_init() )



if __name__ == "__main__":
    parser = ArgumentParser( "This script builds a loop between two places in a structure with the given sequence, and closes the loop."
                    "It is not meant to be the last modeling step, just to create missing density or to prepare for loop modeling.")



    parser.add_argument("--start",
                        help = "Starting resnum.  Ex: 24L",
                        required = True)

    parser.add_argument("--stop",
                        help = "Ending resnum. Ex. 42L. ",
                        required = True)

    parser.add_argument("--sequence",
                        help = "Sequence of the loop",
                        required = True)

    parser.add_argument("--out_prefix",
                        help = "Any prefix to give results. ",
                        default = "loop_built_")

    parser.add_argument("--retain_aligned_roots",
                        help = "Attempt to keep any aligned root residues during the build",
                        default = False,
                        action = "store_true")

    parser.add_argument("--pdb", "-s",
                        help = "Input model",
                        required = True)

    options = parser.parse_args()

    print(options)
    pose = pose_from_pdb(options.pdb)

    refpose_name = "starting_model_refpose"
    pose.reference_pose_from_current(refpose_name, True)

    #core::Size resnum = pose.corresponding_residue_in_current( old_resnum , ref_pose_name)


    rosetta_start=0; rosetta_end = 0


    if options.retain_aligned_roots:

        print("\nAligning roots\n")
        start_search = rosetta.core.pose.parse_resnum(options.start, pose) + 1

        starting_alignment = ""
        for i in range(0, len( options.sequence )):
            model_aa = pose.residue_type( start_search + i ).name1()
            seq_aa = options.sequence[i]

            if model_aa == seq_aa:
                rosetta_start = start_search + i
                starting_alignment+=seq_aa
            else:
                break

        stop_search = rosetta.core.pose.parse_resnum(options.stop, pose)

        ending_alignment = []

        print(options.sequence)
        for i in range(0, len( options.sequence )):
            print(i, len(options.sequence))
            model_aa = pose.residue_type( stop_search - i ).name1()
            seq_aa = options.sequence[ len(options.sequence) - i -1]

            if model_aa == seq_aa:
                rosetta_end= stop_search - i
                ending_alignment.append(seq_aa)
            else:
                break

        ending_alignment.reverse()

        white_spaces = len(options.sequence) - ( len(starting_alignment) + len(ending_alignment) ) - 2
        print("Alignment: \n")

        print(options.sequence)
        print(starting_alignment,"".join(" " for i in range(0, white_spaces) ), "".join(ending_alignment))
        print(pose.pdb_info().pose2pdb(rosetta_start), pose.pdb_info().pose2pdb(rosetta_end))

        if rosetta_start == rosetta_end:
            sys.exit("The sequence is the same as the model!")



    else:
        rosetta_start = rosetta.core.pose.parse_resnum(options.start, pose)
        rosetta_end = rosetta.core.pose.parse_resnum(options.stop, pose)



