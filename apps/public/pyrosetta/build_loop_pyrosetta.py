#!/usr/bin/env python

from __future__ import print_function
import os,re,sys
from argparse import ArgumentParser

from jade.rosetta_jade.flag_util import get_common_flags_string_for_init
from jade.basic.path import get_decoy_name

import pyrosetta
import rosetta

pyrosetta.init( get_common_flags_string_for_init() )


def get_parser():
    parser = ArgumentParser( description="This script builds a loop between two places in a structure with the given sequence, and closes the loop."
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

    parser.add_argument("--kic",
                        help = "Run KIC peruturber after closing the loop?",
                        default = False,
                        action = "store_true")

    parser.add_argument("--dump_midpoints",
                        help = "Dump midpoint PDBs?",
                        default = False,
                        action="store_true")
    return parser

if __name__ == "__main__":


    parser = get_parser()
    options = parser.parse_args()

    print(options)
    pose = pyrosetta.pose_from_pdb(options.pdb)

    refpose_name = "starting_model_refpose"
    pose.reference_pose_from_current(refpose_name, True)

    #core::Size resnum = pose.corresponding_residue_in_current( old_resnum , ref_pose_name)


    rosetta_start=0; rosetta_end = 0

    seq_start = 0; seq_end = len(options.sequence) -1

    decoy_name = get_decoy_name(options.pdb)
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
                seq_start = i
                break

        stop_search = rosetta.core.pose.parse_resnum(options.stop, pose)

        ending_alignment = []

        print(options.sequence)
        for i in range(0, len( options.sequence )):
            #print(i, len(options.sequence))
            model_aa = pose.residue_type( stop_search - i ).name1()
            seq_aa = options.sequence[ len(options.sequence) - i -1]

            if model_aa == seq_aa:
                rosetta_end= stop_search - i
                ending_alignment.append(seq_aa)
            else:
                seq_end = len(options.sequence) - i - 1
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




    #Trim the pose.
    print("Trimming Pose")
    pose.delete_residue_range_slow(rosetta_start+1, rosetta_end-1)

    if options.dump_midpoints:
        pose.dump_pdb(options.out_prefix+"deleted"+"_"+decoy_name+".pdb")

    rsd_set = pose.residue_type_set_for_pose()

    #Build new residues
    print("Building new Residues")
    current_resnum = pose.corresponding_residue_in_current( rosetta_start , refpose_name)
    for i in range(seq_start, seq_end +1):


        #print("Setting", current_resnum, "to 180")
        #print("Adding "+options.sequence[ i ])
        res_type = rsd_set.get_representative_type_name1( options.sequence[ i ])
        res = rosetta.core.conformation.ResidueFactory.create_residue(res_type)
        rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue(pose.conformation(), current_resnum)
        pose.append_polymer_residue_after_seqpos(res, current_resnum, True)

        #print("seq end", seq_end, "i", i)
        #if (i != seq_end - 1):
            #pose.set_omega(current_resnum , 180.0)

        current_resnum+=1

    if options.dump_midpoints:
        pose.dump_pdb(options.out_prefix+"built"+"_"+decoy_name+".pdb")


    #Run CCD
    print("Running CCD to close the loop")

    #rosetta.protocols.grafting.dd_cutpoint_variants_for_ccd(pose, loop)

    new_end_rosettanum = pose.corresponding_residue_in_current( rosetta_end , refpose_name)
    new_start_rosettanum = pose.corresponding_residue_in_current( rosetta_start , refpose_name)

    rosetta.core.pose.add_variant_type_to_pose_residue(pose, rosetta.core.chemical.CUTPOINT_LOWER, new_end_rosettanum - 1)
    rosetta.core.pose.add_variant_type_to_pose_residue(pose, rosetta.core.chemical.CUTPOINT_UPPER, new_end_rosettanum)


    modeling_loop = rosetta.protocols.loops.Loop(new_start_rosettanum, new_end_rosettanum +1, new_end_rosettanum - 1  )
    print(modeling_loop)

    loops = rosetta.protocols.loops.Loops()
    loops.add_loop(modeling_loop)


    tree = rosetta.core.kinematics.FoldTree()
    rosetta.protocols.loops.fold_tree_from_loops( pose, loops, tree, True)

    pose.fold_tree(tree)
    ccd = rosetta.protocols.loops.loop_closure.ccd.CCDLoopClosureMover(modeling_loop)
    ccd.apply(pose)
    #Minimize the new residues to get something that doesn't look like total crap.



    if options.kic:
        print("Running KIC")
        perturber = rosetta.protocols.loops.loop_mover.perturb.LoopMover_Perturb_KIC(loops)


        typeset_swap = rosetta.protocols.simple_moves.SwitchResidueTypeSetMover("centroid")
        return_sidechains = rosetta.protocols.simple_moves.ReturnSidechainMover(pose)

        typeset_swap.apply(pose)

        perturber.apply(pose)

        return_sidechains.apply(pose)

        pose.dump_pdb("modeled.pdb")

    pose.dump_pdb(options.out_prefix+"closed_"+decoy_name+".pdb")