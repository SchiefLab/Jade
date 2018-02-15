#!/usr/bin/env python

from __future__ import print_function
import os,json,re,glob,sys
from argparse import ArgumentParser
from jade.rosetta_jade.score_util import parse_decoy_scores

def get_pdbs(argu):
    if os.path.isdir(argu):
        print("Gathering PDBs: " + argu)
        pdbs = glob.glob(argu+"/*.pdb*")
        return pdbs
    elif os.path.isfile(argu) and not re.search(".pdb", argu) and not re.search(".pdb.gz", argu):
        print("Parsing PDBs: " + argu)
        return [ x.strip() for x in open(argu, 'r').readlines() if not x.startswith('#') and x ]
    else:
        return [argu]


def get_parser():
    parser = ArgumentParser(description="This script creates a Rosetta score file from a set of structures - by parsing the score from them. Pass a directory, a PDBLIST, and/or a list of filenames")


    parser.add_argument("--prefix",
                        help = "Any prefix to use.  ",
                        default = "")

    parser.add_argument("decoys",
                        help = "A directory, a PDBLIST, and/or a list of filenames",
                        default = [],
                        nargs="*")
    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()

    #print(options)

    if len(options.decoys) == 0:
        sys.exit("Please pass decoys to parse score.")

    decoys = []
    for argu in options.decoys:

        pdbs = get_pdbs(argu)
        #print(pdbs)
        decoys.extend(pdbs)

    #print("\n".join(decoys))

    if options.prefix:
        OUTFILE = open(options.prefix+"score.json", 'w')
    else:
        OUTFILE = open(options.prefix + "score.json", 'w')

    scores = []
    decoy_num = 1
    print("Reading",len(decoys), "decoys")
    for decoy in decoys:
        if decoy_num % 50 == 0:
            print("Decoy",decoy_num)
        score_dict = parse_decoy_scores(decoy)
        if not score_dict:
            print("decoy", decoy, "has no score")
        if score_dict:
            OUTFILE.write(json.dumps(score_dict, sort_keys=True)+"\n")

        decoy_num+=1

    #OUTFILE.write("\n".join(json.dumps(scores)))

    OUTFILE.close()
    print("Done")