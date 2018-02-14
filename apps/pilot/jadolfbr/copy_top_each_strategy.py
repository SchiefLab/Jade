#!/usr/bin/env python

#Copies the top N of a particular directory to a new directory.  Used in conjunction with RabD BM.  Got sick of doing this by hand.

from argparse import ArgumentParser
from jade.basic import path
from jade.basic import general
import glob
import os
import sys




def get_parser():
    parser = ArgumentParser()
    parser.add_argument("-n",
                        help = "Number of models to copy. DEFAULT = 2",
                        default = 2)

    parser.add_argument("-i", "--indir",
                        help = "Input directory",
                        required = True)

    parser.add_argument("-o", "--outdir",
                        help = "Output directory",
                        required = True)

    parser.add_argument("-s", "--strategies",
                        help = "The type of strategies we are interested in",
                        default = ["delta_unsats_per_1000_dSASA", "dG_top_Ptotal"],
                        nargs = '*')
    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()


    dirs = path.get_directories_recursively(options.indir)

    dirs.append(options.indir)
    
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    for d in sorted(dirs):
        match, pattern = general.match_patterns(d, options.strategies)
        if not match: continue

        print "\n\n"+"_".join(os.path.basename(d).replace(pattern+"_", "").split(".")[0].split("_")[2:])+" "+pattern
        pdbs = glob.glob(d+"/*.pdb*")
        for pdb in pdbs:

            N = int(os.path.basename(pdb).split("_")[1])
            if N <= int(options.n):
                print os.path.basename(pdb)

                os.system('cp '+pdb+' '+options.outdir+"/"+os.path.basename(pdb))


