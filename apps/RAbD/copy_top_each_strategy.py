#!/usr/bin/env python

#Copies the top N of a particular directory to a new directory.  Used in conjunction with RabD BM.  Got sick of doing this by hand.

from argparse import ArgumentParser
from basic import path
from basic import general
import glob
import os
import sys






if __name__ == "__main__":
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

    options = parser.parse_args()


    dirs = path.get_directories_recursively(options.indir)

    dirs.append(options.indir)
    
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    for d in dirs:
        print d
        if not general.match_patterns(d, options.strategies): continue

        print "\n"+d
        pdbs = glob.glob(d+"/*.pdb*")
        for pdb in pdbs:

            N = int(os.path.basename(pdb).split("_")[1])
            if N <= options.n:
                print os.path.basename(pdb)

                os.system('cp '+pdb+' '+options.outdir+"/"+os.path.basename(pdb))


