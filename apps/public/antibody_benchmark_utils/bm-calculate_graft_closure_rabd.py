#!/usr/bin/env python

import os
import sys
from collections import defaultdict
import re
import glob
from argparse import ArgumentParser
import re
import gzip
from jade.RAbD_BM import tools as bm_tools

cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]


def get_parser():
    parser = ArgumentParser(description="Calculate the frequence of graft closures.")

    ############################
    ## Required Options
    ############################
    parser.add_argument("--dir", "-i",
                        help="Input directory"
                        )

    parser.add_argument("--outfile", "-o",
                        help="Path to outfile"
                        )

    parser.add_argument("--use_ensemble",
                        default=False,
                        action="store_true",
                        help="Use ensembles in calculation")

    parser.add_argument("--match_name",
                        help="Match a subexperiment in the file name such as relax")
    return parser

class CalculateGraftClosure:
    def __init__(self, parse_args = True):
        if parse_args:
            self._parse_args()

        self.totals = ClosureData("totals")
        self.cdr_data = defaultdict()
        for cdr in cdrs:
            self.cdr_data[cdr] = ClosureData(cdr)

    def set_options(self, in_dir, outfile, use_ensemble = False, match_name = ""):
        self.in_dir = in_dir; self.outfile = outfile
        self.use_ensemble = use_ensemble; self.match_name = match_name

    def _parse_args(self):


        parser = get_parser()

        options, args = parser.parse_args()

        self.options = options

        self.outfile = options.outfile
        self.in_dir = options.dir
        self.use_ensemble = options.use_ensemble
        self.match_name = options.match_name

    def run_main(self):

        if not os.path.exists(self.in_dir):
            sys.exit("Inpath does not exist! "+self.in_dir)

        if os.path.exists(self.outfile):
            OUTFILE = open(self.outfile, 'a')
        else:
            outdir = os.path.dirname(self.outfile)
            if not os.path.exists(outdir):
                os.mkdir(outdir)

            OUTFILE = open(self.outfile, 'w')
            line = "#experiment\tclosure\tclosure_freq\ttotal_attempts"
            for cdr in cdrs:
                line = line+"\t"+cdr+"_"+"closure"
            for cdr in cdrs:
                line = line+"\t"+cdr+"_"+"closure_freq"+"\t"+"totals"
            OUTFILE.write(line+"\n")



        exp_name = os.path.basename(self.in_dir)
        files = bm_tools.get_pdb_paths(self.in_dir, exp_name, self.match_name, self.use_ensemble)

        for file in files:
            print file
            self.parse_add_data(file)

        line = exp_name+"\t"+self.totals.get_decimal_format_prob()+"\t"+repr(self.totals.get_closures())+"\t"+repr(self.totals.get_total())

        for cdr in self.cdr_data:
            line = line+"\t"+self.cdr_data[cdr].get_decimal_format_prob()

        for cdr in self.cdr_data:
            line = line + "\t"+repr(self.cdr_data[cdr].get_closures())+"\t"+repr(self.cdr_data[cdr].get_totals())

        OUTFILE.write(line+"\n")

    def parse_add_data(self, file):

        if file.split(".")[-1] =="gz":
            #print "opening gzipped file"
            INFILE = gzip.open(file, 'rb')
        else:
            INFILE = open(file, 'r')

        for line in INFILE.readlines():
            #print line
            if not re.search("GRAFT_CLOSURE", line):
                continue
            #print line
            lineSP = line.split()

            closed = int(lineSP[5])
            cdr = lineSP[7]
            cluster = lineSP[8]
            origin_pdb = lineSP[-1]

            self.totals.total_attempts+=1
            self.cdr_data[cdr].total_attempts+=1
            if closed:
                self.totals.total_closures+=1
                self.cdr_data[cdr].total_closures+=1

        INFILE.close()

class ClosureData:
    def __init__(self, name):
        self.name = name
        self.total_closures = 0
        self.total_attempts = 0

    def get_total(self):
        return self.total_attempts
    def get_totals(self):
        return self.total_attempts

    def get_closures(self):
        return self.total_closures

    def get_prob(self):
        if self.total_attempts == 0:
            return 0

        return self.total_closures/(float(self.total_attempts))

    def get_decimal_format_prob(self):
        return "%.4f"%self.get_prob()

    def get_percent_format_prob(self):
        return "%.2f"%(self.get_prob()*100)






if __name__ == "__main__":
    calculate_closure = CalculateGraftClosure()
    calculate_closure.run_main()





