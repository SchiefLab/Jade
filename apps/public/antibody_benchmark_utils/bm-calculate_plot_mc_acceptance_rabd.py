#!/usr/bin/env python

import os
import sys
from collections import defaultdict
from RAbD_BM import tools as bm_tools
from RAbD_BM.AnalysisInfo import AnalysisInfo

from basic import path as path_tools
import re
import numpy
from basic.plotting.MakeFigure import MakeFigure

import matplotlib.pyplot as plot
import matplotlib as mpl
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, WeekdayLocator


from
import argparse

def main():
    parser = argparse.ArgumentParser(description= "Calculates and plots monte carlo acceptance values for antibody design benchmarking.")


    ############################
    ## Required Options
    ############################
    parser.add_argument("--jsons","-j",
                        help = "Analysis JSONs to use.  See RAbD_MB.AnalysisInfo for more on what is in the JSON."
                               "The JSON allows us to specify the final name, decoy directory, and features db associated with the benchmark as well as all options that went into it.",
                        nargs = "*",
                        required = True)

    parser.add_argument("--data_outdir","-o",
                        help = "Path to outfile. DEFAULT = data",
                        default = "data")

    parser.add_argument("--plot_outdir", "-p",
                        help = "DIR for plots. DEFAULT = plots/mc_benchmarks",
                        default = "plots/mc_benchmarks")

    parser.add_argument("--root_dataset_dir",
                        help = "List of PDBIds to use for individual PDB output. DEFAULT = datasets/pdblists",
                        default = "datasets/pdblists")

    options = parser.parse_args()
    analysis_infos = [AnalysisInfo(json_path) for json_path in options.jsons]

    analyzer = AnalyzeMCAcceptance(analysis_infos, options.root_dataset_dir, options.data_outdir, options.plot_outdir)
    analyzer.analyze()



if __name__ == "__main__":

    main()


