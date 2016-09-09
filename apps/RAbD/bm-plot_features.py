#!/usr/bin/env python

from argparse import ArgumentParser
from RAbD_BM.AnalysisInfo import *
from RAbD_BM.AnalyzeRecovery import *

def main():
    parser = ArgumentParser(description= "Calculates and plots monte carlo acceptance values for antibody design benchmarking.")


    ############################
    ## Required Options
    ############################
    parser.add_argument("--jsons","-j",
                        help = "Analysis JSONs to use.  See RAbD_MB.AnalysisInfo for more on what is in the JSON."
                               "The JSON allows us to specify the final name, decoy directory, and features db associated with the benchmark as well as all options that went into it.",
                        nargs = "*",
                        required = True)

    parser.add_argument("--plot_outdir", "-p",
                        help = "DIR for plots. DEFAULT = plots",
                        default = "plots")

    options = parser.parse_args()
    if not os.path.exists(options.plot_outdir): os.mkdir(options.plot_outdir)
    analysis_infos = [AnalysisInfo(json_file) for json_file in options.jsons]

