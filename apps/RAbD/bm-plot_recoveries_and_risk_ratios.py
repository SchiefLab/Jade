from argparse import ArgumentParser
from RAbD_BM.AnalysisInfo import *
from RAbD_BM.AnalyzeRecovery import *

def main():
    parser = ArgumentParser(description= "Plots recoveries and risk ratios")


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


    options = parser.parse_args()

    analysis_infos = [AnalysisInfo(json_file) for json_file in options.jsons]
    features_db_paths = [analysis_info.get_features_db() for analysis_info in analysis_infos]

    ##Plot this in Python, not R this time.  Control WHICH CDRs are plotted