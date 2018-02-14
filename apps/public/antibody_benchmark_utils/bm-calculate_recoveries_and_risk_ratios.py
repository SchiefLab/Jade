#!/usr/bin/env python

from argparse import ArgumentParser
from jade.RAbD_BM.AnalysisInfo import *
from jade.RAbD_BM.AnalyzeRecovery import *

def get_parser():
    parser = ArgumentParser(description= "Calculates and plots monte carlo acceptance values for antibody design benchmarking.")


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

    return parser

def main():



    #parser.add_argument("--root_dataset_dir",
    #                    help = "List of PDBIds to use for individual PDB output. DEFAULT = datasets/pdblists",
    #                    default = "datasets/pdblists")

    parser = get_parser()
    options = parser.parse_args()

    if not os.path.exists(options.data_outdir): os.mkdir(options.data_outdir)

    analysis_infos = [AnalysisInfo(json_file) for json_file in options.jsons]

    for analysis_info in analysis_infos:
        print "Analyzing "+analysis_info.get_exp()

        ab_db = get_ab_dir(analysis_info.bm_info)
        native_info = NativeInfo(analysis_info.bm_info.get_dataset(), analysis_info.bm_info.get_input_pdb_type())

        #analyze_recovery = AnalyzeRecovery(ab_db, analysis_info, native_info)
        #analyze_recovery.apply(analysis_info.get_features_db(), drop_tables=True)

        #Write the data to another database.  Here, it is written stupidly, so we must recalculate everything.
        analyze_recovery = AnalyzeRecovery(ab_db, analysis_info, native_info)
        analyze_recovery.apply(options.data_outdir+"/all_recovery_and_risk_ratio_data.db", drop_tables=False)




def get_ab_dir(benchmark_info):

    print(benchmark_info.settings)

    if not benchmark_info.settings.has_key("PAPER_DB"):
        return os.getenv("ROSETTA3_DB")+"/sampling/antibodies/antibody_database_rosetta_design.db"

    elif benchmark_info.settings["PAPER_DB"] == 'False':
        return os.getenv("ROSETTA3_DB") +"/sampling/antibodies/antibody_database_rosetta_design_north_paper.db"



if __name__ == "__main__":
    main()