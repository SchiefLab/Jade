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

    parser.add_argument("--native_json",
                        help = "Native JSONS, in the same way as the Analysis Jsons.",
                        required = True)

    parser.add_argument("--data_outdir","-o",
                        help = "Path to outfile. DEFAULT = data",
                        default = "data")

    parser.add_argument("--plot_outdir", "-p",
                        help = "DIR for plots. DEFAULT = plots/mc_benchmarks",
                        default = "plots/mc_benchmarks")

    #parser.add_argument("--root_dataset_dir",
    #                    help = "List of PDBIds to use for individual PDB output. DEFAULT = datasets/pdblists",
    #                    default = "datasets/pdblists")

    options = parser.parse_args()

    analysis_infos = [AnalysisInfo(json_file) for json_file in options.jsons]

    for analysis_info in analysis_infos:
        ab_db = get_ab_dir(analysis_info.bm_info)
        native_info = NativeInfo(analysis_info.bm_info.get_dataset(), analysis_info.bm_info.get_input_pdb_type())
        analyze_recovery = AnalyzeRecovery(ab_db, analysis_info, native_info)
        analyze_recovery.apply(analysis_info.get_features_db(), drop_tables=True)





def get_ab_dir(benchmark_info):

    if not benchmark_info.settings["paper_ab_db"]:
        return os.getenv("ROSETTA3_DB")+"/sampling/antibodies/antibody_database_rosetta_design.db"
    else:
        return os.getenv("ROSETTA3_DB") +"/sampling/antibodies/antibody_database_rosetta_design_north_paper.db"