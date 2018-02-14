#!/usr/bin/env python

import pandas, os
from argparse import ArgumentParser
from jade.RAbD_BM.AnalysisInfo import *

from jade.RAbD_BM import tools_features_db

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

    print "created "+repr(len(analysis_infos))+" "+"Infos"
    dfs = []
    for analysis_info in analysis_infos:

        exp = analysis_info.get_exp()
        print "\nAnalyzing "+ exp

        ab_db = get_ab_dir(analysis_info.bm_info)
        feat_db = analysis_info.get_features_db()
        print "DB AT: "+feat_db

        df = tools_features_db.get_cdr_cluster_df(feat_db)
        df['exp'] = exp

        native_info = NativeInfo(analysis_info.bm_info.get_dataset(), analysis_info.bm_info.get_input_pdb_type())

        for pdbid in native_info.pdbids:
            df_pdbid = df[(df['input_tag'].str.contains(pdbid))]
            df_pdbid['pdbid'] = pdbid
            df_pdbid['decoy'] = df_pdbid['input_tag'].apply(os.path.basename)
            df_pdbid = df_pdbid.drop('input_tag', axis=1)
            df_pdbid = df_pdbid[['exp','pdbid','decoy','CDR','fullcluster','length','normDis_deg','sequence']]
            dfs.append(df_pdbid)

        #dfs.append(df)

    #Output the final data
    final_df = pandas.concat(dfs)
    final_df.to_csv(options.data_outdir+"/"+"all_final_clusters_lengths.csv", doublequote=False,index=False)


def get_ab_dir(benchmark_info):

    if not benchmark_info.settings.has_key("paper_ab_db"):
        return os.getenv("ROSETTA3_DB")+"/sampling/antibodies/antibody_database_rosetta_design.db"
    else:
        return os.getenv("ROSETTA3_DB") +"/sampling/antibodies/antibody_database_rosetta_design_north_paper.db"


if __name__ == "__main__":
    main()