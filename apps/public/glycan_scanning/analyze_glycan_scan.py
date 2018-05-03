#!/usr/bin/env python
from __future__ import print_function
import os,sys
import jade.basic.path as path_utils
import jade.rosetta_jade.ScoreFiles as score_utils
import pandas as pd
import seaborn.apionly as seaborn
import matplotlib as mpl
from matplotlib import cm
import shutil
from argparse import ArgumentParser
from jade.pymol_jade.PyMolScriptWriter import *

def get_position(name):
    return name.split('_')[1]



if __name__ == "__main__":
    parser = ArgumentParser("Analyze data output from a glycan_scanning job created using the JADE create_glycan_scanning_job script")

    parser.add_argument("-s",
                        default= "scores.json",
                        help = "Score file with extra metrics. ")

    parser.add_argument("--outdir_prefix",
                        default="analysis",
                        help = "Prefix for output directories")

    parser.add_argument("--get_top", action = "store_true",
                        default = False,
                        help = "Get pymol sessions and top models")

    parser.add_argument("--get_plots", action = "store_true",
                        default = False,
                        help = "Get plots for analysis")

    parser.add_argument("--native_path",
                        help = "Path to native")

    options = parser.parse_args()

    # Load and organize the Data
    data = score_utils.ScoreFile(options.s)
    df = data.get_Dataframe()
    df['decoy_name'] = df["decoy"].apply(path_utils.get_decoy_name)
    df['position'] = df['decoy_name'].apply(get_position)
    df['selection'] = pd.to_numeric(df['selection'], downcast='unsigned')
    df_origin = df

    if options.get_top:

        #Top Scoring Per-Sequon
        df.reset_index().set_index(['position', 'post-model_total_energy'], inplace=True)
        df.sort_index().reset_index(inplace=True)

        top_scoring = df.groupby('position').head(1)
        top_scoring.reset_index(inplace=True)
        top_scoring.set_index('post-model_total_energy', inplace=True)
        top_scoring.sort_index(inplace=True)
        top_scoring.reset_index(inplace=True)

        paths_top = top_scoring['decoy']
        scores_top = top_scoring['post-model_total_energy']
        names_top = top_scoring['position']
        score_paths_top = zip(scores_top, paths_top, names_top)

        out_name = "all_sequons_top_models"
        maindir = options.outdir_prefix+"_models"
        outdir = maindir+"/"+out_name
        if not os.path.exists(maindir):
            os.mkdir(maindir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        make_pymol_session_on_top_scored(score_paths_top, outdir, outdir, out_name, antibody=False, native_path=options.native_path)

        #Copy Top Models
        for i, tup in enumerate(score_paths_top):
            name = "model_" + repr(i) + "_" + tup[2] + "_%.2f" % (tup[0]) + path_utils.get_decoy_extension(tup[1])
            shutil.copy(tup[1], outdir+'/'+name)


        #Top 15 Scoring Per-Sequon
        top_scoring = df.groupby('position').head(15)
        top_scoring.reset_index(inplace=True)
        top_scoring.set_index('post-model_total_energy', inplace=True)
        top_scoring.sort_index(inplace=True)
        top_scoring.reset_index(inplace=True)

        paths = top_scoring['decoy']
        scores = top_scoring['post-model_total_energy']
        names = top_scoring['position']
        score_paths = zip(scores, paths, names)

        organized_names = defaultdict()
        for score_path in score_paths:
            name = score_path[2]
            if not organized_names.has_key(name):
                organized_names[name] = []
            else:
                organized_names[name].append(name)

        root_out = maindir+"/top_15"
        if not os.path.exists(root_out):
            os.mkdir(root_out)

        for key in sorted(organized_names.keys()):
            tup = organized_names[key]
            outdir = root_out+'/'+key
            if not os.path.exists(outdir):
                os.mkdir(outdir)

            make_pymol_session_on_top_scored(organized_names[key], outdir, outdir, key, antibody=False,
                                             native_path=options.native_path)

            for i, tup in enumerate(organized_names[key]):
                name = "model_" + repr(i) + "_" + tup[2] + "_%.2f" % (tup[0]) + path_utils.get_decoy_extension(tup[1])
                shutil.copy(tup[1], outdir+'/'+name)

        scripter = make_pymol_session_on_top_scored(score_paths, root_out, root_out, "all_top_15", antibody=False,
                                         native_path=options.native_path, run_pymol=False)

        #Group all of them.
        for tup in score_paths_top:
            name = tup[2]
            scripter.add_line("group models_"+name+"_%.2f" % (tup[0])+', models_'+name+'*' )

        pse_path = root_out + "/" + out_name + ".pse"
        scripter.add_save_session(pse_path)
        scripter.write_script("load_align_top.pml")
        run_pymol_script(root_out + "/" + "load_align_top.pml", parellel_process=True)
    print("Done")


