#!/usr/bin/env python
from __future__ import print_function
import os,sys
from jade.basic.path import *
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

def order_by_position(d):
    d.reset_index(inplace=True)
    d.set_index(['selection', 'decoy_name'], inplace=True)
    d.sort_index(inplace=True)
    d.reset_index(inplace=True)


def draw_boxplot_numerical_categorical_cols(df, numerical_col, categorical_col):
    categories = []
    for category in df[categorical_col].unique():
        categories.append({
            'categorical_col': category,
            'numerical_col': df[df[categorical_col] == category][numerical_col]})

    mpl.pyplot.figure(figsize=(11, 8))
    mpl.pyplot.boxplot([e['numerical_col'] for e in categories], labels=[e['categorical_col'] for e in categories],
                       whis=1000)
    mpl.pyplot.xticks(rotation='vertical')
    return mpl.pyplot.gcf()

def get_sds(d, column):
    mean = d[column].mean()
    sd = d[column].std()
    print(column)
    print(repr(mean))
    print(repr(sd))
    #print repr(d[column])
    return (d[column] - mean)/sd

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
    df['decoy_name'] = df["decoy"].apply(get_decoy_name)
    df['position'] = df['decoy_name'].apply(get_position)
    df['selection'] = pd.to_numeric(df['selection'], downcast='unsigned')

    # Here, we calculate diffs.  In order to reduce the amount of changes that need to be done downstream,
    # We are going to keep them the same name.
    df['post-sequon_total_energy_'] = df['post-sequon_total_energy'] - df['post-prerefine_total_energy']
    df['post-model_total_energy_'] = df['post-model_total_energy'] - df['post-sequon_total_energy']

    df['post-sequon_total_energy'] = df['post-sequon_total_energy_']
    df['post-model_total_energy'] = df['post-model_total_energy_']

    df_origin = df

    if options.get_top:

        ## Post-Sequon Total Energy
        df.reset_index().set_index(['position', 'post-sequon_total_energy'], inplace=True)
        df.sort_index().reset_index(inplace=True)

        top_scoring = df.groupby('position').head(1)
        top_scoring.reset_index(inplace=True)
        top_scoring.set_index('post-sequon_total_energy', inplace=True)
        top_scoring.sort_index(inplace=True)
        top_scoring.reset_index(inplace=True)

        paths_top = top_scoring['decoy']
        scores_top = top_scoring['post-sequon_total_energy']
        names_top = top_scoring['position']
        score_paths_top = zip(scores_top, paths_top, names_top)

        out_name = "all_post-sequon_top_model_per_position"
        maindir = options.outdir_prefix+"_models"
        outdir = maindir+"/"+out_name
        if not os.path.exists(maindir):
            os.mkdir(maindir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        make_pymol_session_on_top_scored(score_paths_top, outdir, outdir, out_name, antibody=False, native_path=options.native_path)

        for i, tup in enumerate(score_paths_top):
            print(i)
            print(tup)
            name = "model_" + str(i) + "_" + str(tup[2]) + "_%.2f" % (tup[0]) + get_decoy_extension(get_decoy_path(str(tup[1])))
            shutil.copy(get_decoy_path(str(tup[1])), outdir+'/'+name)

        ## Post-Glycan Total Energy
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

        out_name = "all_post-glycan_top_model_per_position"
        maindir = options.outdir_prefix+"_models"
        outdir = maindir+"/"+out_name
        if not os.path.exists(maindir):
            os.mkdir(maindir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        make_pymol_session_on_top_scored(score_paths_top, outdir, outdir, out_name, antibody=False, native_path=options.native_path)

        #Copy Top Models
        for i, tup in enumerate(score_paths_top):
            print(i)
            print(tup)
            name = "model_" + str(i) + "_" + str(tup[2]) + "_%.2f" % (tup[0]) + get_decoy_extension(get_decoy_path(str(tup[1])))
            shutil.copy(get_decoy_path(str(tup[1])), outdir+'/'+name)



        ## Top 15 Scoring Per-position
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
                organized_names[name].append(score_path)

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
                name = "model_" + str(i) + "_" + str(tup[2]) + "_%.2f" % (tup[0]) + get_decoy_extension(get_decoy_path(str(tup[1])))
                shutil.copy(get_decoy_path(tup[1]), outdir+'/'+name)

        scripter = make_pymol_session_on_top_scored(score_paths, root_out, root_out, "all_top_15", antibody=False,
                                         native_path=options.native_path, run_pymol=False)

        #Group all of them.
        for tup in score_paths_top:
            name = str(tup[2])
            scripter.add_line("group models_"+name+"_%.2f" % (tup[0])+', models_'+name+'*' )

        pse_path = root_out + "/" + out_name + ".pse"
        scripter.add_save_session(pse_path)
        scripter.write_script("load_align_top.pml")
        run_pymol_script(root_out + "/" + "load_align_top.pml", parellel_process=True)

    if options.get_plots:
        df = df_origin

        # Order by position
        df.reset_index(inplace=True)
        df.set_index(['selection', 'decoy_name'], inplace=True)
        df.sort_index(inplace=True)

        # Get stats and output just to have numerical values
        stats = df.groupby(['selection', 'position']).describe()
        stats.to_csv("per_position_stats.csv")

        stats.reset_index(inplace=True)

        outdir = options.outdir_prefix+"_plots"

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-sequon_total_energy', 'mean'), figsize=[11, 8],
                        title="Post-Sequon Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_sequon.pdf", dpi=300)

        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-sequon_total_energy', 'min'), figsize=[11, 8],
                        title="Post-Sequon Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (min)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_sequon_min.pdf", dpi=300)

        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-model_total_energy', 'mean'), figsize=[11, 8],
                        title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_means.pdf", dpi=300)

        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-model_total_energy', 'min'), figsize=[11, 8],
                        title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (min)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_min.pdf", dpi=300)

        stats.sort_values([('post-sequon_total_energy', 'mean')], inplace=True)
        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-sequon_total_energy', 'mean'), figsize=[11, 8],
                        title="Post-Sequon Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_sequon_ordered.pdf", dpi=300)

        stats.sort_values([('post-model_total_energy', 'mean')], inplace=True)
        ax = stats.plot(yerr=stats, kind="bar", x='position', y=('post-model_total_energy', 'mean'), figsize=[11, 8],
                        title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_means_ordered.pdf", dpi=300)

        stats.sort_values([('post-model_total_energy', 'min')], inplace=True)
        ax = stats.plot(kind="bar", x='position', y=('post-model_total_energy', 'min'), figsize=[11, 8],
                        title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (min)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_min_ordered.pdf", dpi=300)

        cmap = cm.get_cmap('inferno')

        ax = stats.plot(kind="scatter", sharex=False, c='selection', colormap=cmap, s=120,
                        x=('post-sequon_total_energy', 'mean'), y=('post-model_total_energy', 'min'), figsize=[11, 8],
                        title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'Glycan $\Delta$ REU (min)')
        ax.set_xlabel(r'Sequon $\Delta$ REU (mean)')
        # ax.legend().set_visible(False)

        fig = ax.get_figure()
        fig.savefig(outdir + "/sequon_vs_glycan_min.pdf", dpi=300)

        df.reset_index(inplace=True)
        df.set_index(['selection', 'post-model_total_energy'], inplace=True)
        df.sort_index(inplace=True)
        df.reset_index(inplace=True)

        df.reset_index(inplace=True)
        df.set_index('selection', inplace=True)
        L = df.groupby('selection').head(10)
        L.reset_index(inplace=True)
        L_stats = L.groupby(['selection', 'position']).describe()
        L_stats.reset_index(inplace=True)
        #L_stats['post-model_total_energy'].tail()

        ax = L_stats.plot(kind="bar", x='position', y=('post-model_total_energy', 'mean'), figsize=[11, 8],
                          title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean_top_10)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_means_top_10.pdf", dpi=300)

        L_stats.sort_values([('post-model_total_energy', 'mean')], inplace=True)
        ax = L_stats.plot(kind="bar", x='position', y=('post-model_total_energy', 'mean'), figsize=[11, 8],
                          title="Post-Glycan Total Energy Change")
        ax.set_ylabel(r'$\Delta$ REU (mean_top_10)')
        ax.set_xlabel('sequon position')
        ax.legend().set_visible(False)
        fig = ax.get_figure()
        fig.savefig(outdir + "/per_position_glycan_means_top_10_ordered.pdf", dpi=300)

        df.reset_index(inplace=True)
        df.set_index('selection', inplace=True)
        L2 = df.groupby('selection').head(int(2500 * .10))

        L2.reset_index(inplace=True)
        L2.set_index('selection', inplace=True)
        L2.sort_index(inplace=True)
        L2.reset_index(inplace=True)

        fig = draw_boxplot_numerical_categorical_cols(L2, 'post-model_total_energy', 'position')
        ax = fig.axes[0]
        ax.set_title("Post-Glycan Total Energy Spread (Top 10%)")
        ax.set_ylabel(r'$\Delta$ REU')
        ax.set_xlabel('sequon position')
        ax.set_ylim([-5, 75])
        ax.legend().set_visible(False)
        fig = ax.get_figure()

        fig.savefig(outdir + "/boxplot_top-10-percent.pdf", dpi=300)

        df.reset_index(inplace=True)
        df.set_index('selection', inplace=True)
        L2 = df.groupby('selection').head(int(2500 * .25))

        L2.reset_index(inplace=True)
        L2.set_index('selection', inplace=True)
        L2.sort_index(inplace=True)
        L2.reset_index(inplace=True)

        fig = draw_boxplot_numerical_categorical_cols(L2, 'post-model_total_energy', 'position')
        ax = fig.axes[0]
        ax.set_title("Post-Glycan Total Energy Spread (Top 25%)")
        ax.set_ylabel(r'$\Delta$ REU')
        ax.set_xlabel('sequon position')
        ax.set_ylim(-5, 80)
        ax.legend().set_visible(False)
        fig = ax.get_figure()

        fig.savefig(outdir + "/boxplot_top-25-percent.pdf", dpi=300)

        df.reset_index(inplace=True)
        df.set_index('selection', inplace=True)
        L2 = df.groupby('selection').head(.01 * 2500)

        L2.reset_index(inplace=True)
        L2.set_index('selection', inplace=True)
        L2.sort_index(inplace=True)
        L2.reset_index(inplace=True)

        fig = draw_boxplot_numerical_categorical_cols(L2, 'post-model_total_energy', 'position')
        ax = fig.axes[0]
        ax.set_title("Post-Glycan Total Energy Spread (Top 1 Percent)")
        ax.set_ylabel(r'$\Delta$ REU')
        ax.set_xlabel('sequon position')
        ax.set_ylim([-5, 75])
        ax.legend().set_visible(False)
        fig = ax.get_figure()

        fig.savefig(outdir + "/boxplot_top-1-percent.pdf", dpi=300)

        # Composite Ranking - Basic

        # Post - Sequon Total Energy
        stats2 = stats.set_index(('post-sequon_total_energy', 'min')).sort_index().reset_index()
        stats2['post_sequon_min_rank'] = stats2.index

        # Post - Glycan Total Energy
        stats2 = stats2.set_index(('post-model_total_energy', 'min')).sort_index().reset_index()
        stats2['post_glycan_min_rank'] = stats2.index

        # Post - Glycan Total Energy (Top 1% Means)
        top = df.reset_index().set_index(['selection', 'post-model_total_energy']).sort_index().groupby(
            'selection').head(25).reset_index()
        test2 = top.groupby('selection').describe()

        stats2 = stats2.set_index('selection').sort_index()
        test2 = test2.sort_index()

        stats2['post-model_total_energy_top_mean'] = test2[('post-model_total_energy', 'mean')].astype(float)

        stats2 = stats2.reset_index().set_index('post-model_total_energy_top_mean').sort_index().reset_index()
        stats2['post_glycan_top-mean_rank'] = stats2.index

        # Rank 1 to N
        stats2['post_sequon_min_rank'] = stats2['post_sequon_min_rank'] + 1
        stats2['post_glycan_min_rank'] = stats2['post_glycan_min_rank'] + 1
        stats2['post_glycan_top-mean_rank'] = stats2['post_glycan_top-mean_rank'] + 1

        stats2['composite_rank'] = (stats2['post_sequon_min_rank'] * (1 / 3.0)) + (
        stats2['post_glycan_top-mean_rank'] * (1 / 3.0)) + (stats2['post_glycan_min_rank'] * (1 / 3.0))

        stats2.set_index('selection', inplace=True)
        stats2.sort_index(inplace=True)
        stats2.reset_index(inplace=True)

        ax = stats2.plot(kind="bar", x='position',
                         y=['post_sequon_min_rank', 'post_glycan_top-mean_rank', 'post_glycan_min_rank'],
                         figsize=[11, 8],
                         title="Composite Rankings")

        ax.set_ylabel('Rank')
        ax.set_xlabel('sequon position')
        ax.set_ylim([0, 57])
        ax.legend(['Post sequon min', 'Post glycan top-mean', 'Post glycan min'], loc="upper right")
        fig = ax.get_figure()
        fig.savefig(outdir + "/composite_ranking_numeric_rank.pdf", dpi=300)

        stats2.set_index('composite_rank', inplace=True)
        stats2.sort_index(inplace=True)
        stats2.reset_index(inplace=True)

        ax = stats2.plot(kind="bar", x='position',
                         y=['post_sequon_min_rank', 'post_glycan_top-mean_rank', 'post_glycan_min_rank',
                            'composite_rank'],
                         figsize=[11, 8],
                         title="Composite Rankings")

        ax.set_ylabel('Rank')
        ax.set_xlabel('sequon position')
        ax.set_ylim([0, 57])
        ax.legend(['Post sequon min', 'Post glycan top-mean', 'Post glycan min', 'Composite', ], loc="upper right")
        fig = ax.get_figure()
        fig.savefig(outdir + "/composite_ranking_numeric_rank_ordered_by_composite.pdf", dpi=300)

        stats2['post_sequon_min_sd'] = get_sds(stats2, ('post-sequon_total_energy', 'min'))
        stats2['post_glycan_min_sd'] = get_sds(stats2, ('post-model_total_energy', 'min'))
        stats2['post_glycan_top-mean_sd'] = get_sds(stats2, 'post-model_total_energy_top_mean')

        stats2['composite_sd'] = (stats2['post_sequon_min_sd'] * (1 / 3.0)) + (
        stats2['post_glycan_min_sd'] * (1 / 3.0)) + (stats2['post_glycan_top-mean_sd'] * (1 / 3.0))

        stats2.set_index('selection', inplace=True)
        stats2.sort_index(inplace=True)
        stats2.reset_index(inplace=True)

        ax = stats2.plot(kind="bar", x='position',
                         y=['post_sequon_min_sd', 'post_glycan_min_sd', 'post_glycan_top-mean_sd'],
                         figsize=[11, 8],
                         title="Composite Rankings")

        ax.set_ylabel('Number of Standard Deviations')
        ax.set_xlabel('sequon position')
        ax.set_ylim([-2.5, 2.5])
        ax.legend(['Post sequon min', 'Post glycan top-mean', 'Post glycan min'], loc="upper right")
        fig = ax.get_figure()
        fig.savefig(outdir + "/composite_ranking_sd_rank.pdf", dpi=300)

        stats2.set_index('composite_sd', inplace=True)
        stats2.sort_index(inplace=True)
        stats2.reset_index(inplace=True)

        ax = stats2.plot(kind="bar", x='position',
                         y=['post_sequon_min_sd', 'post_glycan_min_sd', 'post_glycan_top-mean_sd', 'composite_sd'],
                         figsize=[11, 8],
                         title="Composite Rankings")

        ax.set_ylabel('Number of Standard Deviations')
        ax.set_xlabel('sequon position')
        ax.set_ylim([-2.5, 2.5])
        ax.legend(['Post sequon min', 'Post glycan top-mean', 'Post glycan min', 'Composite'], loc="upper right")
        fig = ax.get_figure()
        fig.savefig(outdir + "/composite_ranking_sd_rank_ordered_by_composite.pdf", dpi=300)


    print("Done")


