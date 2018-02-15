#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Utility script to parse and extract data from score files in JSON format
# @author Luki Goldschmidt <lugo@uw.edu>
# @author Jared Adolf-Bryfogle - Forking.

from argparse import ArgumentParser
from jade.pymol_jade.PyMolScriptWriter import *
from jade.rosetta_jade.ScoreFiles import ScoreFile
from jade.basic.plotting.MakeFigure import *

import shutil

########################################################################

def printVerbose(s):
    print >> sys.stderr, s

def get_parser():

    parser = ArgumentParser(
        description="This utility parses and extracts data from score files in JSON format")

    parser.add_argument("scorefiles", nargs='*', help="A list of scorefiles")

    parser.add_argument("-s", "--scoretypes",
                        default=["dSASA_int", "delta_unsatHbonds", "hbonds_int", "total_score", "dG_separated", "top_n_by_10"],
                        help="List of score terms to extract",
                        nargs='*')

    parser.add_argument("-n", "--top_n",
                        default=-1,
                        type=int,
                        help="Only list Top N when doing top scoring decoys or making pymol sessions"
                             "Default is to print all of them.")

    parser.add_argument('--top_n_by_10',
                        default=10,
                        type=int,
                        help="Top N by 10 percent total score to print out. ")

    parser.add_argument('--top_n_by_10_scoretype',
                        default="dG_separated",
                        help="Scoretype to use for any top N by 10 printing.  If scoretype not present, won't do anything.")

    parser.add_argument("--decoy_names",
                        help="Decoy names to use",
                        default=[],
                        nargs='*')

    parser.add_argument("--list_scoretypes",
                        action="store_true",
                        default=False,
                        help="List score term names")

    parser.add_argument("--pdb_dir", "-d",
                        type=str,
                        help="Directory for PDBs if different than the directory of the scorefile")

    output_opts = parser.add_argument_group("OUTPUT", "General output options.")

    output_opts.add_argument("--summary", "-S",
                             action="store_true",
                             default=False,
                             help="Compute stats summarizing data")

    output_opts.add_argument("--csv", "-c",
                             action="store_true",
                             default=False,
                             help="Output selected columns, top, and decoys as CSV.")

    output_opts.add_argument("--make_pdblist",
                             default=False,
                             action="store_true",
                             help="Output PDBlist file(s)")

    output_opts.add_argument("--pymol_session",
                             help="Make pymol session(s) of the scoretypes specified",
                             default=False,
                             action="store_true")

    output_opts.add_argument("--plot",
                             help="Plot one score type vs another.  Save the plot.  2 or 3 Arguments.  [X, Y, 'Title''] OR [X, 'Title'].  "
                                  "If title has spaces, use quotes. "
                                  "Nothing special, just used for quick info.",

                             nargs='*',
                             default=[])

    output_opts.add_argument("--copy_top_models",
                             help = "Copy the top -n to the output directory for each scorefile passed.",
                             default = False,
                             action = "store_true")


    output_opts.add_argument("--prefix", "-p",
                             default="",
                             help="Prefix to use for any file output. Do not include any _")

    output_opts.add_argument("--outdir", "-o",
                             default=os.getcwd(),
                             type=str,
                             help="Output dir.  Default is current directory.")

    plotting_opts = parser.add_argument_group("PLOTTING", "Options for plot output")

    plotting_opts.add_argument("--plot_type",
                             default = "scatter",
                             choices = ["line", "scatter", "bar", "hist", "box", "kde", "area","pie", "hexbin"],
                             help = "The type of plot we are outputting.")

    plotting_opts.add_argument("--plot_filter",
                             default = 1.0,
                             help = "Filter X to top Percent of this - useful to remove outliers.")

    pymol_opts = parser.add_argument_group("PYMOL", "Options for pymol session output")

    pymol_opts.add_argument("--native",
                            type=str,
                            help="Native structure to use for pymol sessions.")

    pymol_opts.add_argument("--ab_structure",
                            default=False,
                            action="store_true",
                            help="Specify if the module is a renumbered antibody structure.  Will run pymol script for ab-specific selection")

    pymol_opts.add_argument("--super",
                            help="Super this selection instead of align all to.")

    return parser
########################################################################

def main():

    parser = get_parser()
    global options
    options = parser.parse_args()

    if options.decoy_names:
        options.decoy_names = [x.replace(".pdb.gz", "") for x in options.decoy_names]
        options.decoy_names = [x.replace(".pdb", "") for x in options.decoy_names]
        options.decoy_names = [x.replace("pre_model_1__", "pre_model_1_") for x in options.decoy_names]
        options.decoy_names = [os.path.basename(x) for x in options.decoy_names]




    s = 0
    print "Decoy Names:"+ str(len(options.decoy_names))
    for filename in options.scorefiles:
        s += 1
        if filename != "":
            printVerbose("    Scorefile: %s" % filename)

        if filename != "-" and not os.path.isfile(filename):
            print >> sys.stderr, "File not found:", filename
            continue



        sf = ScoreFile(filename)

        # Update Decoys
        if options.decoy_names:
            decoy_dict = defaultdict()
            new_decoy_list = []
            #print sf.decoys

            #Pre-organize data
            for r in sf.decoys:
                decoy_dict[r[sf.decoy_field_name]] = r


            for name in options.decoy_names:
                new_decoy_list.append(decoy_dict[name])

            sf.decoys = new_decoy_list

        df = sf.get_Dataframe()

        printVerbose("  File Decoys: %d" % sf.get_decoy_count())
        printVerbose("  Score terms: %s" % ", ".join(sf.get_scoreterm_names()))
        printVerbose("")

        #Optionally make the ouptut directory.
        if options.outdir and not os.path.exists(options.outdir):
            os.mkdir(options.outdir)

        if not options.pdb_dir:
            pdb_dir = os.path.dirname(filename)
            if not pdb_dir:
                pdb_dir = os.getcwd()

        else:
            pdb_dir = options.pdb_dir

        ### Info handlers

        if options.decoy_names:
            decoy_names = options.decoy_names
        else:
            decoy_names = sf.get_decoy_names()

        scoreterms = sf.get_scoreterm_names()
        if options.list_scoretypes:
            print "\n".join(scoreterms)
            continue

        ### Stats summary
        if options.summary:
            stats = sf.get_stats(options.scoretypes, decoy_names)
            max_width = max([len(x) for x in stats.keys()])
            fmt = "%*s:  %4s  %10s  %10s  %10s  %10s  %10s"
            print "SUMMARY:  " +options.prefix+ "  "+ fmt % (max_width, "TERM", "n", "Min", "Max", "Mean", "Median", "StdDev")

            for column in stats:
                if column == "top_n_by_10": continue
                v = stats[column]
                for f in ['min', 'max', 'mean', 'median', 'stddev']:
                    if not v[f] == None:
                        v[f] = "%.3f" % v[f]
                print "SUMMARY:  " +options.prefix+ "  "+ fmt % (max_width, column, v['n'], v['min'], v['max'], v['mean'], v['median'], v['stddev'])
            continue

        ### Default score list handler
        scores = sf.get_scoreterms(options.scoretypes)
        out = []
        for decoy_name in scores:
            if not decoy_name in decoy_names: continue
            terms = scores[decoy_name]
            decoy_scores = [("decoy", decoy_name)]
            for term in terms:
                decoy_scores.append((term, terms[term]))
            out.append(decoy_scores)

        for term in options.scoretypes:
            if term not in scoreterms: continue
            print "\nBy " + term

            ordered = sf.get_ordered_decoy_list(term, decoy_names=decoy_names, top_n=options.top_n)

            if options.pdb_dir:
                top_decoy_paths = [get_decoy_path(options.pdb_dir + "/" + o[1]) for o in ordered]
            elif os.path.dirname(filename):
                top_decoy_paths = [get_decoy_path(os.path.dirname(filename) + "/" + o[1]) for o in ordered]
            else:
                top_decoy_paths = [get_decoy_path(o[1]) for o in ordered]

            for o in ordered:
                print "%.2f\t" % o[0] + o[1]

            if options.make_pdblist:

                outname = options.outdir + "/PDBLIST_" + options.prefix + term + "_" + repr(s) + ".txt"

                outfile = open(outname, 'w')
                for decoy in top_decoy_paths:
                    outfile.write(decoy + "\n")
                print outname + " created"
                outfile.close()

        ### Top 10 by ten
        if "top_n_by_10" in options.scoretypes and options.top_n_by_10_scoretype in scoreterms:
            top_p = int(len(decoy_names) / 10)
            top_decoys = [o[1] for o in sf.get_ordered_decoy_list("total_score", decoy_names=decoy_names, top_n=top_p)]
            top_by_n_decoys = [o for o in sf.get_ordered_decoy_list(options.top_n_by_10_scoretype, decoy_names=decoy_names) if
                               o[1] in top_decoys][
                              :options.top_n_by_10]

            print "\n\nTop " + options.top_n_by_10_scoretype + " by top 10% Total Score"
            print options.top_n_by_10_scoretype + "\t" + "decoy" + "\t" + "total_score"

            for o in top_by_n_decoys:
                print "%.2f\t" % o[0] + o[1] + "%.2f" % sf.get_score(o[1], "total_score")

        if options.pymol_session:
            print "Making PyMol Session "
            for scoreterm in options.scoretypes:
                print "Scoreterm "+scoreterm
                if scoreterm == "top_n_by_10" and not options.top_n_by_10_scoretype in scoreterms:
                    print "Top N by ten scoreterm not in scoreterms.  Please change the option"
                    continue
                elif scoreterm == "top_n_by_10" and options.top_n_by_10_scoretype in scoreterms:
                    pass
                elif scoreterm in scoreterms:
                    pass
                else:
                    print scoreterm + " specified for pymol session not present. Skipping"
                    continue


                pymol_name = options.prefix + "" + sf.name + scoreterm
                print "Creating: " + pymol_name

                outdir = options.outdir



                print "PDB DIR: " + pdb_dir

                if scoreterm == "top_n_by_10" and "top_n_by_10" in options.scoretypes and options.top_n_by_10_scoretype in scoreterms:
                    top_p = int(len(decoy_names) / 10)
                    top_decoys = [o[1] for o in
                                  sf.get_ordered_decoy_list("total_score", decoy_names=decoy_names, top_n=int(options.top_n))]

                    top_by_n_decoys = [[o[0], pdb_dir + "/" + o[1]] for o in
                                       sf.get_ordered_decoy_list(options.top_n_by_10_scoretype, decoy_names=decoy_names) if
                                       o[1] in top_decoys][
                                      :options.top_n_by_10]


                    if len(top_by_n_decoys) == 0:
                        print "No pdbs found. Skipping"
                    make_pymol_session_on_top_scored(top_by_n_decoys, pdb_dir, outdir, pymol_name,
                                                     int(options.top_n), options.native,
                                                     antibody=options.ab_structure, parellel=False, super=options.super)

                else:
                    ordered = sf.get_ordered_decoy_list(scoreterm, top_n=int(options.top_n), decoy_names=decoy_names)

                    top_decoys = [[o[0], pdb_dir + "/" + o[1]] for o in ordered]
                    #print repr(top_decoys)
                    make_pymol_session_on_top_scored(top_decoys, pdb_dir, outdir, pymol_name,
                                                     int(options.top_n), options.native,
                                                     antibody=options.ab_structure, parellel=False, super=options.super)

        if options.plot:
            if len(options.plot) == 3 and options.plot_type == "scatter":

                x = options.plot[0]
                y = options.plot[1]
                title = options.plot[2]
                outpath = options.outdir+"/"+options.prefix+options.plot_type+"_"+x+"_"+y+"_"+sf.name+".pdf"
                plot_x_vs_y_sea_with_regression(df, title, outpath, x, y, top_p=options.plot_filter)
            else:
                y = None
                z = None
                if len(options.plot) == 3:
                    x = options.plot[0]
                    y = options.plot[1]
                    title = options.plot[2]
                    xy = x+"_"+y
                elif len(options.plot) == 2:
                    x = options.plot[0]
                    title = options.plot[1]
                    xy = x
                else:
                    sys.exit("Plot option must either be a length of 2 or 3.  See help for more information!")

                outpath = options.outdir+"/"+options.prefix+options.plot_type+"_"+xy+"_"+sf.name+".pdf"
                plot_general_pandas(df, title, outpath, options.plot_type, x, y = y, z = z, top_p=options.plot_filter)

            os.system("open "+outpath)

        if options.copy_top_models:
            for scoreterm in options.scoretypes:
                if scoreterm == "top_n_by_10" and "top_n_by_10" in options.scoretypes and options.top_n_by_10_scoretype in scoreterms:
                    top_p = int(len(decoy_names) / 10)
                    top_decoys = [o[1] for o in
                                  sf.get_ordered_decoy_list("total_score", decoy_names=decoy_names, top_n=int(options.top_n))]
                    top_by_n_decoys = [[o[0], pdb_dir + "/" + o[1]] for o in
                                       sf.get_ordered_decoy_list(options.top_n_by_10_scoretype, decoy_names=decoy_names) if
                                       o[1] in top_decoys][
                                      :options.top_n_by_10]


                    if len(top_by_n_decoys) == 0:
                        print "No pdbs found. Skipping"


                    for i in range(0, options.top_n):
                        print get_decoy_path(top_by_n_decoys[i][1])
                        shutil.copy(get_decoy_path(top_by_n_decoys[i][1]), options.outdir)
                else:
                    ordered = sf.get_ordered_decoy_list(scoreterm, top_n=int(options.top_n), decoy_names=decoy_names)
                    top_decoys = [[o[0], pdb_dir + "/" + o[1]] for o in ordered]

                    for i in range(0, options.top_n):
                        print get_decoy_path(top_decoys[i][1])
                        shutil.copy(get_decoy_path(top_decoys[i][1]), options.outdir)

########################################################################

if __name__ == "__main__":
    sys.exit(main())
