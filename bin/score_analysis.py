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
from pymol_gen.PyMolScriptWriter import *
from rosetta_gen.ScoreFiles import ScoreFiles

########################################################################

def output_legacy(out, prefix):
  # Format floats
  for i in range(0, len(out)):
    for j in range(0, len(out[i])):
      t = out[i][j]
      value = t[1]
      if type(value) == float:
        value = "%.03f" % value
        out[i][j] = (t[0], value)

  # Compute column widths
  widths = {}
  for decoy_terms in out:
    for decoy_term in decoy_terms:
      name, value = decoy_term
      if not widths.has_key(name):
        widths[name] = len(name)
      l = len(str(value))
      if widths[name] < l:
        widths[name] = l

  if widths.has_key('decoy'):
    widths['decoy'] = -widths['decoy']

  # Header
  sys.stdout.write(prefix)
  for decoy_term in out[0]:
    name, value = decoy_term
    sys.stdout.write("%*s " % (widths[name], name))
  sys.stdout.write("\n")

  # Scores
  for decoy_terms in out:
    sys.stdout.write(prefix)
    for decoy_term in decoy_terms:
      name, value = decoy_term
      sys.stdout.write("%*s " % (widths[name], value))
    sys.stdout.write("\n")


def output_tab(out, prefix):
  # Header
  names = [decoy_term[0] for decoy_term in out[0]]
  sys.stdout.write(prefix + "\t".join(names) + "\n")

  # Scores
  for decoy_terms in out:
    values = [str(decoy_term[1]) for decoy_term in decoy_terms]
    sys.stdout.write(prefix + "\t".join(values) + "\n")


def output_CSV(out, prefix):
  # Header
  names = [decoy_term[0] for decoy_term in out[0]]
  sys.stdout.write(prefix + ",".join(names) + "\n")

  # Scores
  for decoy_terms in out:
    values = []
    for decoy_term in decoy_terms:
      value = decoy_term[1]
      if type(value) in [str, unicode] and "," in value:
        value = '"%s"' % value
      values.append(str(value))
    sys.stdout.write(prefix + ",".join(values) + "\n")


def printVerbose(s):
  print >> sys.stderr, s


########################################################################

def main(argv):
  parser = ArgumentParser(
    "This utility parses and extracts data from score files in JSON format")

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

  parser.add_argument("-S", "--summary",
                      action="store_true",
                      default=False,
                      help="Compute stats summarizing data", )

  parser.add_argument("--list_scoretypes",
                      action="store_true",
                      default=False,
                      help="List score term names", )

  parser.add_argument("--make_pdblist",
                      default = False,
                      action = "store_true",
                      help = "Output PDBlist file(s)")

  parser.add_argument("--pdblist_prefix",
                      default = "",
                      help = "Prefix to use for PDBLIST outputs")

  parser.add_argument("--pdblist_outdir",
                      default = "pdblists",
                      help = "Output dir for pdblist files")

  pymol_opts = parser.add_argument_group("PyMol", "Options for pymol session output")

  pymol_opts.add_argument("--pymol_session",
                          help="Make pymol session(s) of the scoretypes specified",
                          default = False,
                          action = "store_true")

  pymol_opts.add_argument("--session_prefix",
                          default="",
                          type=str,
                          help="Prefix used for output pymol session")

  pymol_opts.add_argument("--session_outdir",
                          default = "sessions",
                          type=str,
                          help = "Output dir for pymol sessions.")

  pymol_opts.add_argument("--native",
                          type=str,
                          help="Native structure to use for pymol sessions.")

  pymol_opts.add_argument("--top_dir",
                          type=str,
                          help = "Top directory for PDBs if different than the directory of the scorefile")

  pymol_opts.add_argument("--ab_structure",
                          default = False,
                          action="store_true",
                          help = "Specify if the module is a renumbered antibody structure.  Will run pymol script for ab-specific selection")


  global options
  options = parser.parse_args()

  if options.decoy_names:
    options.decoy_names = [x.replace(".pdb.gz", "") for x in options.decoy_names]
    options.decoy_names = [x.replace(".pdb", "") for x in options.decoy_names]
    options.decoy_names = [x.replace("pre_model_1__", "pre_model_1_") for x in options.decoy_names]
    options.decoy_names = [os.path.basename( x ) for x in options.decoy_names]

  s = 0
  for filename in options.scorefiles:
    s +=1
    if filename != "":
      printVerbose("    Scorefile: %s" % filename)

    if filename != "-" and not os.path.isfile(filename):
      print >> sys.stderr, "File not found:", filename
      continue

    sf = ScoreFiles(filename)

    printVerbose("  File Decoys: %d" % sf.getDecoyCount())
    printVerbose("  Score terms: %s" % ", ".join(sf.getScoreTermNames()))
    printVerbose("")

    ### Info handlers
    if options.decoy_names:
      decoy_names = options.decoy_names
    else:
      decoy_names = sf.getDecoyNames()

    scoreterms = sf.getScoreTermNames()
    if options.list_scoretypes:
      print "\n".join(scoreterms)
      continue

    ### Stats summary
    if options.summary:
      stats = sf.getStats(options.scoretypes, decoy_names)
      max_width = max([len(x) for x in stats.keys()])
      fmt = "%*s:  %4s  %10s  %10s  %10s  %10s  %10s"
      print fmt % (max_width, "TERM", "n", "Min", "Max", "Mean", "Median", "StdDev")

      for column in stats:
        if column == "top_n_by_10": continue
        v = stats[column]
        for f in ['min', 'max', 'mean', 'median', 'stddev']:
          if not v[f] == None:
            v[f] = "%.3f" % v[f]
        print fmt % (max_width, column, v['n'], v['min'], v['max'], v['mean'], v['median'], v['stddev'])
      continue

    ### Default score list handler
    scores = sf.getScoreTerms(options.scoretypes)
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

      ordered = sf.getOrdered(term, decoy_names=decoy_names, top_n=options.top_n)

      if options.top_dir:
        top_decoy_paths = [get_decoy_path( options.top_dir+"/"+o[1]  ) for o in ordered]
      elif os.path.dirname(filename):
        top_decoy_paths = [get_decoy_path( os.path.dirname(filename)+"/"+o[1]  ) for o in ordered]
      else:
        top_decoy_paths = [get_decoy_path( o[1]  ) for o in ordered]

      for o in ordered:
        print "%.2f\t" % o[0] + o[1]

      if options.make_pdblist:
        if not os.path.exists(options.pdblist_dir):
          os.mkdir(options.pdblist_dir)

        outname = options.pdblist_dir+"/PDBLIST_"+options.pdblist_prefix+"_"+term+"_"+repr(s)+".txt"

        outfile = open(outname, 'w')
        for decoy in top_decoy_paths:

          outfile.write(decoy+"\n")
        print outname+" created"
        outfile.close()

    ### Top 10 by ten
    if "top_n_by_10" in options.scoretypes and options.top_n_by_10_scoretype in scoreterms:
      top_p = int(len(decoy_names) / 10)
      top_decoys = [o[1] for o in sf.getOrdered("total_score",  decoy_names=decoy_names, top_n=top_p)]
      top_by_n_decoys = [o for o in sf.getOrdered(options.top_n_by_10_scoretype, decoy_names=decoy_names) if o[1] in top_decoys][
                        :options.top_n_by_10]

      print "\n\nTop " + options.top_n_by_10_scoretype + " by top 10% Total Score"
      print options.top_n_by_10_scoretype + "\t" + "decoy" + "\t" + "total_score"

      for o in top_by_n_decoys:
        print "%.2f\t" % o[0] + o[1] + "%.2f" % sf.getScore(o[1], "total_score")

    if options.pymol_session:
      for scoreterm in options.scoretypes:

        if not scoreterm in scoreterms:
          print scoreterm +" specified for pymol session not present. Skipping"
          continue

        if scoreterm == "top_n_by_ten" and options.top_n_by_10_scoretype not in scoreterms:
          print "Top N by ten scoreterm not in scoreterms.  Please change the option"


        pymol_name  = options.session_prefix+""+repr(s)+"_"+scoreterm
        print "Creating: "+pymol_name

        outdir = options.session_outdir
        if not os.path.exists(outdir):
          os.mkdir(outdir)

        if not options.top_dir:
          options.top_dir = os.path.dirname(filename)
          if not options.top_dir:
            options.top_dir = os.getcwd()

        if "top_n_by_10" in options.scoretypes and options.top_n_by_10_scoretype in scoreterms:
          top_p = int(len(decoy_names) / 10)
          top_decoys = [o[1] for o in sf.getOrdered("total_score", decoy_names=decoy_names, top_n=int(options.top_n))]

          top_by_n_decoys = [[o[0], options.top_dir+"/"+o[1] ] for o in sf.getOrdered(options.top_n_by_10_scoretype, decoy_names=decoy_names) if o[1] in top_decoys][
                          :options.top_n_by_10]

          if len(top_by_n_decoys) == 0:
            print "No pdbs found. Skipping"
          make_pymol_session_on_top_scored(top_by_n_decoys, options.top_dir, outdir, pymol_name, int(options.top_n), options.native,
                                           antibody=options.ab_structure, parellel=False)

        else:
          ordered = sf.getOrdered(scoreterm, top_n=int(options.top_n), decoy_names=decoy_names)
          print repr(ordered)
          top_decoys = [[o[0], options.top_dir+"/"+o[1] ] for o in ordered]
          make_pymol_session_on_top_scored(top_decoys, options.top_dir, outdir, pymol_name, int(options.top_n), options.native,
                                           antibody=options.ab_structure, parellel=False)

########################################################################

if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
