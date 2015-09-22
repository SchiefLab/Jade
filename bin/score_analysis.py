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

import sys
import os
import json
from collections import OrderedDict
from argparse import ArgumentParser
from pymol.PyMolScriptWriter import *

class ScoreFiles:
  def __init__(self, filename):
    self.filename = filename
    self.decoys = []
    self.decoy_field_name = "decoy"

    if filename == "-":
      lines = file(sys.stdin).readlines()
    else:
      lines = file(filename).readlines()

    for line in lines:
      try:
        o = json.loads(line)
        # print o[self.decoy_field_name]
        # print repr(o)
        self.decoys.append(o)
      except ValueError:
        print >> sys.stderr, "Failed to parse JSON object; skipping line:\n", line

  def getDecoyCount(self):
    return len(self.decoys)

  def getDecoyNames(self):
    return [str(r[self.decoy_field_name]) for r in self.decoys]

  def getScoreTermNames(self):
    r = []
    for rec in self.decoys:
      for k in rec.keys():
        if k != self.decoy_field_name and not k in r:
          r.append(str(k))
    r.sort()
    return r

  def getScoreTerms(self, scoreterms=""):
    if type(scoreterms) == str:
      if scoreterms in ["", "*"]:
        scoreterms = self.getScoreTermNames()
        scoreterms.sort()
      else:
        scoreterms = scoreterms.split(",")
    r = OrderedDict()
    for rec in self.decoys:
      scores = {}
      for scoreterm in scoreterms:
        scores[scoreterm] = rec.get(scoreterm)
      r[str(rec[self.decoy_field_name])] = scores
    return r

  def getScoreTerm(self, scoreterm):
    r = {}
    for rec in self.decoys:
      r[str(rec[self.decoy_field_name])] = rec.get(scoreterm)
    return r

  def getStats(self, scoreterms="", decoy_names = None):

    print repr(len(decoy_names))
    scores = self.getScoreTerms(scoreterms)
    stats = None

    if not decoy_names:
      decoy_names = self.getDecoyNames()

    # Collect data
    for decoy_name in scores:
      if decoy_name in decoy_names: continue
      decoy = scores[decoy_name]

      if stats == None:
        columns = decoy.keys()
        stats = {k: [] for k in columns}

      for term in decoy:
        score = decoy[term]
        if not score is None and type(score) == float:
          stats[term].append(score)

    # Compute Stats
    calc_stats = OrderedDict()

    def median(s):
      s.sort()
      return s[int(len(s) / 2)]

    def stddev(s):
      mean = sum(s) / len(s)
      r = 0.0
      for x in s:
        r += pow(x - mean, 2)
      return pow(r / len(s), 0.5)

    for column in columns:
      s = stats[column]
      calc_stats[column] = {
        'n': len(s),
        'mean': sum(s) / len(s) if len(s) > 0 else None,
        'median': median(s) if len(s) > 0 else None,
        'stddev': stddev(s) if len(s) > 0 else None,
        'min': min(s) if len(s) > 0 else None,
        'max': max(s) if len(s) > 0 else None
      }

    return calc_stats

  def getOrdered(self, scoreterm, decoy_names = None, top_n=-1, reverse=False, ):
    """
    Get an ordered tuple of [[score, decoy_name], ...]
    Will automatically order some known scoreterms (hbonds_int, dSASA_int)
    """

    if not decoy_names:
      decoy_names = self.getDecoyNames()

    if scoreterm == "hbonds_int" or scoreterm == "dSASA_int": reverse = True
    return sorted([[x[scoreterm], x[self.decoy_field_name]] for x in self.decoys if x[self.decoy_field_name] in decoy_names], reverse=reverse)[:top_n]

  def getScore(self, decoy, scoreterm):
    for o in self.decoys:
      if o[self.decoy_field_name] == decoy:
        return o[scoreterm]


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


  pymol_opts = parser.add_argument_group("PyMol", "Options for pymol session output")

  pymol_opts.add_argument("--pymol_session",
                          help="Make pymol session(s) of the scoretypes specified",
                          default = False,
                          action = "store_true")

  pymol_opts.add_argument("--session_prefix",
                          default="",
                          help="Prefix used for output pymol session")

  pymol_opts.add_argument("--session_outdir",
                          default = "sessions",
                          help = "Output dir for pymol sessions.")

  pymol_opts.add_argument("--native",
                          help="Native structure to use for pymol sessions.")

  pymol_opts.add_argument("--top_dir",
                          help = "Top directory for PDBs if different than the directory of the scorefile")

  global options
  options = parser.parse_args()

  if options.decoy_names:
    options.decoy_names = [x.replace(".pdb.gz", "") for x in options.decoy_names]
    options.decoy_names = [x.replace(".pdb", "") for x in options.decoy_names]


  s = 0
  for filename in options.scorefiles:
    s +=1
    if filename != "":
      printVerbose("    Scorefile: %s" % filename)

    if filename != "-" and not os.path.isfile(filename):
      print >> sys.stderr, "File not found:", filename
      continue

    sf = ScoreFiles(filename)

    printVerbose("       Decoys: %d" % sf.getDecoyCount())
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

      for o in ordered:
        print "%.2f\t" % o[0] + o[1]

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
      scoreterms.append("top_n_by_10")
      for scoreterm in options.scoretypes:

        if not scoreterm in scoreterms:
          sys.exit("Scoreterm specified for pymol session not present.")

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

        if scoreterm == "top_n_by_10":
          top_p = int(len(decoy_names) / 10)
          top_decoys = [o[1] for o in sf.getOrdered("total_score", decoy_names=decoy_names, top_n=int(options.top_n))]
          top_by_n_decoys = [[o[0], options.top_dir+"/"+o[1] ] for o in sf.getOrdered(options.top_n_by_10_scoretype, decoy_names=decoy_names) if o[1] in top_decoys][
                          :options.top_n_by_10]

          make_pymol_session_on_top_scored(top_by_n_decoys, options.top_dir, outdir, pymol_name, int(options.top_n), options.native)

        else:
          ordered = sf.getOrdered(scoreterm, top_n=int(options.top_n), decoy_names=decoy_names)
          print repr(ordered)
          top_decoys = [[o[0], options.top_dir+"/"+o[1] ] for o in ordered]
          make_pymol_session_on_top_scored(top_decoys, options.top_dir, outdir, pymol_name, int(options.top_n), options.native)

########################################################################

if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
