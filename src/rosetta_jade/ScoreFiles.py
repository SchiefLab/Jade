import json

from collections import OrderedDict

import pandas
import re

from pymol_jade.PyMolScriptWriter import *
from basic.pandas.PandasDataFrame import *
from basic.plotting.MakeFigure import *
from basic.path import *


##Original Author: Luki Goldschmidt <lugo@uw.edu>
class ScoreFile:
  def __init__(self, filename):
    self.filename = filename

    if re.search("score_", filename):
        self.name = ".".join(os.path.basename(self.filename).split('.')[0:-1]).replace('score_', "")
    else:
        self.name = ".".join(os.path.basename(self.filename).split('.')[0:-1]).replace('score', "")

    self.decoys = []; #Main holder of data.  Array of dict of scoreterms and 'decoy']

    self.decoy_field_name = "decoy"
    self.decoy_dir = os.path.dirname(filename); #Assume filename is in decoy dir.  This is not nessessarily the case...

    if filename == "-":
      lines = file(sys.stdin).readlines()
    else:
      lines = file(filename).readlines()

    for line in lines:
      try:
        o = json.loads(line)
        # print o[self.decoy_field_name]
        # print repr(o)
        if not re.search("initial_benchmark_perturbation", o[self.decoy_field_name]):
          self.decoys.append(o)
        #self.decoys.append(o)
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

    #print repr(decoy_names)
    # Collect data
    for decoy_name in scores:
      if not decoy_name in [ os.path.basename(x) for x in decoy_names]:
        continue
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
    return sorted([[x[scoreterm], x[self.decoy_field_name]] for x in self.decoys if x[self.decoy_field_name] in decoy_names and scoreterm in x ], reverse=reverse)[:top_n]

  def getScore(self, decoy, scoreterm):
    """
    Get Score of a particular decoy and scoreterm
    :param decoy: str
    :param scoreterm: str
    :rtype: float
    """
    for o in self.decoys:
      if o[self.decoy_field_name] == decoy:
        return o[scoreterm]

  def getScores(self, scoreterm, decoy_names=None, top_n=-1, reverse=False):
    ordered_decoy_names = self.getOrdered(scoreterm, decoy_names, top_n, reverse)
    return [self.getScore( decoy, scoreterm) for decoy in ordered_decoy_names]

  def getDataframe(self, scoreterms=None, order_by="total_score", top_n=-1, reverse=True):
    """
    Get data as a pandas dataframe.  Definitely preferred now.
    :param scoreterms: list
    :param order_by: str
    :param top_n: int
    :param reverse: bool
    :rtype: pandas.DataFrame
    """
    print self.name
    df = pandas.DataFrame.from_dict(self.decoys)
    df = detect_numeric(df)
    df = df.sort(order_by, ascending=reverse)[0:top_n]
    if scoreterms:
      df = get_columns(df, scoreterms)
    df.name = self.name
    df.scorefile = self.filename
    df["experiment"] = self.name
    df["scorefile"] = self.filename
    if os.path.dirname(self.filename):
      df["decoy_path"] = os.path.dirname(self.filename)+"/"+df["decoy"]
      df["decoy_path"] = [get_decoy_path(p) for p in df["decoy_path"]]
    else:
      df["decoy_path"] = df["decoy"]
      df["decoy_path"] = [get_decoy_path(p) for p in df["decoy_path"]]

    return df

def get_scorefiles(indir = os.getcwd()):
    """
    Get Score files from a directory.  Walk through all directories in directory.
    :param indir: str
    :rtype: list
    """
    matches = []
    for root, dirnames, filenames in os.walk(indir):
        for dirname in dirnames:
            matches.extend( glob.glob(root+"/"+dirname+"/*.sc") )

    return matches


def plot_score_vs_rmsd(df, title, outpath, score="total_score", rmsd="looprms", top_p=.95, reverse=True):
  """
  Plot a typical Score VS RMSD using matplotlib, save it somewhere.
  By default, plot the top 95%
  :param df: pandas.DataFrame
  :param outpath: str
  :param score: str
  :param rmsd: str
  :rtype: matplotlib.Axes
  """
  df = detect_numeric(df)
  df = df.sort(score, ascending=reverse)
  slice_top = df[0:int(len(df)*top_p)]
  ax = slice_top.plot(kind="scatter", y=score, x=rmsd, figsize=[11, 8], title = title)
  pad_single_title(ax)
  ax.set_axis_bgcolor('white')
  fig = ax.get_figure()
  fig.savefig(outpath, dpi=300)

  return ax


def pymol_session_on_top_df(df, outdir, decoy_dir,
                         scoreterm="total_score",
                         top_n=10,
                         decoy_column="decoy",
                         native_path = None,
                         out_prefix_override=None,
                         ab_structure = False,
                         superimpose = False,
                         run_pymol = True
                        ):
    """
    Make a PyMol session (or setup a scripter) on top X using a dataframe.
    Return the scripter for extra control.

    df should have an attribute of 'name' or out_prefix_override should be set.

    :param df: pandas.DataFrame
    :param outdir: str
    :param decoy_dir: str
    :param scoreterm: str
    :param top_n: int
    :param decoy_column: str
    :param native_path: str
    :param out_prefix_override: str
    :param ab_structure: boolean
    :param superimpose: boolean
    :rtype: PyMolScriptWriter
    """

    if not os.path.exists(outdir): os.mkdir(outdir)

    if out_prefix_override:
        session_prefix = out_prefix_override
    elif hasattr(df, "name"):
        session_prefix = df.name+"_"
    else:
        session_prefix = "pml_session_"
        print "WARNING: out_prefix_override not set and dataframe has no attribute name!"


    pymol_name  = session_prefix+scoreterm
    print pymol_name

    ascending = True
    if scoreterm == "hbonds_int" or scoreterm == "dSASA_int": reverse = False
    df2 = df.sort(scoreterm, ascending = ascending)
    df2 = df2[0:top_n]
    print df2['total_score'].tail()
    if decoy_dir:
        df2['decoy_path'] = decoy_dir+"/"+df2[decoy_column]
    else:
        df2['decoy_path'] = df2[decoy_column]

    top_decoys = zip(list(df2[scoreterm]), list( df2[ 'decoy_path' ]))
    #print repr(top_decoys)
    scripter = make_pymol_session_on_top_scored(top_decoys, decoy_dir, outdir, pymol_name, int(top_n), native_path,
                                     ab_structure, parellel=False, super = superimpose, run_pymol = run_pymol)

    return scripter