import sys
import json
import re
from collections import OrderedDict

##Original Author: Luki Goldschmidt <lugo@uw.edu>
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
    for o in self.decoys:
      if o[self.decoy_field_name] == decoy:
        return o[scoreterm]
