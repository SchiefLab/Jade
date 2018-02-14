#!/usr/bin/env python2.7

from jade.RAbD_BM.RunBenchmarksRAbD import RunBenchmarksRAbD
from argparse import ArgumentParser

def get_parser():
    parser = ArgumentParser(description="This program runs Rosetta MPI locally or on a cluster using slurm or qsub.  "
                                        "Relative paths are accepted.")
    return parser

if __name__ == "__main__":
    bm = RunBenchmarksRAbD()
    bm.run()