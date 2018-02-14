#!/usr/bin/env python2.7

from jade.RAbD_BM.RunBenchmarksRAbD import RunBenchmarksRAbD

def get_parser():
    bm = RunBenchmarksRAbD()
    return bm.parser

if __name__ == "__main__":
    bm = RunBenchmarksRAbD()
    bm.run()