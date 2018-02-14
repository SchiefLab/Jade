#!/usr/bin/env python

from jade.rosetta_jade.RunRosettaBenchmarks import *


##Main to RunRosetta.

def get_parser():
    run_rosetta = RunRosettaBenchmarks()
    return run_rosetta.parser

if __name__ == "__main__":
    run_rosetta = RunRosettaBenchmarks()
    run_rosetta.run()