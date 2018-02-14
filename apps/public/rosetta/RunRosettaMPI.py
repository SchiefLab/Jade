#!/usr/bin/env python

from jade.rosetta_jade.RunRosetta import *

def get_parser():
    run_rosetta = RunRosetta()
    return run_rosetta.parser

##Main to RunRosetta.

if __name__ == "__main__":
    run_rosetta = RunRosetta()
    run_rosetta.run()