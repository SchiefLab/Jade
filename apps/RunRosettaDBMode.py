#!/usr/bin/env python

from jade.rosetta_jade.RunRosetta import *


##Main to RunRosetta.

if __name__ == "__main__":
    run_rosetta = RunRosetta(db_mode=True)
    run_rosetta.run()