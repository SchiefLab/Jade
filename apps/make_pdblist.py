#!/usr/bin/env python


##This is for when there are too many PDBs in a directory for ls -1 *.pdb* > PDBLIST.txt to work.
#Super basic, as I need it now.

import os
import glob

PDBLIST = open("PDBLIST.txt", 'w')

pdb_files = glob.glob("*.pdb*")
for f in pdb_files:
    PDBLIST.write(f+"\n")
PDBLIST.close()




