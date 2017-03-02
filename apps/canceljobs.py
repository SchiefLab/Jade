#!/usr/bin/env python

import sys
import os

#Cancel jobs for slurm (scancel)

if len(sys.argv) != 3:
    print "Use: canceljobs.py start end"
if sys.argv[1] == "--help":
    print "Use: canceljobs.py start end"


for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    cmd = 'scancel '+repr(i)
    print cmd
    os.system(cmd)
print "Done"
