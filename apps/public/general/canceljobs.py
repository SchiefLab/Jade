#!/usr/bin/env python

import sys
import os
from argparse import ArgumentParser

def get_parser():
    parser = ArgumentParser(description="Call scancel to cancel a consecutive set of cluster job numbers")
    return parser



if __name__ == "__main__":

    # Cancel jobs for slurm (scancel)
    parser = ArgumentParser()

    if len(sys.argv) < 3:
        print "Use: canceljobs.py start end"
        sys.exit()
    if sys.argv[1] == "--help":
        print "Use: canceljobs.py start end"
        sys.exit()

    for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        cmd = 'scancel '+repr(i)
        print cmd
        os.system(cmd)
    print "Done"
