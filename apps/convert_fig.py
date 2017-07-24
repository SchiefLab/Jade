#!/usr/bin/env python
import os
import sys

from __future__ import print_function

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Converts Images to TIFF figures at 300 DPI for publication using sips. "
                 " \n Arguments: Input_Path Export_Path [optional] Non-Tiff_format\n"
                 "Example: convert_to.py in_fig.pdf out_fig.tiff\n"
                 "Example: convert_to.py in_fig.png out_fig.eps eps")

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    outformat = "tiff"
    if len(sys.argv) == 4:
        outformat = sys.argv[3]

    cmd = "sips -s format {format} {in_path} -s dpiHeight 300 -s dpiWidth 300 --out {out_path}".format(
        format = outformat, in_path = in_file, out_path = out_file)

    print(cmd)

    os.system(cmd)
    print("done")



