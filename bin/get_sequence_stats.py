from rosetta import *
from Tkinter import *

from argparse import ArgumentParser
from rosetta_general.DesignBreakdown import DesignBreakdown

if __name__ == '__main__':
    """
    For testing and use outside of GUI.

    """
    Tk()
    rosetta.init()

    parser = ArgumentParser()
    parser.add_argument("--fasta", "-f",
        help = "Path to FASTA file for design comparisons.  Can be full length or of a specific region."
    )
    parser.add_argument("--outpath","-o",
        default="sequence_results",
        help = "Full output directory path.  Default is fasta file path /Results"
    )

    parser.add_argument("--native", "-n",
        default=None,
        help = "Reference pose for numbering and comparison (Required)"
    )

    parser.add_argument("--regions", "-g",
        default=[],
        nargs="*",
        help = "Region(s)  - if none is given in Fasta + Not whole structure used for comparison (region designated as start:end:chain)"
    )

    parser.add_argument("--prefix",
        default = "",
        help = "Prefix to use for output files.  Generally recommended."
    )

    options = parser.parse_args()


    if not options.fasta or not os.path.exists(options.fasta):
        sys.exit("Please specify a FASTA file to use for calculations. You can use get_seq.py or the PyRosetta Toolkit")


    if not options.native:
        sys.exit("Native pdb is required for now for numbering purposes.")

    if not options.prefix:
        sys.exit("A prefix is required for output")

    if options.region:
        for region in options.region:

            prefix = options.prefix+'_'+region.replace(':', '-')

            breakdown = DesignBreakdown(options.fasta, options.native, directory=options.outpath, region=region, prefix=prefix)
            breakdown.run_outputs()
    else:
            breakdown = DesignBreakdown(options.fasta, options.native, directory=options.outpath, prefix=options.prefix)
            breakdown.run_outputs()