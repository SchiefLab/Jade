#!/usr/bin/env python

# Author: Jared Adolf-Bryfogle

#This script will create either cluster features or antibody features json for use in Features R script.
#I am just really sick of doing this by hand.
#Example Cmd-line:  python create_features_json.py --database databases/baseline_comparison.txt --scripts cluster

from jade.rosetta_jade.FeaturesJsonCreator import *
from argparse import ArgumentParser

def get_parser():
    parser = ArgumentParser(description="This script will create either cluster features or antibody features json for use in Features R script.\n"
                            "Example Cmd-line:  python create_features_json.py --database databases/baseline_comparison.txt --scripts cluster")

    parser.add_argument("--databases", "-l",
        help = "List of dbs: db_name,short_name,ref keyword if the reference databaseSeparated by white space. ",
        default = [],
        nargs = "*",


    )

    parser.add_argument("--script", "-s",
        help = "Script type.  Will setup the appropriate output formats and R scripts",
        default = "antibody_minimal",
        choices = ["cluster", "antibody", "interface", "antibody_minimal"]
    )

    parser.add_argument("--db_path", "-p",
        help = "Path to databases.  Default is pwd/databases",
        default=os.getcwd()+"/databases")

    parser.add_argument("--outdir", "-o",
        help = "Where to put the result of the analysis scripts.  Currently unsupported by the features framework.",
        default=os.getcwd()+"/plots")

    parser.add_argument("--outname", "-n",
        help = "Output file name of json file",
        default = "local_json_compare_ss.json" )

    parser.add_argument("--add_comparison_to_this_json", "-a",
        help = "Add all this data to this json as more sample sources.")

    parser.add_argument("--run", "-r",
        help = "Go ahead and run compare_sample_sources.R.  Must be in path!!",
        default = False,
        action = "store_true")

    return parser

if __name__ == "__main__":


    parser = get_parser()

    options = parser.parse_args()

    if not options.databases:
        sys.exit("No databases given!!")

    if not os.path.exists(options.db_path):
        options.db_path = os.getcwd()

    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    json_creator = JsonCreator(options.outdir, options.script)

    for db in options.databases:
        dbSP = db.split(',')

        ref = False
        if len(dbSP) == 3 and dbSP[2] == "ref": ref = True
        db_path = options.db_path+"/"+dbSP[0]; id= dbSP[1]
        json_creator.add_sample_source_info(db_path, id, ref)

    if options.add_comparison_to_this_json:
        json_creator.add_current_sample_sources_to_json_dict(options.add_comparison_to_this_json)

    if not re.search("json", options.outname): options.outname = options.outname+".json"

    json_creator.save_json(options.outdir+"/"+options.outname)

    if options.run:
        json_creator.run_json()

    print "Complete"
