import os
import sys
import re

from jade.basic.general import *
from jade.basic.path import *



def rm_features_dbs(outdir, out_names):

    for out_name in out_names:
        db_name = outdir+'/'+out_name
        if os.path.exists(db_name):
            os.remove(db_name)


def create_features_db(pdb_list,
                       xml_name,
                       compiler,
                       score_weights,
                       out_db_name,
                       out_db_batch,
                       outdir,
                       use_present_dbs,
                       indir = "",
                       mpi = True,
                       np = 5):

    """
    old_db_name = outdir+'/'+out_db_name+'.'+score_weights+".db3"
    new_db_name = outdir+'/'+out_db_name+'.'+xml_name+'.'+score_weights+".db3"
    if os.path.exists(old_db_name):
        os.system('mv '+old_db_name+' '+new_db_name)
        print "Old db name already exists.  Moving."
        return
    """

    features_dir = get_feat_input_path()
    outdir = outdir+"/databases"
    if not os.path.exists(outdir): os.mkdir(outdir)

    if not out_db_batch:
        sys.exit("Must have output db name and db batch to run_features")

    if not use_present_dbs:
        rm_features_dbs(outdir, out_db_name, score_weights, xml_name)

    if mpi:
        feat_dir = "/temp"
    else:
        feat_dir = outdir

    features_command = get_rosetta_program('rosetta_scripts', mpi=mpi, compiler = compiler) +' -parser:protocol '+features_dir+\
                       '/'+xml_name+'.xml @ '+features_dir+'/features.flag -l '+pdb_list+ \
                       ' -parser:script_vars name='+out_db_name+'.'+xml_name+ \
                       ' score='+score_weights+' batch='+out_db_batch+ " -out:path:all "+feat_dir

    if indir:
        features_command = features_command+' in:path:pdb '+indir

    if mpi:
        os.system("mpiexec -np "+np+" "+features_command + " -separate_db_per_mpi_process")
    else:
        os.system(features_command)

    if mpi:
        pass
        #Need to combine MPI databases, move them to outdir, and remove the temp directory.
