#!/usr/bin/env python

import os
import sys
import glob
import re
from collections import defaultdict
from rosetta import *
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.modules.ScoreAnalysis import ScoreAnalysis
from app.pyrosetta_toolkit.modules.Structure import *
from app.pyrosetta_toolkit.modules.definitions import restype_definitions
from app.pyrosetta_toolkit.modules.DesignBreakdown import *

from antibody.analysis.SequenceStats import SequenceStats
from antibody.analysis.PyMolScriptWriter import *
from antibody.analysis.PythonPDB2 import *
from antibody.analysis import create_features_json as json_creator


rosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy false -ex1 -ex2 -use_input_sc")

def main():
    ####################################################################################################################
    ##                                                  OPTIONS
    ####################################################################################################################


    parser = OptionParser()

    ############################
    ## Required Options
    ############################
    parser.add_option("--PDBLIST",
                      help = "Analyze a PDBLIST of the strategy.  Relative or full paths",
                      default=None)

    parser.add_option("--indir",
                      help = "Input directory used for either PDBLIST path (if PDBLIST and this is given or a path full of PDBs to analyze",
                      default = None)

    parser.add_option("--out_name",
                      help = "Name of the output Features db and other out names")


    parser.add_option("--out_db_batch",
                      help="Batch name for databases")

    ###########################
    ## Analysis Types
    ###########################

    parser.add_option("--do_all",
                      help = "Run all available analysis",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_rescore",
                      help = "Rescore the PDBs",
                      default = False,
                      action  = "store_true")

    parser.add_option("--do_run_antibody_features_all",
                      help = "Run the AntibodyFeatures to create databases for the strategy on all PDBs",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_run_antibody_features_top",
                      help = "Run the AntibodyFeatures to create databases for the strategy on the top x PDBs.  If not rescored, will rescore.",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_run_cluster_features_all",
                      help = "Run the ClusterFeatures reporter on all PDBs in list or directory",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_run_cluster_features_top",
                      help = "Run the ClusterFeatures reporter on the top x PDBs of the strategy.  Rescore PDB list if nessessary.",
                      default = False,
                      action = "store_true")

    ############################
    ## Optional
    ############################

    parser.add_option("--outdir",
                      help = "\nOutput directory",
                      default = "antibody_design_analysis_results")

    parser.add_option("--np",
                      help = "Number of processors to use for rescoring PDBs",
                      default = 1)

    parser.add_option("--score_weights",
                      help = "Weights to use during rescore or FeaturesReporters",
                      default = "talaris2013")

    parser.add_option("--rosetta_extension",
                      default = "linuxclangrelease")

    parser.add_option("--use_present_dbs",
                      default = False,
                      help = "Do not attempt to delete features databases present",
                      action = "store_true")

    parser.add_option("--top_score",
                      help = "Number of top scoring models to analyze",
                      default = 10)

    parser.add_option("--use_full_name_for_pymol_sessions",
                      help = "Use the full name for pymol sessions output from rescoring the PDBLIST",
                      default = False,
                      action = "store_true")

    parser.add_option("--native_path",
                      help = "Model of the native.  If passed will add to pymol sessions",
                      default = None)

    (options, args) = parser.parse_args(sys.argv)

    if not options.PDBLIST and not options.indir:
        sys.exit("Cannot analyze strategy without a PDBLIST or an input directory full of PDBs")
    if options.indir and not os.path.exists(options.indir):sys.exit(options.indir+" does not exist - cannot continue")




    ##### Get or make the path to the PDBLIST ######

    if not options.PDBLIST and options.indir:
        pdb_list = make_pdblist(options.indir)
    elif options.PDBLIST and options.indir:
        PDBLIST = open(options.PDBLIST, 'r')
        NEW_PDBLIST = open(options.indir+"/FULLPATH_PDBLIST.txt")

        for line in PDBLIST:
            new_line = options.indir+"/"+line.strip()
            if not os.path.exists(new_line):
                sys.exit(new_line+" does not exist - cannot continue")

            NEW_PDBLIST.write(new_line+"\n")
        NEW_PDBLIST.close()
        PDBLIST.close()
        pdb_list = options.indir+"/FULLPATH_PDBLIST.txt"
    else:
        pdb_list = options.PDBLIST



    if not os.path.exists(pdb_list): sys.exit(pdb_list+" does not exist!")

    os.system('mkdir '+options.outdir)


    do_all = options.do_all




    ####################################################################################################################
    ##                     Analysis Components.  One option for each type of analysis.
    ####################################################################################################################

    top_pdblist = options.outdir+"TOP_"+str(options.top_score)+"_"+options.out_name+"ORDERED_PDBLIST.txt"
    top_dir = options.outdir+"/TOP_"+str(options.top_score)+"_"+options.out_name


    if not options.out_name:
        sys.exit("Output name required")

    #### Rescore ####
    if do_all or options.do_rescore:
        rescore_pdb_list(pdb_list, options.top_score, options.score_weights, options.np, options.outdir, options.out_name,
                         options.use_full_name_for_pymol_sessions, options.native_path)


    #### Create Features Databases ###
    if do_all or options.do_run_antibody_features_all:

        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(pdb_list, fdir,
                       'antibody_features',  options.rosetta_extension,
                       options.score_weights,options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs)

    if do_all or options.do_run_cluster_features_all:

        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(pdb_list, fdir,
                       'cluster_features',  options.rosetta_extension,
                       options.score_weights,options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs)

    if do_all or options.do_run_antibody_features_top:
        if not os.path.exists(top_pdblist):
            rescore_pdb_list(pdb_list, options.top_score, options.score_weights, options.np, options.outdir, options.out_name,
                            options.use_full_name_for_pymol_sessions, options.native_path)


        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(top_pdblist, fdir,
                       'antibody_features',  options.rosetta_extension,
                       options.score_weights,options.out_name+"_TOP_E_"+str(options.top_score),
                       options.out_db_batch, options.outdir, options.use_present_dbs, top_dir)

    if do_all or options.do_run_cluster_features_top:
        if not os.path.exists(top_pdblist):
            rescore_pdb_list(pdb_list, options.top_score, options.score_weights, options.np, options.outdir, options.out_name,
                            options.use_full_name_for_pymol_sessions, options.native_path)

        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(top_pdblist, fdir,
                       'cluster_features',  options.rosetta_extension,
                       options.score_weights,options.out_name+"_TOP_E_"+str(options.top_score),
                       options.out_db_batch, options.outdir, options.use_present_dbs, top_dir)


    print "Complete"


def run_clustal_omega(pdb_list, outdir):
    pass

def create_sequence_fasta(pdb_list, outdir):
    pass

def rm_features_db(outdir, out_name, score_weights, xml_name):
    db_name = outdir+'/'+out_name+'.'+xml_name+'.'+score_weights+".db3"
    if os.path.exists(db_name):
        os.remove(db_name)

def create_features_db(pdb_list, features_dir, xml_name, extension, score_weights, out_db_name, out_db_batch, outdir, use_present_dbs, indir = ""):

    """
    old_db_name = outdir+'/'+out_db_name+'.'+score_weights+".db3"
    new_db_name = outdir+'/'+out_db_name+'.'+xml_name+'.'+score_weights+".db3"
    if os.path.exists(old_db_name):
        os.system('mv '+old_db_name+' '+new_db_name)
        print "Old db name already exists.  Moving."
        return
    """

    outdir = outdir+"/databases"
    if not os.path.exists(outdir): os.mkdir(outdir)

    if not out_db_batch:
        sys.exit("Must have output db name and db batch to run_features")

    if not use_present_dbs:
        rm_features_db(outdir, out_db_name, score_weights, xml_name)

    features_command = 'rosetta_scripts.'+extension +' -parser:protocol '+features_dir+\
                       '/'+xml_name+'.xml @ '+features_dir+'/features.flag -l '+pdb_list+ \
                       ' -parser:script_vars name='+out_db_name+'.'+xml_name+ \
                       ' score='+score_weights+' batch='+out_db_batch

    if indir:
        features_command = features_command+' in:path:pdb '+indir

    os.system(features_command)
    os.system('mv '+out_db_name+'.'+xml_name+'.'+score_weights+'.db3'+' '+outdir)

def get_data_from_features_db(db_path):
    command = """
                SELECT
                    structures.input_tag,
                    interfaces.dG as dG,
                    structure_scores.score_value as total_score,
                    interfaces.dG_cross as dG_cross,
                    interfaces.delta_unsatHbonds as delta_unsatHbonds,
                    interfaces.hbond_E_fraction as hbond_E_fraction,
                    interfaces.dSASA as dSASA,
                    interfaces.interface as interface
                FROM
                    interfaces,
                    score_types,
                    structure_scores,
                    structures
                WHERE
                    score_types.score_type_name='total_score' AND
                    structure_scores.score_type_id = score_types.score_type_id AND
                    structure_scores.struct_id = interfaces.struct_id AND
                    structure_scores.struct_id = structures.struct_id
                ORDER BY total_score
            """

    db = sqlite3.connect(db_path)
    cur = db.cursor()


    for row in cur.execute(command):
        pass

def rescore_pdb_list(pdb_list,top_score, score_weights, np, outdir, out_name, use_full_name_for_pymol_sessions,
                     native_path = None):
    """
    Rescores A PDBLIST,
    creates a Topx PDBLIST,
    copies structures into new output directory for top 10,
    creates a PyMol session for TOP 10 with new names as model:energy for easy analysis.
    """


    top_pdblist = out_name+"_SCORED_PDBLIST_TOP_"+str(top_score)
    scored_pdblist = make_scored_pdblist(pdb_list, score_weights, np, outdir, out_name)
    s = ScoreAnalysis()
    s.read_scores(scored_pdblist)

    top_dir = outdir+"/top_"+str(top_score)+"_"+out_name
    top_scoring_map = s.get_top_scoring_by_number(int(top_score), False)
    s.copy_results2(top_scoring_map, outdir, top_dir, top_pdblist)
    os.system('rm '+top_pdblist); #We don't really need this now that we have ordered PDBLISTS.

    #Add PyMOL Sessions for top scoring in order
    pdb_paths = []
    final_names = []
    load_as = None

    for path in sorted(top_scoring_map, key = top_scoring_map.get):
        new_path = top_dir+"/"+os.path.basename(path)
        pdb_paths.append(new_path)

        if not use_full_name_for_pymol_sessions:
            if not load_as:
                load_as = []
            SP1 = new_path.split('_')
            SP2 = SP1[-1].split('.')
            num = SP2[0]
            new_name = "model_"+num+"_E_"+str(top_scoring_map[path])
            load_as.append(new_name)
            final_names.append(new_name)
            print new_name
        else:
            basename = os.path.basename(new_path)


            final_names.append(basename)

    #If we already have the PSE, do not attempt to overwrite it as we may have already started analysis.
    if os.path.exists(top_dir+"/"+out_name+"_TOP_"+str(top_score)+".pse"):
        return

    make_pymol_session_on_top(top_dir, pdb_paths, load_as, top_dir, str(top_score), out_name, native_path)

def make_pymol_session_on_top(top_dir, pdb_path_list, load_as_list, outdir, out_name, top_num = None, native_path = None):

    if top_num:
        pse_path = outdir+"/"+out_name+"_top_"+str(top_num)+".pse"
    else:
        pse_path = outdir+"/"+out_name+"_all"+".pse"
    if os.path.exists(pse_path):
        print "Not overriding PSE: "+pse_path
        #return

    if len(pdb_path_list) == 0:
        print "PDB list path empty.  Skipping creation of pymol session"
        return
    
    scripter = PyMolScriptWriter(top_dir)

    if native_path:
        scripter.add_load_pdb(native_path, "native_"+os.path.basename(native_path))

    scripter.add_load_pdbs(pdb_path_list, load_as_list)
    scripter.add_align_all_to(scripter.get_final_names()[0])
    scripter.add_show("cartoon")
    scripter.add_line("center")
    scripter.add_save_session(pse_path)
    scripter.write_script("load_align_top.pml")
    run_pymol_script(top_dir+"/"+"load_align_top.pml")

def make_scored_pdblist(pdb_list, score_weights, np, outdir, outname):
    scored_pdblist = output_tools.score_PDBLIST(pdb_list, create_score_function(score_weights), int(np))

    os.system('cp '+scored_pdblist+' '+outdir)
    final_path = outdir+'/'+outname+"_SCORED_PDBLIST.txt"
    os.system('mv '+outdir+'/SCORED_PDBLIST.txt'+' '+final_path)
    print "Copied to "+final_path
    return final_path

def make_pdblist(in_dir):
    pdb_list = glob.glob(in_dir+"/*.pdb*")
    OUTFILE = open(in_dir+"/PDBLIST.txt", 'w')
    for pdb in pdb_list:
        if re.search("initial_benchmark_perturbation", pdb) or re.search("excn", pdb):
            continue
        OUTFILE.write(pdb+"\n")
    OUTFILE.close()

    return in_dir+"/PDBLIST.txt"



if __name__ == "__main__":
    main()