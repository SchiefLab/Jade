#!/usr/bin/env python

#Author: Jared Adolf-Bryfogle

#Purpose:
#  Create a JD3 Job Definition file in to scan compatible surface residues for glycosylation.
#  By default Skips Disulfides, prolines, glycines, makes sure +2 is not buried and not boundary and hydrophobic.
#  Also, by default designs +2 as Threonine as this is shown to increase occupancy of glycan.
#
# This is a very conservative SugarCoat application, but all options can be changed.
#
#  Optionally takes a ResidueSelector to further limit scan.
#
#  The JD3 is then set to run the glycan_scanner.xml script, which does:
#    1) Creates sequon using the CreateSequonMover
#    2) Creates glycans using the SimpleGlycosylateMover
#    3) Models the new glycan using the GlycanTreeMover.
#
# Additionally, the glycan_scanner script will run various SimpleMetrics throughout the run to aid in analysis.
#
#  This is meant to be a stop-gap until we have a proper JD3 SugarCoat/GlycanDesign application.


from __future__ import print_function

from jade.rosetta_jade.flag_util import get_common_flags_string_for_init
from jade.rosetta_jade.jd3 import create_substituted_jd_string
from jade.basic.path import *
import sys,re
from argparse import ArgumentParser
from pyrosetta import *
from pyrosetta.rosetta import *

from pyrosetta.rosetta.core.chemical import *
from pyrosetta.rosetta.core.select.residue_selector import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.core.simple_metrics.metrics import *
from rosetta.utility import *

from collections import defaultdict

init(get_common_flags_string_for_init())



if __name__ == "__main__":
    parser = ArgumentParser('Create a JD3 Job Definition file to scan compatible surface residues and '
                            'do denovo glycosylation and modeling. This is for IDEAL sequons only')

    parser.add_argument('--selections', '-x',
                        help = 'Input XML to pull ResidueSelections to limit scan.',
                        default = 'selections.xml')

    parser.add_argument('--selector', help = 'Name of a residue selector from input xml to limit search to. Loops, regions, etc.  '
                                             'This will be combined as AND logic')

    parser.add_argument('-s', help = "Input PDB we will be scanning.", required=True)

    parser.add_argument('--include_prolines', action='store_true', default=False)

    parser.add_argument('--include_glycines_if_D', action='store_true', default=False,
                        help = "Include glycines?  Glycines can either be D/L conformation.  Set this option to only include glycines if they are D.  Otherwise, skip them at +2")

    parser.add_argument('--design_plus_1', action='store_true', default=False,
                        help = "Design the +1 position?")

    parser.add_argument('--enable_T_and_S', action='store_true', default=False,
                        help = 'Enable Threonine and Serine at +2?  By default we use only threonine as this has been shown to have'
                               'higher occupancy at the glycan site in general. ')

    parser.add_argument('--design_sequon_regardless_of_plus_two_layer', default=False, action = 'store_true',
                        help = 'Enable design of the primary sequon if the +2 residue is part of the core (and not already ser/thr)?')

    parser.add_argument('--include_boundary_buried', default = False, action='store_true',
                        help = "Only use the LayerSelector to define the core for the +2 residues.  Otherwise, we use the SC sasa measured in addition for the boundary residues and do not include these if not designing")

    parser.add_argument('--include_boundary_hydrophobics', default = False,
                        help = "Include hydrophobic boundary residues (for +2 design)?  By default we leave these out to maintain foldability" )

    options = parser.parse_args()

    pose = pose_from_pdb(options.s)

    layer_selector = LayerSelector()

    ###############
    ### Surface ###
    ###############
    #set_layers( bool const pick_core, bool const pick_boundary, bool const pick_surface );
    layer_selector.set_layers(False, False, True)


    selection = layer_selector.apply(pose)
    surface_selection = layer_selector.apply(pose) #avoid shallow copy

    pymol_selection = SelectedResiduesPyMOLMetric(layer_selector)
    sele = pymol_selection.calculate(pose)

    print("\nPyMOL Selection - Surface")
    print(sele.replace('rosetta_sele', 'surface'))


    ############
    ### Core ###
    ############
    layer_selector.set_layers(True, False, False)
    core_selection = layer_selector.apply(pose)

    pymol_selection = SelectedResiduesPyMOLMetric(layer_selector)
    sele = pymol_selection.calculate(pose)

    print("\nPyMOL Selection - Core")
    print(sele.replace('rosetta_sele', 'core'))

    ############
    ### Core extra check - 0 Sasa ###
    ############
    sasa_calc = rosetta.core.scoring.sasa.SasaCalc()
    sasa_calc.calculate(pose)

    sasa_selection = vector1_bool(pose.size())
    res_data = sasa_calc.get_residue_sasa_sc()
    for i in range(1, pose.size() +1 ):
        if res_data[i] > 1:
            sasa_selection[i] = False
        else:
            sasa_selection[i] = True

    return_core = ReturnResidueSubsetSelector(sasa_selection)
    pymol_selection = SelectedResiduesPyMOLMetric(return_core)
    sele = pymol_selection.calculate(pose)

    print("\nPyMOL Selection - Core/SASA")
    print(sele.replace('rosetta_sele', 'core_sasa'))

    ################
    ### Boundary ###
    ################
    layer_selector.set_layers(False, True, False)
    boundary_selection = layer_selector.apply(pose)

    pymol_selection = SelectedResiduesPyMOLMetric(layer_selector)
    sele = pymol_selection.calculate(pose)

    print("\nPyMOL Selection - Boundary")
    print(sele.replace('rosetta_sele', 'boundary'))

    if options.selector:
        if not os.path.exists(options.xml):
            sys.exit("XML file for residue selections does not exist!")

        xml_objects = create_from_file( options.xml)
        user_res_selector = xml_objects.get_residue_selector(options.selector)

        and_selector = AndResidueSelector( layer_selector, user_res_selector)
        selection = and_selector.apply(pose)


    #Limit selection
    disulfides = 0
    prolines = 0


    for i in range(1, pose.size() + 1 ):

        res_type = pose.residue_type(i).aa()
        # Limit to Protein residues
        if not pose.residue_type(i).is_protein():
            selection[i] = False

        # Limit Disulfides
        if pose.residue_type(i).is_disulfide_bonded():
            disulfides+=1
            selection[i] = False

            #Kill -1 position and -2 position for potential glycosylation as the disulfide will be in the sequon (which we are designing).
            if i >= 2:
                selection[ i - 1 ] = False
            if i >= 3:
                selection[ i - 2 ] = False

        # Limit prolines - Remove any sequon that includes a proline
        if (not options.include_prolines) and ( res_type == aa_pro):
            prolines+=1
            selection[i] = False

            if i >= 2 :
                selection[ i - 1 ] = False

            if i >= 3 :
                selection[ i - 2 ] = False

        #Make sure the full sequon is at least one residue away from the end.
        if (not i+3 <= pose.conformation().chain_end(pose.chain(i))):
            selection[ i ] = False
            continue

        ##Ser/Thr - By default we are very conservative here.
        res_type_p_1 = pose.residue_type( i+1 ).aa()
        res_type_p_2 = pose.residue_type( i+2 ).aa()

        if res_type == aa_gly:
            if not options.include_glycines_if_D:
                selection[ i ] = False
            else:
                sys.exit("Contingient inclusion of glycine is not yet implemented.")
        if selection[i] and (res_type_p_2 == aa_gly):
            if not options.include_glycines_if_D:
                selection[i] = False
            else:
                sys.exit("Contingient inclusion of glycine is not yet implemented.")

        ### Plus 2 Contingencies ###
        if surface_selection[i+2]:
            continue

        elif res_type_p_2 == aa_thr:
            continue

        elif res_type_p_2 == aa_ser and options.enable_T_and_S:
            continue

        elif options.design_sequon_regardless_of_plus_two_layer:
            continue

        #Restrict on Core or Boundary depending on options
        elif core_selection[i+2]:
            selection[i] = False
            continue

        #Sasa - fully buried residues that probably shouldn't be designed.
        elif boundary_selection[i+2] and not options.include_boundary_buried and sasa_selection[i+2]:
            selection[i] = False
            continue

        #Hydrophobic -
        # clearly hydrophobic side-chains (FMILYVW).
        # Note that it also excludes glycine, alanine, cysteine, and proline.
        elif boundary_selection[i+2] and not options.include_boundary_hydrophobics and pose.residue_type(i+2).has_property(rosetta.core.chemical.HYDROPHOBIC):
            selection[i] = False
            continue

    #Limit to Contiguous selections. Not Needed as glycosylation is co-translational.
    """
    final_residues = []
    for i in range(1, pose.size() + 1 ):

        if not selection[i]: continue
        #print(i, pose.pdb_info().pose2pdb(i))

        #Make sure Next residue is not end of chain and on
        if i + 1 <= pose.size():
            if not selection[i+1]: continue
        else:
            continue
        #Make sure Third residue in Sequon is not end of chain and on.
        if i + 2 <= pose.size():
            if not selection[i+2]: continue
        else:
            continue


        #We have a contiguous section, add the residue to be the N position of glycosylation
        final_residues.append(i)
    """


    store_subset = ReturnResidueSubsetSelector(selection)
    pymol_selection.set_residue_selector(store_subset)
    final_sele = pymol_selection.calculate(pose)
    print("\nFinal Selection of each start of glycan motif")
    print(final_sele.replace('rosetta_sele', 'final_selection'))

    plus_one = vector1_bool(pose.size())
    plus_two = vector1_bool(pose.size())

    for i in range(1, pose.size()+1):
        if selection[i]:
            plus_one[i+1] = True
            plus_two[i+2] = True

    store_subset.set_residue_subset(plus_one)
    pymol_selection.set_residue_selector(store_subset)
    plus_one_s = pymol_selection.calculate(pose)

    print("\nFinal Selection of Plus ONE")
    print(plus_one_s.replace('rosetta_sele', 'plus_one'))

    store_subset.set_residue_subset(plus_two)
    pymol_selection.set_residue_selector(store_subset)
    plus_two_s = pymol_selection.calculate(pose)

    print("\nFinal Selection of Plus TWO")
    print(plus_two_s.replace('rosetta_sele', 'plus_two'))

    jobs = []

    #Figure out motif we are using
    motif = "N"

    if options.design_plus_1:
        motif+="[^P]"
    else:
        motif+="-"


    if options.enable_T_and_S:
        motif+="[ST]"
    else:
        motif+="T"

    for res in range(1, pose.size() +1):
        if not selection[res]: continue

        job = defaultdict()
        job['fname'] =  options.s
        job['start'] = str(res)
        job['end'] = str( res + 2)
        job['motif'] = motif
        job['start_pdb'] = pose.pdb_info().pose2pdb(res).replace(' ', '')
        jobs.append(job)

    print( "Total Residues:", pose.size())
    print( "Total Jobs: ", len(jobs))
    print( "Skipped",prolines,"prolines")
    print( "Skipped",disulfides, "disulfides")

    jd_scanner = get_xml_scripts_path()+"/glycan_scanner_jd3_jd_base.xml"
    xml_scanner = get_xml_scripts_path()+"/simple_glycan_scanner_manual.xml"

    if not os.path.exists("job_definitions"):
        os.mkdir("job_definitions")

    if not os.path.exists("xmls"):
        os.mkdir("xmls")

    os.system('cp '+xml_scanner+ " xmls")
    print("Copied "+xml_scanner +" to /xmls")

    create_substituted_jd_string(jd_scanner, jobs)