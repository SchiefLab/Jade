#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/protocols/for_design.py
## @brief  design protocols.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.basic.options import get_string_option
from rosetta.basic.options import get_boolean_option
from rosetta.core.pack.task.operation import *
#from rosetta.protocols.forge.remodel import *
from toolbox import mutate_residue

#Python Imports
import time
import os

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox

#Toolkit Imports
from ProtocolBaseClass import ProtocolBaseClass
from app.pyrosetta_toolkit.window_main import global_variables


class DesignProtocols(ProtocolBaseClass):
    def __init__(self, pose, score_class, input_class, output_class):
        ProtocolBaseClass.__init__(self, pose, score_class, input_class, output_class)

    def setupPackDesign(self, main):
        """
        Follows fixbb.cc to allow most user defined options within the GUI.
        Limitations: No symmetry.  Annealers should work fine through the options system.
        UI due to major TKinter bug on my mac.
        """
        if self.pose.total_residue()==0:print "Please load a pose."; return
        def packDesign():
            resfile = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Resfile..")
            if not resfile:return
            global_variables.current_directory = os.path.dirname(resfile)
            mover = DesignWrapper(self.score_class.score, resfile)
            mover.min_sc = min_sc.get()
            mover.min_pack = min_pack.get()
            mover.stochastic_pack = stochastic_pack.get()

            print "Please cite the many references included for Fixed Backbone Design in the Rosetta Manual."
            print "Further options such as annealing use the options system.  Symmetry not supported at this time."
            time.sleep(5)

            self.run_protocol(mover)
            self.main.destroy()

        top_level = Toplevel(main)
        self.main = top_level

        top_level.title("Fixbb Design Setup")
        min_sc = IntVar(); min_sc.set(False)
        min_pack = IntVar(); min_pack.set(False)
        stochastic_pack = IntVar(); stochastic_pack.set(False)

        label_sc = Label(top_level, text = "Do minimization of side chains after rotamer packing")
        label_pack = Label(top_level, text = "Pack and minimize sidechains simultaneously (slower)")
        label_stoch = Label(top_level, text = "Pack using a continuous sidechains rotamer library")

        check_sc = Checkbutton(top_level, text="-minimize_sidechains", variable = min_sc)
        check_pack = Checkbutton(top_level, text = "-min_pack", variable = min_pack)
        check_stoch = Checkbutton(top_level, text = "-stochastic_pack", variable = stochastic_pack)

        button_go = Button(top_level, text = "Run Protocol", command = lambda: packDesign())
        label_sc.grid(row=0, column=1, sticky=W); check_sc.grid(row = 0, column=0, sticky=W)
        label_pack.grid(row=1, column=1, sticky=W); check_pack.grid(row=1, column=0, sticky=W)
        #label_stoch.grid(row=2, column=1, sticky=W); check_stoch.grid(row=2, column=0, sticky=W); stochastic_pack currently does not work for some reason in python.

        button_go.grid(row=3, column=0, columnspan=2, sticky=W+E)

    def packDesign(self, resfile=False):
        """
        #Follows fixbb.cc to allow most user defined options within the GUI.
        #Limitations: No symmetry.  Annealers should work fine through the options system.
        Use setupPackDesign for more options.
        """

        if not resfile:
            resfile = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Resfile..")
            if not resfile:return
            global_variables.current_directory = os.path.dirname(resfile)

        #Can't ask more then 1 dialog in a row on the macbookpro I am using.  So screw it.  Using setupPackDesign instead
        #tkMessageBox.askyesno(title="-min_pack", message="Pack and minimize sidechains simultaneously?")


        print "Please cite the many references included for Fixed Backbone Design in the Rosetta Manual."
        print "Further options such as annealing use the options system.  Symmetry, min_pack and stochastic_min are not supported at this time."
        time.sleep(5)

        self.run_protocol(mover)

    def remodel(self, blueprint=False):
        """
        Ya.  Not yet unfortunately.
        """
        pass

        if not blueprint:
            blueprint = tkFileDialog.askopenfilename(initialdir = global_variables.current_directory, title="Open Blueprint..")
            if not blueprint:return

    def pack_residue(self, res, chain):
        """
        Packs an individual residue.  Scores before and after.
        """

        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        task = TaskFactory.create_packer_task(self.pose)
        task.restrict_to_repacking()
        #task = self._get_set_pack_neighbors(res, task)
        pack_radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=5.0)
        if not pack_radius: return;
        restrict_non_nbrs_from_repacking(self.pose, res, task, pack_radius)
        pack_mover = PackRotamersMover(self.score_class.score, task)
        print self.score_class.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_class.score(self.pose)

    def design_residue(self, res, chain):
        """
        Designs a residue by task not restricting to repacking.
        """
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        pack_radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=5.0)
        if not pack_radius: return;
        task = TaskFactory.create_packer_task(self.pose)
        for i in range(1, self.pose.total_residue()+1):
            if i != res:
                task.nonconst_residue_task(i).restrict_to_repacking()
        restrict_non_nbrs_from_repacking(self.pose, res, task, pack_radius)
        pack_mover = PackRotamersMover(self.score_class.score, task)
        print self.score_class.score(self.pose)
        pack_mover.apply(self.pose)
        print self.score_class.score(self.pose)

    def mutateRes(self, res, chain, new_res):
        new_res = new_res.split(":")
        new_res = new_res[2]
        res = self.pose.pdb_info().pdb2pose(chain, int(res))
        pack_radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=5.0)
        if pack_radius == None: return;
        #print radius
        print self.score_class.score(self.pose)
        print "Mutating to " + new_res
        mutate_residue(self.pose, res, new_res, pack_radius, self.score_class.score);
        print self.score_class.score(self.pose)
        print "Mutagenesis Complete."
        return self.pose

    def _get_set_pack_neighbors(self, res, task):
        """
        Asks for packing radius.  Sets task to prevent repacking for all but within that radius.
        Original Author: Evan H. Baugh, Johns Hopkins University. from mutants.py.
        Somehow not working with task
        """
        radius = tkSimpleDialog.askfloat(title="Pack radius", prompt = "Please enter the desired neighbor packing radius (A)", initialvalue=0.0)

        center = self.pose.residue( res ).nbr_atom_xyz()
        for i in range( 1 , self.pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
            if not i == res or center.distance_squared(self.pose.residue( i ).nbr_atom_xyz() ) > radius**2:
                task.nonconst_residue_task( i ).prevent_repacking()

        return task

class DesignWrapper:
    """
    Wrapper to re-init tasks for resfile if multiple rounds are given.
    """
    def __init__(self, score, resfile):
        self.score = score
        self.resfile = resfile
        self.min_pack = False
        self.min_sc = False
        self.stochastic_pack = False

        """
        try:
            self.min_pack = get_boolean_option('fixbb:minimize_sidechains')
            self.min_sc = get_boolean_option('fixbb:min_pack')
            self.stochastic_pack = get_boolean_option('fixbb:stochastic_pack')
        except PyRosettaException:
            pass
        """
    def set_minimize_sidechains(self):
        self.min_sc = True

    def apply(self, pose):
        task = TaskFactory()
        task.push_back(InitializeFromCommandline())
        task.push_back(ReadResfile(self.resfile))

        s_mover = SequenceMover()

        if self.min_pack or self.stochastic_pack:
            design_mover = MinPackMover()
            design_mover.task_factory(task)
            design_mover.score_function(self.score)
            #Stochastic Packing is resulting in segfault.  Use dun10 anyway for this.
            #if self.stochastic_pack:
                #design_mover.stochastic_pack(True)
            s_mover.add_mover(design_mover)
        else:
            design_mover = PackRotamersMover()
            design_mover.task_factory(task)
            design_mover.score_function(self.score)
            s_mover.add_mover(design_mover)


        if self.min_sc:
            mm = MoveMap(); #Empty movemap
            min_mover = MinMover(mm, self.score, get_string_option('run:min_type'), 0.01, True)
            task_min_mover = TaskAwareMinMover(min_mover, task)
            s_mover.add_mover(task_min_mover)

        s_mover.apply(pose)
