from __future__ import print_function

from collections import namedtuple
from pymol import cmd


def select_neighbors(distance=5.0, selection_name='sele'):
    """
    Select neighbors around the selection or passed in selections at the given distance

    :param distance: float
    :param selection_name: A selection object, or a specific selection

    :rtype: None
    """
    cmd.select('root_sele', selection_name)
    cmd.select('sele', 'br. {} around {}'.format(str(selection_name), distance))
    cmd.enable('sele')


def show_polar_contacts(sel1='all', sel2="all"):
    """
    Show polar contacts between the two contacts
    Super stupid-simple, but a much better name then distance.  Who would have thought I could do this?

    USE: show_polar_contacts antigen, ab

    :param selection_name:
    :return:
    """
    cmd.distance(name='{}_{}_polars'.format(str(sel1), str(sel2)), selection1=str(sel1), selection2=str(sel2), cutoff=3.2, mode=2)



# Copyright (c) 2010 Robert L. Campbell
# Adapted, and modified original list_hb function.
# This function is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
# (http://creativecommons.org/licenses/by-nc-sa/4.0/)
def show_hbonds(selection='all',selection2=None,cutoff=3.2,angle=55,hb_list_name='hbonds'):
    """
      USAGE

      list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                         [angle (default=55)], [mode (default=1)],
                         [hb_list_name (default='hbonds')]

      The script automatically adds a requirement that atoms in the
      selection (and selection2 if used) must be either of the elements N or
      O.

      e.g.
      To get a list of all H-bonds within chain A of an object
        list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds

      To get a list of H-bonds between chain B and everything else:
        list_hb 1tl9 & c. b, 1tl9 &! c. b

    """


    cutoff = float(cutoff)
    angle = float(angle)
    # ensure only N and O atoms are in the selection
    selection = str(selection) + " & e. n+o"
    if not selection2:
        hb = cmd.find_pairs(selection, selection, mode=1, cutoff=cutoff, angle=angle)
    else:
        selection2 = str(selection2) + " & e. n+o"
        hb = cmd.find_pairs(selection, selection2, mode=1, cutoff=cutoff, angle=angle)

    # sort the list for easier reading
    hb.sort(lambda x, y: (cmp(x[0][1], y[0][1])))

    for pairs in hb:
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'print "%1s/%3s`%s/%-4s " % (chain,resn,resi,name),')
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'print "%1s/%3s`%s/%-4s " % (chain,resn,resi,name),')
        print("%.2f" % cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]),
                                    "%s and index %s" % (pairs[1][0], pairs[1][1])))

    cmd.h_fix()

### Hydrogen Adding and Showing ###
def add_hydrogens(show_polar_only=True):
    cmd.h_add()
    if show_polar_only:
        hide_hydrogens()


def hide_hydrogens(keep_polar_only=True):
    representations = ['lines', 'sticks']
    for rep in representations:
        if keep_polar_only:
            cmd.hide(rep, 'elem H and neighbor elem C')
        else:
            cmd.hide(rep, 'hydro')


def show_hydrogens(polar_only=True, representation='sticks'):
    if polar_only:
        cmd.show(representation, 'elem H and neighbor elem O')
        cmd.show(representation, 'elem H and neighbor elem N')
    else:
        cmd.show(representation, 'elem H')

def generate_symmetry_mates(selection = 'all', distance=6.0):
    cmd.symexp("symm", selection, '('+selection+')', distance)

### Some glycan extensions ###



### Extend the command interface ###
cmd.extend('select_neighbors', select_neighbors)
cmd.extend('show_polar_contacts', show_polar_contacts)
cmd.extend('show_hbonds', show_hbonds)
cmd.extend('add_hydrogens', add_hydrogens)
cmd.extend('show_hydrogens', show_hydrogens)
cmd.extend('hide_hydrogens', hide_hydrogens)
cmd.extend('generate_symmetry_mates', generate_symmetry_mates)