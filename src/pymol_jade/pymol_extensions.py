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

cmd.extend('select_neighbors', select_neighbors)