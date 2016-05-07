import sys
import re
import itertools

def get_rosetta_program(program, mpi = True, compiler = 'gcc'):
    """
    Get the set program
    """

    if get_platform() == 'macos':
        compiler = 'clang'

    if mpi:
        return program +".mpi."+get_platform() + compiler+"release"
    else:
        return program +get_platform() + compiler+"release"

def get_platform():
    """
    Get OS of the particular platform the toolkit is being run on.
    """

    plat = sys.platform
    if re.search("darwin", plat):
        return "macos"
    elif re.search("linux", plat):
        return "linux"
    elif re.search("win", plat):
        return "windows"
    else:
        print "Platform Not Found"
        return "error"

def get_all_combos(list_of_lists):
    """
    Get all the position-specific combos of a list of lists.

    This is taken directly from Stack Overflow:
       http://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists

    :param list_of_lists: A list of lists we would like combos of.
    :rtype: list[list]
    """

    return list(itertools.product(*list_of_lists))