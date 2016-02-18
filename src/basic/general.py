import sys
import re

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
