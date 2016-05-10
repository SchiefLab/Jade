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

def fix_input_args():
    """
    Enables options to be passed to ArgumentParser with dashes, but not single charactor ones.
    Example:
      --rosetta_args "-out:prefix test -out:path:all my/dir/"

    Normally, this would fail if you had declared an -o option to the ArgumentParser.
      This happens because although the quotes are being parsed correctly, the system is looking or options using the
      starting '-' charactor.  If you give a quote and then a space, you will recieve no error.

    This code essentially checks for single dashes and puts a space in front of them.  Note that this does not work with single
     charactor options you are hoping to pass with a quote.  Because there is no way to grab the input string from the system and fix it myself,
     for these it will have to have a space after the quotes.  This at least fixes the most common use cases (Mostly for use with Rosetta.).


    """
    new_argv = []
    for arg in sys.argv:

        if arg[0:2] == "--":
            new_argv.append(arg)
        elif arg[0] == '-' and len(arg) > 2 and arg[2] != ' ':

            new_arg = ' '+arg
            new_argv.append(new_arg)
        else:
            new_argv.append(arg)

    sys.argv = new_argv