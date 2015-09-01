import sys
import re

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
