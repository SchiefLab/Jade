from jade.basic.path import *


def get_common_flags_string_for_init(flags_name = "common_flags.flags"):
    """
    Get a string of common flags as specified in the database.
    :return: str
    """

    return " ".join([ line.strip() for line in open(get_rosetta_flags_path()+'/'+flags_name, 'r') if line and not line.startswith('#')])
