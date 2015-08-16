import os
import re
import

#A Collection of pathing functions


def make_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_file_names(pattern, dir):
    pass
