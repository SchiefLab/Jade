import ast



def deduce_str_type(s):
    """
    Deduce the type of a string.  Either return the string as the literal, or as the string if not possible.
    http://stackoverflow.com/questions/13582142/deduce-the-type-of-data-in-a-string

    :param s: str
    :return:
    """
    try:
        return ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return s