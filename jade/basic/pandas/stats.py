import math
from collections import defaultdict
import pandas




def calculate_stddev(df, x, y, hue=None):
    """

    Calcuates standard deviations for a normal distribution (Numerical data) over X and Hue categories.

    If hue is given, the hue column will be added, and the overall will be of 'ALL'

    Example DataFrame output (x='exp', y= 'length_recovery_freq', hue = 'cdr':

                  SD cdr                     exp                     y
        20  6.739596  H2                     ALL  length_recovery_freq
        21  7.373650  H2    min.remove_antigen-F  length_recovery_freq
        22  6.400637  ALL   min.remove_antigen-T  length_recovery_freq


    :param df: pandas.DataFrame
    :param x: str
    :param y: str
    :param total_column: str
    :param hue: str
    :rtype: pandas.DataFrame
    """

    flat_dict = defaultdict(list)

    for x_name in df[x].unique():
        local = df[df[x] == x_name]
        flat_dict[x].append(x_name)
        flat_dict['SD'].append(local[y].std())
        flat_dict['y'].append(y)

        if hue:
            flat_dict[hue].append('ALL')

            for hue_name in df[hue].unique():
                local2 = local[df[hue] == hue_name]

                flat_dict[x].append(x_name)
                flat_dict[hue].append(hue_name)
                flat_dict['SD'].append(local2[y].std())
                flat_dict['y'].append(y)

    # Calculate Hue overall SDs.
    if hue:
        for x_name in df[hue].unique():
            local = df[df[hue] == x_name]

            flat_dict[hue].append(x_name)
            flat_dict['SD'].append(local[y].std())
            flat_dict['y'].append(y)
            if hue:
                flat_dict[x].append('ALL')


    # print repr(flat_dict)
    stddev_df = pandas.DataFrame.from_dict(flat_dict)
    return stddev_df

def calculate_stddev_binomial_distribution(df, x, y, total_column, y_mean_column, hue=None):
    """

    Calculates standard deviations for a binomial distribution (like experiment True/False values) over X and Hue categories..

    Typically used for bar-plot.

    If hue is given the hue column will be added, and the overall will be of 'ALL', plus that of Hue

    Example DataFrame output (x='exp', y= 'length_recovery_freq', hue = 'cdr':

                  SD cdr                     exp                     y
        20  6.739596  H2                     ALL  length_recovery_freq
        21  7.373650  H2    min.remove_antigen-F  length_recovery_freq
        22  6.400637  ALL   min.remove_antigen-T  length_recovery_freq


    :param df: pandas.DataFrame
    :param x: str
    :param y: str
    :param total_column: str
    :param hue: str
    :rtype: pandas.DataFrame
    """

    flat_dict = defaultdict(list)

    for x_name in df[x].unique():
        local = df[df[x] == x_name]
        mean = local[y].mean()
        total = local[total_column].sum()
        dev = math.sqrt(mean * (1 - mean / total * 1.0))

        flat_dict[x].append(x_name)
        flat_dict['SD'].append(dev)
        flat_dict['y'].append(y_mean_column)
        if hue:
            flat_dict[hue].append('ALL')

            for hue_name in df[hue].unique():
                # print x_name+" "+hue_name
                local2 = local[df[hue] == hue_name]
                mean = local2[y].mean()
                total = local2[total_column].sum()
                # print x_name+" "+hue_name+" "+repr(mean)+" "+repr(total)
                dev = math.sqrt(mean * (1 - mean / float(total)))

                flat_dict[x].append(x_name)
                flat_dict[hue].append(hue_name)
                flat_dict['SD'].append(dev)
                flat_dict['y'].append(y_mean_column)

    # Calculate Hue overall SDs.
    if hue:
        for x_name in df[hue].unique():
            local = df[df[hue] == x_name]
            mean = local[y].mean()
            total = local[total_column].sum()
            dev = math.sqrt(mean * (1 - mean / total * 1.0))

            flat_dict[hue].append(x_name)
            flat_dict['SD'].append(dev)
            flat_dict['y'].append(y_mean_column)
            if hue:
                flat_dict[x].append('ALL')

    # print repr(flat_dict)
    stddev_df = pandas.DataFrame.from_dict(flat_dict)
    print stddev_df.tail()
    return stddev_df


