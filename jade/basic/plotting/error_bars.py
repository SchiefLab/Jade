import seaborn as sea
import pandas
import matplotlib as mpl
import numpy

from jade.basic.pandas.stats import calculate_stddev_binomial_distribution
from jade.basic.pandas.stats import calculate_stddev



def calculate_set_errorbars_hist(ax, data, x, y,
                                 binomial_distro = True, total_column = 'total_entries', y_freq_column = None,
                                 x_order = None, hue_order = None,
                                 hue=None, caps=True, color='k', linewidth=.75, base_columnwidth=.8, full = True):
    """
    Calculates the standard deviation of the data, sets erorr bars for a bar chart.
    Default base_columnwidth for seaborn plots is .8

    Optionally give x_order and/or hue_order for the ordering of the columns.  Make sure to pass this while plotting.
    Note:
     If Hue is enabled, this base is divided by the number of hue_names for the final width used for plotting.

    :param ax: mpl.Axes
    :param data: pandas.DataFrame
    :param x: str
    :param y: str
    :param binomial_distro: bool
    :param total_column: str
    :param y_freq_column: str
    :param x_order: list
    :param hue_order: list
    :param hue: str
    :param caps: bool
    :param color: str
    :param linewidth: float
    :param base_columnwidth:  float
    :param full: bool
    :rtype: None
    """

    def get_sd(errors, x_name, hue_name = None):
        if hue:
            return errors[errors[x] == x_name][errors[hue] == hue_name][errors['y'] == y].iloc[0]['SD']
        else:
            return errors[errors[x] == x_name][errors['y'] == y].iloc[0]['SD']

    def get_mean(x_name, hue_name = None):
        if hue:
            # print "WTF?" + repr(data[data[x] == x_name][data[hue] == hue_name][y])
            f = data[data[x] == x_name][data[hue] == hue_name][y]
            return sum(float(embedding) for embedding in f) / len(f)
            # return data[data[x] == x_name][data[hue] == hue_name][y], dtype=float).mean()
            # return data[data[x] == x_name][data[hue] == hue_name][y].mean()
        else:
            return data[data[x] == x_name][y].mean()

    # This makes it easier for frequencies of x/100, instead of passing two columns
    #   - one to calc mean (y), and one for freq.


    if not y_freq_column:
        y_freq_column = y

    if binomial_distro:
        error_dfs = calculate_stddev_binomial_distribution(data, x, y_freq_column, total_column, y, hue)
    else:
        error_dfs = calculate_stddev(data, x, y, hue)

    if not hue and hue in error_dfs.columns:

        error_dfs = error_dfs[error_dfs[hue] == 'ALL']

    # Make sure ALL is not plotted.

    # Need X columns, and y data to plot

    # Need to only plot upper bars
    # (_, caps1, _) = ax.errorbar(x = [(.2)], y=[.25], yerr = [(0,), (.5,)], ls = 'None', capsize=5, color ='k', lw=1)


    x_indexes = []
    y_means = []
    yerr = []
    zeros = []

    x_names = sea.utils.categorical_order(data[x], x_order)

    # Start the coordinates at Zero, then minus i+base/2 from everything to get it centered at each i.
    if hue:
        hue_names = sea.utils.categorical_order(data[hue], hue_order)

        w = ( base_columnwidth /float(len(hue_names)))
        base_w = base_columnwidth

        ##Check Ordering here, should be good?
        for i, x_name in enumerate(x_names):

            for z, hue_name in enumerate(hue_names):
                print x_name + " "+ hue_name

                index = (w * z) + (w / 2) - base_w / float(2) + i
                x_indexes.append(index)
                # print "Index: "+repr(index)
                print "SD: " + repr(get_sd(error_dfs, x_name, hue_name))
                # print "mean: " + repr(get_mean(x_name, hue_name))
                yerr.append(get_sd(error_dfs, x_name, hue_name))
                y_means.append(get_mean(x_name, hue_name))


    else:
        for i, x_name in enumerate(x_names):
            x_indexes.append(i)
            yerr.append(get_sd(error_dfs, x_name))
            y_means.append(get_mean(x_name))

    zeros = list(numpy.zeros(len(x_indexes)))

    if full:
        (_, caps_list, _) = ax.errorbar(x=x_indexes, y=y_means, yerr=yerr,
                                        ls='None', capsize=5, color=color, lw=linewidth)
    else:
        (_, caps_list, _) = ax.errorbar(x=x_indexes, y=y_means, yerr=[zeros, yerr],
                                        ls='None', capsize=5, color=color, lw=linewidth)
    if caps:
        for cap in caps_list:
            cap.set_linewidth(linewidth)
            cap.set_markeredgewidth(linewidth)


def calculate_set_errorbars_scatter(ax, data, x, y,
                                    binomial_distro=False, total_column='total_entries',
                                    caps=False, color='k', lw=1.5):
    """
    (Untested) - Calculates the standard deviation of the data, sets error bars for a typical scatter plot
    """

    if binomial_distro:
        error_dfs = calculate_stddev_binomial_distribution(data, x, y, total_column, hue=None)
    else:
        error_dfs = calculate_stddev(data, x, y, hue=None)