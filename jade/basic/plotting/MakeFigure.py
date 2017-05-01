from collections import defaultdict
import os
import sys

from jade.basic.pandas.PandasDataFrame import *
import seaborn.apionly as sea

import pandas
import matplotlib.pyplot as plot
import matplotlib as mpl
import matplotlib.axes

from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, WeekdayLocator

# Jared Adolf-Bryfogle




'''
NOTES:

Matplotlib class organization generally makes no overall sense, so here is some notes (to myself) to help:

1) Figures are made of axes.
2) subplot is just an axes on a grid!  Hahahahahah
 -> 'Axes are very similar to subplots but allow placement of plots at any location in the figure.
     So if we want to put a smaller plot inside a bigger one we do so with axes.'

3) axes are the actual fucking PLOTS. They should be called PLOTS. Why it is called this is beyond me.

You can get the figure of an axes, and the [axes] of a figure.

fig.axes
axes.figure

Boom.  Like an access pointer in C++ I guess.

Almost everything you need can and should be controlled by individual axes (USING ACTUAL SETTERS!).
In order to give common x or y labels, or titles, see the functions below.

Passing x= or y= to labels like titles, will control the offset.
x=0 and y=1.0 are the default.

y=1.05 will raise it up a bit - just enough to make the titles look great.
Title size of 18 and x/y label size of 16 seems to look ok.

'''



def plot_x_vs_y_sea_with_regression(df, title, outpath, x, y, top_p = .95, reverse = True):
    """
    Plot X vs Y using a Pandas Dataframe and Seaborn, with regression line., save the figure, and return the Axes.

    If you are doing this multiple times in a Notebook:
        Don't forget to call (matplotlib.pyplot)
             plot.show()
             plot.close()

    :param df: pandas.DataFrame
    :param title: str
    :param outpath: str
    :param x: str
    :param y: str
    :param top_p: float
    :param reverse:  bool
    :rtype: matplotlib.Axes
    """
    df = detect_numeric(df)
    df = df.sort(x, ascending=reverse)
    slice_top = df[0:int(len(df)*float(top_p))]
    ax = sea.regplot(x=x, y=y, data=slice_top)
    ax.set_title(title)
    pad_single_title(ax)
    fig = ax.get_figure()
    fig.savefig(outpath, dpi=300)

    return ax

def plot_general_pandas(df, title, outpath, plot_type, x, y = None, z = None, top_p = .95, reverse = True):
    """
    Plot anything in pandas.  Make it look descent.  Save the figure.

    If you are doing this multiple times in a Notebook:
        Don't forget to call (matplotlib.pyplot)
             plot.show()
             plot.close()


    :param df: pandas.DataFrame
    :param title: str
    :param outpath: str
    :param plot_type: str
    :param x: str
    :param y: str
    :param z: str
    :param top_p: float
    :param reverse: bool
    :rtype: matplotlib.Axes
    """
    print "Plotting "+plot_type+" For X: "+x +" and y: "+repr(y)

    df = detect_numeric(df)
    df = df.sort(x, ascending=reverse)
    slice_top = df[0:int(len(df)*float(top_p))]

    if y:
        ax = slice_top.plot(x, y,  kind=plot_type, figsize=[11, 8], title = title)
    elif plot_type=="kde":
        #Y=X is not an error.  For KDE plots, pandas uses y instead of X.  Super stupid.
        ax = slice_top.plot(y=x, kind=plot_type, figsize=[11, 8], title = title)
    else:
        ax = slice_top.plot(x=x, kind=plot_type, figsize=[11, 8], title = title)
    pad_single_title(ax)
    ax.set_axis_bgcolor('white')
    fig = ax.get_figure()
    fig.savefig(outpath, dpi=300)
    return ax


def pad_single_title(ax, x=.5, y=1.05):
    """
    Move the Title up in reference to the plot, essentially adding padding.
    SINGLE AXES
    :param ax:Axes
    :param x:
    :param y:
    :return:
    """
    ttl = ax.title
    ttl.set_position([x, y])

def set_common_title(fig, title, size=16, x=0, y=1.05):
    """
    for FACETED plots, add a common title.

    :param fig: Figure
    :param title: str
    :param x: int
    :param y: int
    :return:
    """

    fig.suptitle(title, size=size, x = x, y= y)

def set_common_x_y_label(fig, x_text, y_text):
    """
    For FACETED plots, add a common X or Y.

    :param fig: Figure
    :param x_text: str
    :param y_text: str
    :return:
    """
    for ax in fig.axes:
        ax.set_xlabel("")
        ax.set_ylabel("")

    fig.text(0.5, -.04, x_text, ha='center', fontsize=18)
    fig.text(-.04, 0.5, y_text, va='center', rotation='vertical', fontsize=18)





class MakeFigure:
    """
    Deprecated.
    NOW - GO Checkout SEABORN instead of this class!
    Essentially, this is an interface to a facet grid.  Seaborn does this awesomely.

    My take on a plotting interface.  Because I think matplotlib's interface sucks.

    I wrote this before I knew of pandas.

    You need to know the number of plots ahead of time by passing the grid.
        1x1 will make one plot.
        2x2 will make a grid of 4 plots.
        1x3 is 3 columns of grids horizontally
        3x1 is a list of figures.

        share_x and share_y tell the full sublplot to share the axis.

    """
    def __init__(self, rows = 1, columns = 1, share_x = True, share_y = True):

        self.reset(rows, columns, share_x, share_y)

        self.labels = []
        self.data = defaultdict(dict)
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        self.linestyles = ['_', '-', '--', ':']

    def reset(self, rows = 1, columns = 1, share_x = True, share_y = True):
        p, axes = plot.subplots(rows, columns, share_x, share_y)

        try:
            len(axes)
        except TypeError:
            axes = [axes]

        self.p = p
        self.axes = axes
        self.current_plot = 0

    def set_data(self, data_dict):
        self.data = data_dict

    def merge_data(self, data_dict, replace_current_labels = False):
        for label in data_dict:
            if self.data.has_key(label):

                if replace_current_labels:
                    self.data[label] = data_dict[label]
                else:
                    self.data[label].extend(data_dict[label]['x'])
                    self.data[label].extend(data_dict[label]['y'])

            if not label in self.labels:
                self.labels.append(label)

    def add_data(self, x, y, label):
        print label
        self.data[label]['x'] = x
        self.data[label]['y'] = y
        if not self.labels.count(label):
            self.labels.append(label)

    def get_x_data(self, label):
        return self.data[label]['x']

    def get_y_data(self, label):
        return self.data[label]['y']

    def get_data(self, label):
        return self.data[label]

    def get_labels(self):
        return self.labels

    def get_y_as_list(self, labels):
        return [self.data[label]['y'] for label in labels]

    ################################

    def get_plot(self, n =0):
        """

        :param n: int
        :return: mpl.axes.SubplotBase
        """
        return self.axes[n]

    def set_x_limits(self, min, max, plot_num = None):
        if plot_num:
            self.axes[plot_num].set_xlim(min, max)
        else:
            for p in self.axes:
                p.set_xlim(min, max)

    def set_y_limits(self, min, max, plot_num = None):
        if plot_num:
            self.axes[plot_num].set_ylim(min, max)
        else:
            for p in self.axes:
                p.set_ylim(min, max)

    def set_x_scale(self, scale = 'log', plot_num = None):

        #Here, I think we could return functions to generalize everything.
        if plot_num:
            self.axes[plot_num].set_xscale(scale)
        else:
            for p in self.axes:
                p.set_sclale(scale)

    def set_y_scale(self, scale = 'log', plot_num = None):
        if plot_num:
            self.axes[plot_num].set_yscale(scale)
        else:
            for p in self.axes:
                p.set_yscale(scale)

    def fill_subplot(self, title, labels, x_axis_label = None, y_axis_label = None,
                     index = None, grid = None,
                     add_legend = False, linestyle="--", marker="^", colors = None):
        """
        This will add data to a particular subplot/plot.

        : title:
        : labels:
        : x_axis_label:
        : y_axis_label:
        : specify_index:
        : add_legend:
        : linestyle:
        : marker:
        : colors:
        :return:
        """


        if colors and len(colors) != len(labels): sys.exit("Labels must match colors length")


        if index:
            p = self.get_plot(index)

        else:
            p = self.get_plot(self.current_plot)
            self.current_plot+=1

        assert isinstance(p, mpl.axes.SubplotBase)

        p.set_title(title)
        i = 0
        for label in labels:

            if colors:
                color = colors[i]
            else:
                color = self.colors[i]
            #i +=1
            #linestyle=linestyle, marker=marker,
            p.plot(self.get_x_data(label), self.get_y_data(label), label = label)


        if add_legend:
            p.legend(labels, loc="upper left",  fancybox=True)

        if x_axis_label:
            p.set_xlabel(x_axis_label, fontweight="bold")

        if y_axis_label:
            p.set_ylabel(y_axis_label, fontweight="bold")

    def add_grid(self, x_grid = True, y_grid = True):
        for p in self.axes:
            if x_grid:
                p.xaxis.grid()
            if y_grid:
                p.yaxis.grid()

    def save_plot(self, outpath, tight = True):

        print "Saving: "+outpath
        self.p.subplots_adjust(wspace=0.2, hspace=0.2)
        #self.p.autoscale_view()

        if tight:
            self.p.savefig(outpath, dpi=300, bbox_inches='tight')
        else:
            self.p.savefig(outpath, dpi=300, bbox_inches='tight')

