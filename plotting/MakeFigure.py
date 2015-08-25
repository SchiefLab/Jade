from collections import defaultdict
import os
import sys

import matplotlib.pyplot as plot
import matplotlib as mpl
import matplotlib.axes

from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, WeekdayLocator



class MakeFigure:
    """
    My take on a plotting interface.  Because I think matplotlib's interface sucks.
    Probably should have just used ggplot now that I see its been released...oh well.
    You need to know the number of plots ahead of time by passing the grid.
        1x1 will make one plot.
        2x2 will make a 4x4 grid of plots.
        1x3 is 3 columns of grids horizontally
        3x1 is a list of figures.

        share_x and share_y tell the full sublplot to share the axis.

    """
    def __init__(self, rows = 1, columns = 1, share_x = True, share_y = True):

        self.reset(rows, columns, share_x, share_y)

        self.labels = []
        self.data = defaultdict(dict)
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

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
        self.data[label]['x'] = x
        self.data[label]['y'] = y
        self.labels.append(label)

    def get_x_data(self, label):
        return self.data[label]['x']

    def get_y_data(self, label):
        return self.data[label]['y']

    def get_data(self, label):
        return self.data[label]

    def get_labels(self):
        return self.labels

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

    def fill_subplot(self, title, labels, x_axis_label, y_axis_label,
                     index = None, grid = None,
                     add_legend = True, linestyle="--", marker="^", colors = None):
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

            #linestyle=linestyle, marker=marker,
            p.plot(self.get_x_data(label), self.get_y_data(label), color=color, label = label)


        if add_legend:
            p.legend(labels)

        p.set_xlabel(x_axis_label, fontweight="bold")
        p.set_ylabel(y_axis_label, fontweight="bold")

    def save_plot(self, outpath, tight = True):

        self.p.grid( True )

        print "Saving: "+outpath
        self.p.subplots_adjust(wspace=0.3, hspace=0.3)
        #self.p.autoscale_view()

        if tight:
            self.p.savefig(outpath, dpi=300, bbox_inches='tight')
        else:
            self.p.savefig(outpath, dpi=300, bbox_inches='tight')