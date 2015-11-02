from collections import defaultdict
import os
import sys

import matplotlib.pyplot as plot
import matplotlib as mpl
import matplotlib.axes

from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, WeekdayLocator

# Jared Adolf-Bryfogle


class MakeFigure:
    """
    My take on a plotting interface.  Because I think matplotlib's interface sucks.
    Probably should have just used ggplot now that I see its been released...oh well.
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