import pandas
import sys,os,re, sqlite3, copy
from jade.RAbD_BM.AnalyzeRecovery import *
from jade.basic.pandas.PandasDataFrame import *
from jade.basic.pandas.stats import *
import scipy
import jade.basic.plotting.MakeFigure as plotting

import matplotlib.pyplot as plt
import seaborn.apionly as sea
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict

#Help in calculating risk ratios and plotting them for antibody design paper.

def load_precomputed_recoveries(db_path="data/all_recovery_and_risk_ratio_data.db", table="full_data"):
    """
    Reads recovery data from a database created via script.

    rtype: pandas.Dataframe
    """
    con = sqlite3.connect(db_path)
    # df = pandas.read_sql_table(table, con) Not supported!!! Bah!
    df = pandas.read_sql("select * from " + table, con)
    return df

def order_by_row_group(df, column, groups):
    """
    Order a dataframe by groups.  Return the dataframe.
    Probably a better way to do this already, but I don't know what it is.
    """
    dfs = []
    for group in groups:
        dfs.append(df[df[column] == group])
    out_df = pandas.concat(dfs)
    out_df = out_df.reset_index()
    return out_df

def remove_pdb_and_cdr(df, pdbid, cdr):
    """
    Removes a particular pdbid and cdr from the db. Returns the new df.
    """

    df_new = df[~((df['cdr'] == cdr) & (df['pdbid'] == pdbid.lower()))]
    return df_new


def calculate_geometric_means_rr(df, x, y, hue=None):
    """
    Example use:
    rr_data_lengths = calculate_geometric_means_rr(df_all, x='cdr', y='length_rr', hue='exp')
    rr_data_clusters = calculate_geometric_means_rr(df_all, x='cdr', y='cluster_rr', hue='exp')
    """
    flat_dict = defaultdict(list)

    for x_name in df[x].unique():
        local = df[df[x] == x_name]
        # print local
        logs = []
        for index, row in local.iterrows():
            # print row
            rr = row[y]
            if pandas.isnull(rr):
                continue
            rtype = y.split("_")[0]

            rec = row[rtype + "_recovery"]
            obs = row['native_' + rtype + 's_observed'] / float(row['total_grafts'])
            # print "Rec: "+repr(rec)
            # print "Obs: "+repr(obs)


            if obs == 0:
                continue
            elif rec == 0:
                l = math.log((.01 / obs + .01))
            else:
                l = math.log(rr)

            logs.append(l)
            # print "RR: "+repr(rr)+" lnRR: "+repr(l)

        #print repr(logs)
        m = numpy.array(logs).mean()
        # print "Mean "+repr(m)
        geometric_mean = math.exp(m)

        flat_dict[x].append(x_name)
        flat_dict['N'].append(len(logs))
        flat_dict[y].append(geometric_mean)
        flat_dict['raw_rr'].append(str(logs))
        if hue:
            flat_dict[hue].append('ALL')

            for hue_name in df[hue].unique():
                # print x_name+" "+hue_name
                local2 = local[df[hue] == hue_name]

                # print type(local2)
                logs = []

                # print local2
                for index, row in local2.iterrows():
                    # print row
                    rr = row[y]
                    if pandas.isnull(rr):
                        continue

                    rec = row[rtype + "_recovery"]
                    obs = row['native_' + rtype + 's_observed'] / float(row['total_grafts'])
                    # print "Rec: "+repr(rec)
                    # print "Obs: "+repr(obs)


                    if obs == 0:
                        continue;
                    elif rec == 0:
                        l = math.log((0.01 / obs + 0.01))
                    else:
                        l = math.log(rr)

                    logs.append(l)
                    # print "RR: "+repr(rr)+" lnRR: "+repr(l)

                # print repr(logs)
                m = numpy.array(logs).mean()
                # print "Mean "+repr(m)
                geometric_mean = math.exp(m)

                flat_dict['N'].append(len(logs))
                flat_dict[x].append(x_name)
                flat_dict[hue].append(hue_name)
                flat_dict[y].append(geometric_mean)
                flat_dict['raw_rr'].append(str(logs))

    # Calculate Hue overall SDs.
    if hue:
        for x_name in df[hue].unique():
            local = df[df[hue] == x_name]

            logs = []
            for index, row in local.iterrows():

                rr = row[y]
                if pandas.isnull(rr):
                    continue

                rec = row[rtype + "_recovery"]
                obs = row['native_' + rtype + 's_observed'] / float(row['total_grafts'])
                # print "Rec: "+repr(rec)
                # print "Obs: "+repr(obs)


                if obs == 0:
                    continue
                elif rec == 0:
                    l = math.log((.01 / obs + .01))
                else:
                    l = math.log(rr)

                # print "RR: "+repr(rr)+" lnRR: "+repr(l)
                logs.append(l)
            # print repr(logs)
            m = numpy.array(logs).mean()
            # print "Mean "+repr(m)
            geometric_mean = math.exp(m)

            flat_dict[hue].append(x_name)
            flat_dict[y].append(geometric_mean)
            flat_dict['N'].append(len(logs))
            flat_dict['raw_rr'].append(str(logs))

            if hue:
                flat_dict[x].append('ALL')

    # print repr(flat_dict)
    means = pandas.DataFrame.from_dict(flat_dict)
    # print stddev_df.tail()
    return means

def calculate_stddev_binomial_distribution2(df, x, y, total_column, y_mean_column, hue=None, percent=True):
    """
    Calcuates stddeviations for a binomial distribution.  Returns a dataframe of stddevs
    If percent=True, we dived by the total to normalize the standard deviation.
    SD of 'mean' = SQRT(n*p*q) where p is probability of success and q is probability of failure.
    """

    # Because these are percent, and we don't have 100 total_grafts, we need to devide to get the ratio for the stddevs.
    # Right?

    flat_dict = defaultdict(list)

    for x_name in df[x].unique():
        local = df[df[x] == x_name]
        total = local[total_column].sum()
        p = local[y].sum() / float(total)

        if percent:
            dev = math.sqrt(total * p * (1.0 - p)) / total * 100
        else:
            dev = math.sqrt(total * p * (1.0 - p))

        flat_dict[x].append(x_name)
        flat_dict['SD'].append(dev)
        flat_dict['y'].append(y_mean_column)
        flat_dict['total'].append(total)
        flat_dict['p'].append(p)
        if hue:
            flat_dict[hue].append('ALL')

            for hue_name in df[hue].unique():
                # print x_name+" "+hue_name
                local2 = local[df[hue] == hue_name]
                mean = local2[y].mean()
                total = local2[total_column].sum()
                p = local2[y].sum() / float(total)
                # print x_name+" "+hue_name+" "+repr(mean)+" "+repr(total)
                if percent:
                    dev = math.sqrt(total * p * (1.0 - p)) / total * 100
                else:
                    dev = math.sqrt(total * p * (1.0 - p))

                flat_dict[x].append(x_name)
                flat_dict[hue].append(hue_name)
                flat_dict['SD'].append(dev)
                flat_dict['y'].append(y_mean_column)
                flat_dict['total'].append(total)
                flat_dict['p'].append(p)

    # Calculate Hue overall SDs.
    if hue:
        for x_name in df[hue].unique():
            local = df[df[hue] == x_name]

            total = local[total_column].sum()
            p = local[y].sum() / float(total)

            if percent:
                dev = math.sqrt(total * p * (1.0 - p)) / total * 100
            else:
                dev = math.sqrt(total * p * (1.0 - p))

            flat_dict[hue].append(x_name)
            flat_dict['SD'].append(dev)
            flat_dict['y'].append(y_mean_column)
            flat_dict['total'].append(total)
            flat_dict['p'].append(p)
            if hue:
                flat_dict[x].append('ALL')

    # print repr(flat_dict)
    stddev_df = pandas.DataFrame.from_dict(flat_dict)
    # print stddev_df.tail()
    return stddev_df

def calculate_rr_errors(df_all_errors):
    """
    Calculates the risk ratio errors for cluster and lengths using propagation error equations calculated for
    the recovery itself.  Which is the same for percent as it would be raw data, as the N cancels out in the equations.
    http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm
    """
    df_length_recovered = df_all_errors[df_all_errors['y'] == 'length_recovery_freq'].reset_index()
    df_cluster_recovered = df_all_errors[df_all_errors['y'] == 'cluster_recovery_freq'].reset_index()
    df_cluster_observed = df_all_errors[df_all_errors['y'] == 'cluster_observed_perc'].reset_index()
    df_length_observed = df_all_errors[df_all_errors['y'] == 'length_observed_perc'].reset_index()

    # print repr(df_length_recovered)
    # print df_length_recovered
    # print df_cluster_observed
    length_rr = df_length_recovered['p'] / \
                df_length_observed['p']

    cluster_rr = df_cluster_recovered['p'] / \
                 df_cluster_observed['p']

    # print length_rr
    # print repr(length_rr)


    # Clusters
    a = ((1 - df_cluster_recovered['p']) / (df_cluster_recovered['total'] * df_cluster_recovered['p']))
    b = ((1 - df_cluster_observed['p']) / (df_cluster_observed['total'] * df_cluster_observed['p']))

    se_log = numpy.sqrt(a + b)
    df_cluster_recovered['se_log'] = se_log

    # Lengths
    a = ((1 - df_length_recovered['p']) / (df_length_recovered['total'] * df_length_recovered['p']))
    b = ((1 - df_length_observed['p']) / (df_length_observed['total'] * df_length_observed['p']))

    se_log = numpy.sqrt(a + b)
    df_length_recovered['se_log'] = se_log

    # i = ((df_cluster_recovered['SD']/100)/df_cluster_recovered['p'])**2
    # ii = ((df_cluster_observed['SD']/100)/df_cluster_observed['p'])**2
    # df_cluster_recovered['SD'] = cluster_rr * numpy.sqrt( i + ii)

    # x = ((df_length_recovered['SD']/100)/df_length_recovered['p'])**2
    # xx = ((df_length_observed['SD']/100)/df_length_observed['p'])**2
    # df_length_recovered['SD'] = length_rr * numpy.sqrt( x + xx)


    # Mutate to a new dataframe.
    df_errors_length = df_length_recovered
    df_errors_cluster = df_cluster_recovered

    df_errors_length['y'] = 'length_rr'
    df_errors_cluster['y'] = 'cluster_rr'

    df_errors_length['p'] = length_rr
    df_errors_cluster['p'] = cluster_rr

    df_errors = pandas.concat([df_errors_length, df_errors_cluster])
    return df_errors

def set_errorbars_bar(ax, data, x, y, error_dfs,
                      x_order=None, hue_order=None,
                      hue=None, caps=False, color='k', linewidth=.75, base_columnwidth=.8, full=True):
    """
    Sets erorr bars for a bar chart.

    Default base_columnwidth for seaborn plots is .8

    Optionally give x_order and/or hue_order for the ordering of the columns.  Make sure to pass this while plotting.

    Notes:
     1) If Hue is enabled, this base is divided by the number of hue_names for the final width used for plotting.

     2) Caps are the line horizontal lines in the errorbar.

     3) 'full' means error bars on both vertical sides of the histogram bar.

    Warning:
      linewidth of .5 does not show up in all PDFs for all bars.

    """
    print x + " " + y + " " + repr(hue)

    def get_sd(errors, x_name, hue_name=None):
        if hue:
            return errors[errors[x] == x_name][errors[hue] == hue_name][errors['y'] == y].iloc[0]['SD']
        else:
            return errors[errors[x] == x_name][errors['y'] == y].iloc[0]['SD']

    def get_mean(x_name, hue_name=None):
        if hue:
            # print "WTF?" + repr(data[data[x] == x_name][data[hue] == hue_name][y])
            f = data[data[x] == x_name][data[hue] == hue_name][y]
            m = sum(float(embedding) for embedding in f) / len(f)
            p = error_dfs[error_dfs[x] == x_name][error_dfs[hue] == hue_name][error_dfs['y'] == y].iloc[0]['p']
            # return data[data[x] == x_name][data[hue] == hue_name][y], dtype=float).mean()
            # return data[data[x] == x_name][data[hue] == hue_name][y].mean()

            print 'MEAN: ' + repr(m) + " p: " + repr(p)
            return m
        else:
            m = data[data[x] == x_name][y].mean()
            p = error_dfs[error_dfs[x] == x_name][error_dfs['y'] == y].iloc[0]['p']
            print 'MEAN:' + repr(m) + " p: " + repr(p)
            return m

    x_indexes = []
    y_means = []
    yerr = []
    zeros = []

    x_names = sea.utils.categorical_order(data[x], x_order)

    # Start the coordinates at Zero, then minus i+base/2 from everything to get it centered at each i.
    if hue:
        hue_names = sea.utils.categorical_order(data[hue], hue_order)

        w = (base_columnwidth / float(len(hue_names)))
        base_w = base_columnwidth

        ##Check Ordering here, should be good?
        for i, x_name in enumerate(x_names):

            for z, hue_name in enumerate(hue_names):
                print x_name + " " + hue_name

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

    print repr(y_means)
    print repr(yerr)
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

def set_errorbars_bar_rr(ax, data, x, y, error_dfs,
                         x_order=None, hue_order=None,
                         hue=None, caps=False, color='k', linewidth=.75, base_columnwidth=.8, full=True):
    """
    Sets erorr bars for a bar chart.

    Default base_columnwidth for seaborn plots is .8

    Optionally give x_order and/or hue_order for the ordering of the columns.  Make sure to pass this while plotting.

    Notes:
     1) If Hue is enabled, this base is divided by the number of hue_names for the final width used for plotting.

     2) Caps are the line horizontal lines in the errorbar.

     3) 'full' means error bars on both vertical sides of the histogram bar.

    Warning:
      linewidth of .5 does not show up in all PDFs for all bars.

    """
    print x + " " + y + " " + repr(hue)

    def get_sd(errors, x_name, hue_name=None):
        if hue:

            se_log = errors[errors[x] == x_name][errors[hue] == hue_name][errors['y'] == y].iloc[0]['se_log']
            m = get_mean(x_name, hue_name)
            ci_log_max = math.log(m) + (1.96 * se_log)
            ci_max = math.exp(ci_log_max)
            error = m - ci_max
            return error

        else:
            se_log = errors[errors[x] == x_name][errors['y'] == y].iloc[0]['se_log']
            m = get_mean(x_name, hue_name)
            ci_log_max = math.log(m) + (1.96 * se_log)
            ci_max = math.exp(ci_log_max)
            error = m - ci_max
            return error

    def get_mean(x_name, hue_name=None):
        if hue:
            # print "WTF?" + repr(data[data[x] == x_name][data[hue] == hue_name][y])
            f = data[data[x] == x_name][data[hue] == hue_name][y].iloc[0]

            print "MEAN: " + repr(f)

            return f
        else:
            # m = data[data[x] == x_name][y].mean()
            f = data[data[x] == x_name][y].iloc[0]
            print "MEAN: " + repr(f)
            return f

    x_indexes = []
    y_means = []
    yerr = []
    zeros = []

    x_names = sea.utils.categorical_order(data[x], x_order)

    # Start the coordinates at Zero, then minus i+base/2 from everything to get it centered at each i.
    if hue:
        hue_names = sea.utils.categorical_order(data[hue], hue_order)

        w = (base_columnwidth / float(len(hue_names)))
        base_w = base_columnwidth

        ##Check Ordering here, should be good?
        for i, x_name in enumerate(x_names):

            for z, hue_name in enumerate(hue_names):
                print x_name + " " + hue_name

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

    print repr(y_means)
    print repr(yerr)
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

def plot_rr(data, x, y, hue=None, ci=None):
    if not hue:

        if x == 'exp':
            h = 'cdr'
        else:
            h = 'exp'

        if h in data.columns:
            data2 = data[data[h] == 'ALL']

    else:

        data2 = data[data[hue] != 'ALL']
        data2 = data2[data[x] != 'ALL']

    print data2
    ax = sea.barplot(data=data2, x=x, y=y, hue=hue, ci=None)
    set_errorbars_bar_rr(ax, data2, x, y, df_stddev_rr, hue=hue)
    return ax

def calculate_set_errorbars_hist(ax, data, x, y,
                                 binomial_distro=True, total_column='total_entries', y_freq_column=None,
                                 x_order=None, hue_order=None,
                                 hue=None, caps=False, color='k', linewidth=.75, base_columnwidth=.8, full=True):
    """
    Calculates the standard deviation of the data, sets erorr bars for a histogram.
    Default base_columnwidth for seaborn plots is .8

    Optionally give x_order and/or hue_order for the ordering of the columns.  Make sure to pass this while plotting.

    Notes:
     1) If Hue is enabled, this base is divided by the number of hue_names for the final width used for plotting.

     2) Caps are the line horizontal lines in the errorbar.

     3) 'full' means error bars on both vertical sides of the histogram bar.

    Warning:
      linewidth of .5 does not show up in all PDFs for all bars.

    """

    # This makes it easier for frequencies of x/100, instead of passing two columns
    #   - one to calc mean (y), and one for freq.
    if not y_freq_column:
        y_freq_column = y

    if binomial_distro:
        error_dfs = calculate_stddev_binomial_distribution2(data, x, y_freq_column, total_column, y, hue)
    else:
        error_dfs = calculate_stddev(data, x, y, hue)

    if not hue and hue in error_dfs.columns:
        error_dfs = error_dfs[error_dfs[hue] == 'ALL']

    set_errorbars_bar(ax, data, x, y, error_dfs, x_order=x_order, hue_order=hue_order,
                      hue=hue, caps=caps, color=color, linewidth=linewidth, base_columnwidth=base_columnwidth,
                      full=full)

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


