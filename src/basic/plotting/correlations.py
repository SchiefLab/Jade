import numpy
import pandas
import matplotlib
from scipy import stats


def annotate_r_value(data, x, y, ax, func=stats.pearsonr, template=None, stat=None, loc="best"):
    """
    Forked from seaborn JointPlot for use with regplot, scatter, etc.  Woot.  Needs to actually go into Seaborn Now!

    Annotate the plot with a statistic about the relationship.

    Parameters
    ----------

    data: pandas.DataFrame
    x: str
    y: str
    ax: matplotlib.Axes

    func : callable
        Statistical function that maps the x, y vectors either to (val, p)
        or to val.
    template : string format template, optional
        The template must have the format keys "stat" and "val";
        if `func` returns a p value, it should also have the key "p".
    stat : string, optional
        Name to use for the statistic in the annotation, by default it
        uses the name of `func`.
    loc : string or int, optional
        Matplotlib legend location code; used to place the annotation.


    """
    # Possibly extract the variables from a DataFrame
    if data is not None:
        if x in data:
            x = data[x]
        if y in data:
            y = data[y]

    # Convert the x and y data to arrays for plotting
    x = numpy.asarray(x)
    y = numpy.asarray(y)

    default_template = "{stat} = {val:.2g}; p = {p:.2g}"

    # Call the function and determine the form of the return value(s)
    out = func(x, y)
    try:
        val, p = out
    except TypeError:
        val, p = out, None
        default_template, _ = default_template.split(";")

    # Set the default template
    if template is None:
        template = default_template

    # Default to name of the function
    if stat is None:
        stat = func.__name__

    # Format the annotation
    if p is None:
        annotation = template.format(stat=stat, val=val)
    else:
        annotation = template.format(stat=stat, val=val, p=p)

    # Draw an invisible plot and use the legend to draw the annotation
    # This is a bit of a hack, but `loc=best` works nicely and is not
    # easily abstracted.
    phantom, = ax.plot(x, y, linestyle="", alpha=0)
    ax.legend([phantom], [annotation], loc=loc)
    phantom.remove()