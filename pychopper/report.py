import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Report:

    # Maybe it would be a good idea to convert these utility methods to object
    # oriented matplotlib style.

    def __init__(self, pdf):
        """Class for plotting utilities on the top of matplotlib. Plots are saved in the specified file through the PDF backend.

        :param self: object.
        :param pdf: Output pdf.
        :returns: The report object.
        :rtype: Report

        """
        self.pdf = pdf
        self.plt = plt
        self.pages = PdfPages(pdf)

    def _set_properties_and_close(self, fig, title, xlab, ylab):
        """Utility method to set title, axis labels and close the figure.

        :param self: object.
        :param fig: The current figure.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :returns: None
        :rtype: object
        """
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.close(fig)

    def save_close(self):
        """Utility method to save and close figure."""
        self.pages.savefig()
        plt.close()
        plt.clf()

    def plot_arrays(self, data_map, title="", xlab="", ylab="", marker='.', legend_loc='best', legend=True, vlines=None, vlcolor='green', vlwitdh=0.5):
        """Plot multiple pairs of data arrays.

        :param self: object.
        :param data_map: A dictionary with labels as keys and tupples of data arrays (x,y) as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param marker: Marker passed to the plot function.
        :param legend_loc: Location of legend.
        :param legend: Plot legend if True
        :param vlines: Dictionary with labels and positions of vertical lines to draw.
        :param vlcolor: Color of vertical lines drawn.
        :param vlwidth: Width of vertical lines drawn.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, data_arrays in data_map.items():
            plt.plot(data_arrays[0], data_arrays[1], marker, label=label)
        if vlines is not None:
            for label, pos in vlines.items():
                plt.axvline(x=pos, label=label, color=vlcolor, lw=vlwitdh)
        if legend:
            plt.legend(loc=legend_loc)
        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_histograms(self, data_map, title="", xlab="", ylab="", bins=50, alpha=0.7, legend_loc='best', legend=True, vlines=None):
        """Plot histograms of multiple data arrays.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data arrays as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param bins: Number of bins.
        :param alpha: Transparency value for histograms.
        :param legend_loc: Location of legend.
        :param legend: Plot legend if True.
        :param vlines: Dictionary with labels and positions of vertical lines to draw.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, data in data_map.items():
            if len(data) > 0:
                plt.hist(data, bins=bins, label=label, alpha=alpha)
        if vlines is not None:
            for label, pos in vlines.items():
                plt.axvline(x=pos, label=label)
        if legend:
            plt.legend(loc=legend_loc)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_bars_simple(self, data_map, title="", xlab="", ylab="", alpha=0.6, xticks_rotation=0, auto_limit=False):
        """Plot simple bar chart from input dictionary.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param alpha: Alpha value.
        :param xticks_rotation: Rotation value for x tick labels.
        :param auto_limit: Set y axis limits automatically.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        labels = list(data_map.keys())
        data = list(data_map.values())
        positions = np.arange(len(labels))
        plt.bar(positions, data, align='center', alpha=alpha)
        plt.xticks(positions, labels, rotation=xticks_rotation)

        if auto_limit:
            low, high = min(data), max(data)
            plt.ylim([(low - 0.5 * (high - low)), (high + 0.5 * (high - low))])

        self._set_properties_and_close(fig, title, xlab, ylab)

    def close(self):
        """Close PDF backend. Do not forget to call this at the end of your script or your output will be damaged!

        :param self: object
        :returns: None
        :rtype: object
        """
        self.pages.close()
