import matplotlib.pyplot as plt
from matplotlib.widgets import Lasso
from matplotlib import path
import matplotlib
import numpy as np

class Scatter(object):

    """
    Scatter is an optional viewer that plugs into a SelectionHub.
    It displays catalog properties in a scatter plot. 
    Users can select scatter points directly by clicking and dragging
    a lasso around points of interest. These selected points' 
    corresponding structures will then be highlighted in all other 
    viewers.

    Example use:

        >>> from astrodendro.scatter import Scatter
        # ... code to create a dendrogram (d) and catalog ...
        >>> dv = d.viewer()
        >>> ds = Scatter(d, dv.hub, catalog, 'radius', 'v_rms')
        >>> dv.show()

    For more information on using Scatter, see the online
    documentation.

    """

    def __init__(self, dendrogram, hub, catalog, xaxis, yaxis):

        self.hub = hub
        self.hub.add_callback(self.update_selection)
        self.dendrogram = dendrogram
        self.structures = list(self.dendrogram.all_structures)

        self.fig = plt.figure()
        self.axes = plt.subplot(1,1,1)

        self.catalog = catalog
        self.xdata = catalog[xaxis]
        self.ydata = catalog[yaxis]

        self.xys = np.column_stack((self.xdata, self.ydata))

        self.x_column_name = xaxis
        self.y_column_name = yaxis

        self.lines2d = {} # selection_id -> matplotlib.lines.Line2D

        # This is a workaround for a (likely) bug in matplotlib.widgets. Lasso crashes without this fix.
        if matplotlib.get_backend() == 'MacOSX':
            self.fig.canvas.supports_blit = False

        self._draw_plot()
        self.hub.add_callback(self.update_selection)

        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onpress)

    def _draw_plot(self):

        self.axes.plot(self.xdata, self.ydata, 'o', color='w', mec='k', zorder=-5)

        self.axes.set_xlabel(self.x_column_name)
        self.axes.set_ylabel(self.y_column_name)

        self.fig.canvas.draw()

    # This is a closure - we have to pass the input key to callback somehow.
    def callback_generator(self, event):

        input_key = event.button

        def callback(verts):
            p = path.Path(verts)

            indices = np.where(p.contains_points(self.xys))[0]
            selected_structures = [self.dendrogram[i] for i in indices]

            if len(selected_structures) == 0:
                selected_structures = [None]

            self.hub.select(input_key, selected_structures, subtree=False)

            self.fig.canvas.draw_idle()
            del self.lasso

        return callback

    def onpress(self, event):
        if event.canvas.toolbar.mode != '':
            return
        if event.inaxes is None: 
            return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback_generator(event))

    def update_selection(self, selection_id):
        """Highlight seleted structures"""
        
        if selection_id in self.lines2d:
            if self.lines2d[selection_id] is not None:
                self.lines2d[selection_id].remove()
                del self.lines2d[selection_id]

        structures = self.hub.selections[selection_id]
        struct = structures[0]

        if struct is None:
            self.fig.canvas.draw()
            return
        if self.hub.select_subtree[selection_id]:
            selected_indices = [leaf.idx for leaf in struct.descendants + [struct]]
        else:
            selected_indices = [leaf.idx for leaf in structures]

        self.lines2d[selection_id] = self.axes.plot(
            self.xdata[selected_indices], 
            self.ydata[selected_indices], 
            'o', color=self.hub.colors[selection_id], zorder=struct.height)[0]

        self.fig.canvas.draw()
