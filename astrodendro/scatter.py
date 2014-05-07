import matplotlib.pyplot as plt
from matplotlib.widgets import Lasso
from matplotlib import path
import matplotlib
import numpy as np

class Scatter(object):
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

        self.xys = [(x, y) for x, y in zip(self.xdata, self.ydata)]

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

    # This might be what they call a 'closure'? (We have to pass the input key to callback somehow.)
    def callback_generator(self, event):

        input_key = event.button

        def callback(verts):
            p = path.Path(verts)
            bool_array = p.contains_points(self.xys)

            indices = np.arange(len(bool_array))[bool_array]

            selected_structures = [structure for structure in self.structures if structure.idx in indices]

            self.hub.select(input_key, selected_structures, subtree=False)

            self.fig.canvas.draw_idle()
            self.fig.canvas.widgetlock.release(self.lasso)
            del self.lasso

        return callback

    def onpress(self, event):
        if self.fig.canvas.widgetlock.locked(): return
        if event.inaxes is None: return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback_generator(event))
        # acquire a lock on the widget drawing
        self.fig.canvas.widgetlock(self.lasso)

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
