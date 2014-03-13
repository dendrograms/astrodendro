import matplotlib.pyplot as plt

class DendroScatter(object):
    def __init__(self, dendrogram, hub, catalog, xaxis, yaxis):

        self.hub = hub
        self.hub.add_callback(self.update_selection)
        self.dendrogram = dendrogram

        self.fig = plt.figure()
        self.axes = plt.subplot(1,1,1)

        self.catalog = catalog
        self.xdata = catalog[xaxis]
        self.ydata = catalog[yaxis]

        self.x_column_name = xaxis
        self.y_column_name = yaxis

        self.lines2d = {} # selection_id -> matplotlib.lines.Line2D

        self._draw_plot()
        self.hub.add_callback(self.update_selection)

    def _draw_plot(self):

        self.axes.plot(self.xdata, self.ydata, 'o', color='w', mec='k', zorder=-5)

        self.axes.set_xlabel(self.x_column_name)
        self.axes.set_ylabel(self.y_column_name)

        self.fig.canvas.draw()

    def update_selection(self, selection_id):
        """Highlight seleted structures"""
        
        if selection_id in self.lines2d:
            if self.lines2d[selection_id] is not None:
                self.lines2d[selection_id].remove()
                del self.lines2d[selection_id]

        struct = self.hub.selections[selection_id][0]
        if struct is None:
            self.fig.canvas.draw()
            return
        selected_indices = [leaf.idx for leaf in struct.descendants + [struct]]

        self.lines2d[selection_id] = self.axes.plot(
            self.xdata[selected_indices], 
            self.ydata[selected_indices], 
            'o', color=self.hub.colors[selection_id], zorder=struct.height)[0]

        self.fig.canvas.draw()

    def select(self, structure, index):
        "Select a given structure and assign it to index-th selection"
        raise NotImplementedError()