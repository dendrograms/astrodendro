import matplotlib.pyplot as plt

class Selection(object):
    def __init__(self):
        self.selected = {} # selection id -> list of dendrogram structure ids
        self.colors = {}  # selection id -> a color


class DendroScatter(object):

    def __init__(self, dendrogram, catalog, xaxis, yaxis):
        self.fig = plt.figure()
        self.axes = plt.subplot(1,1,1)
        self.catalog = catalog
        self.xdata = catalog[xaxis]
        self.ydata = catalog[yaxis]

        self.x_column_name = xaxis
        self.y_column_name = yaxis

    def _draw_plot(self):

        self.axes.plot(self.xdata, self.ydata, 'o', color='w', mec='k')

        self.axes.set_xlabel(self.x_column_name)
        self.axes.set_ylabel(self.y_column_name)

        self.fig.canvas.draw()

    def update_selection(self, selection):
        """Highlight seleted structures"""
        
        self.fig.cla()
        self._draw_plot()

        for key in selection.selected.keys():
            self.axes.plot(
                self.xdata[selection.selected[key]], 
                self.ydata[selection.selected[key]], 
                'o', color=selection.colors[key])


    def select(self, structure, index):
        "Select a given structure and assign it to index-th selection"
        raise NotImplementedError()