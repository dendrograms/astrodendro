import numpy as np
import matplotlib.pyplot as plt
from astrodendro.plot import DendrogramPlotter


class BasicDendrogramViewer(object):

    def __init__(self, array, dendrogram):

        if array.ndim != 2:
            raise ValueError("Only 2-dimensional arrays are supported at this time")

        self.array = array
        self.dendrogram = dendrogram
        self.plotter = DendrogramPlotter(dendrogram)
        self.plotter.sort(reverse=True)

        # Get the lines as individual elements, and the mapping from line to structure
        self.lines = self.plotter.get_lines()

        # Define the currently selected subtree
        self.selected = None
        self.selected_contour = None

        # Initiate plot
        self.fig = plt.figure(figsize=(14, 8))
        self.ax1 = self.fig.add_axes([0.1, 0.1, 0.4, 0.8])
        self.ax1.imshow(array, origin='lower')
        self.ax2 = self.fig.add_axes([0.6, 0.3, 0.35, 0.4])
        self.ax2.add_collection(self.lines)
        self.selected_label = self.fig.text(0.6, 0.75, "No structure selected", fontsize=18)
        x = [p.vertices[:, 0] for p in self.lines.get_paths()]
        y = [p.vertices[:, 1] for p in self.lines.get_paths()]
        xmin = np.min(x)
        xmax = np.max(x)
        ymin = np.min(y)
        ymax = np.max(y)
        self.lines.set_picker(2.)
        dx = xmax - xmin
        self.ax2.set_xlim(xmin - dx * 0.1, xmax + dx * 0.1)
        self.ax2.set_ylim(ymin * 0.5, ymax * 2.0)
        self.ax2.set_yscale('log')

        self.fig.canvas.mpl_connect('pick_event', self.line_picker)
        self.fig.canvas.mpl_connect('button_press_event', self.select_from_map)

        plt.show()

    def select_from_map(self, event):

        if event.inaxes is self.ax1:

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iy = int(round(event.ydata))

            # Select the structure
            structure = self.dendrogram.node_at((iy, ix))
            self.select(structure)

            # Re-draw
            event.canvas.draw()

    def line_picker(self, event):

        # event.ind gives the indices of the paths that have been selected

        # Find levels of selected paths
        peaks = [event.artist.structures[i].get_peak(subtree=True)[1] for i in event.ind]

        # Find position of minimum level (may be duplicates, let Numpy decide)
        ind = event.ind[np.argmax(peaks)]

        # Select the structure
        self.select(event.artist.structures[ind])

        # Re-draw
        event.canvas.draw()

    def select(self, structure):

        # Remove previously selected collection
        if self.selected is not None:
            self.ax2.collections.remove(self.selected)
            self.selected = None

        # Remove previously selected contour
        if self.selected_contour is not None:
            for collection in self.selected_contour.collections:
                self.ax1.collections.remove(collection)
            self.selected_contour = None

        if structure is None:
            self.selected_label.set_text("No structure selected")
            return

        self.selected_label.set_text("Selected structure: {0}".format(structure.idx))

        # Get collection for this substructure
        self.selected = self.plotter.get_lines(structure=structure)
        self.selected.set_color('red')
        self.selected.set_linewidth(2)
        self.selected.set_alpha(0.5)

        # Add to axes
        self.ax2.add_collection(self.selected)

        # Draw contour
        mask = structure.get_mask(self.array.shape, subtree=True)
        self.selected_contour = self.ax1.contour(mask, colors='red', linewidths=2, levels=[0.5], alpha=0.5)
