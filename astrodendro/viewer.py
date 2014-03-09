# Licensed under an MIT open source license - see LICENSE

from collections import defaultdict

import numpy as np
import matplotlib.backend_bases
from .plot import DendrogramPlotter


class BasicDendrogramViewer(object):

    def __init__(self, dendrogram):

        if dendrogram.data.ndim not in [2, 3]:
            raise ValueError("Only 2- and 3-dimensional arrays are supported at this time")

        self.array = dendrogram.data
        self.dendrogram = dendrogram
        self.plotter = DendrogramPlotter(dendrogram)
        self.plotter.sort(reverse=True)

        # Get the lines as individual elements, and the mapping from line to structure
        self.lines = self.plotter.get_lines()

        # Define the currently selected subtree
        self.selected = {}
        self.selected_lines = {}
        self.selected_contour = {}

        self.colordict = defaultdict(lambda: 'red')
        self.colordict[1] = 'red'
        self.colordict[2] = 'green'
        self.colordict[3] = 'magenta'

        # Initiate plot
        import matplotlib.pyplot as plt
        self.fig = plt.figure(figsize=(14, 8))

        self.ax1 = self.fig.add_axes([0.1, 0.1, 0.4, 0.7])

        from matplotlib.widgets import Slider

        self._clim = (np.min(self.array[~np.isnan(self.array) & ~np.isinf(self.array)]),
                      np.max(self.array[~np.isnan(self.array) & ~np.isinf(self.array)]))

        if self.array.ndim == 2:

            self.slice = None
            self.image = self.ax1.imshow(self.array, origin='lower', interpolation='nearest', vmin=self._clim[0], vmax=self._clim[1], cmap=plt.cm.gray)

            self.slice_slider = None

        else:

            if self.array.shape[0] > 1:

                self.slice = int(round(self.array.shape[0] / 2.))

                self.slice_slider_ax = self.fig.add_axes([0.1, 0.95, 0.4, 0.03])
                self.slice_slider_ax.set_xticklabels("")
                self.slice_slider_ax.set_yticklabels("")
                self.slice_slider = Slider(self.slice_slider_ax, "3-d slice", 0, self.array.shape[0], valinit=self.slice, valfmt="%i")
                self.slice_slider.on_changed(self.update_slice)
                self.slice_slider.drawon = False

            else:

                self.slice = 0
                self.slice_slider = None

            self.image = self.ax1.imshow(self.array[self.slice, :,:], origin='lower', interpolation='nearest', vmin=self._clim[0], vmax=self._clim[1], cmap=plt.cm.gray)

        self.vmin_slider_ax = self.fig.add_axes([0.1, 0.90, 0.4, 0.03])
        self.vmin_slider_ax.set_xticklabels("")
        self.vmin_slider_ax.set_yticklabels("")
        self.vmin_slider = Slider(self.vmin_slider_ax, "vmin", self._clim[0], self._clim[1], valinit=self._clim[0])
        self.vmin_slider.on_changed(self.update_vmin)
        self.vmin_slider.drawon = False

        self.vmax_slider_ax = self.fig.add_axes([0.1, 0.85, 0.4, 0.03])
        self.vmax_slider_ax.set_xticklabels("")
        self.vmax_slider_ax.set_yticklabels("")
        self.vmax_slider = Slider(self.vmax_slider_ax, "vmax", self._clim[0], self._clim[1], valinit=self._clim[1])
        self.vmax_slider.on_changed(self.update_vmax)
        self.vmax_slider.drawon = False

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
        self.fig.canvas.mpl_connect('key_press_event', self.select_from_map)

        plt.show()

    def update_slice(self, pos=None):
        if self.array.ndim == 2:
            self.image.set_array(self.array)
        else:
            self.slice = int(round(pos))
            self.image.set_array(self.array[self.slice, :,:])

        # This is gonna have to be a little more nuanced, because we'll want to
        # update contours for ALL of the currently-defined Input Keys.
        self.remove_all_contours()
        self.update_contour()

        self.fig.canvas.draw()

    def update_vmin(self, vmin):
        if vmin > self._clim[1]:
            self._clim = (self._clim[1], self._clim[1])
        else:
            self._clim = (vmin, self._clim[1])
        self.image.set_clim(*self._clim)
        self.fig.canvas.draw()

    def update_vmax(self, vmax):
        if vmax < self._clim[0]:
            self._clim = (self._clim[0], self._clim[0])
        else:
            self._clim = (self._clim[0], vmax)
        self.image.set_clim(*self._clim)
        self.fig.canvas.draw()

    def select_from_map(self, event):

        # Only do this if no tools are currently selected
        if event.canvas.toolbar.mode != '':
            return

        if event.inaxes is self.ax1:

            if isinstance(event, matplotlib.backend_bases.MouseEvent):
                print "Mouse button pressed: ", event.button
                input_key = event.button
            elif isinstance(event, matplotlib.backend_bases.KeyEvent):
                print "Key pressed: ", event.key
                input_key = event.key

            # Find pixel co-ordinates of click
            ix = int(round(event.xdata))
            iy = int(round(event.ydata))

            if self.array.ndim == 2:
                indices = (iy, ix)
            else:
                indices = (self.slice, iy, ix)

            # Select the structure
            structure = self.dendrogram.structure_at(indices)
            self.select(structure, input_key=input_key)

            # Re-draw
            event.canvas.draw()

    def line_picker(self, event):

        # Only do this if no tools are currently selected
        if event.canvas.toolbar.mode != '':
            return

        print "Mouse button pressed: ", event.mouseevent.button
        input_key = event.mouseevent.button

        # event.ind gives the indices of the paths that have been selected

        # Find levels of selected paths
        peaks = [event.artist.structures[i].get_peak(subtree=True)[1] for i in event.ind]

        # Find position of minimum level (may be duplicates, let Numpy decide)
        ind = event.ind[np.argmax(peaks)]

        # Extract structure
        structure = event.artist.structures[ind]

        # If 3-d, select the slice
        if self.slice_slider is not None:
            peak_index = structure.get_peak(subtree=True)
            self.slice_slider.set_val(peak_index[0][0])

        # Select the structure
        self.select(structure, input_key=input_key)

        # Re-draw
        event.canvas.draw()

    def select(self, structure, input_key=1):

        # Remove previously selected collection
        if input_key in self.selected_lines:
            self.ax2.collections.remove(self.selected_lines[input_key])
            del self.selected_lines[input_key]

        self.remove_contour(input_key)

        if structure is None:
            self.selected_label.set_text("No structure selected")
            self.fig.canvas.draw()
            return

        self.selected[input_key] = structure

        self.selected_label.set_text("Selected structure: {0}".format(structure.idx))

        # Get collection for this substructure
        self.selected_lines[input_key] = self.plotter.get_lines(structure=structure)
        self.selected_lines[input_key].set_color(self.colordict[input_key])
        self.selected_lines[input_key].set_linewidth(2)
        self.selected_lines[input_key].set_alpha(0.5)

        # Add to axes
        self.ax2.add_collection(self.selected_lines[input_key])

        self.update_contour()

    def remove_contour(self, input_key=1):

        if input_key in self.selected_contour:
            for collection in self.selected_contour[input_key].collections:
                self.ax1.collections.remove(collection)
            del self.selected_contour[input_key]

    def remove_all_contours(self):
        """ Remove all selected contours. """
        for key in self.selected_contour.keys():
            self.remove_contour(key)


# Perhaps this shouldn't take input_key as a keyword, but instead should go through 
# ALL the input keys available?
    def update_contour(self):

        keys = self.selected.keys()

        for input_key in keys:
            if input_key in self.selected:
                mask = self.selected[input_key].get_mask(subtree=True)
                if self.array.ndim == 3:
                    mask = mask[self.slice, :,:]
                self.selected_contour[input_key] = self.ax1.contour(mask, 
                    colors=self.colordict[input_key], linewidths=2, levels=[0.5], alpha=0.5)
