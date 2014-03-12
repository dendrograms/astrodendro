

import matplotlib.pyplot as plt

class Selection(object):
    selected = {} # selection id -> list of dendrogram structure ids
    colors = {}  # selection id -> a color


class DendroScatter(object):

    def __init__(self, dendrogram, catalog, xaxis, yaxis):
        self.fig = plt.figure()
        pass

    def _draw_plot(self):
        pass

    def update_selection(self, selection):
        """Highlight seleted structures"""
        pass

    def select(self, structure, index):
        "Select a given structure and assign it to index-th selection"
        raise NotImplementedError()