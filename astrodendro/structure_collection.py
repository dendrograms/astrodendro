# Licensed under an MIT open source license - see LICENSE

from matplotlib.collections import LineCollection


class StructureCollection(LineCollection):

    @property
    def structures(self):
        return self._structures

    @structures.setter
    def structures(self, values):
        self._structures = values
