# Licensed under an MIT open source license - see LICENSE

from matplotlib.collections import LineCollection

__all__ = ['StructureCollection']


class StructureCollection(LineCollection):

    @property
    def structures(self):
        return self._structures

    @structures.setter
    def structures(self, values):
        self._structures = values
