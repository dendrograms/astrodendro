Astronmical Dendrograms
=======================

The aim of this module is to provide an easy way to compute dendrograms of
observed or simulated Astronomical data in Python. The easiest way to think of
a dendrogram is to think of a tree that represents the hiearchy of the
structures in your data. For an example of use of dendrograms on real data,
see `Goodman, A. et al (2009) <http://adsabs.harvard.edu/abs/2009Natur.457...63G>`_.

Requirements
------------

This package has the following depdenencies:

* `Python <http://www.python.org>`_ 2.6 or later (Python 3.x is supported)
* `Numpy <http://www.numpy.org>`_ 1.4.1 or later
* `Astropy <http://www.astropy.org>`_ 0.2.0 or later, optional (needed for reading/writing FITS files)
* `h5py <http://www.h5py.org>`_ 0.2.0 or later, optional (needed for reading/writing HDF5 files)

Installation
------------

At this time, there are no stable releases of the core dendrogram code, so you
will need to install the package from the git repository::

    git clone https://github.com/dendrograms/dendro-core.git
    cd dendro-core
    python setup.py install

You may need to add the ``--user`` option to the last line if you do not have
root access.

Terminology
-----------

In a dendrogram, there are two types of structures: *branches*, which are
structures which split into multiple sub-structures, and branches, which are
structures that have no sub-structure. Branches can split up into branches and
leaves, resulting in a hiearchical structure.

Usage
-----

Computing a Dendrogram
^^^^^^^^^^^^^^^^^^^^^^

Dendrograms can be computed from an n-dimensional array using:

    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram.compute(array)

where ``array`` is a Numpy array and ``d`` is then an instance of the
``Dendrogram`` class, which can be used to access the computed dendrogram (see
`Exploring the Dendrogram`_ below). Where the ``array`` comes from is not
important - for example it can be read in from a FITS file, from an HDF5 file,
or it can be generated in memory. If you are interested in making a dendrogram
from data in a FITS file, you can do:

    >>> from astropy.io import fits
    >>> array = fits.getdata('observations.fits')
    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram.compute(array)

The computation may take anywhere between less than a second to several
minutes depending on the size and complexity of the data. By default, the
above command will compute a dendrogram where there are as many levels in the
dendrograms as pixels, which is likely not what you are interested in. There
are several options to control the computation of the dendrogram and can be
passed to the ``compute`` method:

* ``min_intensity``: the minimum value to consider in the dataset - any value
  lower than this will not be considered in the dendrogram. If you are working
  with observations, it is likely that you will want to set this to the
  `detection` level, for example 3- or 5-sigma, so that only significant
  values are included in the dendrogram. By default, all values are used.

* ``min_npix``: the minimum number of pixels/values needed for a leaf to be
  considered an indepdendent entity. When the dendrogram is being computed,
  and when a leaf is about to be joined onto a branch or another leaf, if the
  leaf has fewer than this number of pixels, then it is combined with the
  branch or leaf it is being merged with and is no longer considered a
  separate entity. By default, this parameter is set to zero, so there is no
  minimum number of pixels required for leaves to remain indepdendent entities.

* ``min_delta``: the minimum `height` a leaf has to have in order to be
  considered an independent entity. The `height` of the leaf is the difference
  between its peak flux and the value at which it is being merged into the
  tree at. If you are working with observational data, then you could set this
  to e.g. 1-sigma, which means that any leaf that is less than 1-sigma tall is
  combined with its neighboring leaf or branch and is no longer considered a
  separate entity.

For example, if you have an observational dataset with values in mJy/beam, and
a noise level of 0.1 mJy/beam, you could use::

   >>> sigma = 0.1
   >>> d = Dendrogram.compute(array, min_intensity=3 * sigma,
                             min_npix=10, min_delta=sigma)

which will compute a dendrogram using only values above 3-sigma, and in which
all leaves will have 10 or more pixels and will have a height of at least
1-sigma.

Exploring the Dendrogram
^^^^^^^^^^^^^^^^^^^^^^^^

Once the dendrogram has been computed, you will want to explore/visualize it.
At this time, there are no visualization tools included by default, but you
can access the full tree from the computed dendrogram. Assuming that you have
computed a dendrogram with::

    >>> d = Dendrogram.compute(array, ...)

you can now access the full tree from the ``d`` variable.

The first place to start is the *trunk* of the tree (the ``trunk`` attribute),
which is a list of all the structures at the lowest level. Unless you left
``min_intensity`` to the default setting which means that all values in the
dataset are used, it's possible that not all structures are connected. So the
``trunk`` is a collection of items at the lowest level, each of which could be
a leaf or a branch (itself having leaves and branches). This is accessed with::

    >>> d.trunk
    [<astrodendro.components.Branch object at 0x10279b250>]

In the above case, the trunk only contains a single branch. Since ``trunk`` is
just a list, you can access items in it with e.g.::

    >>> d.trunk[0]
    <astrodendro.components.Branch object at 0x10279b250>

Branches have an ``items`` attribute which returns a list of all
sub-structures, which can include branches and leaves. Thus, we can return the
sub-structures of the above branch with::

    >>> d.trunk[0].items
    [<astrodendro.components.Branch object at 0x10279b1d0>, <astrodendro.components.Leaf object at 0x10279b210>]

which shows that the branch is composed of another branch, and a leaf. We can
therefore access the sub-structures of this branch with::

    >>> d.trunk[0].items[0].items
    [<astrodendro.components.Leaf object at 0x10279b150>, <astrodendro.components.Branch object at 0x10279b190>]

which again shows the branch splitting into a leaf and a branch.

We can access the properties of leaves as follows::

    >>> leaf = d.trunk[0].items[1]
    >>> leaf.coords
    [(230, 50, 75),
    (230, 50, 74),
    (229, 50, 74),
    (229, 51, 74)]
    >>> leaf.f
    [1.4287809,
    1.4096074,
    1.4536692,
    1.4319911]

The following attributes are available for leaves:

* ``coords``: a list of tuples giving the n-dimensional co-ordinates of the
  pixels in the leaf.

* ``f``: a list of the pixel values in the leaf, in the same order as
  ``coords``.

* ``fmin`` and ``fmax``: the minimum and maximum flux of values in the leaf
  (this is equivalent to ``f.max()`` and ``f.min()``)

* ``height``: if the leaf is not attached to the tree, then this is simply
  ``fmax - fmin``. If the leaf is attached to a tree, then it is the
  difference between the leaf and the value at which the leaf was merged into
  the tree (which will be the next value that would have been included in the
  leaf had the leaf not been merged).

* ``parent``: the structure directly containing the leaf.

* ``ancestor``: the largest structure containing the leaf.

* ``level``: the level of the leaf in the tree, i.e. how many structures and
  sub-structures need to be traversed to reach the leaf.

Branches include the same attributes, with the following addition:

* ``descendents``: a flattened list of all leaves and branches that are
  sub-structures of the present branch.

Saving the dendrogram
^^^^^^^^^^^^^^^^^^^^^

A ``Dendrogram`` object can be exported to an HDF5 file (requires h5py) and
loaded at a later time (FITS support is currently planned). To export the
dendrogram to an HDF5 file, use::

    >>> d.save_to('my_dendrogram.hdf5')

and to load and existing dendrogram::

    >>> d = Dendrogram.load_from('my_other_dendrogram.hdf5')

Reporting issues
----------------

Please help us improve this package by reporting issues via `GitHub
<https://github.com/dendrograms/dendro-core/issues>`_.

Developers
----------

This package was developed by:

* Braden MacDonald
* Chris Beaumont
* Thomas Robitaille
* Erik Rosolowsky

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

