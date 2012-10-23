About
=====

The aim of this Python module is to make it easy to compute dendrograms
of Astronomical data. The use of dendrograms to represent Astronomical
data is described in detail in [Goodman, A. (2009,
Nature)](http://adsabs.harvard.edu/abs/2009Natur.457...63G).

**DISCLAIMER**: The code has not yet been thoroughly tested.

Screenshot
==========

In this screenshot, you can see what the generated dendrograms look like.

![Screenshot](http://i.imgur.com/QDePB.png)

Installing
==========

To install the ``astrodendro`` module, simply do:

    pip install git+git://github.com/dendrograms/dendro-core.git#egg=dendro-core

Using
=====

Dendrograms can be computed from any 2- or 3-D Numpy array:

    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram.compute(array)

where ``array`` can be read in from a FITS file for example:

    >>> import pyfits
    >>> array = pyfits.getdata('observations.fits')

(but ``array`` can also be generated in memory, or read in from other
files, e.g. HDF5)

The tree can be explored from the trunk:

    >>> d.trunk
    [<astrodendro.components.Branch object at 0x10279b250>]

This is a list of the lowest structures. We can select one of these:

    >>> d.trunk[0]
    <astrodendro.components.Branch object at 0x10279b250>

Sub-structures of branches can then be retrieved as a list with the ``children`` attribute:

    >>> d.trunk[0].children
    [<astrodendro.components.Branch object at 0x10279b1d0>, <astrodendro.components.Leaf object at 0x10279b210>]

    >>> d.trunk[0].children[0].children
    [<astrodendro.components.Leaf object at 0x10279b150>, <astrodendro.components.Branch object at 0x10279b190>]

The pixel positions and intensities of all the pixels in a leaf can then be retrieved with the ``coords`` and ``f`` attributes:

    >>> In [7]: d.trunk[0].children[0].coords
    [(230, 50, 75),
    (230, 50, 74),
    (229, 50, 74),
    (229, 51, 74)]
    >>> d.trunk[0].children[0].children[0].f
    [1.4287809,
    1.4096074,
    1.4536692,
    1.4319911]

It is also possible to directly retrieve all the leaves:

    >>> leaves = d.get_leaves()

Computation Options
===================

There are several options that can be used when computing a
``Dendrogram`` object:

* ``min_intensity`` - the minimum intensity of pixels to be used in the
  dendrogram (default is -infinity)
* ``min_npix`` - the minimum number of pixels necessary to create a
  new leaf (default is 0)
* ``min_delta`` - the minimum intensity difference for two structures to
  be treated as being separate (minimum is 0)
* ``verbose`` - set `True` to display a progress bar during computation

For example:

    d = Dendrogram.compute(array, min_intensity=1.2, min_npix=10, min_delta=0.1)
 
Import/Export
=============

Dendrograms can be written out and read in from various format:

    # Create dendrogram and write it out as HDF5
    d = Dendrogram.compute(array)
    d.save_to('observations_dendrogram.hdf5')

    # Read a dendrogram from a FITS file
    d2 = Dendrogram.load_from('observations_dendrogram.fits')

Unit Tests and Benchmarks
=========================

Several unit tests are included, and are also installed with the package.
To run the unit tests, simply run the command

    python setup.py test

A benchmark is also included that uses realistic data to determine how fast the
dendrograms are being generated. To run the benchmark, use the command

    python -m astrodendro.test.benchmark
