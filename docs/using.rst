Using ``astrodendro``
=====================

Computing a Dendrogram
----------------------

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

* ``min_value``: the minimum value to consider in the dataset - any value
  lower than this will not be considered in the dendrogram. If you are working
  with observations, it is likely that you will want to set this to the
  `detection` level, for example 3- or 5-sigma, so that only significant
  values are included in the dendrogram. By default, all values are used.

* ``min_delta``: the minimum `height` a leaf has to have in order to be
  considered an independent entity. The `height` of the leaf is the difference
  between its peak flux and the value at which it is being merged into the
  tree at. If you are working with observational data, then you could set this
  to e.g. 1-sigma, which means that any leaf that is less than 1-sigma tall is
  combined with its neighboring leaf or branch and is no longer considered a
  separate entity.

* ``min_npix``: the minimum number of pixels/values needed for a leaf to be
  considered an independent entity. When the dendrogram is being computed,
  and when a leaf is about to be joined onto a branch or another leaf, if the
  leaf has fewer than this number of pixels, then it is combined with the
  branch or leaf it is being merged with and is no longer considered a
  separate entity. By default, this parameter is set to zero, so there is no
  minimum number of pixels required for leaves to remain independent entities.

For example, if you have an observational dataset with values in mJy/beam, and
a noise level of 0.1 mJy/beam, you could use::

   >>> sigma = 0.1
   >>> d = Dendrogram.compute(array, min_value=3 * sigma,
                              min_delta=sigma, min_npix=10)

which will compute a dendrogram using only values above 3-sigma, and in which
all leaves will have 10 or more pixels and will have a height of at least
1-sigma.

Exploring the Dendrogram
------------------------

Once the dendrogram has been computed, you will want to explore/visualize it.
At this time, there are no visualization tools included by default, but you
can access the full tree from the computed dendrogram. Assuming that you have
computed a dendrogram with::

    >>> d = Dendrogram.compute(array, ...)

you can now access the full tree from the ``d`` variable.

The first place to start is the *trunk* of the tree (the ``trunk`` attribute),
which is a list of all the structures at the lowest level. Unless you left
``min_value`` to the default setting which means that all values in the
dataset are used, it's possible that not all structures are connected. So the
``trunk`` is a collection of items at the lowest level, each of which could be
a leaf or a branch (itself having leaves and branches). This is accessed with::

    >>> d.trunk
    [<Branch>]

In the above case, the trunk only contains a single branch. Since ``trunk`` is
just a list, you can access items in it with e.g.::

    >>> d.trunk[0]
    [<Structure type=branch>]

Branches have an ``children`` attribute which returns a list of all
sub-structures, which can include branches and leaves. Thus, we can return the
sub-structures of the above branch with::

    >>> d.trunk[0].items
    [<Structure type=branch>, <Structure type=leaf>]

which shows that the branch is composed of another branch, and a leaf. We can
therefore access the sub-structures of this branch with::

    >>> d.trunk[0].items[0].items
    [<Structure type=leaf>, <Structure type=branch>]

which again shows the branch splitting into a leaf and a branch.

We can access the properties of leaves as follows::

    >>> leaf = d.trunk[0].children[1]
    >>> leaf.indices
    [(230, 50, 75),
    (230, 50, 74),
    (229, 50, 74),
    (229, 51, 74)]
    >>> leaf.values
    [1.4287809,
    1.4096074,
    1.4536692,
    1.4319911]

The following attributes are available for leaves and branches:

* ``indices``: a list of tuples giving the n-dimensional indices of the pixels
  in the structure (excluding sub-structures).

* ``values``: a list of the pixel values in the leaf, in the same order as
  ``indices`` (excluding sub-structures)

* ``indices_all``: as for ``indices`` but including sub-structures

* ``values_all``: as for ``values`` but including sub-structures

* ``vmin`` and ``vmax``: the minimum and maximum flux of values in the structure
  (this is equivalent to ``values.max()`` and ``values.min()``)

* ``height``: if the structure is not attached to a tree, then this is simply
  ``vmax - vmin``. If the leaf is attached to a tree, then it is the difference
  between the leaf and the value at which the structure was merged into the
  tree (which will be the next value that would have been included in the leaf
  had the structure not been merged).

* ``children``: all direct sub-structures to the present structure.

* ``parent``: the structure directly containing the structure.

* ``ancestor``: the largest structure containing the structure.

* ``level``: the level of the structure in the tree, i.e. how many structures
  and sub-structures need to be traversed to reach the present structure.

* ``descendents``: a flattened list of all sub-structures of the present
  structure.

Computing Dendrogram Statistics
-------------------------------
For 2D (PP) and 3D (PPV) observational data, you can use the ``pp_catalog``
and ``ppv_catalog`` functions to compute basic properties for each
Dendrogram structure::

   >>> import numpy as np
   >>> from astrodendro import Dendrogram, ppv_catalog
   >>> d = Dendrogram.compute(np.random.random((10, 10, 10)))
   >>> metadata = {}
   >>> cat = ppv_catalog(d, metadata)

   WARNING: Missing Metadata:
    bmaj (Beam major axis, sigma)
     Defaulting to bmaj=0 [astrodendro.analysis]
   WARNING: Missing Metadata:
    bmin (Beam minor axis, sigma)
     Defaulting to bmin=0 [astrodendro.analysis]
   WARNING: Missing Metadata:
    bunit (Unit of intensity)
     Defaulting to bunit=1 [astrodendro.analysis]
   WARNING: Missing Metadata:
    dist (Distance)
     Defaulting to dist=1 [astrodendro.analysis]
   WARNING: Missing Metadata:
    dv (Velocity channel width)
     Defaulting to dv=1 [astrodendro.analysis]
   WARNING: Missing Metadata:
    dx (Angular length of a pixel)
     Defaulting to dx=1 [astrodendro.analysis]
   WARNING: Missing Metadata:
    vaxis (Index of velocity axis (numpy convention))
     Defaulting to vaxis=1 [astrodendro.analysis]

   >>> print cat[:3]
   _idx      flux         luminosity    ...    sky_rad         vrms
   ---- ------------- ----------------- ... ------------- -------------
    191 64.1306480569   0.0195353125403 ... 2.85334306153 2.96246166695
     12 4.63582743919   0.0014121537931 ...  3.2987034401  3.5720567466
    >>> print cat.columns
       <TableColumns names=('_idx','flux','luminosity','sky_deconvolved_rad','sky_maj','sky_min','sky_pa','sky_rad','vrms')>

The catalog functions return an astropy `Table` object.

The ``metadata`` dictionary provides information about how to convert pixel-level quantities to meaningful units. By default, ``ppv_catalog`` generates warnings about missing metadata items (these can be suppressed by setting ``verbose=False`` in the call to ``ppv_catalog``).

Here's a sensible looking metadata dictionary::

    >>> import astropy.units as u
    >>> md = dict(dv=0.5 * u.km / u.s,
    >>>           vaxis=0,
    >>>           dx=.002 * u.deg,
    >>>           dist=100 * u.pc,
    >>>           bunit=u.K,
    >>>           bmaj=.004 * u.deg,
    >>>           bmin=.004 * u.deg)
    >>> cat = ppv_catalog(d, md)
    >>> for c in cat.columns:
    >>>     print c, cat[c].units
   _idx None
   flux deg2 K km / (s)
   luminosity K km pc2 / (s)
   sky_deconvolved_rad deg
   sky_maj deg
   sky_min deg
   sky_pa None
   sky_rad deg
   vrms km / (s)

Here's a brief description of each quantity computed in the catalog functions:

* ``_idx`` : The structure ``.idx`` that this row describes
* ``flux`` : The integrated intensity of each structure
* ``luminosity`` : ``flux * d^2``
* ``sky_mag`` : The intensity-weighted second moment of emission, along the major axis of the structure projected onto the sky
* ``sky_min`` : The intensity-weighted second moment of emission, perpendicular to the major axis of the structure projected onto the sky
* ``sky_pa`` : The position angle of the structure projected onto the sky. Given in radians CCW from the +x axis (note that this is the +x axis in pixel coordinates, which is the ``-x`` axis for conventional astronomy images)
* ``sky_rad`` : The geometric mean of ``sky_maj`` and ``sky_min``
* ``vrms`` : The intensity-weighted second moment of emission, along the velocity axis. The velocity axis is given by the ``vaxis`` metadata item. This axis is in Numpy convention, which is the reverse of FITS convention (that is, if an array is read from a FITS file where ``AXIS3`` is the velocity axis, then ``vaxis=0``).
* ``sky_deconvolved_rad``: The size of the structure, corrected for the effects of beam-smearing.

For more information on these quantities, consult the paper on `Bias Free Measurements of Molecular Cloud Properties <http://adsabs.harvard.edu/abs/2006PASP..118..590R>`_ or `the original dendrogram paper <http://adsabs.harvard.edu/abs/2008ApJ...679.1338R>`_. In the terminology of the dendrogram paper, the quantities in ``ppv_catalog`` and ``pp_catalog`` adopt the "bijection" paradigm.

Saving the dendrogram
---------------------

A ``Dendrogram`` object can be exported to an HDF5 file (requires h5py) and
loaded at a later time (FITS support is currently planned). To export the
dendrogram to an HDF5 file, use::

    >>> d.save_to('my_dendrogram.hdf5')

and to load and existing dendrogram::

    >>> d = Dendrogram.load_from('my_other_dendrogram.hdf5')
