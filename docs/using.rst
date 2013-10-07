Computing and exploring dendrograms
===================================

For a graphical description of the actual algorithm used to compute
dendrograms, see :doc:`algorithm`.

Computing a Dendrogram
----------------------

Dendrograms can be computed from an n-dimensional array using:

    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram.compute(array)

where ``array`` is a Numpy array and ``d`` is then an instance of the
:class:`~astrodendro.dendrogram.Dendrogram` class, which can be used to access
the computed dendrogram (see `Exploring the Dendrogram`_ below). Where the
``array`` comes from is not important - for example it can be read in from a
FITS file, from an HDF5 file, or it can be generated in memory. If you are
interested in making a dendrogram from data in a FITS file, you can do:

    >>> from astropy.io import fits
    >>> array = fits.getdata('observations.fits')
    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram.compute(array)

.. There should probably be a stronger/bolder warning against doing this
   example blindly because it will create a LOT of dendro branches.
   The computation may take anywhere between less than a second to several
   minutes depending on the size and complexity of the data. By default, the
   above command will compute a dendrogram where there are as many levels in the
   dendrograms as pixels, which is likely not what you are interested in. There
   are several options to control the computation of the dendrogram and can be
   passed to the :meth:`~astrodendro.dendrogram.Dendrogram.compute` method:

* ``min_value``: the minimum value to consider in the dataset - any value
  lower than this will not be considered in the dendrogram. If you are working
  with observations, it is likely that you will want to set this to the
  "detection level," for example 3- or 5-sigma, so that only significant
  values are included in the dendrogram. By default, all values are used.

* ``min_delta``: how significant a leaf has to be in order to be considered an
  independent entity. The significance is measured from the difference between
  its peak flux and the value at which it is being merged into the tree. If
  you are working with observational data, then you could set this to, e.g.,
  1-sigma, which means that any leaf that is locally less than 1-sigma tall is
  combined with its neighboring leaf or branch and is no longer considered a
  separate entity.

* ``min_npix``: the minimum number of pixels/values needed for a leaf to be
  considered an independent entity. When the dendrogram is being computed,
  and when a leaf is about to be joined onto a branch or another leaf, if the
  leaf has fewer than this number of pixels, then it is combined with the
  branch or leaf it is being merged with and is no longer considered a
  separate entity. By default, this parameter is set to zero, so there is no
  minimum number of pixels required for leaves to remain independent entities.

These options are illustrated graphically in :doc:`algorithm`.

As an example, we can use a publicly available extinction map of the Perseus
star-formation region from the The COordinated Molecular Probe Line Extinction
Thermal Emission (COMPLETE) Survey of Star Forming Regions
(:download:`PerA_Extn2MASS_F_Gal.fits`, originally obtained from
`<http://hdl.handle.net/10904/10080>`_). The units of the map are magnitudes of
extinction, and we want to make a dendrogram of all structures above a minimum
value of 2 magnitudes, and we only consider leaves with at least 10 pixels and
which have a peak to base difference larger than one magnitude of extinction::

    >>> from astrodendro import Dendrogram
    >>> from astropy.io import fits
    >>> image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')
    >>> d = Dendrogram.compute(image, min_value=2.0, min_delta=1., min_npix=10)

By default, the computation will be silent, but for large dendrograms, it can
be useful to have an idea of how long the computation will take::

    >>> d = Dendrogram.compute(image, min_value=2.0, min_delta=1., min_npix=10,
                               verbose=True)
    Generating dendrogram using 6,386 of 67,921 pixels (9% of data)
    [=========================>               ] 64%

The '9% of data' indicates that only 9% of the data are over the `min_value`
threshold.

Exploring the Dendrogram
------------------------

Once the dendrogram has been computed, you will want to explore/visualize it.
You can access the full tree from the computed dendrogram. Assuming that you
have computed a dendrogram with::

    >>> d = Dendrogram.compute(array, ...)

you can now access the full tree from the ``d`` variable.

The first place to start is the *trunk* of the tree (the ``trunk`` attribute),
which is a `list` of all the structures at the lowest level. Unless you left
``min_value`` to the default setting, which would mean that all values in the
dataset are used, it's likely that not all structures are connected. So the
``trunk`` is a collection of items at the lowest level, each of which could be
a leaf or a branch (itself having leaves and branches). In the case of the
Perseus extinction map, we get::

    >>> d.trunk
    [<Structure type=leaf idx=101>,
     <Structure type=branch idx=2152>,
     <Structure type=leaf idx=733>,
     <Structure type=branch idx=303>]

In the above case, the trunk contains two leaves and two branches. Since
``trunk`` is just a list, you can access items in it with e.g.::

    >>> d.trunk[1]
    <Structure type=branch idx=2152>

Branches have a ``children`` attribute that returns a list of all
sub-structures, which can include branches and leaves. Thus, we can return the
sub-structures of the above branch with::

    >>> d.trunk[1].children
    [<Structure type=branch idx=1680>,
     <Structure type=branch idx=5771>]

which shows that the child branch is composed of two more branches. We can
therefore access the sub-structures of these branch with e.g.::

    >>> d.trunk[1].children[0].children
    [<Structure type=leaf idx=1748>,
     <Structure type=leaf idx=1842>]

which shows this branch splitting into two leaves.

We can access the properties of leaves as follows::

    >>> leaf = d.trunk[1].children[0].children[0]
    >>> leaf.indices
    (array([143, 142, 142, 142, 139, 141, 141, 141, 143, 140, 140]),
     array([116, 114, 115, 116, 115, 114, 115, 116, 115, 115, 114]))
    >>> leaf.values
    array([ 2.7043395 ,  2.57071948,  3.4551146 ,  3.29953575,  2.53844047,
            2.59633183,  3.11309052,  2.70936489,  2.81024122,  2.76864815,
            2.52840114], dtype=float32)

A full list of attributes and methods for leaves and branches (i.e. structures)
is available from the :class:`~astrodendro.structure.Structure` page, while a
list of attributes and methods for the dendrogram itself is available from the
:class:`~astrodendro.dendrogram.Dendrogram` page.

Saving and loading the dendrogram
---------------------------------

A :class:`~astrodendro.dendrogram.Dendrogram` object can be exported to an HDF5 file (requires h5py) or FITS file (requires astropy). To export the
dendrogram to a file, use::

    >>> d.save_to('my_dendrogram.hdf5')

or::

    >>> d.save_to('my_dendrogram.fits')

and to load and existing dendrogram::

    >>> d = Dendrogram.load_from('my_other_dendrogram.hdf5')

or::

    >>> d = Dendrogram.load_from('my_other_dendrogram.fits')

If you wish, you can use this to separate the computation and analysis of the
dendrogram into two scripts, to ensure that the dendrogram is only computed
once. For example, you could have a script ``compute.py`` that contains::

    from astropy.io import fits
    from astrodendro import Dendrogram

    array = fits.getdata('observations.fits')
    d = Dendrogram.compute(array)
    d.save_to('dendrogram.fits')

and a second file containing::

    from astrodendro import Dendrogram
    d = Dendrogram.load_from('dendrogram.fits')

    # any analysis code here
