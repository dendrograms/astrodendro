Advanced topics
===============

Specifying a custom structure merging strategy
----------------------------------------------

By default, the decision about whether a leaf remains independent when merged
is made based on the ``min_delta`` and ``min_npix`` parameters, but in some
cases, you may want to use more specialized criteria. For example, you may want
only leaves overlapping with a certain position, or you may want leaves with a
certain spatial or velocity extent, or a minimum peak value, to be considered
independent structures.

In order to accomodate this, the
:meth:`~astrodendro.dendrogram.Dendrogram.compute` method can optionally take
an ``is_independent`` argument which should be a function with the following
call signature::

    def is_independent(structure, index=None, value=None):
        ...

where ``structure`` is the :class:`~astrodendro.structure.Structure` object
that is being considered, and ``index`` and ``value`` are the pixel index and
value of the pixel that is linking the structure to the rest of the tree. These
last two values are only set when calling the ``is_independent`` function
during the tree computation, but the ``is_independent`` function is also used
at the end of the computation to prune leaves that are not attached to the
tree, and in this case ``index`` and ``value`` are not set.

The following example compares the dendrogram obtained without and with a custom ``is_independent`` function:

.. plot::
   :include-source:

    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astrodendro import Dendrogram

    image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')

    fig = plt.figure(figsize=(10,5))

    # Default merging strategy

    d1 = Dendrogram.compute(image, min_value=2.0)
    p1 = d1.plotter()

    ax = fig.add_subplot(1, 2, 1)
    p1.plot_tree(ax, color='black')
    ax.set_xlabel("Structure")
    ax.set_ylabel("Flux")
    ax.set_title("Default merging")

    # Require minimum peak value

    def custom_independent(structure, index=None, value=None):
        peak_index, peak_value = structure.get_peak()
        return peak_value > 3.5

    d2 = Dendrogram.compute(image, min_value=2.0, is_independent=custom_independent)
    p2 = d2.plotter()

    ax = fig.add_subplot(1, 2, 2)
    p2.plot_tree(ax, color='black')
    ax.set_xlabel("Structure")
    ax.set_ylabel("Flux")
    ax.set_title("Custom merging")
