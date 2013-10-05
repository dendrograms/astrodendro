Migration guide for previous users of ``astrodendro``
=====================================================

The ``astrodendro`` package has been in development for a couple of years, and
we have recently undertaken an effort to prepare the package for a first
release, which involved tidying up the programming interface to the package,
and re-writing large sections. This means that the present version of
``astrodendro`` will likely not work with scripts you had if you were using the
original astrodendro packages from @astrofrog and @brandenmacdonald's
repositories.

This page summarizes the main changes in the new code, and how to adapt your
code to ensure that it will work correctly. This only covers changes that will
*break* your code, but you are encouraged to look through the rest of the
documentation to read about new features! Also, only the main
backward-incompatible changes are mentioned, but for any questions on changes
not mentioned here, please open an issue on `GitHub
<https://github.com/dendrograms/astrodendro/issues>`_.

Computing a dendrogram
----------------------

Rather than computing a dendrogram using::

    d = Dendrogram(data)
    d.compute(...)

you should now use::

    d = Dendrogram.compute(data)

In addition, the following options for ``compute`` have been renamed:

* ``minimum_flux`` is now ``min_value`` (since we expect dendrograms to be used
  not only for images, but also e.g. density fields).

* ``minimum_delta`` is now ``min_delta``

* ``minimum_npix`` is now ``min_npix``

Dendrogram methods and attributes
---------------------------------

The following dendrogram methods have changed:

* ``get_leaves()`` has now been replaced by a ``leaves`` attribute (it is no
  longer a method.)

* the ``to_hdf5()`` and ``from_hdf5()`` methods have been replaced by
  :meth:`~astrodendro.dendrogram.Dendrogram.save_to`

``Leaf`` and ``Branch`` classes
-------------------------------

The ``Leaf`` and ``Branch`` classes no longer exist, and have been replaced by
a single :class:`~astrodendro.structure.Structure` class that instead has
``is_leaf`` and ``is_branch`` attributes. Thus, if you were checking if
something was a leaf by doing e.g.::

    if type(s) == Leaf:
        # code here

or::

    if isinstance(s, Leaf):
         # code here

then you will instead need to use::

    if s.is_leaf:
         # code here

Leaf and branch attributes
--------------------------

The following leaf and branch attributes have changed:

* ``f`` has been replaced by a method called :meth:`~astrodendro.structure.Structure.values` that can take a
  ``subtree=`` option that indicates whether pixels in sub-structures should be
  included.

* ``coords`` has been replaced by a method called :meth:`~astrodendro.structure.Structure.indices` that can take a
  ``subtree=`` option that indicates whether pixels in sub-structures should be
  included.

* ``height`` now has a different definition - it is ``vmax`` for a leaf, or the
  smallest ``vmin`` of the children for a branch - this is used when plotting
  the dendrogram, to know at what height to plot the structure.

Interactive visualization
-------------------------

Visualizing the results of the dendrogram is now much easier, and does not
require the additional ``astrocube`` package. To launch the interactive viewer
(which requires only Matplotlib), once the dendrogram has been computed, you can do:

    >>> d.viewer()

and the interactive viewer will launch. It will however no longer have the
option to re-compute the dendrogram from the window, and will also no longer
have an IPython terminal. For the latter, we recommend you consider using the
`Glue <http://www.glue-viz.org>`_ package.

  