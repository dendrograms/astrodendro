Plotting Dendrograms
====================

Once you have computed a dendrogram, you will likely want to plot it as well as
over-plot the structures on your original image.

Interactive Visualization
-------------------------

One you have computed your dendrogram, the easiest way to view it interactively
is to use the :meth:`~astrodendro.dendrogram.Dendrogram.viewer` method::

    d = Dendrogram.compute(...)
    d.viewer()

This will launch an interactive window showing the original data, and the
dendrogram itself. Note that the viewer is only available for 2 or 3-d
datasets. The main window will look like this:

.. image:: viewer_screenshot.png
   :width: 100%

Within the viewer, you can:

**Highlight structures:** either click on structures in the dendrogram to
highlight them, which will also show them in the image view on the left, or
click on pixels in the image and have the corresponding structure be
highlighted in the dendrogram plot. Clicking on a branch in the dendrogram plot
or in the image will highlight that branch and all sub-structures.

**Change the image stretch:** use the ``vmin`` and ``vmax`` sliders above the
image to change the lower and upper level of the image stretch.

**Change slice in a 3-d cube:** if you select a structure in the dendrogram for
a 3-d cube, the cube will automatically change to the slice containing the peak
pixel of the structure (including sub-structures). However, you can also change
slice manually by using the ``slice`` slider.

**View the selected structure ID:** in a computed dendrogram, every structure
has a unique integer ID (the ``.idx`` attribute) that can be used to recognize
the identify the structure when computing catalogs or making plots manually
(see below).

Making plots for publications
-----------------------------

While the viewer is useful for exploring the dendrogram, it does not allow one
to produce publication-quality plots. For this, you can use the non-interactive
plotting interface. To do this, you can first use the
:meth:`~astrodendro.dendrogram.Dendrogram.plotter` method to provide a plotting
tool::

    d = Dendrogram.compute(...)
    p = d.plotter()

and then use this to make the plot you need. The following complete example
shows how to make a plot of the dendrogram of an extinction map of the Perseus
region using the :meth:`~astrodendro.plot.DendrogramPlotter.plot_tree`,
highlighting two of the main branches:

.. plot::
   :include-source:

    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astrodendro import Dendrogram

    image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')
    d = Dendrogram.compute(image, min_value=2.0, min_delta=1., min_npix=10)
    p = d.plotter()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Plot the whole tree
    p.plot_tree(ax, color='black')

    # Highlight two branches
    p.plot_tree(ax, structure=2077, color='red', lw=2, alpha=0.5)
    p.plot_tree(ax, structure=3262, color='orange', lw=2, alpha=0.5)

    # Add axis labels
    ax.set_xlabel("Structure")
    ax.set_ylabel("Flux")

You can find out the structure ID you need either from the interactive viewer
presented above, or programmatically by accessing the ``idx`` attribute of a
Structure.

A :meth:`~astrodendro.plot.DendrogramPlotter.plot_contour` method is also
provided to outline the contours of structures. Calling
:meth:`~astrodendro.plot.DendrogramPlotter.plot_contour` without any arguments
results in a contour corresponding to the value of ``min_value`` used being
shown.

.. plot::
   :include-source:

    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astrodendro import Dendrogram

    image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')
    d = Dendrogram.compute(image, min_value=2.0, min_delta=1., min_npix=10)
    p = d.plotter()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(image, origin='lower', interpolation='nearest', cmap=plt.cm.Blues, vmax=4.0)

    # Show contour for ``min_value``
    p.plot_contour(ax, color='black')

    # Highlight two branches
    p.plot_contour(ax, structure=2077, lw=3, colors='red')
    p.plot_contour(ax, structure=3262, lw=3, colors='orange')


