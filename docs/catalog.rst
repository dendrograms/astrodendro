Computing Dendrogram Statistics
===============================

For 2D position-position (PP) and 3D position-position-velocity (PPV)
observational data, the :doc:`api/astrodendro.analysis` module can be used to
compute basic properties for each Dendrogram structure. There are two ways to
compute statistics - on a structure-by-structure basis and as a catalog, both
of which are described below.

Deriving statistics for individual structures
---------------------------------------------

In order to derive statistics for a given structure, you will need to use the
:class:`~astrodendro.analysis.PPStatistic` or the
:class:`~astrodendro.analysis.PPVStatistic` classes from the
:doc:`api/astrodendro.analysis` module, e.g.::

   >>> from astrodendro.analysis import PPStatistic
   >>> stat = PPStatistic(structure)

where ``structure`` is a :class:`~astrodendro.structure.Structure` instance
from a dendrogram. The resulting object then has methods to compute various
statistics. Using the example data from :doc:`using`::

    >>> from astrodendro import Dendrogram
    >>> from astropy.io import fits
    >>> image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')
    >>> d = Dendrogram.compute(image, min_value=2.0, min_delta=1., min_npix=10)

we can get statistics for the first structure in the trunk, which is a leaf::

    >>> from astrodendro.analysis import PPStatistic
    >>> d.trunk[0]
    <Structure type=leaf idx=101>
    >>> stat = PPStatistic(d.trunk[0])
    >>> stat.major_sigma  # length of major axis on the sky
    <Quantity 1.882980574564531 pix>
    >>> stat.minor_sigma  # length of minor axis on the sky
    <Quantity 1.4639300383020182 pix>
    >>> stat.position_angle  # position angle on the sky
    <Quantity 134.61988014787443 deg>

.. note:: The objects returned are Astropy :class:`~astropy.units.quantity.Quantity`
          objects, which are Numpy scalars or arrays with units attached. For more
          information, see the `Astropy Documentation
          <http://docs.astropy.org/en/stable/units/index.html>`_. In most cases,
          you should be able to use these objects directly, but if for any
          reason you need to access the underlying value, then you can do so
          with the ``value`` and ``unit`` attributes::

              >>> q = 1.882980574564531 * u.pix
              >>> q.unit
              Unit("pix")
              >>> q.value
              1.882980574564531

Specifying meta-data when computing statistics
----------------------------------------------

In some cases, meta-data can or should be specified. To demonstrate this, we
will use a different data set which is a small section
(:download:`L1551_scuba_850mu.fits`) of a SCUBA 850 micron map from the `SCUBA
Legacy Catalog
<http://www3.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/scubalegacy/>`_. This
map has a pixel scale of 6 arcseconds per pixel and a circular beam with a
full-width at half maximum (FWHM) of 22.9 arcseconds. First, we compute the
dendrogram as usual::

    >>> from astropy.io import fits
    >>> from astrodendro import Dendrogram
    >>> image = fits.getdata('L1551_scuba_850mu.fits')
    >>> d = Dendrogram.compute(image, min_value=0.1, min_delta=0.02)

then we set up a Python dictionary containing the required meta-data::

    >>> from astropy import units as u
    >>> metadata = {}
    >>> metadata['data_unit'] = u.Jy / u.beam
    >>> metadata['spatial_scale'] =  6 * u.arcsec
    >>> metadata['beam_major'] =  22.9 * u.arcsec # FWHM
    >>> metadata['beam_minor'] =  22.9 * u.arcsec # FWHM

Note that the meta-data required depends on the units of the data and whether
you are working in position-position or position-position-velocity (see
`Required metadata`_).

Finally, as before, we use the :class:`~astrodendro.analysis.PPStatistic` class
to extract properties for the first structure::

    >>> from astrodendro.analysis import PPStatistic
    >>> stat = PPStatistic(d.trunk[0], metadata=metadata)
    >>> stat.major_sigma
    <Quantity 20.34630778380526 arcsec>
    >>> stat.minor_sigma
    <Quantity 8.15504176035544 arcsec>
    >>> stat.position_angle
    <Quantity 85.14309012311242 deg>
    >>> stat.flux
    <Quantity 0.24119688679751278 Jy>

Note that the major and minor sigma on the sky of the structures are now in
arcseconds since the spatial scale was specified, and the flux (density) has
been converted from Jy/beam to Jy. Note also that the flux does not include any
kind of background subtraction, and is just a plain sum of the values in the
structure, converted to Jy

Making a catalog
----------------

In order to produce a catalog of properties for all structures, it is also
possible to make use of the :func:`~astrodendro.analysis.pp_catalog` and
:func:`~astrodendro.analysis.ppv_catalog` functions. We demonstrate this using
the same SCUBA data as used above::

    >>> from astropy.io import fits
    >>> from astrodendro import Dendrogram, pp_catalog
    >>> image = fits.getdata('L1551_scuba_850mu.fits')
    >>> d = Dendrogram.compute(image, min_value=0.1, min_delta=0.02)

    >>> from astropy import units as u
    >>> metadata = {}
    >>> metadata['data_unit'] = u.Jy / u.beam
    >>> metadata['spatial_scale'] =  6 * u.arcsec
    >>> metadata['beam_major'] =  22.9 * u.arcsec
    >>> metadata['beam_minor'] =  22.9 * u.arcsec

    >>> cat = pp_catalog(d, metadata)
    >>> cat.pprint(show_unit=True, max_lines=10)
    _idx       flux       major_sigma   minor_sigma  ...     radius        x_cen         y_cen
                Jy           arcsec        arcsec    ...     arcsec         pix           pix
    ---- --------------- ------------- ------------- ... ------------- ------------- -------------
       7  0.241196886798 20.3463077838 8.15504176036 ... 12.8811874315 168.053017504 3.98809714744
      51  0.132470059814 14.2778133293 4.81100492125 ...  8.2879810685  163.25495657 9.13394216473
      60 0.0799106574322 9.66298008473 3.47364264736 ... 5.79359471511 169.278409915 15.1884110291
     ...             ...           ...           ... ...           ...           ...           ...
    1203  0.183438198239 22.7202518034 4.04690367115 ... 9.58888264776 15.3760934458 100.136384362
    1384   2.06217635837 38.1060171889  19.766115194 ... 27.4446338168 136.429313911 107.190835447
    1504   1.90767291972 8.64476839751 8.09070477357 ... 8.36314946298  68.818705665 120.246719845

The catalog function returns an Astropy :class:`~astropy.table.table.Table` object.

Note that :func:`~astrodendro.analysis.pp_catalog` and
:func:`~astrodendro.analysis.ppv_catalog` generate warnings if required
meta-data is missing and sensible defaults can be assumed. If no sensible
defaults can be assumed (e.g. for ``data_unit``) then an exception is raised.

Unlike clumpfind-style algorithms, there is not a one-to-one mapping between
identifiers and pixels in the map: each pixel may belong to multiple nested
branches in the catalog.

Available statistics
--------------------

For a full list of available statistics for each type of statistic class, see
:class:`~astrodendro.analysis.PPStatistic` and
:class:`~astrodendro.analysis.PPVStatistic`. For more information on the
quantities listed in these pages, consult the paper on `Bias Free Measurements
of Molecular Cloud Properties
<http://adsabs.harvard.edu/abs/2006PASP..118..590R>`_ or `the original
dendrogram paper <http://adsabs.harvard.edu/abs/2008ApJ...679.1338R>`_. In the
terminology of the dendrogram paper, the quantities in
:class:`~astrodendro.analysis.PPStatistic` and
:class:`~astrodendro.analysis.PPVStatistic` adopt the "bijection" paradigm.

Required metadata
-----------------

As described above, the metadata needed by the statistic routines depends on
what statistics are required and on the units of the data. With the exception
of ``wcs``, all meta-data should be specified as `Astropy Quantity
<http://docs.astropy.org/en/stable/units/index.html>`_ objects (e.g., ``3 *
u.arcsec``):

* ``data_unit`` is **required** in order to compute the flux, so it is needed
  for both the :func:`~astrodendro.analysis.pp_catalog` and
  :func:`~astrodendro.analysis.ppv_catalog` functions, as well as for the
  ``flux`` attribute of the :class:`~astrodendro.analysis.PPStatistic` and
  :class:`~astrodendro.analysis.PPVStatistic` classes.
  Note: if ``data_unit`` is given as ``K``, it is interpreted as units of 
  main beam brightness temperature, following the conventions in the astropy 
  `Brightness Temperature / Flux Density equivalency <http://docs.astropy.org/en/stable/units/equivalencies.html#brightness-temperature-flux-density-equivalency>`_ .

* ``spatial_scale`` is **required** if the data are in units of surface
  brightness (e.g. ``MJy/sr`` or ``Jy/beam``) so as to be able to convert the
  surface brightness to the flux in each pixel. Even if the data are not in
  units of surface brightness, the ``spatial_scale`` can **optionally** be
  specified, causing any derived size (e.g. ``major_sigma``) to be in the
  correct units instead of in pixels.

* ``velocity_scale`` can **optionally** be specified for PPV data, causing
  ``v_rms`` to be in the correct units instead of in pixels.

* ``beam_major`` and ``beam_minor`` are **required** if the data units depend
  on the beam (e.g. ``Jy/beam`` or ``K``).

* ``vaxis`` can **optionally** be specified when using 3-dimensional data to
  indicate which dimension corresponds to the velocity. By default, this is
  ``0``, which corresponds to the third axis in a FITS file (because the
  dimensions are reversed in numpy).

* ``wavelength`` is **required** if the data are in monochromatic flux
  densities per unit wavelength because the fluxes need to be converted to
  monochromatic flux densities per unit frequency. It is also required if
  the data are in brightness temperature units of ``K``.

* ``wcs`` can **optionally** be specified and should be a
  :class:`~astropy.wcs.WCS` instance. If specified, it allows ``x_cen``,
  ``y_cen``, and ``v_cen`` to be in the correct world coordinates rather than
  in pixel coordinates.

Example
-------

The following example shows how to combine the plotting functionality in
:doc:`plotting` and the analysis tools shown above, to overlay ellipses
approximating the structures on top of the structures themselves:

.. plot::
   :include-source:

    from astropy.io import fits

    from astrodendro import Dendrogram
    from astrodendro.analysis import PPStatistic

    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    hdu = fits.open('PerA_Extn2MASS_F_Gal.fits')[0]

    d = Dendrogram.compute(hdu.data, min_value=2.0, min_delta=1., min_npix=10)
    p = d.plotter()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.imshow(hdu.data, origin='lower', interpolation='nearest',
              cmap=plt.cm.Blues, vmax=6.0)

    for leaf in d.leaves:

        p.plot_contour(ax, structure=leaf, lw=3, colors='red')

        s = PPStatistic(leaf)
        ellipse = s.to_mpl_ellipse(edgecolor='orange', facecolor='none')

        ax.add_patch(ellipse)

    ax.set_xlim(75., 170.)
    ax.set_ylim(120., 260.)

As shown above, the :class:`~astrodendro.analysis.PPStatistic` and
:class:`~astrodendro.analysis.PPVStatistic` classes have a
:meth:`~astrodendro.analysis.PPStatistic.to_mpl_ellipse` method to convert the
first and second moments of the structures into schematic ellipses.


