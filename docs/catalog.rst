Computing Dendrogram Statistics
===============================

For 2D position-position (PP) and 3D position-position-velocity (PPV)
observational data, you can use the :func:`~astrodendro.analysis.pp_catalog` and
:func:`~astrodendro.analysis.ppv_catalog` functions to compute basic properties
for each Dendrogram structure::

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
   _idx      flux         luminosity    ...  sky_radius        vrms
   ---- ------------- ----------------- ... ------------- -------------
    191 64.1306480569   0.0195353125403 ... 2.85334306153 2.96246166695
     12 4.63582743919   0.0014121537931 ...  3.2987034401  3.5720567466
    >>> print cat.columns
       <TableColumns names=('_idx','flux','luminosity','sky_deconvolved_rad','sky_maj','sky_min','sky_pa','sky_radius','vrms')>

The catalog functions return an Astropy :class:`~astropy.table.table.Table` object.

The ``metadata`` dictionary provides information about how to convert
pixel-level quantities to meaningful units. By default,
:func:`astrodendro.analysis.ppv_catalog` generates warnings about missing
metadata items (these can be suppressed by setting ``verbose=False`` in the
call to :func:`astrodendro.analysis.ppv_catalog`).

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
   sky_radius deg
   vcen None
   vrms km / (s)
   xcen None
   ycen None

Here's a brief description of each quantity computed in the catalog functions:

* ``_idx`` : The structure ``.idx`` that this row describes
* ``flux`` : The integrated intensity of each structure
* ``luminosity`` : ``flux * d^2``
* ``sky_mag`` : The intensity-weighted second moment of emission, along the major axis of the structure projected onto the sky
* ``sky_min`` : The intensity-weighted second moment of emission, perpendicular to the major axis of the structure projected onto the sky
* ``sky_pa`` : The position angle of the structure projected onto the sky. Given in radians CCW from the +x axis (note that this is the +x axis in pixel coordinates, which is the ``-x`` axis for conventional astronomy images)
* ``sky_radius`` : The geometric mean of ``sky_maj`` and ``sky_min``
* ``vrms`` : The intensity-weighted second moment of emission, along the velocity axis. The velocity axis is given by the ``vaxis`` metadata item. This axis is in Numpy convention, which is the reverse of FITS convention (that is, if an array is read from a FITS file where ``AXIS3`` is the velocity axis, then ``vaxis=0``).
* ``sky_deconvolved_rad``: The size of the structure, corrected for the effects of beam-smearing.
* ``xcen`` : X-position of intensity-weighted centroid (in world units if a ``WCS`` object is stored in ``metadta['wcs']``
* ``ycen`` : Y-position of intensity-weighted centroid (see above)
* ``vcen`` : V-position of intensity-weighted centroid (see above)

For more information on these quantities, consult the paper on `Bias Free
Measurements of Molecular Cloud Properties
<http://adsabs.harvard.edu/abs/2006PASP..118..590R>`_ or `the original
dendrogram paper <http://adsabs.harvard.edu/abs/2008ApJ...679.1338R>`_. In the
terminology of the dendrogram paper, the quantities in
:func:`astrodendro.analysis.pp_catalog` and
:func:`astrodendro.analysis.ppv_catalog` adopt the "bijection" paradigm.

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
        ax.add_patch(Ellipse((s.xcen(), s.ycen()),
                              s.sky_maj(), s.sky_min(), angle=s.sky_pa(),
                              edgecolor='orange', facecolor='none'))

    ax.set_xlim(75., 170.)
    ax.set_ylim(120., 260.)
