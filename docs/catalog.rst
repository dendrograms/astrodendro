Computing Dendrogram Statistics
===============================

For 2D position-position (PP) and 3D position-position-velocity (PPV)
observational data, the :doc:`api/astrodendro.analysis` module can be used to
compute basic properties for each Dendrogram structure. There are two ways to
compute statistics - on a structure-by-structure basis, and as a catalog, both
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
    >>> stat.sky_major_sigma  # length of major axis on the sky
    3.7659611491290619
    >>> stat.sky_minor_sigma  # length of minor axis on the sky
    2.9278600766040364
    >>> stat.sky_pa  # position angle on the sky
    134.61988014787443
    >>> stat.flux  # total flux contained in structure
    88.605692148208618

Making a catalog
----------------

In order to produce a catalog of properties for all structures, it is also
possible to make use of the :func:`~astrodendro.analysis.pp_catalog` and
:func:`~astrodendro.analysis.ppv_catalog` functions::

   >>> import numpy as np
   >>> from astrodendro import Dendrogram, ppv_catalog
   >>> d = Dendrogram.compute(np.random.random((10, 10, 10)))
   >>> metadata = {}
   >>> cat = ppv_catalog(d, metadata)

    WARNING: bmaj (Beam major axis, sigma) missing, defaulting to 0 [astrodendro.analysis]
    WARNING: bmin (Beam minor axis, sigma) missing, defaulting to 0 [astrodendro.analysis]
    WARNING: bunit (Unit of intensity) missing, defaulting to 1 [astrodendro.analysis]
    WARNING: dist (Distance) missing, defaulting to 1 [astrodendro.analysis]
    WARNING: velocity_scale (Velocity channel width) missing, defaulting to 1 [astrodendro.analysis]
    WARNING: spatial_scale (Angular length of a pixel) missing, defaulting to 1 [astrodendro.analysis]
    WARNING: vaxis (Index of velocity axis (numpy convention)) missing, defaulting to 1 [astrodendro.analysis]
    WARNING: wcs (WCS object) missing, defaulting to None [astrodendro.analysis]

   >>> print cat[:3]
   _idx      flux         luminosity    ...  sky_radius        vrms
   ---- ------------- ----------------- ... ------------- -------------
    191 64.1306480569   0.0195353125403 ... 2.85334306153 2.96246166695
     12 4.63582743919   0.0014121537931 ...  3.2987034401  3.5720567466
    >>> print cat.columns
       <TableColumns names=('_idx','flux','luminosity','sky_deconvolved_radius','sky_major_sigma','sky_minor_sigma','sky_pa','sky_radius','vrms')>

The catalog functions return an Astropy :class:`~astropy.table.table.Table` object.

The ``metadata`` dictionary provides information about how to convert
pixel-level quantities to meaningful units. By default,
:func:`~astrodendro.analysis.ppv_catalog` generates warnings about missing
metadata items (these can be suppressed by setting ``verbose=False`` in the
call to :func:`~astrodendro.analysis.ppv_catalog`).

Here's a sensible looking metadata dictionary::

    >>> import astropy.units as u
    >>> md = dict(velocity_scale=0.5 * u.km / u.s,
    >>>           vaxis=0,
    >>>           spatial_scale=.002 * u.deg,
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
   sky_deconvolved_radius deg
   sky_major_sigma deg
   sky_minor_sigma deg
   sky_pa None
   sky_radius deg
   vcen None
   vrms km / (s)
   xcen None
   ycen None

Available statistics
--------------------

For a full list of available statistics for each type of statistic class, see
:class:`~astrodendro.analysis.PPStatistic` and
:class:`~astrodendro.analysis.PPVStatistic`.

Here's a more detailed description of the available quantities:

* ``_idx`` : The structure ``.idx`` that this row describes
* ``flux`` : The integrated intensity of each structure
* ``luminosity`` : ``flux * d^2``
* ``sky_mag`` : The intensity-weighted second moment of emission, along the major axis of the structure projected onto the sky
* ``sky_minor_sigma`` : The intensity-weighted second moment of emission, perpendicular to the major axis of the structure projected onto the sky
* ``sky_pa`` : The position angle of the structure projected onto the sky. Given in radians CCW from the +x axis (note that this is the +x axis in pixel coordinates, which is the ``-x`` axis for conventional astronomy images)
* ``sky_radius`` : The geometric mean of ``sky_major_sigma`` and ``sky_minor_sigma``
* ``vrms`` : The intensity-weighted second moment of emission, along the velocity axis. The velocity axis is given by the ``vaxis`` metadata item. This axis is in Numpy convention, which is the reverse of FITS convention (that is, if an array is read from a FITS file where ``AXIS3`` is the velocity axis, then ``vaxis=0``).
* ``sky_deconvolved_radius``: The size of the structure, corrected for the effects of beam-smearing.
* ``xcen`` : X-position of intensity-weighted centroid (in world units if a ``WCS`` object is stored in ``metadta['wcs']``
* ``ycen`` : Y-position of intensity-weighted centroid (see above)
* ``vcen`` : V-position of intensity-weighted centroid (see above)

For more information on these quantities, consult the paper on `Bias Free
Measurements of Molecular Cloud Properties
<http://adsabs.harvard.edu/abs/2006PASP..118..590R>`_ or `the original
dendrogram paper <http://adsabs.harvard.edu/abs/2008ApJ...679.1338R>`_. In the
terminology of the dendrogram paper, the quantities in
:func:`~astrodendro.analysis.pp_catalog` and
:func:`~astrodendro.analysis.ppv_catalog` adopt the "bijection" paradigm.

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
        ax.add_patch(Ellipse((s.xcen, s.ycen),
                              s.sky_major_sigma * 2.3548,
                              s.sky_minor_sigma * 2.3548,
                              angle=s.sky_pa,
                              edgecolor='orange', facecolor='none'))

    ax.set_xlim(75., 170.)
    ax.set_ylim(120., 260.)
