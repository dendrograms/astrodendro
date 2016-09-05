Example Uses of astrodendro
===========================

All of these examples use the Perseus data set

Example 1
---------

Extract the dendrogram structure using a minimum-peak cutoff, then make some
intricate plots.

These plots show two derived quantities - the "flux" and the effective radius -
plotted against each other in two different ways.  In order to visualize the
structure as plotted, the lines are color- and marker- coded to match each
other.


.. plot::
   :include-source:

    from astrodendro import Dendrogram
    from astropy.io import fits
    from astrodendro.pruning import min_peak
    from astrodendro.analysis import PPStatistic
    import itertools
    import pylab as pl
    import numpy as np

    image = fits.getdata('PerA_Extn2MASS_F_Gal.fits')

    from astropy import units as u
    # this metadata isn't really appropriate for this data set, but it is
    # useful for the example
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam
    metadata['spatial_scale'] =  6 * u.arcsec
    metadata['beam_major'] =  22.9 * u.arcsec # FWHM
    metadata['beam_minor'] =  22.9 * u.arcsec # FWHM

    d = Dendrogram.compute(image, min_value=1.0, min_delta=1., min_npix=10,
                           is_independent=min_peak(3), verbose=True)

    f1 = pl.figure(1,figsize=(15,5))
    f1.clear()
    ax1,ax2,ax3 = (pl.subplot(1,3,ii) for ii in xrange(1,4))

    p = d.plotter()
    p.plot_tree(ax3, color='k', alpha=0.5)

    # quantities to plot:
    flux = {}
    radius = {}

    markers = itertools.cycle(['o','h','s','^','v','<','>','p','D','*'])

    for leaf in d.leaves:
        id = leaf.idx
        print("Plotting ID ",id)
        # need one suite of plotted things for each leaf
        flux[id] = []
        radius[id] = []

        while hasattr(leaf,'parent'):
            stat = PPStatistic(leaf, metadata=metadata)
            flux[id].append(stat.flux.value)
            radius[id].append(stat.radius.value)
            leaf = leaf.parent
        
        L, = ax1.loglog(np.array(radius[id])**2 * np.pi,
                   flux[id],
                   '-',marker=markers.next(),
                   alpha=0.5)
        ax2.loglog(np.array(radius[id])**2 * np.pi,
                   np.array(flux[id]) / (4/3. * np.pi * np.array(radius[id])**2),
                   'o-', color=L.get_color(), marker=L.get_marker(),
                   alpha=0.5)
        p.plot_tree(ax3, structure=id, color=L.get_color(), linewidth=2, alpha=0.5)
        lines = p.get_lines(structure=id, color=L.get_color(), linewidth=2, alpha=0.5)
        ax3.plot(*lines.get_segments()[0].T,marker=L.get_marker())



    ax1.set_xlabel('Area (arcsec$^2$)')
    ax1.set_ylabel('Flux (Jy)')
    ax2.set_xlabel('Area (arcsec$^2$)')
    ax2.set_ylabel('"Column Density" (Jy/arcsec$^3$)')
    ax3.set_xlabel("Structure")
    ax3.set_ylabel("Flux")
