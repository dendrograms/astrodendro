import os

import numpy as np
from astropy.io import fits

from astrodendro import Dendrogram

BENCHMARKS = {'2d1.fits': {'min_value': 2.5, 'min_npix': 20},
              '2d2.fits': {'min_value': 2.5},
              '2d3.fits': {'min_value': 2.5, 'min_delta': 1},
              '3d1.fits': {'min_value': 4.5, 'min_npix': 20},
              '3d2.fits': {'min_value': 4.5},
              '3d3.fits': {'min_value': 4.5, 'min_delta':1}}

def main():

    data = fits.getdata(os.path.join('benchmark_data', '2d.fits'))
    for outfile in '2d1.fits 2d2.fits 2d3.fits'.split():
        d = Dendrogram.compute(data, verbose=True, **BENCHMARKS[outfile])
        d.save_to(os.path.join('benchmark_data', outfile))

    data = fits.getdata(os.path.join('benchmark_data', '3d.fits'))
    for outfile in '3d1.fits 3d2.fits 3d3.fits'.split():
        d = Dendrogram.compute(data, verbose=True, **BENCHMARKS[outfile])
        d.save_to(os.path.join('benchmark_data', outfile))


if __name__ == "__main__":
    main()
