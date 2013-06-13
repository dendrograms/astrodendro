Installing ``astrodendro``
==========================

Requirements
------------

This package has the following depdenencies:

* `Python <http://www.python.org>`_ 2.6 or later (Python 3.x is supported)
* `Numpy <http://www.numpy.org>`_ 1.4.1 or later
* `Astropy <http://www.astropy.org>`_ 0.2.0 or later, optional (needed for reading/writing FITS files)
* `h5py <http://www.h5py.org>`_ 0.2.0 or later, optional (needed for reading/writing HDF5 files)

Installation
------------

At this time, there are no stable releases of the core dendrogram code, so you
will need to install the package from the git repository::

    git clone https://github.com/dendrograms/dendro-core.git
    cd dendro-core
    python setup.py install

You may need to add the ``--user`` option to the last line if you do not have
root access.