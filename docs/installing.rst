Installing ``astrodendro``
==========================

Requirements
------------

This package has the following depdenencies:

* `Python <http://www.python.org>`_ 2.6 or later (Python 3.x is supported)
* `Numpy <http://www.numpy.org>`_ 1.4.1 or later
* `Astropy <http://www.astropy.org>`_ 0.2.0 or later, optional (needed for reading/writing FITS files and for analysis code)
* `h5py <http://www.h5py.org>`_ 0.2.0 or later, optional (needed for reading/writing HDF5 files)

Installation
------------

To install the latest stable release, you can type::

    pip install astrodendro

or you can download the latest tar file from
`PyPI <https://pypi.python.org/pypi/astrodendro>`_ and install it using::

    python setup.py install

Developer version
-----------------

If you want to install the latest developer version of the dendrogram code, you
can do so from the git repository::

    git clone https://github.com/dendrograms/astrodendro.git
    cd astrodendro
    python setup.py install

You may need to add the ``--user`` option to the last line if you do not have
root access.