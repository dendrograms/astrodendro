import os

from astropy.io import fits

path = os.path.dirname(__file__)

def load_perseus():
    """ Load the Primary HDU for the Perseus extinction map
    used throughout the documentation.

    For more details about the data, see http://hdl.handle.net/10904/10080
    """
    pth = os.path.join(path, 'PerA_Extn2MASS_F_Gal.fits')
    return fits.open(pth)[0]
