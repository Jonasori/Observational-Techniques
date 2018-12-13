"""
Dark Current Lab
ASTR 522
Jonas Powell

Copied/pasted from Roy's histimage.py:
Python program to load a FITS image and display it
"""


# import needed extensions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits


# Read in the file
h = fits.open("/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/bias/bias_-5deg.fits")


def get_img_stats(img, plot_field=False, plot_hist=True):
    """ Print some statistics about the image.

    img: a 1D image array to plot. I think this has to be 1D?
    """
    print('Image minimum = ', np.nanmin(img))
    print('Image maximum = ', np.nanmax(img))
    print('Image mean = ', np.mean(img))
    print('Image standard deviation = ', np.std(img))

    if plot_field is True:
        colmap = plt.get_cmap('gray') # load gray colormap
        plt.figure(1)
        plt.imshow(img, cmap=colmap)

    if plot_hist is True:
        plt.figure(2)
        plt.hist(img, bins=100, histtype='stepfilled')
        plt.show()


# DO DATA STUFF
img = h[0].data
# Get array dimensions
nx, ny = img.shape
# Change the array to 1D for histogram plotting.
imgh = np.reshape(img, nx*ny)

# Cut out bad data (dead/hot pixels)
plow = 1861.0
phi = 2000.0
q = np.where((imgh >= plow) & (imgh <= phi))
imghcut = imgh[q]


get_img_stats(imghcut)





# The End
