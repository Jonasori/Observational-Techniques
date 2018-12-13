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
from scipy import stats

colmap = plt.get_cmap('gray') # load gray colormap


def get_img_stats(img, plow=1861, phi=2000, plot_field=False, plot_hist=True, save_hist=False):
    """ Print some statistics about the image.

    img: a 2D image array to plot. I think this has to be 1D?
    """
    plt.close()

    nx, ny = img.shape
    # Change the array to 1D for histogram plotting.
    imgh = np.reshape(img, nx*ny)

    # Cut out bad data (dead/hot pixels)
    q = np.where((imgh >= plow) & (imgh <= phi))
    imghcut = imgh[q]

    print('Image minimum = ', np.nanmin(imghcut))
    print('Image maximum = ', np.nanmax(imghcut))
    print('Image mean = ', np.mean(imghcut))
    print('Image standard deviation = ', np.std(imghcut))

    if plot_field is True:
        plt.imshow(img, cmap=colmap)
        plt.show()

    if plot_hist is True or save_hist is True:
        plt.figure(2)
        plt.hist(imghcut, bins=100, histtype='stepfilled')
        if save_hist is True:
            # name = raw_input('Please enter ouput name: ')
            plt.savefig('img_hist.pdf')
        if plot_hist is True:
            plt.show()





# DO DATA STUFF
# Read in the file
img = fits.open("/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/bias/bias_-5deg.fits")[0].data

get_img_stats(img, plow=1861, phi=2000, plot_field=True)


# DIFFIMAGE Stuff
img1 = fits.open("/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/bias/bias1_0deg.fits")[0].data
img2 = fits.open("/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/bias/bias2_0deg.fits")[0].data
img_diff = img1 - img2

# Shouldn't some of these be negative?
# Maybe should be removing outliers before diff'ing
get_img_stats(img_diff, plow=-100, phi=100, plot_field=True)


"""
Edit the program to produce a histogram that covers only the difference values of interest
in the bias frame. You histogram should extend out to +/- about 3 to 5 times the standard
deviation. Save the histogram plot. Record the mean and standard deviation of the
differences. The standard deviation is a good estimate of the read noise.
"""

diff_val_of_interest_lo = 1
diff_val_of_interest_hi = 1


gain = 'something'





# DARK CURRENT vs. TIME
def darktime():

    basepath = "/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/"

    # Read in the bias file
    bias = fits.open(basepath + 'darks/dark_0sec_0deg.fits',)[0].data

    # list of dark files
    darkfile = [basepath + 'darks/dark_1sec_0deg.fits',
                basepath + 'darks/dark_5sec_0deg.fits',
                basepath + 'darks/dark_20sec_0deg.fits',
                basepath + 'darks/dark_60sec_0deg.fits',
                basepath + 'darks/dark_180sec_0deg.fits',
                basepath + 'darks/dark_600sec_0deg.fits',]
    # array of corresponding times
    time = np.array([1, 5, 20, 60, 180, 600])

    # arrays to hold results
    c_mean = 0.0*time
    c_median = 0.0*time
    c_rms = 0.0*time

    # wait after each set of plots
    plt.ioff()

    # process the files
    for i in range(len(darkfile)):
      # Read in the file
      # i = 0
      img = fits.open(darkfile[i])[0].data
      # print 'Read ', darkfile[i]

      diff = img - bias
      nx, ny = diff.shape # find the size of the array
      diffh = np.reshape(diff, nx*ny) # change the shape to be 1d

      # choose selection region
      f = 0.003 # ignore first and last fraction f of points
      s = np.sort(diffh)
      vmin = s[int(f*len(s))]
      vmax = s[int((1-f)*len(s))]
      print 'Excluding lowest and highest pixels, fraction = ', f
      print 'bounds = ', vmin, vmax

      plt.figure(1)

      plt.subplot(121)
      plt.imshow(diff, cmap=colmap, vmin=vmin, vmax=vmax)

      plt.subplot(122)
      ht = plt.hist(diffh, bins=vmax-vmin, range=[vmin,vmax], histtype='stepfilled', color='k')

      # Select only values within ranges
      q = np.where((diffh >= vmin) & (diffh <= vmax))
      diffhcut = diffh[q]

      # Find and save statistics
      nc = max(ht[0]) # maximum value in plotted histogram
      c_mean[i] = np.mean(diffhcut)
      c_median[i] = np.median(diffhcut)
      c_rms[i] = np.std(diffhcut)

      print 'Mean, median, RMS = ', c_mean[i], c_median[i], c_rms[i]

      plt.plot([c_mean[i], c_mean[i]], [0, nc], '-g')
      plt.plot([c_median[i], c_median[i]], [0, nc], '--g')
      plt.show()
      # raw_input("Press Enter to continue...")

    # Plot median/mean versus time
    m = c_mean
    plt.figure(2)
    plt.xlabel('Time (s)')
    plt.ylabel('Mean Counts')
    plt.plot(time, m, 'd')

    # Do a linear fit to the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(time, m)
    print('Slope = ', slope)
    print('Intercept = ', intercept)
    print('Correlation coefficient r =', r_value)

    # Plot the fit
    plt.plot(time, intercept + time*slope, '--r')
    plt.show()



darktime()





# DARK TEMP
# Same as above but with different input files





# The End
