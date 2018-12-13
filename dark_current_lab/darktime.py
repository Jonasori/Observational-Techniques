# Python program to load a FITS image and display it

# import needed extensions
from numpy import *
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
from astropy.io import fits
from scipy import stats

# name of file containing the bias frame
biasfile = 'dark_20c_0s.fit' # 'dark_25c_1s.fit' 
# read in the bias file
h = fits.open(biasfile)
# copy the image data into a numpy (numerical python) array
bias = h[0].data

# list of dark files
darkfile = ['dark_20c_1s.fit', 'dark_20c_2s.fit', 'dark_20c_5s.fit', 
            'dark_20c_20s.fit', 'dark_20c_60s.fit', 
            'dark_20c_180s.fit', 'dark_20c_600s.fit']
# array of corresponding times
time = array([1.0, 2, 5, 20, 60, 180, 600])

# arrays to hold results
c_mean = 0.0*time
c_median = 0.0*time
c_rms = 0.0*time

# wait after each set of plots
plt.ioff()
colmap = plt.get_cmap('gray') # load gray colormap

# process the files
for i in range(len(darkfile)):
  # read in the file
  h = fits.open(darkfile[i])
  print 'Read ', darkfile[i]
  # copy the image data into a numpy (numerical python) array
  img = h[0].data
  # find the difference fromthe bias frames
  diff = img - bias
  # diff is a 2-d array, need to change to 1-d to make a histogram
  nx, ny = diff.shape # find the size of the array
  diffh = reshape(diff, nx*ny) # change the shape to be 1d
  # choose selection region
  f = 0.003 # ignore first and last fraction f of points
  s = sort(diffh)
  vmin = s[f*len(s)]
  vmax = s[(1-f)*len(s)]
  print 'Excluding lowest and highest pixels, fraction = ', f
  print 'bounds = ', vmin, vmax
  # one figure with two plots
  plt.figure(1)
  # plot the first image
  plt.subplot(121)
  plt.imshow(diff, cmap=colmap, vmin=vmin, vmax=vmax) # plot image using gray colorbar
  # plot the second image
  plt.subplot(122)
  # plot a histogram of the image values
  ht = plt.hist(diffh, bins=vmax-vmin, range=[vmin,vmax], histtype='stepfilled', color='k')
  nc = max(ht[0]) # maximum value in plotted histogram
  # select only values within ranges
  q = where((diffh >= vmin) & (diffh <= vmax))
  diffhcut = diffh[q]
  # find statistics
  c_mean[i] = mean(diffhcut)
  c_median[i] = median(diffhcut)
  c_rms[i] = std(diffhcut)
  print 'mean, median, rms = ', c_mean[i], c_median[i], c_rms[i]
  plt.plot([c_mean[i], c_mean[i]], [0, nc], '-g')
  plt.plot([c_median[i], c_median[i]], [0, nc], '--g')
  plt.show() # display the plots
  print
  raw_input("Press Enter to continue...")

# plot median versus time
m = c_mean
#m = c_median
plt.figure(2)
plt.xlabel('Time (s)')
plt.ylabel('Mean Counts')
#plt.ylabel('Median Counts')
plt.plot(time, m, 'd')
# do a linear fit to the data
slope, intercept, r_value, p_value, std_err = stats.linregress(time, m)
print('slope = ', slope)
print('intercept = ', intercept)
print('correlation coefficient r =', r_value)
# plot the fit
plt.plot(time, intercept + time*slope, '--r')
plt.show()
