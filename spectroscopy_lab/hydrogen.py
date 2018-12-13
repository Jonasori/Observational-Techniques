# Python program to load a FITS spectrum and display it

# import needed extensions
from numpy import *
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
from matplotlib.colors import LogNorm
from astropy.io import fits

# read in the files
# change the file names as appropriate
basepath = "/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/spectroscopy_lab/data/"

s1raw = fits.open(basepath + 'H300s.FIT')[0].data
dark1 = fits.open(basepath + 'H300s_dark.FIT')[0].data

# do the dark subtraction
s1sub = s1raw-dark1

plt.ion() # do plots in interactive mode
colmap = plt.get_cmap('gray') # load gray colormap

# plot the difference image
plt.figure(1)
# plot image using gray colorbar on log scale
# adjust vmin and vmax based on spectrum
plt.imshow(s1sub, cmap=colmap, norm=LogNorm(vmin=5, vmax=6E4))
plt.show() # display the image


# calculate a spectrum
# y0 is the center of the band over which the spectrum is extracted
y0 = 220
# dy sets the width of the band (y0-dy to y0+dy)
dy = 5
# figure out dimensions of spectrum image
(ny, nx) = shape(s1sub)
# find a 1-d spectrum by integrating across the band
s1 = zeros(nx)
for i in range(nx):
  s1[i] = sum(s1sub[(y0-dy):(y0+dy+1), i])

# plot the spectrum versus pixel number
p = 1+arange(len(s1)) # pixel numbers
plt.figure(2)
plt.xlabel('Pixel number')
plt.ylabel('Counts')
#plt.yscale('log')
#plt.axis([0, 750, 0, 1E5])
plt.plot(p, s1, '-b')
plt.show() # display the plot

print('Line centroids in pixels')
# calculate centroids for each line
linec = array([    172,     485]) # center of interval (pixel)
lined = array([      4,       4]) # width of interval (pixel)
linew = array([486.133, 656.285]) # wavelength of line (nm)
# maximum value in spectrum (only for plotting)
smax = max(s1)
centroid = 0.0*linec
for i in range(len(linec)):
  # array elements included in this line
  k = range(linec[i]-lined[i], linec[i]+lined[i]+1)
  # find statistics for this line
  centroid[i] = sum(p[k]*s1[k])/sum(s1[k])
  plt.plot([centroid[i], centroid[i]], [0, smax], '-g')
  plt.plot([linec[i]-lined[i], linec[i]-lined[i]], [0, smax], '--g')
  plt.plot([linec[i]+lined[i], linec[i]+lined[i]], [0, smax], '--g')
  print('center, range, mean = ', linec[i], lined[i], centroid[i])
plt.show() # display the plot

print('Calibration')
# calibration
centralp = mean(centroid)
centralw = mean(linew)
slope = (linew[1]-linew[0])/(centroid[1]-centroid[0])
print('central pixel = ', centralp)
print('central wavelength = ', centralw)
print('slope = ', slope)

# plot the spectrum versus wavelength
w = centralw + slope*(p-centralp)
plt.figure(3)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Counts')
#plt.yscale('log')
#plt.axis([0, 750, 0, 1E5])
plt.plot(w, s1, '-b')
plt.show() # display the plot

print('Line centroids in wavelength (nm)')
# calculate centroids for each line
# list of lines of interest
linew = array([486.133, 656.285]) # wavelength of line (nm)
lined = array([    2.0,     2.0]) # width of interval (nm)
# maximum value in spectrum (only for plotting)
smax = max(s1)
centroid = 0.0*linew
for i in range(len(linew)):
  # array elements included in this line
  wlow = linew[i]-lined[i]
  k0 = int(floor((wlow-centralw)/slope+centralp))
  whi = linew[i]+lined[i]
  k1 = int(ceil((whi-centralw)/slope+centralp))
  k = range(k0, k1)
  # find statistics for this line
  centroid[i] = sum(w[k]*s1[k])/sum(s1[k])
  plt.plot([centroid[i], centroid[i]], [0, smax], '-g')
  print('wavelength = ', centroid[i])
plt.show() # display the plot
