# Python program to load a FITS image and display it

# import needed extensions
from numpy import *
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
from astropy.io import fits


# read in the file
# change input.fits to the name of your file
h = fits.open("/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/dark_current_lab/bias/bias_-5deg.fits")

# copy the image data into a numpy (numerical python) array
img = h[0].data

plt.ion() # do plots in interactive mode
colmap = plt.get_cmap('magma') # load gray colormap

# plot the image on the screen
plt.figure(1)
plt.imshow(img, cmap=colmap) # plot image using gray colorbar

# img is a 2-d array, need to change to 1-d to make a histogram
#imgh = 1.0*img # make a copy
nx, ny = img.shape # find the size of the array
imgh = reshape(img, nx*ny) # change the shape to be 1d

# print some statistics about the image
print('Image minimum = ', min(imgh))
print('Image maximum = ', max(imgh))
print('Image mean = ', mean(imgh))
print('Image standard deviation = ', std(imgh))

# now plot a histogram of the image values
plt.figure(2)
plt.hist(imgh, bins=100, histtype='stepfilled')

plt.show() # display the plots

plow = 5000.0
phi = 14000.0
q = where((imgh >= plow) & (imgh <= phi))
imghcut = imgh[q]

print('Image minimum = ', min(imghcut))
print('Image maximum = ', max(imghcut))
print('Image mean = ', mean(imghcut))
print('Image standard deviation = ', std(imghcut))
