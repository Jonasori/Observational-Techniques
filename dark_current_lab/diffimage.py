# Python program to load a FITS image and display it

# import needed extensions
from numpy import *
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
from astropy.io import fits

# read in the files
# change the file names as appropriate
h1 = fits.open('testpattern1.fit')
h2 = fits.open('testpattern2.fit')

# copy the image data into a numpy (numerical python) array
img1 = h1[0].data
img2 = h2[0].data

plt.ion() # do plots in interactive mode
colmap = plt.get_cmap('gray') # load gray colormap

# plot the first image
plt.figure(1)  
plt.imshow(img1, cmap=colmap) # plot image using gray colorbar
#plt.show() # display the image

# plot the second image in another window
plt.figure(2)  
plt.imshow(img2, cmap=colmap) # plot image using gray colorbar
#plt.show() # display the image

# find the difference in the images
diff = img2-img1

# plot the difference image
plt.figure(3)
plt.imshow(diff, cmap=colmap) # plot image using gray colorbar

plt.show() # display the images


