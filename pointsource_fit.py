#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 15:22:17 2022

@author: dexter
"""
import numpy as np
from astropy.modeling import models, fitting
from astropy.io import fits
import matplotlib.pyplot as plt


img_location = "/home/dexter/Downloads/ngc4689_vcc2058res1.fits"
hdul = fits.open(img_location)
pic = hdul[0].data
print(hdul[0].data)

#plt.imshow(hdul[0].data, cmap="hot",clim=(0, 10))

radius = 150
x0 = 2333
y0 = 5000-2669


#pic_m = np.fliplr(pic)
pic_m2 = np.flipud(pic)
pic2 = pic_m2[x0-radius:x0+radius,y0-radius:y0+radius]

#print(pic2.shape)
plt.imshow(pic2, cmap="hot",clim=(0, 75))
plt.hlines(pic2.shape[0]/2, pic2.shape[0]/2-20-10,pic2.shape[0]/2-10,color='lime',lw=3)
plt.vlines(pic2.shape[0]/2, pic2.shape[0]/2-20-10,pic2.shape[0]/2-10,color='lime',lw=3)
#plt.imshow(pic_m2, cmap="hot",clim=(0, 75))

x,y = pic2
z = pic2[:]

p_init = models.Moffat2D()
fit_p = fitting.LevMarLSQFitter()

print(p_init)    

p = fit_p(p_init,x,y,z)

