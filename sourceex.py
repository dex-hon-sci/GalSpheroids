#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 02:05:55 2020

@author: dexter
"""

import numpy as np
import sep


# additional setup for reading the test image and displaying plots
import fitsio
import matplotlib.pyplot as plt
from matplotlib import rcParams

#%matplotlib inline
rcParams['figure.figsize'] = [10., 8.]


data = fitsio.read("/home/dexter/result/Exclude_case/Bin3/NGC5375.fits")
objects = sep.extract(data, 0.5 )


from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data), np.std(data)
im = ax.imshow(data, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=27*objects['a'][i],
                height=27*objects['b'][i],
                angle=(objects['theta'][i] * 180. / np.pi))
    
    width=20*objects['a'][i]
    height=20*objects['b'][i]
    angle= (objects['theta'][i] * 180. / np.pi) +90.
    r = (width*height)/np.sqrt((width**2)*(np.sin(angle)**2)+(height**2)*(np.cos(angle)**2))
    
    print('ellipse(',objects['x'][i],',',objects['y'][i],',',r,',',width/height,',',angle,')')
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

