#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 05:22:56 2020

@author: dexter
"""



img ="/home/dexter/result/image_plot/fit_example/NGC2872.fits"
md = "/home/dexter/result/image_plot/fit_example/md1_NGC2872.fits"
res= "/home/dexter/result/image_plot/fit_example/res1_NGC2872.fits"

Plot2D.plot_galaxy_3plot(img,md,
                         res, (252, 432),r_max=207, a =15)
plt.show() #2872


#Plot2D.plot_galaxy_3plot(img,md,
  #                       res, (1088, 789),r_max=255, a =15)
#plt.show()4045
#
#Plot2D.plot_galaxy_3plot(img,md,
#                         res, (271, 1080),r_max=250, a =15)
#plt.show() #3675
