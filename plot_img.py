#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 05:22:56 2020

@author: dexter
"""
import matplotlib.pyplot as plt

import SphPlot as SPlot

img ="/home/dexter/result/image_plot/fit_example/NGC2872.fits"
md = "/home/dexter/result/image_plot/fit_example/md1_NGC2872.fits"
res= "/home/dexter/result/image_plot/fit_example/res1_NGC2872.fits"


img1 ="/home/dexter/result/image_plot/fit_example/NGC4045.fits"
md1 = "/home/dexter/result/image_plot/fit_example/md1_NGC4045.fits"
res1= "/home/dexter/result/image_plot/fit_example/res1_NGC4045.fits"

img2 ="/home/dexter/result/image_plot/fit_example/NGC3675.fits"
md2 = "/home/dexter/result/image_plot/fit_example/md1_NGC3675.fits"
res2= "/home/dexter/result/image_plot/fit_example/res1_NGC3675.fits"

SPlot.Plot2D.plot_galaxy_3plot(img,md,
                         res, (252, 432),r_max=207, a = 3)
#plt.rcParams["font.family"] = "Times New Roman"
plt.show() #2872


SPlot.Plot2D.plot_galaxy_3plot(img1,md1,
                        res1, (1088, 789),r_max=255, a =3)
plt.show()
#4045
SPlot.Plot2D.plot_galaxy_3plot(img2,md2,
                         res2, (271, 1080),r_max=250, a =3)
plt.show() 

#3675
