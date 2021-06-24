#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 07:51:49 2021

@author: dexter

This script is meant to plot 2D galaxy images

Loop through the list of image, model, residual. and produce 
all the three plots. 

"""
# import 
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAna
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl
import cmasher as cmr
import matplotlib.gridspec as gridspec

"""
Input area
"""
# Read the name of the ISOFIT output from a list
outlist = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')

# Read the geometry file for the sample
geom_file = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat")
geom_file_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')

bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
bundle = SRead.read_list(bundle_name)



image_list = "/media/dexter/My Passport/Ozstar_back/NGC3665.fits"
md_list= "/media/dexter/My Passport/Ozstar_back/md5_NGC3665.fits"
res_list = "/media/dexter/My Passport/Ozstar_back/res5_NGC3665.fits"

B = SPlot.Plot2D.read_fits_img(image_list)
B2 = SPlot.Plot2D.read_fits_img(md_list)
B3 = SPlot.Plot2D.read_fits_img(res_list)

avg_data0, std_data0 = np.average(B), np.std(B)
val_min0, val_max0 =  avg_data0 , avg_data0+1.5*std_data0


img ="/home/dexter/result/image_plot/fit_example/NGC2872.fits"
md = "/home/dexter/result/image_plot/fit_example/md1_NGC2872.fits"
res= "/home/dexter/result/image_plot/fit_example/res1_NGC2872.fits"


SPlot.Plot2D.plot_galaxy_3plot(img,md, res, (252, 432), r_max=350, 
                               name = r"$\rm NGC~2872$",alp =15)
plt.show() #2872

centre = (179, 776)


SPlot.Plot2D.plot_galaxy_3plot(image_list,md_list, res_list, centre, r_max=650, 
                               name = r"$\rm NGC~3665$",alp =10)
plt.show() #2872
#SPlot.Plot2D.plot_galaxy_3plot(image_list,md_list, res_list, 
#                               centre,r_max=400, alp=15)

"""
Plotting
"""

def plot_all_gal_3plot():
    """
    A function to plot all galaxy image,model,residual 3 plots
    
    """
    return None