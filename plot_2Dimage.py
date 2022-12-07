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
img_list = SRead.read_table("/home/dexter/result/stat/completeness/gal_img_list.txt",
                               dtype='str')
md_list = SRead.read_table("/home/dexter/result/stat/completeness/gal_md_list.txt",
                               dtype='str')
res_list = SRead.read_table("/home/dexter/result/stat/completeness/gal_res_list.txt",
                               dtype='str')

# Read the geometry file for the sample
geom_file = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat")
geom_file_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')



#SPlot.Plot2D.plot_galaxy_3plot(image_list,md_list, res_list, 
#                               centre,r_max=400, alp=15)


print(img_list)
"""
Plotting
"""

def plot_all_gal_3plot():
    """
    A function to plot all galaxy image,model,residual 3 plots
    
    """
    name = geom_file_n[:,0]
    centre_x, centre_y = geom_file[:,3], geom_file[:,4]
    R_max = geom_file[:,1]
    
    for i in range(len(name)):
        centre = (centre_x[i],centre_y[i])    
        r_max = int(R_max[i]/0.4+50)
        
        img = str(img_list[i])
        md = str(md_list[i])
        res = str(res_list[i])       
        
        nam = name[i]
        
        SPlot.Plot2D.plot_galaxy_3plot(img,md, res, 
                               centre,r_max=r_max, name= nam, alp=15)
        plt.ioff()
    return None

#plot_all_gal_3plot()

# list all CF3 distance 

import pycf3

table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_2.txt")
table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_2.txt",
    dtype=str)
name = table_n[:,0]
ra, dec = table[:,1], table[:,2]
vel = table[:,-1]#helio
dist = table[:,24]


def cal_CF3_distance_list():
    cf3 = pycf3.CF3() 
    NAM = pycf3.NAM()
    
    for i in range(len(name)):
      ## or nam = pycf3.NAM()
      if vel[i] > 2400:
          result = cf3.calculate_distance(velocity=vel[i], ra=ra[i], dec=dec[i])
          print(name[i], dist[i], result.observed_distance_)

      elif vel[i] <2400:
          result = NAM.calculate_distance(velocity=vel[i], ra=ra[i], dec=dec[i])
          print(name[i], dist[i], result.observed_distance_)
    return None


cal_CF3_distance_list()

