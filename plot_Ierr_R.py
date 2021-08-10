#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 23:43:07 2021

@author: dexter
"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

import matplotlib.gridspec as gridspec



# read the list

filelist0 = "/home/dexter/result/stat/completeness/gal_output_list_all.dat"

geomlist = "/home/dexter/result/stat/completeness/gal_geom_all.dat"

Rmax = SRead.read_table(geomlist)[:,1] 
#Re = SRead.read_table(geomlist)[:,2] 

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0

def read_Ierr_R(filelist):
    """
    Read the output file from ISOFIT.
    Store R(pixel), I2, I_err (sigma/Npix) in to an matrix

    Returns
    -------
    item_bundle:
        [
        [[r1,r2..],[I1/Ierr1,I2/Ierr2,...]] # galaxy a
        [[r1,r2..],[I1/Ierr1,I2/Ierr2,...]] # galaxy b     
        ...
        ]
    """
    item_bundle = []
    #open the list of output
    filelist = SRead.read_table(filelist,dtype='str')
    
    for i in filelist:
        item_list, I_Ierr = [], []
        
        #loop through the list, read individual file
        item = SRead.read_table(i)
        R_arcsec = item[:,0]*0.4 #convert from pix to arcsec
        I2, pix_var, NPIX = item[:,1], item[:,3], item[:,25]
        
        I_err = np.sqrt(pix_var/NPIX)
        
        I_Ierr.append(I_err/I2)
        
        item_list.append(R_arcsec)
        item_list.append(I_Ierr)
        
        item_bundle.append(item_list)
        
    return item_bundle


def plot_Ierr_R(item_bundle, Rnorm):
    """
    Plotting the Ierr/I- R/Rnorm.

    Parameters
    ----------
    item_bundle : list
        The matrix that store the information of .
    Rnorm : list
        A 1D list that contain the normalization radius Rnorm.

    Returns
    -------
    plot.

    """
    
    fig = plt.figure(figsize=(6.4, 4.8))

    I_dump = []
    for i in range(len(item_bundle)):
        #print(type(item_bundle[i][0][0]))
#        print(Rnorm[i])
        
        R = np.array(item_bundle[i][0])/ np.full((1,len(item_bundle[i][0])),Rnorm[i])
        I = np.array(item_bundle[i][1])
        
        # dumping all Ierr/I into one list to calculate sigma
        #I_dump = I_dump + list(item_bundle[i][1])
        
        plt.plot(R,I,'bo',ms=0.3)
    s_I = np.std(I_dump)
    
    plt.hlines(0, 0, 1) #centre line at y=0
    
    plt.hlines(s_I, 0, 1)
    
    plt.ylim(0, 0.05)
    plt.xlim(0, 1.0)

    plt.xlabel(r"$R/\rm R_{max}$", fontsize=16)
    plt.ylabel(r"$|I_\mathrm{err}|/I$", fontsize=16)    
    
    plt.tight_layout()
    return None

#print(read_Ierr_R(filelist0))

plot_Ierr_R(read_Ierr_R(filelist0), Rmax)

