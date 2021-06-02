#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 16:54:05 2021

@author: dexter

This script is dedicated to 
analyzing the surface brightness (SB) profile of galaxies

This script will produce the following data:
1) The various radius for the spheroids

This script also plot:
1) The illustration of percentage covered by the calculation


"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAna
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

import matplotlib.gridspec as gridspec

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

#scan, trim, info

NGC4772SB ='/home/dexter/result/NGC4772/out1_NGC4772.dat'

NGC3646SB = '/home/dexter/result/NGC3646/out1_NGC3646.dat'

NGC5382SB = '/home/dexter/result/NGC5382/out1_NGC5382.dat'


pix = SRead.read_table(NGC4772SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R = SRead.read_table(NGC4772SB)[:,0]
R_s = SRead.read_table(NGC4772SB)[:,0]*0.4
e = SRead.read_table(NGC4772SB)[:,5]
Rmax= max(R)

R_e = SAna.Isophote.circularized(R_s, e)

Rmax_e = SAna.Isophote.circularized(Rmax,e[-1])

mu = SAna.Isophote.pix_val_to_mu(pix)

print(Rmax_e)
#print("totalmag",cal_SB_mag(R,pix,Rmax,e)) #10.3

print("NGC4772 totalmag",SAna.Isophote.cal_SB_mag(R,mu,max(R),e),'right answer',10.3) #10.3 #no conv 10.299

pix_1 = SRead.read_table(NGC3646SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R_1 = SRead.read_table(NGC3646SB)[:,0]
R_s_1 = SRead.read_table(NGC3646SB)[:,0]*0.4
e_1 = SRead.read_table(NGC3646SB)[:,5]
Rmax_1= max(R_1)

R_e_1 = SAna.Isophote.circularized(R_s_1, e_1)

Rmax_e_1 = SAna.Isophote.circularized(Rmax_1,e_1[-1])

mu_1 = SAna.Isophote.pix_val_to_mu(pix_1)

#29
print("NGC3646 totalmag",SAna.Isophote.cal_SB_mag(R_1,mu_1,38/0.4,e_1),'right answer',11.07) #11.07 # no cov

pix_2 = SRead.read_table(NGC5382SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R_2 = SRead.read_table(NGC5382SB)[:,0]
R_s_2 = SRead.read_table(NGC5382SB)[:,0]*0.4
e_2 = SRead.read_table(NGC5382SB)[:,5]
Rmax_2 = 46

R_e_2 = SAna.Isophote.circularized(R_s_2, e_2)

Rmax_e_2 = SAna.Isophote.circularized(Rmax_2,e_2[-1])

mu_2 = SAna.Isophote.pix_val_to_mu(pix_2)

#42
print("NGC5382 totalmag",SAna.Isophote.cal_SB_mag(R_2,mu_2,42.0/0.4,e_2),'right answer', 11.54) #11.54 # no cov 11.4695

def scan_basic1D(input_array_x, input_array_y, 
                 percentage=1.0, start_point=0.5):
    # Find max
    
    # start with start_point, find area, 
    return None #x,y coordinate for the percentage