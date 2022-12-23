#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 19:20:11 2020

@author: dexter
"""
from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt",
    dtype = 'str')

D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")


mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]


#################################
#plot percentage

A ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
Y = D0_all_table[:,29]

SPlot.ShowcaseIndi.plot_hist_percentage(A, Y)

plt.show()

A1 ="/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt"
Y1 = D0_Bin1_table[:,29]

SPlot.ShowcaseIndi.plot_hist_percentage(A1, Y1)

plt.show()

A2 ="/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt"
Y2 = D0_Bin2_table[:,29]

SPlot.ShowcaseIndi.plot_hist_percentage(A2, Y2)

plt.show()

A3 ="/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt"
Y3 = D0_Bin3_table[:,29]

SPlot.ShowcaseIndi.plot_hist_percentage(A3, Y3)

plt.show()