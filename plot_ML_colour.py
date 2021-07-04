#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:45:13 2021

@author: dexter

This script is for plotting the M/L ratio and (g-i) colour vs mass 
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


"""
Input file

"""

master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"

#master_file="/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt"

broad_cut_file="/home/dexter/result/stat/completeness/diagonal_selection_bag3_2.dat"

master = SRead.read_table(master_file)
master_n = SRead.read_table(master_file,dtype='str')

broad_cut = SRead.read_table(broad_cut_file)
broad_cut_str = SRead.read_table(broad_cut_file,dtype="str")


dc = master[:,3]
initial_mass = master[:,4]

name_b = broad_cut_str[:,4]
dc_b = broad_cut[:,2]
initial_mass_b = broad_cut[:,3]

# stich up M2000 and z-indepent distance
D = SRead.grab_dist(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",                    
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_2.dat",
              name_index = 4, vel_index = 6, vel_err_index = 7, scale_index = 8,
              dir_name_index = 0, dir_dist_index = 3, dir_dist_err_index =4, 
              dir_method_index = 5, dir_method_flag_index = 7)

name1 = D["Gal_name"]
DD1 = D["Dist"]
DD1_err = D["Dist_err"]
scale1 = D["Scale"]


mag_g, mag_i = broad_cut[:,10],broad_cut[:,9]

ML_select_T11= SPlot.MLRelationIband(mag_g,mag_i).Taylor11_MassRatio
ML_select_Z09= SPlot.MLRelationIband(mag_g,mag_i).Zibetti09_MassRatio
ML_select_RC15= SPlot.MLRelationIband(mag_g,mag_i).Roediger15BC03_MassRatio
ML_select_IP13= SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio


M = SPlot.MassCalculation(mag_i, DD1, 4.53,mag_g,mag_i)

E_IP13 = M.cal_Mass(ML_select_IP13)
E_T11 = M.cal_Mass(ML_select_T11)
E_Z09 = M.cal_Mass(ML_select_Z09)
E_RC15 = M.cal_Mass(ML_select_RC15)

from scipy.optimize import curve_fit

    
E_IP13n, ML_select_IP13n = np.nan_to_num(np.log10(E_IP13)), np.nan_to_num(np.log10(ML_select_IP13))
E_T11n, ML_select_T11n = np.nan_to_num(np.log10(E_T11)), np.nan_to_num(np.log10(ML_select_T11))
E_Z09n, ML_select_Z09n = np.nan_to_num(np.log10(E_Z09)), np.nan_to_num(np.log10(ML_select_Z09))
E_RC15n, ML_select_RC15n = np.nan_to_num(np.log10(E_RC15)), np.nan_to_num(np.log10(ML_select_RC15))

popt_T11,pcov_T11 = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                E_T11n, ML_select_T11n)
popt_Z09,pcov_Z09 = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                E_Z09n, ML_select_Z09n)
popt_RC15,pcov_RC15 = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                E_RC15n, ML_select_RC15n)
popt_IP13,pcov_IP13 = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                E_IP13n, ML_select_IP13n)

print(E_IP13n,ML_select_IP13n)
print(*popt_IP13, *pcov_IP13)
"""
Plotting

"""

x = np.linspace(0.5e9, 2e13,num=50)

def plot_ML_mass(x,y):
    fig = plt.figure()
    ax0 = plt.subplot()
    
    # data points
    
    ax0.plot(E_T11,ML_select_T11,'o',ms=8, color = "r", alpha = 0.4,label='$\mathrm{T11}$')
    ax0.plot(E_Z09,ML_select_Z09,'o',ms=8, color = "b", alpha = 0.4,label='$\mathrm{Z09}$')
    ax0.plot(E_RC15,ML_select_RC15,'o',ms=8, color = "k", alpha = 0.4, label='$\mathrm{RC15}$')
    ax0.plot(E_IP13,ML_select_IP13,'o',ms=8, color = "g", alpha = 0.4, label='$\mathrm{IP13}$')
    
    
    # the line

    ax0.set_ylabel("$ M_*/L$", fontsize=16)
    ax0.set_xlabel(r"$ M_*/\rm M_{\odot}$", fontsize=16)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    
    ax0.set_xlim(1e9,1e13)
    ax0.set_ylim(0.5,6)
    
    line_T11 = SAna.AnalyticFunctions.linear_func_1D(x,*list(popt_T11))
    line_Z09 = SAna.AnalyticFunctions.linear_func_1D(x,*list(popt_Z09))
    line_RC15 = SAna.AnalyticFunctions.linear_func_1D(x,*list(popt_RC15))
    line_IP13 = SAna.AnalyticFunctions.linear_func_1D(x,*list(popt_IP13))
    
    ax0.plot(x,line_T11,'-', color = "r")
    ax0.plot(x,line_Z09,'-', color = "b")
    ax0.plot(x,line_RC15,'-', color = "k")
    ax0.plot(x,line_IP13,'-', color = "g")   

    ax0.legend(loc="lower right")
    
    plt.show()
    return fig

plot_ML_mass(E_T11,ML_select_T11)