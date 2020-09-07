#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 30 01:30:40 2020

@author: dexter
"""

from astropy.cosmology import FlatLambdaCDM

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)




#######################################################
################reading the new set of data############
#######################################################


#Plot the distances showcase
#############################
D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt")


D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt",
    dtype = 'str')


# name for the galaxies
d_name1 = D0_Bin1_table_n[:,0]
d_name2 = D0_Bin2_table_n[:,0]
d_name3 = D0_Bin3_table_n[:,0]

# The corrected distance, Mould 2000
DDD1 = D0_Bin1_table[:,7] / 68.0
DDD2 = D0_Bin2_table[:,7] / 68.0
DDD3 = D0_Bin3_table[:,7] / 68.0

# The corrected distance error, Mould 2000

DDD1_err = D0_Bin1_table[:,8] / 68.0
DDD2_err = D0_Bin2_table[:,8] / 68.0
DDD3_err = D0_Bin3_table[:,8] / 68.0

#dc## Willick2007 redshift , h68

ddc1 = D0_Bin1_table[:,3]
ddc2 = D0_Bin2_table[:,3]
ddc3 = D0_Bin3_table[:,3]

ddc1_err = D0_Bin1_table[:,28]
ddc2_err = D0_Bin2_table[:,28]
ddc3_err = D0_Bin3_table[:,28]

corr_dist1, corr_E1 = D0_Bin1_table[:,12],  D0_Bin1_table[:,13]
corr_dist2, corr_E2 = D0_Bin2_table[:,12],  D0_Bin2_table[:,13]
corr_dist3, corr_E3 = D0_Bin3_table[:,12],  D0_Bin3_table[:,13]



# Cosmicflow-3 distance
d1 = D0_Bin1_table[:,24] * (75.0/75)
d2 = D0_Bin2_table[:,24] * (75.0/75)
d3 = D0_Bin3_table[:,24] * (75.0/75)

d1_uerr = (D0_Bin1_table[:,25] - D0_Bin1_table[:,24])*(75.0/75)
d2_uerr = (D0_Bin2_table[:,25] - D0_Bin2_table[:,24])*(75.0/75)
d3_uerr = (D0_Bin3_table[:,25] - D0_Bin3_table[:,24])*(75.0/75)

d1_lerr = (D0_Bin1_table[:,24]- D0_Bin1_table[:,23]) * (75.0/75)
d2_lerr = (D0_Bin2_table[:,24]- D0_Bin2_table[:,23]) * (75.0/75)
d3_lerr = (D0_Bin3_table[:,24]- D0_Bin3_table[:,23]) * (75.0/75)



# special distance, z-indepnednet
d_spec_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",
    dtype="float")
d_spec_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",
    dtype="str")

d_spec = d_spec_table[:,3] 
d_spec_name = d_spec_table_n[:,0] 
d_spec_err = d_spec_table[:,4]


d_spec_index, d_spec_d, d_spec_d_err = [],[],[]

decision_table1 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin1.dat", dtype='str')
decision_table2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin2.dat", dtype='str')
decision_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin3.dat", dtype='str')

decision1 = decision_table1[:,1]
decision2 = decision_table2[:,1]
decision3 = decision_table3[:,1]


def plot_dist_difference2():
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD1,ddc1,d1,d_name1,110, 
                                                 DD_err=DDD1_err,
                                                 dc_err=ddc1_err, 
                                                 d_err= [d1_lerr,d1_uerr],                                                 
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD2,ddc2,d2,d_name2,75, 
                                                 DD_err=DDD2_err,
                                                 dc_err=ddc2_err, 
                                                 d_err=[d2_lerr,d2_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD3,ddc3,d3,d_name3,45, 
                                                 DD_err=DDD3_err,
                                                 dc_err=ddc3_err, 
                                                 d_err=[d3_lerr,d3_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name , 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    plt.show()
    
plot_dist_difference2()
########################################################



###############################################
# plot data selection

from astropy.io import fits


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

ML_select= SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
M = SPlot.MassCalculation(mag_i, DD1, 4.53,mag_g,mag_i)
E = M.cal_Mass(ML_select)


import statistics
import matplotlib.gridspec as gridspec

fig = plt.figure()

gs = gridspec.GridSpec(2,1) 

axs1 = plt.subplot(gs[1])
axs0 = plt.subplot(gs[0])


#SPlot.ShowcaseIndi.show_name(DD1,E,name_b)


axs1.hlines(2e12,0,110,linestyle="solid")
axs1.hlines(5e11,0,110,linestyle="solid")
axs1.hlines(2e11,0,75,linestyle="solid")
axs1.hlines(1e11,0,45,linestyle="solid")

axs1.vlines(110,5e11,2e12,linestyle="solid")
axs1.vlines(75,2e11,5e11,linestyle="solid")
axs1.vlines(45,1e11,2e11,linestyle="solid")

x_edge1,y_edge1= [0,0,110,110], [5e11,2e12,2e12,5e11]
x_edge2,y_edge2= [0,0,75,75], [2e11,5e11,5e11,2e11]
x_edge3,y_edge3= [0,0,45,45], [1e11,2e11,2e11,1e11]

axs0.plot(dc,initial_mass,'x',color='grey', alpha=0.7, label='All')

axs0.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)

axs0.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)
#axs0.plot(dc,E,'x',color='green', alpha=0.4, label='New Distance')
axs0.plot(dc_b,initial_mass_b,'x',color='green', alpha=0.4, label='Broad Selection')


axs1.plot(DD1,E,'x',color='#a51a74', alpha=0.4, label='New Distance')

axs1.plot(corr_dist1,corr_E1,'ro',label='Bin1 samples', ms=12,alpha=0.5)
axs1.plot(corr_dist2,corr_E2,'bo',label='Bin2 samples', ms=12,alpha=0.5)
axs1.plot(corr_dist3,corr_E3,'ko',label='Bin3 samples', ms=12,alpha=0.5)

axs1.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)
axs1.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)

axs1.fill(x_edge1,y_edge1, alpha=0.3, color='green')
axs1.fill(x_edge2,y_edge2, alpha=0.2, color='green')
axs1.fill(x_edge3,y_edge3, alpha=0.1, color='green')
axs0.set_ylim(top =5e12 , bottom = 1e9)
axs0.set_xlim(left=0, right = 120)


axs0.set_ylabel("$M_*/M_\odot(Into Mass)$",fontsize=24)
axs0.set_xlabel("$Distance(H_0 = 68, \Omega_m =0.3)/Mpc$",fontsize=24)
axs0.set_yscale( 'log' )

axs1.set_ylabel("$M_*/M_\odot(Into Mass)$",fontsize=24)
axs1.set_xlabel("$Distance(H_0 = 68, \Omega_m =0.3)/Mpc$",fontsize=24)
axs1.set_yscale( 'log' )

axs1.set_ylim(top =5e13 , bottom = 1e10)
axs1.set_xlim(left=0, right = 120)

axs0.legend()
axs1.legend()
plt.show()
