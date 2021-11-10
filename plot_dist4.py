#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 04:56:31 2021

@author: dexter
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Vertical Bins, Kcorrection+EXTINCTION

from astropy.cosmology import FlatLambdaCDM

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)

#plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

Bin1_lowerlimit = 4.75e11
Bin2_lowerlimit = 1.9e11
Bin3_lowerlimit = 9.95e10

#######################################################
################reading the new set of data############
#######################################################
#read K correction extinction files
K_table1 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1V_Kcorr_EXT.dat")
K_table2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2V_Kcorr_EXT.dat")
K_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr_EXT.dat")

K_table1_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1V_Kcorr_EXT.dat", 
    dtype='str')
K_table2_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2V_Kcorr_EXT.dat",
    dtype='str')
K_table3_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr_EXT.dat", 
    dtype='str')

K_name1, K_name2, K_name3 = K_table1_n[:,4], K_table2_n[:,4], K_table3_n[:,4]

mag_g1, mag_i1 = K_table1[:,10], K_table1[:,9]
mag_g2, mag_i2 = K_table2[:,10], K_table2[:,9]
mag_g3, mag_i3 = K_table3[:,10], K_table3[:,9]

mag_g1_kcorr, mag_i1_kcorr = K_table1[:,19], K_table1[:,18]
mag_g2_kcorr, mag_i2_kcorr = K_table2[:,19], K_table2[:,18]
mag_g3_kcorr, mag_i3_kcorr = K_table3[:,19], K_table3[:,18]

g1_EXT, i1_EXT = K_table1[:,23], K_table1[:,24]
g2_EXT, i2_EXT = K_table2[:,23], K_table2[:,24]
g3_EXT, i3_EXT = K_table3[:,23], K_table3[:,24]

g1_kcorr, i1_kcorr = K_table1[:,25], K_table1[:,26]
g2_kcorr, i2_kcorr = K_table2[:,25], K_table2[:,26]
g3_kcorr, i3_kcorr = K_table3[:,25], K_table3[:,26]

# the corrected mag g and i, Kcorrection+EXTINCTIOn
mag_g1_corr, mag_i1_corr = mag_g1-g1_kcorr, mag_i1-i1_kcorr
mag_g2_corr, mag_i2_corr = mag_g2-g2_kcorr, mag_i2-i2_kcorr
mag_g3_corr, mag_i3_corr = mag_g3-g3_kcorr, mag_i3-i3_kcorr



#Plot the distances showcase
#############################
D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_3.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_3.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_3.txt")


D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_3.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_3.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_3.txt",
    dtype = 'str')


# name for the galaxies
d_name1 = D0_Bin1_table_n[:,0]
d_name2 = D0_Bin2_table_n[:,0]
d_name3 = D0_Bin3_table_n[:,0]

# The corrected distance, Mould 2000 (vel/68)
DDD1 = D0_Bin1_table[:,7] / 68.0
DDD2 = D0_Bin2_table[:,7] / 68.0
DDD3 = D0_Bin3_table[:,7] / 68.0

# The corrected distance error, Mould 2000 (vel/68)
DDD1_err = D0_Bin1_table[:,8] / 68.0
DDD2_err = D0_Bin2_table[:,8] / 68.0
DDD3_err = D0_Bin3_table[:,8] / 68.0

#dc## Willick2007 redshift , h68 (h68_dist)
ddc1 = D0_Bin1_table[:,3]
ddc2 = D0_Bin2_table[:,3]
ddc3 = D0_Bin3_table[:,3]

ddc1_err = D0_Bin1_table[:,28]
ddc2_err = D0_Bin2_table[:,28]
ddc3_err = D0_Bin3_table[:,28]

g1 , i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
g2 , i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
g3 , i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]


#SPlot.ShowcaseCompare2.plot_compare_generic(g1, i1,para_name="Bin1",label=["g","i"])
#SPlot.ShowcaseCompare2.plot_compare_generic(g2, i2,para_name="Bin1",label=["g","i"])
#SPlot.ShowcaseCompare2.plot_compare_generic(g3, i3,para_name="Bin1",label=["g","i"])

corr_dist1, corr_E1 = D0_Bin1_table[:,12],  D0_Bin1_table[:,13]
corr_dist2, corr_E2 = D0_Bin2_table[:,12],  D0_Bin2_table[:,13]
corr_dist3, corr_E3 = D0_Bin3_table[:,12],  D0_Bin3_table[:,13]

# Cosmicflow-3 distance (dist in the table)
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
    "/home/dexter/result/stat/completeness/dist_judgement_Bin1V.dat", dtype='str')
decision_table2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin2V.dat", dtype='str')
decision_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin3V.dat", dtype='str')

decision1 = decision_table1[:,1]
decision2 = decision_table2[:,1]
decision3 = decision_table3[:,1]


# The final distance we used
D1 = D0_Bin1_table[:,29]
D2 = D0_Bin2_table[:,29]
D3 = D0_Bin3_table[:,29]

Mag_i1_kcorr_cDis = mag_i1_corr-25-5*np.log10(corr_dist1) 
Mag_i2_kcorr_cDis = mag_i2_corr-25-5*np.log10(corr_dist2) 
Mag_i3_kcorr_cDis = mag_i3_corr-25-5*np.log10(corr_dist3) 

Mag_i1_kcorr = mag_i1_corr-25-5*np.log10(D1) 
Mag_i2_kcorr = mag_i2_corr-25-5*np.log10(D2) 
Mag_i3_kcorr = mag_i3_corr-25-5*np.log10(D3) 


########################################################

def plot_dist_difference2():
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD1,ddc1,d1,d_name1,
                                                 75, 110, 
                                                 DD_err=DDD1_err,
                                                 dc_err=ddc1_err, 
                                                 d_err= [d1_lerr,d1_uerr],                                                 
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD2,ddc2,d2,d_name2,
                                                 45, 75, 
                                                 DD_err=DDD2_err,
                                                 dc_err=ddc2_err, 
                                                 d_err=[d2_lerr,d2_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD3,ddc3,d3,d_name3,
                                                 0, 45, 
                                                 DD_err=DDD3_err,
                                                 dc_err=ddc3_err, 
                                                 d_err=[d3_lerr,d3_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name , 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])
    plt.show()
    
#plot_dist_difference2()
########################################################
##plot all dist dist in one single big graph


# Stich all array together 

DDD = np.concatenate((DDD1, DDD2, DDD3))
ddc = np.concatenate((ddc1, ddc2, ddc3))
ddc_err= np.concatenate((ddc1_err, ddc2_err, ddc3_err)) 
d = np.concatenate((d1, d2, d3))
d_name = np.concatenate((d_name1, d_name2, d_name3))
DDD_err = np.concatenate((DDD1_err, DDD2_err, DDD3_err))
d_lerr = np.concatenate((d1_lerr, d2_lerr, d3_lerr))
d_uerr = np.concatenate((d1_uerr, d2_uerr, d3_uerr))


SPlot.ShowcaseCompare2.plot_dist_difference3_summary(DDD,ddc,d,d_name,
                                                 DD_err=DDD_err,
                                                 dc_err=ddc_err, 
                                                 d_err= [d_lerr,d_uerr],                                                 
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = [])

###############################################
# plot data selection


master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"

#master_file="/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt"


broad_cut_file="/home/dexter/result/stat/completeness/diagonal_selection_bag3_4.dat"


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


# calculate the corrected magnitude from the 708 broadcut objects
mag_g_corr = broad_cut[:,10]-broad_cut[:,13]#-broad_cut[:,11]
mag_i_corr = broad_cut[:,9]-broad_cut[:,14]#-broad_cut[:,12]

# calculate their mass-to-light ratio
ML_select_T11= SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Taylor11_MassRatio
ML_select_Z09= SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Zibetti09_MassRatio
ML_select_RC15= SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Roediger15BC03_MassRatio
ML_select_IP13= SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Into13_MassRatio


M = SPlot.MassCalculation(mag_i_corr, DD1, 4.53,mag_g_corr,mag_i_corr)

# calculate thier mass
E = M.cal_Mass(ML_select_IP13)
E_T11 = M.cal_Mass(ML_select_T11)
E_Z09 = M.cal_Mass(ML_select_Z09)
E_RC15 = M.cal_Mass(ML_select_RC15)

#selects for the galaxies in the three different Bins

#for Bin1
QQQ1 = SSort.zone_in_2D([E,DD1],low_bound=[5e11,75.0],up_bound=[130e11,110])
QQQ2 = SSort.zone_in_2D([E,DD1],low_bound=[2e11,45.0],up_bound=[130e11,75.0])
QQQ3 = SSort.zone_in_2D([E,DD1],low_bound=[1e11,0.0],up_bound=[130e11,45.0])

def list_selection(Q):
    for i in range(len(Q['index'])):
        print(name_b[Q['index']][i])
        
list_selection(QQQ1)
print("=========================")
list_selection(QQQ2)
print("=========================")
list_selection(QQQ3)
print("=========================")

# same mass calculation but seperated by bin
ML_select1_T11= SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Taylor11_MassRatio
ML_select1_Z09= SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Zibetti09_MassRatio
ML_select1_RC15= SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Roediger15BC03_MassRatio
ML_select1_IP13= SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Into13_MassRatio

ML_select2_T11= SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Taylor11_MassRatio
ML_select2_Z09= SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Zibetti09_MassRatio
ML_select2_RC15= SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Roediger15BC03_MassRatio
ML_select2_IP13= SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Into13_MassRatio

ML_select3_T11= SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Taylor11_MassRatio
ML_select3_Z09= SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Zibetti09_MassRatio
ML_select3_RC15= SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Roediger15BC03_MassRatio
ML_select3_IP13= SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Into13_MassRatio


M1 = SPlot.MassCalculation(mag_i1_corr, corr_dist1, 4.53,mag_g1_corr, mag_i1_corr)
M2 = SPlot.MassCalculation(mag_i2_corr, corr_dist2, 4.53,mag_g2_corr, mag_i2_corr)
M3 = SPlot.MassCalculation(mag_i3_corr, corr_dist3, 4.53,mag_g3_corr, mag_i3_corr)

# use Cosmicflow-3 distance
M1_d = SPlot.MassCalculation(mag_i1_corr, d1, 4.53,mag_g1_corr, mag_i1_corr)
M2_d = SPlot.MassCalculation(mag_i2_corr, d2, 4.53,mag_g2_corr, mag_i2_corr)
M3_d = SPlot.MassCalculation(mag_i3_corr, d3, 4.53,mag_g3_corr, mag_i3_corr)

E1, E2, E3 = M1.cal_Mass(ML_select1_IP13), M2.cal_Mass(ML_select2_IP13), M3.cal_Mass(ML_select3_IP13)
E1_T11, E2_T11, E3_T11 = M1.cal_Mass(ML_select1_T11), M2.cal_Mass(ML_select2_T11), M3.cal_Mass(ML_select3_T11)
E1_Z09, E2_Z09, E3_Z09 = M1.cal_Mass(ML_select1_Z09), M2.cal_Mass(ML_select2_Z09), M3.cal_Mass(ML_select3_Z09)
E1_RC15, E2_RC15, E3_RC15  = M1.cal_Mass(ML_select1_RC15), M2.cal_Mass(ML_select2_RC15), M3.cal_Mass(ML_select3_RC15)

E1, E2, E3 = np.log10(E1),np.log10(E2),np.log10(E3)
E1_T11, E2_T11, E3_T11 = np.log10(E1_T11),np.log10(E2_T11),np.log10(E3_T11)
E1_Z09, E2_Z09, E3_Z09 = np.log10(E1_Z09),np.log10(E2_Z09),np.log10(E3_Z09)
E1_RC15, E2_RC15, E3_RC15  = np.log10(E1_RC15), np.log10(E2_RC15), np.log10(E3_RC15)

E1_d, E2_d, E3_d = M1_d.cal_Mass(ML_select1_IP13), M2_d.cal_Mass(ML_select2_IP13), M3_d.cal_Mass(ML_select3_IP13)
E1_T11_d, E2_T11_d, E3_T11_d = M1_d.cal_Mass(ML_select1_T11), M2_d.cal_Mass(ML_select2_T11), M3_d.cal_Mass(ML_select3_T11)
E1_Z09_d, E2_Z09_d, E3_Z09_d = M1_d.cal_Mass(ML_select1_Z09), M2_d.cal_Mass(ML_select2_Z09), M3_d.cal_Mass(ML_select3_Z09)
E1_RC15_d, E2_RC15_d, E3_RC15_d  = M1_d.cal_Mass(ML_select1_RC15), M2_d.cal_Mass(ML_select2_RC15), M3_d.cal_Mass(ML_select3_RC15)

E1_d, E2_d, E3_d = np.log10(E1_d),np.log10(E2_d),np.log10(E3_d)
E1_T11_d, E2_T11_d, E3_T11_d = np.log10(E1_T11_d),np.log10(E2_T11_d),np.log10(E3_T11_d)
E1_Z09_d, E2_Z09_d, E3_Z09_d = np.log10(E1_Z09_d),np.log10(E2_Z09_d),np.log10(E3_Z09_d)
E1_RC15_d, E2_RC15_d, E3_RC15_d  = np.log10(E1_RC15_d), np.log10(E2_RC15_d), np.log10(E3_RC15_d)


print('min',min(E1-E1_T11),min(E1-E1_Z09), min(E1-E1_RC15))
print('max',max(E1-E1_T11),max(E1-E1_Z09), max(E1-E1_RC15))
print('delta',max(E1-E1_T11) - min(E1-E1_T11),max(E1-E1_Z09) -min(E1-E1_Z09), max(E1-E1_RC15) - min(E1-E1_RC15))
        # look for ultra-massive sample
def look_for_ultramassive():
    print('-----------------------')
    for i in range(len(name1)):
        if DD1[i] <110 and E[i] >2e12:
            print(name1[i],DD1[i],E[i], mag_g[i],mag_i[i])
    print('-----------------------')


import statistics
import matplotlib.gridspec as gridspec

Bin1_limit = np.repeat(Bin1_lowerlimit,len(E1))
Bin2_limit = np.repeat(Bin2_lowerlimit,len(E2))
Bin3_limit = np.repeat(Bin3_lowerlimit,len(E3))

key_delta_Bin1 = np.log10((10**E1)-Bin1_limit)
key_delta_Bin2 = np.log10((10**E2)-Bin2_limit)
key_delta_Bin3 = np.log10((10**E3)-Bin3_limit)

delta_E_ET11, s_E_ET11 = np.median(np.nan_to_num(abs(E-E_T11))),np.std(np.nan_to_num(abs(E-E_T11)))
delta_E_EZ09, s_E_EZ09 =  np.median(np.nan_to_num(abs(E-E_Z09))),np.std(np.nan_to_num(abs(E-E_Z09)))
delta_E_ERC15, s_E_ERC15 = np.median(np.nan_to_num(abs(E-E_RC15))),np.std(np.nan_to_num(abs(E-E_RC15)))

print('============')
print(E1_T11-key_delta_Bin1)
print(E2_T11-key_delta_Bin2)
print(E3_T11-key_delta_Bin3)
print('============')
print(np.std(E1_T11-key_delta_Bin1))
print(np.std(E2_T11-key_delta_Bin2))
print(np.std(E3_T11-key_delta_Bin3))
print('============')
print(np.average(E1_T11-key_delta_Bin1))
print(np.average(E2_T11-key_delta_Bin2))
print(np.average(E3_T11-key_delta_Bin3))
print('============')

#print("delta(E-E_T11)", np.log10(delta_E_ET11), np.log10(s_E_ET11))
#print("delta(E-E_TZ09)", np.log10(delta_E_EZ09), np.log10(s_E_EZ09))
#print("delta(E-E_TRC15)", np.log10(delta_E_ERC15), np.log10(s_E_ERC15))


#print("Bin1:", "T11",5e11-delta_E_ET11, "Z09", 5e11-delta_E_EZ09,"IP13", 5e11-delta_E_ERC15)
#print("Bin2:", "T11",2e11-delta_E_ET11, "Z09", 2e11-delta_E_EZ09,"IP13", 2e11-delta_E_ERC15)
#print("Bin3:", "T11",1e11-delta_E_ET11,"Z09", 1e11-delta_E_EZ09,"IP13", 1e11-delta_E_ERC15)

delta_E1_ET11, s_E1_ET11 = np.median(np.nan_to_num(abs(E1-E1_T11))),np.std(np.nan_to_num(abs(E1-E1_T11)))
delta_E1_EZ09, s_E1_EZ09 =  np.median(np.nan_to_num(abs(E1-E1_Z09))),np.std(np.nan_to_num(abs(E1-E1_Z09)))
delta_E1_ERC15, s_E1_ERC15 = np.median(np.nan_to_num(abs(E1-E1_RC15))),np.std(np.nan_to_num(abs(E1-E1_RC15)))

delta_E2_ET11, s_E2_ET11 = np.median(np.nan_to_num(abs(E2-E2_T11))),np.std(np.nan_to_num(abs(E2-E2_T11)))
delta_E2_EZ09, s_E2_EZ09 =  np.median(np.nan_to_num(abs(E2-E2_Z09))),np.std(np.nan_to_num(abs(E2-E2_Z09)))
delta_E2_ERC15, s_E2_ERC15 = np.median(np.nan_to_num(abs(E2-E2_RC15))),np.std(np.nan_to_num(abs(E2-E2_RC15)))

delta_E3_ET11, s_E3_ET11 = np.median(np.nan_to_num(abs(E3-E3_T11))),np.std(np.nan_to_num(abs(E3-E3_T11)))
delta_E3_EZ09, s_E3_EZ09 =  np.median(np.nan_to_num(abs(E3-E3_Z09))),np.std(np.nan_to_num(abs(E3-E3_Z09)))
delta_E3_ERC15, s_E3_ERC15 = np.median(np.nan_to_num(abs(E3-E3_RC15))),np.std(np.nan_to_num(abs(E3-E3_RC15)))


max_E1_ET11, min_E1_ET11 = max(E1-E1_T11),min(E1-E1_T11)
max_E1_EZ09, min_E1_EZ09 =  max(E1-E1_Z09),min(E1-E1_Z09)
max_E1_ERC15, min_E1_ERC15 = max(E1-E1_RC15),min(E1-E1_RC15)

max_E2_ET11, min_E2_ET11 = max(E2-E2_T11),min(E2-E2_T11)
max_E2_EZ09, min_E2_EZ09 =  max(E2-E2_Z09),min(E2-E2_Z09)
max_E2_ERC15, min_E2_ERC15 = max(E2-E2_RC15),min(E2-E2_RC15)

max_E3_ET11, min_E3_ET11 = max(E3-E3_T11),min(E3-E3_T11)
max_E3_EZ09, min_E3_EZ09 =  max(E3-E3_Z09),min(E3-E3_Z09)
max_E3_ERC15, min_E3_ERC15 = max(E3-E3_RC15),min(E3-E3_RC15)

print('----------------------------')
print("delta(E-E_T11)_1", (delta_E1_ET11), (s_E1_ET11))
print("delta(E-E_TZ09)_1", (delta_E1_EZ09), (s_E1_EZ09))
print("delta(E-E_TRC15)_1", (delta_E1_ERC15), (s_E1_ERC15))
print('----------------------------')
print("delta(E-E_T11)_2", (delta_E2_ET11), (s_E2_ET11))
print("delta(E-E_TZ09)_2", (delta_E2_EZ09), (s_E2_EZ09))
print("delta(E-E_TRC15)_2", (delta_E2_ERC15), (s_E2_ERC15))
print('----------------------------')
print("delta(E-E_T11)_3", (delta_E3_ET11), (s_E3_ET11))
print("delta(E-E_TZ09)_3", (delta_E3_EZ09), (s_E3_EZ09))
print("delta(E-E_TRC15)_3", (delta_E3_ERC15), (s_E3_ERC15))
print('----------------------------')

#print("Bin1:", "T11",5e11-delta_E1_ET11, "Z09", 5e11-delta_E1_EZ09,"IP13", 5e11-delta_E1_ERC15)
#print("Bin2:", "T11",2e11-delta_E2_ET11, "Z09", 2e11-delta_E2_EZ09,"IP13", 2e11-delta_E2_ERC15)
#print("Bin3:", "T11",1e11-delta_E3_ET11,"Z09", 1e11-delta_E3_EZ09,"IP13", 1e11-delta_E3_ERC15)

print("Bin1:", "T11",np.log10(4.39e11)-delta_E1_ET11, "Z09", 4.39e11-delta_E1_EZ09,"IP13", 4.39e11-delta_E1_ERC15)
print("Bin2:", "T11",1.8e11-delta_E2_ET11, "Z09", 1.8e11-delta_E2_EZ09,"IP13", 1.8e11-delta_E2_ERC15)
print("Bin3:", "T11",9.3e10-delta_E3_ET11,"Z09", 9.3e10-delta_E3_EZ09,"IP13", 9.3e10-delta_E3_ERC15)

####################
# find outlier in Bin 1, 2, 3

Bin2_out = SSort.selection_generic(corr_dist2,10**E2, np.repeat(75.0,len(E2)),direction="high",axis='x')
Bin3_out = SSort.selection_generic(corr_dist3,10**E3, np.repeat(45.0,len(E3)),direction="high",axis='x')

print('Bin2_out',Bin2_out)
print('Bin3_out',Bin3_out)

#####################
#SPlot.ShowcaseIndi.show_name(DD1,E,name_b,A=axs1)


def mass_distance_2plots():
    fig = plt.figure(figsize=(6.4, 9.6))
    gs = gridspec.GridSpec(2,1) 

    axs1 = plt.subplot(gs[1])
    axs0 = plt.subplot(gs[0])
    
    max_mass = 5e13

    #axs1.hlines(2e12,0,110,linestyle="solid")
    axs1.hlines(Bin1_lowerlimit,75,110,linestyle="solid")
    axs1.hlines(Bin2_lowerlimit,45,75,linestyle="solid")
    axs1.hlines(Bin3_lowerlimit,0,45,linestyle="solid")

    axs1.vlines(110,Bin1_lowerlimit,max_mass,linestyle="solid")
    axs1.vlines(75,Bin2_lowerlimit,max_mass,linestyle="solid")
    axs1.vlines(45,Bin3_lowerlimit,max_mass,linestyle="solid")

    #x_edge1,y_edge1= [75,75,110,110], [5e11,max_mass,max_mass,5e11]
    #x_edge2,y_edge2= [45,45,75,75], [2e11,max_mass,max_mass,2e11]
    #x_edge3,y_edge3= [0,0,45,45], [1e11,max_mass,max_mass,1e11]
    
    x_edge1,y_edge1= [75,75,110,110], [Bin1_lowerlimit,max_mass,max_mass,Bin1_lowerlimit]
    x_edge2,y_edge2= [45,45,75,75], [Bin2_lowerlimit,max_mass,max_mass,Bin2_lowerlimit]
    x_edge3,y_edge3= [0,0,45,45], [Bin3_lowerlimit,max_mass,max_mass,Bin3_lowerlimit]

    axs0.plot(dc,initial_mass,'x',color='grey', alpha=0.8, label='All')

    axs0.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)

    axs0.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)
    #axs0.plot(dc,E,'x',color='green', alpha=0.4, label='New Distance')
    axs0.plot(dc_b,initial_mass_b,'x',color='green', alpha=0.8, label='Broad Selection')

    axs1.plot(QQQ1["data"][1],QQQ1["data"][0],'rX',label='Bin1 samples', ms=16,alpha=0.7)
    axs1.plot(QQQ2["data"][1],QQQ2["data"][0],'bX',label='Bin2 samples', ms=16,alpha=0.7)
    axs1.plot(QQQ3["data"][1],QQQ3["data"][0],'kX',label='Bin3 samples', ms=16,alpha=0.7)  
    
    axs1.plot(DD1,E,'x',color='#a51a74', alpha=0.8, label='New Distance')
    #topcaaxs1.plot(DD1,E_T11,'x',color='g', alpha=0.8, label='T11')
    #axs1.plot(DD1,E_Z09,'x',color='b', alpha=0.8, label='Z09')
    #axs1.plot(DD1,E_RC15,'x',color='k', alpha=0.8, label='RC15')

    #axs1.plot(corr_dist1,corr_E1,'ro',label='Bin1 samples', ms=12,alpha=0.8)
    #axs1.plot(corr_dist2,corr_E2,'bo',label='Bin2 samples', ms=12,alpha=0.8)
    #axs1.plot(corr_dist3,corr_E3,'ko',label='Bin3 samples', ms=12,alpha=0.8)
    
    #axs1.plot(d1,10**E1_d,'ro',label='Bin1 samples', ms=12,alpha=0.8)
    #axs1.plot(d2,10**E2_d,'bo',label='Bin2 samples', ms=12,alpha=0.8)
    #axs1.plot(d3,10**E3_d,'ko',label='Bin3 samples', ms=12,alpha=0.8)    
    
    axs1.plot(corr_dist1,10**E1,'ro',label='Bin1 samples', ms=12,alpha=1.0)
    axs1.plot(corr_dist2,10**E2,'bo',label='Bin2 samples', ms=12,alpha=1.0)
    axs1.plot(corr_dist3,10**E3,'ko',label='Bin3 samples', ms=12,alpha=1.0)
    
    
    QQQ1
    
    #plot also the outliers, overlay on the old points
    axs1.plot(Bin2_out['bag_x'],Bin2_out['bag_y'],'ro', ms=12,alpha=1.0)
    axs1.plot(Bin3_out['bag_x'],Bin3_out['bag_y'],'bo', ms=12,alpha=1.0)
    
  
    
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

    axs0.set_ylabel(r"$  M_*/\rm M_\odot(IP13)$",fontsize=16)
    axs0.set_xlabel(r"$ \rm Distance/Mpc$",fontsize=16)
    axs0.set_yscale( 'log' )

    axs1.set_ylabel(r"$ M_*/ \rm M_\odot(IP13)$",fontsize=16)
    axs1.set_xlabel(r"$ \rm Distance/Mpc$",fontsize=16)
    axs1.set_yscale( 'log' )

    axs1.set_ylim(top =2e12 , bottom = 1e10)
    axs1.set_xlim(left=0, right = 120)

    axs0.legend(numpoints=1,scatterpoints=1,loc=4,fontsize=10)
    axs1.legend(numpoints=1,scatterpoints=1,loc=4, fontsize=10)
    plt.tight_layout()

    plt.show()

def mass_distance_bin_compare():
    
    fig = plt.figure(figsize=(6.4, 4.8))

    gs = gridspec.GridSpec(1,1) 

    axs0 = plt.subplot(gs[0])
    max_mass = 5e13

    #axs1.hlines(2e12,0,110,linestyle="solid")
    axs0.hlines(5e11,75,110,linestyle="solid")
    axs0.hlines(2e11,45,75,linestyle="solid")
    axs0.hlines(1e11,0,45,linestyle="solid")

    #q1 = np.log10(5e11)-np.average(E1_T11-key_delta_Bin1)-0.5*np.std(E1_T11-key_delta_Bin1)
    #q2 = np.log10(2e11)-np.average(E2_T11-key_delta_Bin2)-0.5*np.std(E2_T11-key_delta_Bin2)
    #q3 = np.log10(1e11)-np.average(E3_T11-key_delta_Bin3)-0.5*np.std(E3_T11-key_delta_Bin3)
    
    #q1 = np.log10(5e11)-np.average(E1_Z09-key_delta_Bin1)-0*np.std(E1_Z09-key_delta_Bin1)
    #q2 = np.log10(2e11)-np.average(E2_Z09-key_delta_Bin2)-0*np.std(E2_Z09-key_delta_Bin2)
    #q3 = np.log10(1e11)-np.average(E3_Z09-key_delta_Bin3)-0*np.std(E3_Z09-key_delta_Bin3)
    
    #q1 = np.log10(5e11)-0.45*np.average(E1_RC15-key_delta_Bin1)-0*np.std(E1_Z09-key_delta_Bin1)
    #q2 = np.log10(2e11)-0.45*np.average(E2_RC15-key_delta_Bin2)-0*np.std(E2_Z09-key_delta_Bin2)
    #q3 = np.log10(1e11)-0.45*np.average(E3_RC15-key_delta_Bin3)-0*np.std(E3_Z09-key_delta_Bin3)
    
    #q1 = np.log10(5e11)-np.average(E1_T11-key_delta_Bin1)
    #q2 = np.log10(2e11)-np.average(E2_T11-key_delta_Bin2)
    #q3 = np.log10(1e11)-np.average(E3_T11-key_delta_Bin3)
    
    #q1 = np.log10(5e11)-max_E1_ERC15
    #q2 = np.log10(2e11)-max_E2_ERC15
    #q3 = np.log10(1e11)-max_E3_ERC15
    
    #q1 = np.log10(5e11)-max_E1_EZ09
    #q2 = np.log10(2e11)-max_E2_EZ09
    #q3 = np.log10(1e11)-max_E3_EZ09
    
    Bin1_border =5e11
    Bin2_border =2e11
    Bin3_border =1e11
    
    q1 = np.log10(Bin1_border)-max_E1_ET11
    q2 = np.log10(Bin2_border)-max_E2_ET11
    q3 = np.log10(Bin3_border)-max_E3_ET11   
    
    print("q1 , q2, q3", q1 , q2, q3)
    axs0.hlines(10**q1,75,110,linestyle="solid")
    axs0.hlines(10**q2,45,75,linestyle="solid")
    axs0.hlines(10**q3,0,45,linestyle="solid")

    axs0.vlines(110,5e11,max_mass,linestyle="solid")
    axs0.vlines(75,2e11,max_mass,linestyle="solid")
    axs0.vlines(45,1e11,max_mass,linestyle="solid")

    x_edge1,y_edge1= [75,75,110,110], [5e11,max_mass,max_mass,5e11]
    x_edge2,y_edge2= [45,45,75,75], [2e11,max_mass,max_mass,2e11]
    x_edge3,y_edge3= [0,0,45,45], [1e11,max_mass,max_mass,1e11]

    axs0.plot(DD1,E,'x',color='#a51a74', alpha=0.8, label='New Distance IP13')
    axs0.plot(DD1,E_T11,'x',color='g', alpha=0.8, label='T11')
    #axs1.plot(DD1,E_Z09,'x',color='b', alpha=0.8, label='Z09')
    #axs1.plot(DD1,E_RC15,'x',color='k', alpha=0.8, label='RC15')

    axs0.plot(corr_dist1,E1,'ro',label='Bin1 samples', ms=12,alpha=0.8)
    axs0.plot(corr_dist2,E2,'bo',label='Bin2 samples', ms=12,alpha=0.8)
    axs0.plot(corr_dist3,E3,'ko',label='Bin3 samples', ms=12,alpha=0.8)
        
    axs0.plot(corr_dist1,10**E1_T11,'rs',label='Bin1 samples', ms=12,alpha=0.8)
    axs0.plot(corr_dist2,10**E2_T11,'bs',label='Bin2 samples', ms=12,alpha=0.8)
    axs0.plot(corr_dist3,10**E3_T11,'ks',label='Bin3 samples', ms=12,alpha=0.8)

    #axs0.plot(corr_dist1,10**E1_Z09,'rs',label='Bin1 samples', ms=12,alpha=0.8)
    #axs0.plot(corr_dist2,10**E2_Z09,'bs',label='Bin2 samples', ms=12,alpha=0.8)
    #axs0.plot(corr_dist3,10**E3_Z09,'ks',label='Bin3 samples', ms=12,alpha=0.8)
    
    #axs0.plot(corr_dist1,10**E1_RC15,'rs',label='Bin1 samples', ms=12,alpha=0.8)
    #axs0.plot(corr_dist2,10**E2_RC15,'bs',label='Bin2 samples', ms=12,alpha=0.8)
    #axs0.plot(corr_dist3,10**E3_RC15,'ks',label='Bin3 samples', ms=12,alpha=0.8)
    
    axs0.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)
    axs0.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)

    axs0.fill(x_edge1,y_edge1, alpha=0.3, color='green')
    axs0.fill(x_edge2,y_edge2, alpha=0.2, color='green')
    axs0.fill(x_edge3,y_edge3, alpha=0.1, color='green')

    axs0.set_ylim(top =5e12 , bottom = 1e9)
    axs0.set_xlim(left=0, right = 120)

    axs0.set_ylabel(r"$ M_*/ \rm M_\odot(IP13)$",fontsize=16)
    axs0.set_xlabel(r"$ \rm Distance/Mpc$",fontsize=16)
    axs0.set_yscale( 'log' )

    axs0.set_ylim(top =5e13 , bottom = 1e10)
    axs0.set_xlim(left=0, right = 120)

    axs0.legend(numpoints=1,scatterpoints=1,loc=4,fontsize=10)
    plt.tight_layout()

    plt.show()
    
mass_distance_2plots()

mass_distance_bin_compare()