#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 12:14:11 2020

@author: dexter
"""
from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

## style setting
import matplotlib.style
import matplotlib as mpl
plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0
##Define the selection volume
D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180))-np.cos(np.pi/2))

#V1,V2,V3=voll[2],voll[1],voll[0]

V1_V = voll[2]-voll[1]
V2_V = voll[1]-voll[0]
V3_V = voll[0]

V1,V2,V3 = V1_V,V2_V,V3_V

#voll[0],voll[1],voll[2]

print("volume", V1, V2, V3)
print("1/V", 1/V1, 1/V2, 1/V3)
print("2/V", 2/V1, 2/V2, 2/V3)
print("5/V", 5/V1, 5/V2, 5/V3)

M_sun = 4.53

#Read data#############################################
D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_3.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_3.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_3.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_3.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_3.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_3.txt",
    dtype = 'str')

D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt")

#RA,DEC 
RA_1,DEC_1 = D0_Bin1_table[:,1], D0_Bin1_table[:,2]
RA_2,DEC_2 = D0_Bin2_table[:,1], D0_Bin2_table[:,2]
RA_3,DEC_3 = D0_Bin3_table[:,1], D0_Bin3_table[:,2]

mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

corr_D1 = D0_Bin1_table[:,12]
corr_D2 = D0_Bin2_table[:,12]
corr_D3 = D0_Bin3_table[:,12]

# calculate the ellipticity of the Sersic dist2D fit
b_a_1 = D0_Bin1_table[:,34]
b_a_2 = D0_Bin2_table[:,34]
b_a_3 = D0_Bin3_table[:,34]

#extended disk ellipicity
elle1 = D0_Bin1_table[:,-1] 
elle2 = D0_Bin2_table[:,-1] 
elle3 = D0_Bin3_table[:,-1] 

seeing1 = D0_Bin1_table[:,15]
seeing2 = D0_Bin2_table[:,15]
seeing3 =  D0_Bin3_table[:,15]

#Get the name of the galaxies
name_D1 = D0_Bin1_table_n[:,0]
name_D2 = D0_Bin2_table_n[:,0]
name_D3 = D0_Bin3_table_n[:,0]

#Get the morphology of the galaxies
morph1 = D0_Bin1_table_n[:,17]
morph2 = D0_Bin2_table_n[:,17]
morph3 = D0_Bin3_table_n[:,17]

morph1_new = D0_Bin1_table_n[:,-3]
morph2_new = D0_Bin2_table_n[:,-3]
morph3_new = D0_Bin3_table_n[:,-3]

Rmax1 = D0_Bin1_table_n[:,-1]
Rmax2 = D0_Bin1_table_n[:,-1]
Rmax3 = D0_Bin1_table_n[:,-1]

corr_mass1 = D0_Bin1_table[:,13]
corr_mass2 = D0_Bin2_table[:,13]
corr_mass3 = D0_Bin3_table[:,13]

#####################################################
# red ReRmax
geom_file = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat")

geom_file1 = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_Bin1.dat")
geom_file2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_Bin2.dat")
geom_file3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_Bin3.dat")

Rmax = geom_file[:,2]
Rmax1, Rmax2, Rmax3 = geom_file1[:,2], geom_file2[:,2], geom_file3[:,2]

Rmax1_spc = SSort.morph_str_selection(Rmax1,morph1_new)
Rmax2_spc = SSort.morph_str_selection(Rmax2,morph2_new)
Rmax3_spc = SSort.morph_str_selection(Rmax3,morph3_new)

############### reading result files###############
master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"

name1 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
name2 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])


mu_e_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 0) 
mu_e_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 0) 
mu_e_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 0) 

Re_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

n1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 2)
n2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 2)
n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 2)

ars = (4.84814e-6)*1e3 #1 arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale1 = D1* ars
scale2 = D2* ars
scale3 = D3* ars

scale1_lerr, scale1_uerr = (D1-D1_lerr)*ars, (D1+D1_uerr)*ars
scale2_lerr, scale2_uerr =  (D2-D2_lerr)*ars, (D2+D2_uerr)*ars
scale3_lerr, scale3_uerr =  (D3-D3_lerr)*ars, (D3+D3_uerr)*ars

Re_1_kpc = Re_1* scale1
Re_2_kpc = Re_2* scale2
Re_3_kpc = Re_3* scale3

Re_1_kpc_lerr, Re_1_kpc_uerr = abs(Re_1_kpc - Re_1* scale1_lerr) , abs(Re_1* scale1_uerr - Re_1_kpc)
Re_2_kpc_lerr, Re_2_kpc_uerr = abs(Re_2* scale2_lerr - Re_2_kpc) , abs(Re_2* scale2_uerr - Re_2_kpc)
Re_3_kpc_lerr, Re_3_kpc_uerr = abs(Re_3* scale3_lerr - Re_3_kpc), abs(Re_3* scale3_uerr - Re_3_kpc)

Re_1_kpc_err =[Re_1_kpc_lerr, Re_1_kpc_uerr]
Re_2_kpc_err =[Re_2_kpc_lerr, Re_2_kpc_uerr]
Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]

############### reading result files###############
master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"

name1 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
name2 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])


Sersic_n1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 2) #get n
Sersic_n2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 2) #get n
Sersic_n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 2) #get n


mu_e_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 0) 
mu_e_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 0) 
mu_e_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 0) 

Re_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

#Get Re major axis
Re_1_major = SRead.grab_parameter("F_Gal_bundle_major_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2_major = SRead.grab_parameter("F_Gal_bundle_major_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3_major = SRead.grab_parameter("F_Gal_bundle_major_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

n1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 2)
n2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 2)
n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 2)

ars = (4.84814e-6)*1e3 #1 arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale1 = D1* ars
scale2 = D2* ars
scale3 = D3* ars

scale1_lerr, scale1_uerr = (D1-D1_lerr)*ars, (D1+D1_uerr)*ars
scale2_lerr, scale2_uerr =  (D2-D2_lerr)*ars, (D2+D2_uerr)*ars
scale3_lerr, scale3_uerr =  (D3-D3_lerr)*ars, (D3+D3_uerr)*ars

Re_1_kpc = Re_1* scale1
Re_2_kpc = Re_2* scale2
Re_3_kpc = Re_3* scale3

Re_1_kpc_major = Re_1_major* scale1
Re_2_kpc_major = Re_2_major* scale2
Re_3_kpc_major = Re_3_major* scale3

Re_1_kpc_lerr, Re_1_kpc_uerr = abs(Re_1_kpc - Re_1* scale1_lerr) , abs(Re_1* scale1_uerr - Re_1_kpc)
Re_2_kpc_lerr, Re_2_kpc_uerr = abs(Re_2* scale2_lerr - Re_2_kpc) , abs(Re_2* scale2_uerr - Re_2_kpc)
Re_3_kpc_lerr, Re_3_kpc_uerr = abs(Re_3* scale3_lerr - Re_3_kpc), abs(Re_3* scale3_uerr - Re_3_kpc)

Re_1_kpc_err =[Re_1_kpc_lerr, Re_1_kpc_uerr]
Re_2_kpc_err =[Re_2_kpc_lerr, Re_2_kpc_uerr]
Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]


############# calculating spheroid mass ########
ML_select1_IP13 = SPlot.MLRelationIband(mag_g1,mag_i1).Into13_MassRatio
ML_select1_R15BC = SPlot.MLRelationIband(mag_g1,mag_i1).Roediger15BC03_MassRatio
ML_select1_Z09 = SPlot.MLRelationIband(mag_g1,mag_i1).Zibetti09_MassRatio
ML_select1_T11 = SPlot.MLRelationIband(mag_g1,mag_i1).Taylor11_MassRatio

M1 = SPlot.MassCalculation(sph_mag1, D1, 4.53,mag_g1,mag_i1)

E1_IP13 = M1.cal_Mass(ML_select1_IP13)
E1_R15BC = M1.cal_Mass(ML_select1_R15BC)
E1_Z09 = M1.cal_Mass(ML_select1_Z09)
E1_T11 = M1.cal_Mass(ML_select1_T11)


ML_select2_IP13 = SPlot.MLRelationIband(mag_g2,mag_i2).Into13_MassRatio
ML_select2_R15BC = SPlot.MLRelationIband(mag_g2,mag_i2).Roediger15BC03_MassRatio
ML_select2_Z09 = SPlot.MLRelationIband(mag_g2,mag_i2).Zibetti09_MassRatio
ML_select2_T11 = SPlot.MLRelationIband(mag_g2,mag_i2).Taylor11_MassRatio

M2 = SPlot.MassCalculation(sph_mag2, D2, 4.53,mag_g2,mag_i2)

E2_IP13 = M2.cal_Mass(ML_select2_IP13)
E2_R15BC = M2.cal_Mass(ML_select2_R15BC)
E2_Z09 = M2.cal_Mass(ML_select2_Z09)
E2_T11 = M2.cal_Mass(ML_select2_T11)


ML_select3_IP13 = SPlot.MLRelationIband(mag_g3,mag_i3).Into13_MassRatio
ML_select3_R15BC = SPlot.MLRelationIband(mag_g3,mag_i3).Roediger15BC03_MassRatio
ML_select3_Z09 = SPlot.MLRelationIband(mag_g3,mag_i3).Zibetti09_MassRatio
ML_select3_T11 = SPlot.MLRelationIband(mag_g3,mag_i3).Taylor11_MassRatio

M3 = SPlot.MassCalculation(sph_mag3, D3, 4.53,mag_g3,mag_i3)

E3_IP13 = M3.cal_Mass(ML_select3_IP13)
E3_R15BC = M3.cal_Mass(ML_select3_R15BC)
E3_Z09 = M3.cal_Mass(ML_select3_Z09)
E3_T11 = M3.cal_Mass(ML_select3_T11)

################################
#Calculate mass with K-correction

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
mag_g1_corr, mag_i1_corr = mag_g1-g1_kcorr-g1_EXT, mag_i1-i1_kcorr-i1_EXT
mag_g2_corr, mag_i2_corr = mag_g2-g2_kcorr-g2_EXT, mag_i2-i2_kcorr-i2_EXT
mag_g3_corr, mag_i3_corr = mag_g3-g3_kcorr-g3_EXT, mag_i3-i3_kcorr-i3_EXT

Mag_i1_kcorr_cDis = mag_i1_corr-25-5*np.log10(corr_D1) 
Mag_i2_kcorr_cDis = mag_i2_corr-25-5*np.log10(corr_D2) 
Mag_i3_kcorr_cDis = mag_i3_corr-25-5*np.log10(corr_D3) 

Mag_i1_kcorr = mag_i1_corr-25-5*np.log10(D1) 
Mag_i2_kcorr = mag_i2_corr-25-5*np.log10(D2) 
Mag_i3_kcorr = mag_i3_corr-25-5*np.log10(D3) 


###################
#calculate the stellar mass after extinction and Kcorrection with corr_Dist
ML_select1_IP13_corr = SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Into13_MassRatio
M1_corr = SPlot.MassCalculation(mag_i1_corr, corr_D1, 4.53,mag_g1_corr,mag_i1_corr)

ML_select2_IP13_corr = SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Into13_MassRatio
M2_corr = SPlot.MassCalculation(mag_i2_corr, corr_D2, 4.53,mag_g2_corr,mag_i2_corr)

ML_select3_IP13_corr = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Into13_MassRatio
M3_corr = SPlot.MassCalculation(mag_i3_corr, corr_D3, 4.53,mag_g3_corr,mag_i3_corr)

E1_gal_cDis = M1_corr.cal_Mass(ML_select1_IP13_corr)
E2_gal_cDis = M2_corr.cal_Mass(ML_select2_IP13_corr)
E3_gal_cDis = M3_corr.cal_Mass(ML_select3_IP13_corr)


def list_gi_corr():
    for i in range(len(mag_g3)):
        print(name3[i], mag_i3_corr[i],Mag_i3_kcorr_cDis[i],Mag_i3_kcorr[i], E3_gal_cDis[i],corr_mass3[i])
    print('--------------------------------------')

#list_gi_corr()


#####################
# calculate the M/L ratio
ML_select1_IP13_K = SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Into13_MassRatio
ML_select1_R15BC_K = SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Roediger15BC03_MassRatio
ML_select1_Z09_K = SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Zibetti09_MassRatio
ML_select1_T11_K = SPlot.MLRelationIband(mag_g1_corr,mag_i1_corr).Taylor11_MassRatio

ML_select2_IP13_K = SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Into13_MassRatio
ML_select2_R15BC_K = SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Roediger15BC03_MassRatio
ML_select2_Z09_K = SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Zibetti09_MassRatio
ML_select2_T11_K = SPlot.MLRelationIband(mag_g2_corr,mag_i2_corr).Taylor11_MassRatio

ML_select3_IP13_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Into13_MassRatio
ML_select3_R15BC_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Roediger15BC03_MassRatio
ML_select3_Z09_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Zibetti09_MassRatio
ML_select3_T11_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Taylor11_MassRatio

############# calculating  stellar mass

##########################################
#K correction and extinction for sph mag
sph_mag1 = sph_mag1 - i1_EXT - i1_kcorr
sph_mag2 = sph_mag2 - i2_EXT - i2_kcorr
sph_mag3 = sph_mag3 - i3_EXT - i3_kcorr

sph_abs_mag1 = SPlot.MassCalculation(mag_g1_corr, mag_i1_corr, sph_mag1, D1, M_sun).cal_abs_mag(sph_mag1, D1)
sph_abs_mag2 = SPlot.MassCalculation(mag_g2_corr, mag_i2_corr, sph_mag2, D2, M_sun).cal_abs_mag(sph_mag2, D2)
sph_abs_mag3 = SPlot.MassCalculation(mag_g3_corr, mag_i3_corr, sph_mag3, D3, M_sun).cal_abs_mag(sph_mag3, D3)


sph_abs_mag1 = SPlot.MassCalculation(mag_g1_corr, mag_i1_corr,sph_mag1, D1, M_sun).cal_abs_mag(sph_mag1, D1)
sph_abs_mag2 = SPlot.MassCalculation(mag_g2_corr, mag_i2_corr,sph_mag2, D2, M_sun).cal_abs_mag(sph_mag2, D2)
sph_abs_mag3 = SPlot.MassCalculation(mag_g3_corr, mag_i3_corr,sph_mag3, D3, M_sun).cal_abs_mag(sph_mag3, D3)


#######################
M1_K = SPlot.MassCalculation(sph_mag1, D1, 4.53,mag_g1_corr,mag_i1_corr)

E1_IP13_K = M1_K.cal_Mass(ML_select1_IP13_K)
E1_R15BC_K = M1_K.cal_Mass(ML_select1_R15BC_K)
E1_Z09_K = M1_K.cal_Mass(ML_select1_Z09_K)
E1_T11_K = M1_K.cal_Mass(ML_select1_T11_K)


M2_K = SPlot.MassCalculation(sph_mag2, D2, 4.53,mag_g2_corr,mag_i2_corr)

E2_IP13_K = M2_K.cal_Mass(ML_select2_IP13_K)
E2_R15BC_K = M2_K.cal_Mass(ML_select2_R15BC_K)
E2_Z09_K = M2_K.cal_Mass(ML_select2_Z09_K)
E2_T11_K = M2_K.cal_Mass(ML_select2_T11_K)


M3_K = SPlot.MassCalculation(sph_mag3, D3, 4.53, mag_g3_corr,mag_i3_corr)

E3_IP13_K = M3_K.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K = M3_K.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K = M3_K.cal_Mass(ML_select3_Z09_K)
E3_T11_K = M3_K.cal_Mass(ML_select3_T11_K)


#####################
E_R15BC_K = np.concatenate((E1_R15BC_K,E2_R15BC_K, E3_R15BC_K))
E_T11_K = np.concatenate((E1_T11_K,E2_T11_K, E3_T11_K))


############# calculating spheroid absoulte magnitude
Abs_sph_mag1 = M1.cal_abs_mag(sph_mag1,D1)
Abs_sph_mag2 = M2.cal_abs_mag(sph_mag2,D2)
Abs_sph_mag3 = M3.cal_abs_mag(sph_mag3,D3)

################################
#calculate the mass error

mag_e = 0.28 #magnitude error
MLR_e_Z09 = 10**0.125
MLR_e_T11 = 10**0.1
MLR_e_IP13 = 10**0.14
MLR_e_R15BC = 10**0.13


def cal_mass_error(mag_e, D, D_err, MLR ,MLR_e):
    return np.sqrt(((mag_e/2.5)**2)+
                   ((2*D_err/(D*np.log(10)))**2)+
                   ((MLR_e/(MLR*np.log(10)))**2))

mass_Z09_uerr1 = E1_Z09_K* cal_mass_error(mag_e, D1, D1_uerr, ML_select1_Z09_K ,MLR_e_Z09)
mass_T11_uerr1 = E1_T11_K * cal_mass_error(mag_e, D1, D1_uerr, ML_select1_T11_K ,MLR_e_T11)
mass_IP13_uerr1 = E1_IP13_K* cal_mass_error(mag_e, D1, D1_uerr, ML_select1_IP13_K ,MLR_e_IP13)
mass_R15BC_uerr1 = E1_R15BC_K* cal_mass_error(mag_e, D1, D1_uerr, ML_select1_R15BC_K ,MLR_e_R15BC)

mass_Z09_uerr2 = E2_Z09_K* cal_mass_error(mag_e, D2, D2_uerr, ML_select2_Z09_K ,MLR_e_Z09)
mass_T11_uerr2 = E2_T11_K * cal_mass_error(mag_e, D2, D2_uerr, ML_select2_T11_K ,MLR_e_T11)
mass_IP13_uerr2 = E2_IP13_K* cal_mass_error(mag_e, D2, D2_uerr, ML_select2_IP13_K ,MLR_e_IP13)
mass_R15BC_uerr2 = E2_R15BC_K* cal_mass_error(mag_e, D2, D2_uerr, ML_select2_R15BC_K ,MLR_e_R15BC)

mass_Z09_uerr3 = E3_Z09_K* cal_mass_error(mag_e, D3, D3_uerr, ML_select3_Z09_K ,MLR_e_Z09)
mass_T11_uerr3 = E3_T11_K *cal_mass_error(mag_e, D3, D3_uerr, ML_select3_T11_K ,MLR_e_T11)
mass_IP13_uerr3 = E3_IP13_K*cal_mass_error(mag_e, D3, D3_uerr, ML_select3_IP13_K ,MLR_e_IP13)
mass_R15BC_uerr3 = E3_R15BC_K*cal_mass_error(mag_e, D3, D3_uerr, ML_select3_R15BC_K ,MLR_e_R15BC)

mass_Z09_lerr1 = E1_Z09_K*cal_mass_error(mag_e, D1, D1_lerr, ML_select1_Z09_K ,MLR_e_Z09)
mass_T11_lerr1 = E1_T11_K *cal_mass_error(mag_e, D1, D1_lerr, ML_select1_T11_K ,MLR_e_T11)
mass_IP13_lerr1 = E1_IP13_K*cal_mass_error(mag_e, D1, D1_lerr, ML_select1_IP13_K ,MLR_e_IP13)
mass_R15BC_lerr1 = E1_R15BC_K*cal_mass_error(mag_e, D1, D1_lerr, ML_select1_R15BC_K ,MLR_e_R15BC)

mass_Z09_lerr2 = E2_Z09_K*cal_mass_error(mag_e, D2, D2_lerr, ML_select2_Z09_K ,MLR_e_Z09)
mass_T11_lerr2 = E2_T11_K *cal_mass_error(mag_e, D2, D2_lerr, ML_select2_T11_K ,MLR_e_T11)
mass_IP13_lerr2 = E2_IP13_K*cal_mass_error(mag_e, D2, D2_lerr, ML_select2_IP13_K ,MLR_e_IP13)
mass_R15BC_lerr2 = E2_R15BC_K*cal_mass_error(mag_e, D2, D2_lerr, ML_select2_R15BC_K ,MLR_e_R15BC)

mass_Z09_lerr3 = E3_Z09_K*cal_mass_error(mag_e, D3, D3_lerr, ML_select3_Z09_K ,MLR_e_Z09)
mass_T11_lerr3 = E3_T11_K *cal_mass_error(mag_e, D3, D3_lerr, ML_select3_T11_K ,MLR_e_T11)
mass_IP13_lerr3 = E3_IP13_K*cal_mass_error(mag_e, D3, D3_lerr, ML_select3_IP13_K ,MLR_e_IP13)
mass_R15BC_lerr3 = E3_R15BC_K*cal_mass_error(mag_e, D3, D3_lerr, ML_select3_R15BC_K ,MLR_e_R15BC)

################################
#calculate the mass error

MLR1 = ML_select1_R15BC_K
MLR_e1 = 10**0.13

MLR2 = ML_select2_R15BC_K
MLR_e2 = 10**0.13

MLR3 = ML_select3_R15BC_K
MLR_e3 = 10**0.13

mag_e = 0.3 #magnitude error

mass_uerr1 = np.sqrt(((mag_e/2.5)**2)+((2*D1_uerr/(D1*np.log(10)))**2)+((MLR_e1/(MLR1*np.log(10)))**2))
mass_uerr2 = np.sqrt(((mag_e/2.5)**2)+((2*D2_uerr/(D2*np.log(10)))**2)+((MLR_e2/(MLR2*np.log(10)))**2))
mass_uerr3 = np.sqrt(((mag_e/2.5)**2)+((2*D3_uerr/(D3*np.log(10)))**2)+((MLR_e3/(MLR3*np.log(10)))**2))

mass_err1 = mass_uerr1
mass_err2 = mass_uerr2
mass_err3 = mass_uerr3

##########################################
# read NSA-Sloan catalog
nsa = SRead.read_table('/home/dexter/result/stat/completeness/nsa_sizemass.dat')

nsa_z = nsa[:,0]

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)
nsa_d = cosmo.comoving_distance(nsa_z).value

nsa_scale = nsa_d * ars

nsa_Re = nsa[:,1]*nsa_scale
nsa_mass_T11 = nsa[:,2]
nsa_mass_Z09 = nsa[:,3]
nsa_mass_R15BC = nsa[:,4]
nsa_mass_IP13 = nsa[:,5]


# read Benzanson catalog

Benzanson2015 =SRead.read_table("/home/dexter/result/Bezanson2015table.txt")
Zahid2015 = SRead.read_table("/home/dexter/result/Zahid2015table.txt")

Benzanson2015_Re = np.array(Benzanson2015[:,4])
Benzanson2015_mass = np.array(10**Benzanson2015[:,7])

Zahid2015_Re  = np.array(Zahid2015[:,2])
Zahid2015_mass  = np.array(10**Zahid2015[:,4])



#################################
#T11 for mass
#mass1,mass2,mass3 = E1_T11_K, E2_T11_K, E3_T11_K

#Z09 for mass
#mass1,mass2,mass3 = E1_Z09_K, E2_Z09_K, E3_Z09_K

#RC15 for mass
mass1, mass2, mass3 = E1_R15BC_K, E2_R15BC_K, E3_R15BC_K

#IP13 for mass
#mass1, mass2, mass3 = E1_IP13_K, E2_IP13_K, E3_IP13_K


##################################
# Dust Correction and morph seperation

# correct for galaxies g and i in bulge+Disk

def dust_correction_new_sizemass(mass1,mass2,mass3,mass_err1,mass_err2,mass_err3):
    # Calculate the absoulte magnitude and the stellar mass
    Abs_sph_mag1 = M1_K.cal_abs_mag(sph_mag1,D1)
    Abs_sph_mag2 = M2_K.cal_abs_mag(sph_mag2,D2)
    Abs_sph_mag3 = M3_K.cal_abs_mag(sph_mag3,D3)
    
    # calculate the dust corrected version of abs mag for ALL sample 
    # E and S0 does not require that but I calculate them nonethesless
    Abs_sph_mag1_dustCorr = M1_K.dust_correction_Driver08(Abs_sph_mag1,elle1)
    Abs_sph_mag2_dustCorr = M2_K.dust_correction_Driver08(Abs_sph_mag2,elle2)
    Abs_sph_mag3_dustCorr = M3_K.dust_correction_Driver08(Abs_sph_mag3,elle3)
    
    # calculate the dust corrected version of stellar mass for ALL sample 
    mass1_K_dustCorr = MLR1*(10**((4.53-Abs_sph_mag1_dustCorr)/2.5))
    mass2_K_dustCorr = MLR2*(10**((4.53-Abs_sph_mag2_dustCorr)/2.5))
    mass3_K_dustCorr = MLR3*(10**((4.53-Abs_sph_mag3_dustCorr)/2.5))
    
    #Seperate the sample into E, S0, S
    index1 = np.arange(len(Abs_sph_mag1))
    index2 = np.arange(len(Abs_sph_mag2))
    index3 = np.arange(len(Abs_sph_mag3))
    
    morph_dict1 = SSort.morph_str_selection(index1, morph1_new)
    morph_dict2 = SSort.morph_str_selection(index2, morph2_new)
    morph_dict3 = SSort.morph_str_selection(index3, morph3_new)
    
    #cherry pick the distance in each bin by morphology
    D1_E = SSort.cherry_pick(morph_dict1['E'], D1)
    D1_S0 = SSort.cherry_pick(morph_dict1['S0'], D1)
    D1_S = SSort.cherry_pick(morph_dict1['S'], D1)
    
    D2_E = SSort.cherry_pick(morph_dict2['E'], D2)
    D2_S0 = SSort.cherry_pick(morph_dict2['S0'], D2)
    D2_S = SSort.cherry_pick(morph_dict2['S'], D2)
    
    D3_E = SSort.cherry_pick(morph_dict3['E'], D3)
    D3_S0 = SSort.cherry_pick(morph_dict3['S0'], D3)
    D3_S = SSort.cherry_pick(morph_dict3['S'], D3)

    # Stich the distance together
    D1_reshp = np.concatenate((D1_E, D1_S0, D1_S))
    D2_reshp = np.concatenate((D2_E, D2_S0, D2_S))
    D3_reshp = np.concatenate((D3_E, D3_S0, D3_S))
    
    #cherry pick the sph mag in each bin by morphology
    mass1_E = SSort.cherry_pick(morph_dict1['E'], mass1)
    mass1_S0 = SSort.cherry_pick(morph_dict1['S0'], mass1)
    mass1_dustCorr_S = SSort.cherry_pick(morph_dict1['S'], mass1_K_dustCorr)
    
    mass2_E = SSort.cherry_pick(morph_dict2['E'], mass2)
    mass2_S0 = SSort.cherry_pick(morph_dict2['S0'], mass2)
    mass2_dustCorr_S = SSort.cherry_pick(morph_dict2['S'], mass2_K_dustCorr)

    mass3_E = SSort.cherry_pick(morph_dict3['E'], mass3)
    mass3_S0 = SSort.cherry_pick(morph_dict3['S0'], mass3)
    mass3_dustCorr_S = SSort.cherry_pick(morph_dict3['S'], mass3_K_dustCorr)

    mass_err1_E = SSort.cherry_pick(morph_dict1['E'], mass_err1)
    mass_err1_S0 = SSort.cherry_pick(morph_dict1['S0'], mass_err1)
    mass_err1_S = SSort.cherry_pick(morph_dict1['S'], mass_err1)
    
    mass_err2_E = SSort.cherry_pick(morph_dict2['E'], mass_err2)
    mass_err2_S0 = SSort.cherry_pick(morph_dict2['S0'], mass_err2)
    mass_err2_S = SSort.cherry_pick(morph_dict2['S'], mass_err2)

    mass_err3_E = SSort.cherry_pick(morph_dict3['E'], mass_err3)
    mass_err3_S0 = SSort.cherry_pick(morph_dict3['S0'], mass_err3)
    mass_err3_S = SSort.cherry_pick(morph_dict3['S'], mass_err3) 

    # Stich the mass together, E,S0 no dust correction, S with dust correction
    # stellar mass reshape
    mass1_dustCorr_reshp = np.concatenate((mass1_E, mass1_S0, mass1_dustCorr_S))
    mass2_dustCorr_reshp = np.concatenate((mass2_E, mass2_S0, mass2_dustCorr_S))
    mass3_dustCorr_reshp = np.concatenate((mass3_E, mass3_S0, mass3_dustCorr_S))

    mass_err1_reshp = np.concatenate((mass_err1_E, mass_err1_S0, mass_err1_S))
    mass_err2_reshp = np.concatenate((mass_err2_E, mass_err2_S0, mass_err2_S))
    mass_err3_reshp = np.concatenate((mass_err3_E, mass_err3_S0, mass_err3_S))

    #cherry pick the Re in equvi in each bin by morphology
    Re_1_kpc_E = SSort.cherry_pick(morph_dict1['E'], Re_1_kpc)
    Re_1_kpc_S0 = SSort.cherry_pick(morph_dict1['S0'], Re_1_kpc)
    Re_1_kpc_S = SSort.cherry_pick(morph_dict1['S'], Re_1_kpc)

    Re_2_kpc_E = SSort.cherry_pick(morph_dict2['E'], Re_2_kpc)
    Re_2_kpc_S0 = SSort.cherry_pick(morph_dict2['S0'], Re_2_kpc)
    Re_2_kpc_S = SSort.cherry_pick(morph_dict2['S'], Re_2_kpc)

    Re_3_kpc_E = SSort.cherry_pick(morph_dict3['E'], Re_3_kpc)
    Re_3_kpc_S0 = SSort.cherry_pick(morph_dict3['S0'], Re_3_kpc)
    Re_3_kpc_S = SSort.cherry_pick(morph_dict3['S'], Re_3_kpc)

    # Stich the Re together
    # Re reshape
    Re_1_kpc_reshp = np.concatenate((Re_1_kpc_E, Re_1_kpc_S0, Re_1_kpc_S))
    Re_2_kpc_reshp = np.concatenate((Re_2_kpc_E, Re_2_kpc_S0, Re_2_kpc_S))
    Re_3_kpc_reshp = np.concatenate((Re_3_kpc_E, Re_3_kpc_S0, Re_3_kpc_S))

    #cherry pick the Re in major in each bin by morphology
    Re_1_kpc_major_E = SSort.cherry_pick(morph_dict1['E'], Re_1_kpc_major)
    Re_1_kpc_major_S0 = SSort.cherry_pick(morph_dict1['S0'], Re_1_kpc_major)
    Re_1_kpc_major_S = SSort.cherry_pick(morph_dict1['S'], Re_1_kpc_major)

    Re_2_kpc_major_E = SSort.cherry_pick(morph_dict2['E'], Re_2_kpc_major)
    Re_2_kpc_major_S0 = SSort.cherry_pick(morph_dict2['S0'], Re_2_kpc_major)
    Re_2_kpc_major_S = SSort.cherry_pick(morph_dict2['S'], Re_2_kpc_major)

    Re_3_kpc_major_E = SSort.cherry_pick(morph_dict3['E'], Re_3_kpc_major)
    Re_3_kpc_major_S0 = SSort.cherry_pick(morph_dict3['S0'], Re_3_kpc_major)
    Re_3_kpc_major_S = SSort.cherry_pick(morph_dict3['S'], Re_3_kpc_major)
    
    # Stich the Re in major together
    # Re reshape in major
    Re_1_kpc_major_reshp = np.concatenate((Re_1_kpc_major_E, Re_1_kpc_major_S0, Re_1_kpc_major_S))
    Re_2_kpc_major_reshp = np.concatenate((Re_2_kpc_major_E, Re_2_kpc_major_S0, Re_2_kpc_major_S))
    Re_3_kpc_major_reshp = np.concatenate((Re_3_kpc_major_E, Re_3_kpc_major_S0, Re_3_kpc_major_S))
    
    Re_1_kpc_lerr_E,  Re_1_kpc_uerr_E = SSort.cherry_pick(morph_dict1['E'], Re_1_kpc_lerr), SSort.cherry_pick(morph_dict1['E'], Re_1_kpc_uerr)
    Re_1_kpc_lerr_S0,  Re_1_kpc_uerr_S0 = SSort.cherry_pick(morph_dict1['S0'], Re_1_kpc_lerr), SSort.cherry_pick(morph_dict1['S0'], Re_1_kpc_uerr)
    Re_1_kpc_lerr_S,  Re_1_kpc_uerr_S = SSort.cherry_pick(morph_dict1['S'], Re_1_kpc_lerr), SSort.cherry_pick(morph_dict1['S'], Re_1_kpc_uerr)
    
    Re_2_kpc_lerr_E,  Re_2_kpc_uerr_E = SSort.cherry_pick(morph_dict2['E'], Re_2_kpc_lerr), SSort.cherry_pick(morph_dict2['E'], Re_2_kpc_uerr)
    Re_2_kpc_lerr_S0,  Re_2_kpc_uerr_S0 = SSort.cherry_pick(morph_dict2['S0'], Re_2_kpc_lerr), SSort.cherry_pick(morph_dict2['S0'], Re_2_kpc_uerr)
    Re_2_kpc_lerr_S,  Re_2_kpc_uerr_S = SSort.cherry_pick(morph_dict2['S'], Re_2_kpc_lerr), SSort.cherry_pick(morph_dict2['S'], Re_2_kpc_uerr)
    
    Re_3_kpc_lerr_E,  Re_3_kpc_uerr_E = SSort.cherry_pick(morph_dict3['E'], Re_3_kpc_lerr), SSort.cherry_pick(morph_dict3['E'], Re_3_kpc_uerr)
    Re_3_kpc_lerr_S0,  Re_3_kpc_uerr_S0 = SSort.cherry_pick(morph_dict3['S0'], Re_3_kpc_lerr), SSort.cherry_pick(morph_dict3['S0'], Re_3_kpc_uerr)
    Re_3_kpc_lerr_S,  Re_3_kpc_uerr_S = SSort.cherry_pick(morph_dict3['S'], Re_3_kpc_lerr), SSort.cherry_pick(morph_dict3['S'], Re_3_kpc_uerr)


    Re_1_kpc_lerr_reshp = np.concatenate((Re_1_kpc_lerr_E, Re_1_kpc_lerr_S0, Re_1_kpc_lerr_S))
    Re_2_kpc_lerr_reshp = np.concatenate((Re_2_kpc_lerr_E, Re_2_kpc_lerr_S0, Re_2_kpc_lerr_S))
    Re_3_kpc_lerr_reshp = np.concatenate((Re_3_kpc_lerr_E, Re_3_kpc_lerr_S0, Re_3_kpc_lerr_S))

    Re_1_kpc_uerr_reshp = np.concatenate((Re_1_kpc_uerr_E, Re_1_kpc_uerr_S0, Re_1_kpc_uerr_S))
    Re_2_kpc_uerr_reshp = np.concatenate((Re_2_kpc_uerr_E, Re_2_kpc_uerr_S0, Re_2_kpc_uerr_S))
    Re_3_kpc_uerr_reshp = np.concatenate((Re_3_kpc_uerr_E, Re_3_kpc_uerr_S0, Re_3_kpc_uerr_S))
    
    print(len(Re_1_kpc_lerr),len(Re_1_kpc_uerr))
    print(len(Re_2_kpc_lerr),len(Re_2_kpc_uerr))
    print(len(Re_3_kpc_lerr),len(Re_3_kpc_uerr))

    Re_1_kpc_err =[Re_1_kpc_lerr, Re_1_kpc_uerr]
    Re_2_kpc_err =[Re_2_kpc_lerr, Re_2_kpc_uerr]
    Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]

    # The mass, and size of E galaxies
    mass_E = np.concatenate((mass1_E,mass2_E,mass3_E))

    #E_R15BC_K_E = np.concatenate((E1_R15BC_K_E,E2_R15BC_K_E,E3_R15BC_K_E))
    Re_kpc_E = np.concatenate((Re_1_kpc_E,Re_2_kpc_E,Re_3_kpc_E))
    Re_kpc_major_E = np.concatenate((Re_1_kpc_major_E,Re_2_kpc_major_E,Re_3_kpc_major_E))
    return mass_E,Re_kpc_E, Re_kpc_major_E
    
mass_E,Re_kpc_E,Re_kpc_major_E = dust_correction_new_sizemass(mass1,mass2,mass3,mass_err1,mass_err2,mass_err3)

################################
# The three ES are at 7 and 17 in Bin1 and at 4 in Bin2

mass_ES = np.array([mass1[17],mass2[4]])

E_R15BC_K_ES = np.array([E1_R15BC_K[17],E2_R15BC_K[4]])
Re_kpc_ES = np.array([Re_1_kpc[17],Re_2_kpc[4]])
Re_kpc_major_ES = np.array([Re_1_kpc_major[17],Re_2_kpc_major[4]])

ES_index = [33,93]

################################
# Hightlight outliers

# define outliers as g-i > 1.4

# find outliers mass
gi_out_1_m = SSort.selection_generic(mass1, mag_g1- mag_i1, np.repeat(1.5,len(mass1)),direction='high',axis='y')
gi_out_2_m = SSort.selection_generic(mass2, mag_g2- mag_i2, np.repeat(1.5,len(mass2)),direction='high',axis='y')
gi_out_3_m = SSort.selection_generic(mass3, mag_g3- mag_i3, np.repeat(1.5,len(mass3)),direction='high',axis='y')

# find outliers radius
gi_out_1_R = SSort.selection_generic(Re_1_kpc, mag_g1- mag_i1, np.repeat(1.5,len(mass1)),direction='high',axis='y')
gi_out_2_R = SSort.selection_generic(Re_2_kpc, mag_g2- mag_i2, np.repeat(1.5,len(mass2)),direction='high',axis='y')
gi_out_3_R = SSort.selection_generic(Re_3_kpc, mag_g3- mag_i3, np.repeat(1.5,len(mass3)),direction='high',axis='y')

gi_out_1_Rm = SSort.selection_generic(Re_1_kpc_major, mag_g1- mag_i1, np.repeat(1.5,len(mass1)),direction='high',axis='y')
gi_out_2_Rm = SSort.selection_generic(Re_2_kpc_major, mag_g2- mag_i2, np.repeat(1.5,len(mass2)),direction='high',axis='y')
gi_out_3_Rm = SSort.selection_generic(Re_3_kpc_major, mag_g3- mag_i3, np.repeat(1.5,len(mass3)),direction='high',axis='y')


mass1_out = gi_out_1_m['bag_x']
mass2_out = gi_out_2_m['bag_x']
mass3_out = gi_out_3_m['bag_x']
 
R1_out = gi_out_1_R['bag_x']
R2_out = gi_out_2_R['bag_x']
R3_out = gi_out_3_R['bag_x']

R1_out_m = gi_out_1_Rm['bag_x']
R2_out_m = gi_out_2_Rm['bag_x']
R3_out_m = gi_out_3_Rm['bag_x']

name1_out = name1[gi_out_1_R['index']]
name2_out = name2[gi_out_2_R['index']]
name3_out = name3[gi_out_3_R['index']]

morph3_out = morph3[gi_out_3_R['index']]

print('name_out',name1_out,name2_out,name3_out)
print('morph3_out',morph3_out)

#plotting function 
def plot_gi_outliers(A,scale='log',alpha=0.65):
    A.plot(mass1_out, R1_out,marker='*',color='k',label='outliers', 
              ms =15, alpha=1.0,linestyle="None")
    A.plot(mass2_out, R2_out,marker='*',color='k',label='', 
              ms =15, alpha=1.0,linestyle="None")    
    A.plot(mass3_out, R3_out,marker='*',color='k',label='', 
              ms =15, alpha=1.0,linestyle="None")
    
    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )
    
def plot_gi_outliers_major(A,scale='log',alpha=0.65):
    A.plot(mass1_out, R1_out_m,marker='*',color='k',label='outliers', 
              ms =15, alpha=1.0,linestyle="None")
    A.plot(mass2_out, R2_out_m,marker='*',color='k',label='', 
              ms =15, alpha=1.0,linestyle="None")    
    A.plot(mass3_out, R3_out_m,marker='*',color='k',label='', 
              ms =15, alpha=1.0,linestyle="None")
    
    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )
    
################################
# Define ploting in RDJ15
def plot_dexter_sample_Bin():
    
    
    SPlot.ShowcaseIndi.Mass_Re_plot(mass1, Re_1_kpc, yerr = Re_1_kpc_err, 
                                xerr = mass_err1*mass1,
                                colour = '#a5200b',
                                legend='Bin1',ms=8,alpha0 = 0.2, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(mass2, Re_2_kpc, yerr = Re_2_kpc_err, 
                                xerr = mass_err2*mass2,
                                colour = '#0b5786',
                                legend='Bin2',ms=8,alpha0 = 0.2, lw=3)
    SPlot.ShowcaseIndi.Mass_Re_plot(mass3, Re_3_kpc, yerr = Re_3_kpc_err, 
                                xerr = mass_err3*mass3,
                                colour='#2a3236',
                                legend='Bin3',ms=8,alpha0 = 0.2, lw=3)
    
def plot_dexter_sample_all():
    
    SPlot.ShowcaseIndi.Mass_Re_plot(mass1, Re_1_kpc, yerr = Re_1_kpc_err,
                                xerr = mass_err1*mass1,
                                colour = '#a5200b',
                                legend='spheroids',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(mass2, Re_2_kpc, yerr = Re_2_kpc_err,
                                xerr = mass_err2*mass2,
                                colour = '#a5200b',
                                legend='',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(mass3, Re_3_kpc, yerr = Re_3_kpc_err,
                                xerr = mass_err3*mass3,
                                colour='#a5200b',
                                legend='',ms=8,alpha0 = 0.4,lw=3 )
    



txsep1 = 40#10*5.5
txsep2 = 65#10*6.0
txsep3 = 80#10*6.5

x,y = 8.6e8,100

circle_size = 310

# The cutting board
def plot_Barro_cut_all(AX):
    
    Bcut1 = SPlot.SelectionCut(mass1, D1).Barro13_cut()
    Bcut2 = SPlot.SelectionCut(mass2, D2).Barro13_cut()
    Bcut3 = SPlot.SelectionCut(mass3, D3).Barro13_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc, Bcut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc, Bcut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc, Bcut3)
    
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    AX.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
def plot_vdWel_cut_all(AX):

    vdWcut1 = SPlot.SelectionCut(mass1, D1).vdWel14_cut()
    vdWcut2 = SPlot.SelectionCut(mass2, D2).vdWel14_cut()
    vdWcut3 = SPlot.SelectionCut(mass3, D3).vdWel14_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc_major, vdWcut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc_major, vdWcut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc_major, vdWcut3)
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)
    
    AX.scatter(S1["bag_x"],S1["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
   
    
def plot_vDokkum_cut_all(AX):
    
    vDcut1 = SPlot.SelectionCut(mass1, D1).vDokkum15_cut()
    vDcut2 = SPlot.SelectionCut(mass2, D2).vDokkum15_cut()
    vDcut3 = SPlot.SelectionCut(mass3, D3).vDokkum15_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc, vDcut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc, vDcut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc, vDcut3)
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)
    
    AX.scatter(S1["bag_x"],S1["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)


def plot_Damjanov_cut_all(AX):
    

    Damcut1 = SPlot.SelectionCut(mass1, D1).Damjanov14_cut()
    Damcut2 = SPlot.SelectionCut(mass2, D2).Damjanov14_cut()
    Damcut3 = SPlot.SelectionCut(mass3, D3).Damjanov14_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc, Damcut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc, Damcut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc, Damcut3)

    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    AX.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
def plot_Cassata_cut_all(AX):

    Cascut1 = SPlot.SelectionCut(mass1, D1).Cassata11_cut()
    Cascut2 = SPlot.SelectionCut(mass2, D2).Cassata11_cut()
    Cascut3 = SPlot.SelectionCut(mass3, D3).Cassata11_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc, Cascut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc, Cascut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc, Cascut3)

    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    AX.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
def plot_Graham_cut_all(AX):
    

    Gcut1 = SPlot.SelectionCut(mass1, D1).Graham15_broad_cut()
    Gcut2 = SPlot.SelectionCut(mass2, D2).Graham15_broad_cut()
    Gcut3 = SPlot.SelectionCut(mass3, D3).Graham15_broad_cut()

    S1 = SSort.selection_generic(mass1, Re_1_kpc, Gcut1)
    S2 = SSort.selection_generic(mass2, Re_2_kpc, Gcut2)
    S3 = SSort.selection_generic(mass3, Re_3_kpc, Gcut3)

    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    AX.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    AX.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    AX.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    AX.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    AX.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    AX.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)


def plot_SDSS_hist():
    fig, ax = plt.subplots()        

    counts,xbins,ybins,image = plt.hist2d(np.log10(nsa_mass_IP13),np.log10(nsa_Re),bins=(1000,1000)
                                      , norm=LogNorm()
                                      , cmap = plt.cm.gray, alpha=0.3)

    plt.plot(np.log10(E1_IP13), np.log10(Re_1_kpc),'ro',ms=14)
    plt.plot(np.log10(E2_IP13), np.log10(Re_2_kpc),'bo',ms=14)
    plt.plot(np.log10(E3_IP13), np.log10(Re_3_kpc),'ko',ms=14)

    plt.show()

# Plotting size-mass diagram
def plot_sizemass_SDSS_mine():
    fig1, ax1 = plt.subplots()        


    SPlot.ShowcaseIndi.Mass_Re_plot(nsa_mass_IP13, nsa_Re, marker='x', colour = 'grey',
                                ms = 5, legend='SDSS galaxies',alpha0 = 0.2)

    #SPlot.SelectionCut(mass0,Dist0).plot_cut()
    plot_dexter_sample_Bin()
    plt.show()

def plot_sizemass_CasCut_mine(ax):

    plot_Cassata_cut_all(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Cassata",
                                                      "Cassata et al. 2011",
                                                      alpha0=0, AX=ax)
    plot_dexter_sample_all2(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Cassata","", AX=ax)

    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()

def plot_sizemass_DamCut_mine(ax):

    plot_Damjanov_cut_all(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Damjanov",
                                                      "Damjanov et al. 2014",
                                                      alpha0=0, AX=ax)
    plot_dexter_sample_all2(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Damjanov","", AX=ax)

    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()

def plot_sizemass_BarroCut_mine(ax):
    plot_Barro_cut_all(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      "Barro et al. 2013", 
                                                      alpha0=0,AX=ax)
    plot_dexter_sample_all2(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro","", AX=ax)

    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()

def plot_sizemass_vDokkumCut_mine(ax):
    
    plot_vDokkum_cut_all(ax)

    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum",
                                                      "van Dokkum et al. 2015",
                                                      alpha0=0, AX=ax)
    plot_dexter_sample_all2(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", "", AX=ax)

    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()

def plot_sizemass_vdWelCut_mine(ax):
    plot_vdWel_cut_all(ax)

    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      "van der Wel et al. 2014",
                                                      alpha0=0, AX=ax)
    plot_dexter_sample_all2_major(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", "", AX=ax)

    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()


def plot_sizemass_GrahamCut_mine(ax):
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Graham", 
                                                      "Graham et al. 2015",
                                                      alpha0=0, AX=ax)
    plot_dexter_sample_all2(ax)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Graham", "", AX=ax)

    plot_Graham_cut_all(ax)
    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],xlim[1])
    plt.show()


def plot_dexter_sample_Bin2(A):
    #SDSS
    A.scatter(nsa_mass_R15BC, nsa_Re, marker='x', color = 'grey', 
              label = "SDSS galaxies",
                                s = 16,alpha = 0.65)
    #Bin1
    A.scatter(mass1, Re_1_kpc,marker='o',c='#a5200b',label='Bin1', 
              s =70, alpha=0.7)
    A.errorbar(mass1, Re_1_kpc, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*mass1, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=0.65, marker='o')
    #Bin2
    A.scatter(mass2, Re_2_kpc,marker='o',c='#0b5786',label='Bin2', 
              s =70,alpha=0.71)
    A.errorbar(mass2, Re_2_kpc, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*mass2,ls='none',linewidth=4,
                  color = '#0b5786',
                  ecolor='#0b5786',capsize=0, alpha=0.65, marker='o')
    #Bin3
    A.scatter(mass3, Re_3_kpc,marker='o',c='#2a3236',label='Bin3', 
              s =70, alpha=0.7)
    A.errorbar(mass3, Re_3_kpc, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*mass3,ls='none',linewidth=4,
                  color = '#2a3236',
                  ecolor='#2a3236',capsize=0,
                  alpha=0.65, marker='o')   
    
    
    A.set_ylim(ylim[0],ylim[1])
    A.set_xlim(xlim[0],xlim[1])
    
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )    

def plot_dexter_sample_all2(A,alpha=0.65):
    
    # Plot E galaxies
    A.plot(mass_E, Re_kpc_E ,marker='s',color='#1b872a',label='E', 
              ms =10, alpha=1.0,linestyle="None")   
    # Plot ES galaxies
    A.plot(mass_ES, Re_kpc_ES,marker='s',color='#ccab05',label='ES', 
              ms =10, alpha=1.0,linestyle="None")
    
    #Bin1
    A.scatter(mass1, Re_1_kpc,marker='o',c='#a5200b',label='Spheroids', 
              s =70, alpha=alpha)
    A.errorbar(mass1, Re_1_kpc, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*mass1, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=alpha, marker='o')
    #Bin2
    A.scatter(mass2, Re_2_kpc,marker='o',c='#a5200b',label='', 
              s =70,alpha=alpha)
    A.errorbar(mass2, Re_2_kpc, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*mass2,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0, alpha=alpha, marker='o')
    #Bin3
    A.scatter(mass3, Re_3_kpc,marker='o',c='#a5200b',label='', 
              s =70, alpha=alpha)
    A.errorbar(mass3, Re_3_kpc, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*mass3,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0,
                  alpha=alpha, marker='o')   
    
    #plot outliers
    plot_gi_outliers(A)
    
    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )
    

def plot_dexter_sample_all2_major(A,scale='log',alpha=0.65):
    # Plot E galaxies
    A.plot(mass_E, Re_kpc_major_E ,marker='s',color='#1b872a',label='E', 
              ms =10, alpha=1.0,linestyle="None")   
    # Plot ES galaxies
    A.plot(mass_ES, Re_kpc_major_ES,marker='s',color='#ccab05',label='ES', 
              ms =10, alpha=1.0,linestyle="None")
    
    #Bin1
    A.scatter(mass1, Re_1_kpc_major, marker='o',c='#a5200b',label='Spheroids', 
              s =70, alpha=alpha)
    A.errorbar(mass1, Re_1_kpc_major, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*mass1, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=alpha, marker='o')
    #Bin2
    A.scatter(mass2, Re_2_kpc_major,marker='o',c='#a5200b',label='', 
              s =70,alpha=alpha)
    A.errorbar(mass2, Re_2_kpc_major, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*mass2,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0, alpha=alpha, marker='o')
    #Bin3
    A.scatter(mass3, Re_3_kpc_major,marker='o',c='#a5200b',label='', 
              s =70, alpha=alpha)
    A.errorbar(mass3, Re_3_kpc_major, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*mass3,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0,
                  alpha=alpha, marker='o')   
    #plot outliers
    plot_gi_outliers_major(A)
    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale(scale)
    A.set_yscale(scale)

###### Ploting ##################

xlim = [3e8,1.3e12]
ylim = [0.08,167]

mass0 = np.linspace(2e8,0.5e13,2000)
Dist0 = np.linspace(0,120,2000)


###### Ploting ##################



#6plots################
import matplotlib.gridspec as gridspec


def plot_sizemass_6plot():
    
    fig = plt.figure(figsize=(12.8, 18.4))
    gs = gridspec.GridSpec(ncols=2, nrows=3,
                               hspace=0, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])
    plot_dexter_sample_Bin2(axs0)
    axs0.legend(fontsize = 10, loc=2)
    axs0.set_ylabel(r"$R_\mathrm{e,Sph} \rm~(kpc)$",fontsize=16)
    #axs0.grid(True)

    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharey=axs0)
    plot_sizemass_DamCut_mine(axs1)
    axs1.legend(fontsize = 10, loc=4)
    plt.setp(axs1.get_yticklabels(), visible=False)
    #axs1.set_yticks([])
    #axs1.grid(True)
    
    #plot Panel (2)
    #axs1 = plt.subplot(gs[1],sharey=axs0)
    #plot_sizemass_CasCut_mine(axs1)
    #axs1.legend(loc=4)
    #plt.setp(axs1.get_yticklabels(), visible=False)
    ##axs1.set_yticks([])
    ##axs1.grid(True)
    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2])
    plot_sizemass_BarroCut_mine(axs2)
    axs2.legend(fontsize = 10, loc=4)
    axs2.set_ylabel(r"$R_\mathrm{e,Sph} \rm~(kpc)$",fontsize=16)
    #axs2.grid(True)
    
    #plot Panel (4)
    axs3 = plt.subplot(gs[3],sharey=axs2)
    plot_sizemass_vDokkumCut_mine(axs3)
    axs3.legend(fontsize = 10, loc=4)
    #axs3.set_yticks([])
    plt.setp(axs3.get_yticklabels(), visible=False)
    #axs3.grid(True)

    #plot Panel (5)
    axs4 = plt.subplot(gs[4])
    plot_sizemass_vdWelCut_mine(axs4)
    axs4.legend(fontsize = 10, loc=4)
    axs4.set_ylabel(r"$ R_\mathrm{e,Sph} \rm~(kpc)$", fontsize=16)
    axs4.set_xlabel(r"$ M_\mathrm{*,Sph} / \rm M_{\odot} (RC15)$", fontsize=16)
    #axs4.grid(True)
   
    #plot Panel (6)
    axs5 = plt.subplot(gs[5],sharey=axs4)
    plot_sizemass_GrahamCut_mine(axs5)
    axs5.legend(fontsize = 10, loc=4)
    plt.setp(axs5.get_yticklabels(), visible=False)
    #axs5.set_yticks([])
    axs5.set_xlabel(r"$ M_\mathrm{*,Sph} / \rm M_{\odot} (RC15)$", fontsize=16)
    #axs5.grid(True)

    plt.show()


plot_sizemass_6plot()

########################
#comparison with Sahu 2019

#read data
Savorgnan_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Savorgnan_sizemass.dat")
Savorgnan_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Savorgnan_sizemass.dat",dtype='str')

Davis_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Davis_sizemass.dat")
Davis_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Davis_sizemass.dat",dtype='str')

Sahu_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Sahu_sizemass2.dat")
Sahu_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Sahu_sizemass2.dat",dtype='str')

Savorgnan_name = Savorgnan_data_n[:,0]
Savorgnan_size_eq_kpc = Savorgnan_data[:,8]
Savorgnan_mass_36 = Savorgnan_data[:,10]

Savorgnan_mass_T11 = 10**(0.88 * Savorgnan_mass_36+1.02)

Davis_name = Davis_data_n[:,0]
Davis_size_eq_kpc = Davis_data[:,8]
Davis_mass_36 = Davis_data[:,10]

Davis_mass_T11 = 10**(0.88 * Davis_mass_36+1.02)

Sahu_name = Sahu_data_n[:,0]
Sahu_size_eq_kpc = Sahu_data[:,8]
Sahu_mass_36 = Sahu_data[:,10]

Sahu_mass_T11 = 10**(0.88 * Sahu_mass_36+1.02)

import SphAnalysis as SAnalysis
#Shen 2003 size-mass relation

# Lange 2016 size-mass relation
Lange2016_E = SAnalysis.AnalyticFunctions.size_mass_powerlaw(mass0,2.114,0.329)
Lange2016_E_M1e10 = SAnalysis.AnalyticFunctions.size_mass_powerlaw(mass0,1.382,0.643)
Lange2016_ETG_bulge = SAnalysis.AnalyticFunctions.size_mass_powerlaw(mass0,1.836,0.267)

s = SAnalysis.AnalyticFunctions.size_mass_powerlaw

E_T11 = np.concatenate((E1_T11, E2_T11, E3_T11,Savorgnan_mass_T11,Davis_mass_T11,Sahu_mass_T11))
Re_kpc = np.concatenate((Re_1_kpc, Re_2_kpc, Re_3_kpc,Savorgnan_size_eq_kpc,Davis_size_eq_kpc,Sahu_size_eq_kpc))

E_T11_mine = np.concatenate((E1_T11, E2_T11, E3_T11))
Re_kpc_mine = np.concatenate((Re_1_kpc, Re_2_kpc, Re_3_kpc))

Sersic_n_mine = np.concatenate((Sersic_n1, Sersic_n2, Sersic_n3))
Abs_sph_mag_mine = np.concatenate((sph_abs_mag1, sph_abs_mag2, sph_abs_mag3))
mu_mine = np.concatenate((mu_e_1, mu_e_2, mu_e_3))

########################
#Curve fit for our sample
from scipy.optimize import curve_fit

# Before graham equation size mass fit
# plot Mag_i-n graph to find the linear relatio 

def linear_func_1D(x,A,B):
    
    return A*x+B

from scipy.special import gamma
from scipy.special import gammainc

def b_value(n,x):
   bn=[]
   diff=[]
   b= np.linspace(0.1, 50.0, 50000) # works for most z
#   x=1/z=10 # 10, 2, 10/9    
#   b= np.linspace(0.1, 25.0, 2500)  #Re  2
#   b= np.linspace(0.1, 50.0, 25000)  #R10  10
#   b= np.linspace(0.1, 100.0, 50000)  #R90  10/9
   for j in b:        
        g2= x*gammainc(2*n,j)
        g22=np.round(g2,2)
        if g22==1.00: ## gammainc(2*n,j) is normalized by dividing with gamma(2*n)
            k=np.round(j,4)
            bn.append(k)
            dif=abs(1-g2)
            diff.append(dif)
   
   diff_min= min(diff)
   for s, d in zip(bn, diff):
       if d==diff_min:
           b_final=s
       else:
           continue     
   
   return b_final

def plot_Magi_to_n(Sersic_n,Abs_mag):
    fig, ax = plt.subplots()

    n_line = np.linspace(-5, 20, 30)
    
    ax.plot(Sersic_n,Abs_mag, 'o', label=r"This work")
    
    popt_lin,pcov_lin = curve_fit(linear_func_1D, np.log10(Sersic_n_mine), Abs_sph_mag_mine)
    print("Linear fit",*popt_lin)
    ax.plot(n_line, linear_func_1D(np.log10(n_line),*popt_lin), label=r"linear fit")
    
    #ax.text(10,-16, "$\rm Mag = {:.2e}+ {:.2e}log_{10}(n)$".format(popt_lin[0], popt_lin[1]))
    
    ax.legend(fontsize = 12,loc=2)
    #plt.grid(True)
    plt.ylabel(r"$\rm Mag_{sph,i-band}$",fontsize=16)
    plt.xlabel(r"$\rm n$",fontsize=16)
    plt.gca().invert_yaxis()
    ax.set_xscale( 'log' )
   
    
    plt.show()
    
def plot_Magi_to_mu0(mu,Abs_mag,Sersic_n):
    fig, ax = plt.subplots()

    mu_line = np.linspace(5, 30, 30)
    
    #calculate the corresponding b_n from sersic indices 
    bn=[]
    for i in range(len(mu)):
        n = Sersic_n[i]
        b = SAnalysis.b_value(n, 1.0/0.5)
        bn.append(b)
        
    bn = np.array(bn)
    
    mu0 = mu-2.5*(bn/np.log(10))
    
    ax.plot(mu0,Abs_mag, 'o', label=r"This work")
    
    popt_mu,pcov_mu = curve_fit(linear_func_1D, mu0, Abs_sph_mag_mine)
    print("Linear fit",*popt_mu)
    ax.plot(mu_line, linear_func_1D(mu_line,*popt_mu), label=r"linear fit")
    
    #ax.text(10,-16, "$\rm Mag = {:.2e}+ {:.2e}log_{10}(n)$".format(popt_lin[0], popt_lin[1]))
    
    ax.legend(fontsize = 12,loc=2)
    #plt.grid(True)
    plt.ylabel(r"$\rm Mag_{sph,i-band}$",fontsize=16)
    plt.xlabel(r"$\rm \mu_{0,sph,i-band}$",fontsize=16)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()

    #ax.set_xscale( 'log' )    
    plt.show()

#plot_Magi_to_n(Sersic_n_mine,Abs_sph_mag_mine)
#plot_Magi_to_mu0(mu_mine,Abs_sph_mag_mine,Sersic_n_mine)

## Nandini's script insert(temp) #####################


#l=55
z=0.5 #z=0.1 for 10% radius, z=0.5 for half-light radius, z=0.9 for 90% radius
Rz=[]  # radius
#mag=[] # galaxy absolute mag
Mass_gal=[] #galaxy mass
for i in range(70):
    mag1= 0.0-10-0.1*i   #range of magnitude -11.9 to -22.9 (better take -12 to -23)
    #mag.append(mag1)
    mass=(1.6*10**(0.4*(4.65-mag1))) #1.6 is a typical T11 M/L ratio, 4.65 is M_sun
    Mass_gal.append(mass)
    n1=10**(0.0 - (mag1 +14.3)/9.4)   #equation 16 from Graham(2019)
    bz=b_value(n1,1/z)  #1/0.9        #exact value of bn, defined on the top
    fnz= (z*2*n1*np.exp(bz)*gamma(2*n1))/(bz**(2*n1)) #equation 20 from Graham(2019)
    #equation 25 (combines eq.17 and eq.24) from Graham(2019) 
    R_z= (mag1/10) + 0.5*(np.log10(z)-np.log10(fnz))+0.217*bz+1.2874 #in log10 #  
    Rz.append(R_z)
    
xdata4=np.array(Mass_gal)
ydata4=10**np.array(Rz)


def graham_equ(mass, B):
    """
    Graham 2019  Equation 25
    
    R_e = (mag/10) + 0.5*(np.log10(0.5)-np.log10(fnz))+0.217*bz+1.2874
    
    """
    mag = 4.65 - (2.5)*np.log10(mass/1.6)

    n1=10**(0.0 - (mag -19.6)/-3.8)   #equation 16 from Graham(2019)
    
    
    bn=[]
    for i in range(len(Sersic_n_mine)):
        n = Sersic_n_mine[i]
        b = SAnalysis.b_value(n, 1.0/0.5)
        bn.append(b)
        
    #bz=b_value(n1,1.0/0.5)  #1/0.9   #exact value of bn, defined on the top
    fnz= (0.5*2*n1*np.exp(bz)*gamma(2*n1))/(bn**(2*n1)) #equation 20 from Graham(2019)
    #fnzz.append(fnz)    
    
    R_e = (mag/10) + 0.5*(np.log10(0.5) - np.log10(fnz)) +  bz/(2*np.log(10)) + B


    return R_e
   
## End Nadini'script  ################################

ss = graham_equ

#print(len(E_T11),len(Re_kpc))

##fitting
#popt,pcov = curve_fit(s, E_T11, Re_kpc)
#popt2,pcov2 = curve_fit(s,E_T11_mine, Re_kpc_mine)
#popt3,pcov3 = curve_fit(s,Savorgnan_mass_T11, Savorgnan_size_eq_kpc)
#popt4,pcov4 = curve_fit(s,Davis_mass_T11, Davis_size_eq_kpc)
#popt5,pcov5 = curve_fit(s,Sahu_mass_T11, Sahu_size_eq_kpc)
#
#popt_g,pcov_g = curve_fit(ss, E_T11, Re_kpc)

#print("graham_equ_para", *popt_g )
#plotting

def plot_sizemass_z0comparison():
    fig, ax = plt.subplots()
    
    plot_dexter_sample_all2_T11(ax,alpha = 0.2)

    ax.plot(mass0,Lange2016_E, '--', label=r"E in Lange et al. 2016")
    ax.plot(mass0,Lange2016_ETG_bulge, '--', label=r"early-type bulge in Lange et al. 2016")
    ax.plot(mass0,Lange2016_E_M1e10, '--', label=r"E ($M_*<10^{10}$)in Lange et al. 2016")
  
    ax.plot(mass0,s(mass0, *popt), "-",linewidth=6, label=r"our fit")
    ax.plot(mass0,s(mass0, *popt2), "-",linewidth=6, label=r"my fit", color='#a5200b')
    ax.plot(mass0,ss(mass0, *popt_g),"--",linewidth=10, label=r"Graham equ.")

   # ax.plot(mass0,s(mass0, *popt3), "-", linewidth=6, label=r"Savorgnan et al. 2016 fit", color='#b940c8')
   # ax.plot(mass0,s(mass0, *popt4), "-", linewidth=6, label=r"Davis et al. 2019 fit", color='#2e417b')
   # ax.plot(mass0,s(mass0, *popt5), "-", linewidth=6, label=r"Sahu et al. 2019 fit", color='#e1a000')
    
    
    SPlot.ShowcaseIndi.Mass_Re_plot(Savorgnan_mass_T11, Savorgnan_size_eq_kpc, 
                                    yerr = None,
                                xerr = None,
                                colour='#b940c8',
                                name=None,legend='Savorgnan et al. 2016',
                                ms=10,alpha0 = 0.4,lw=3)
        
    SPlot.ShowcaseIndi.Mass_Re_plot(Davis_mass_T11, Davis_size_eq_kpc, 
                                    yerr = None,
                                xerr = None,
                                colour='#2e417b',
                                name=None,legend='Davis et al. 2019',
                                ms=10,alpha0 = 0.4,lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(Sahu_mass_T11, Sahu_size_eq_kpc, 
                                    yerr = None,
                                xerr = None,
                                colour='#e1a000',
                                name=None,legend='Sahu et al. 2019',
                                ms=10,alpha0 = 0.4,lw=3)


    ax.set_xlim(left = xlim[0], right = xlim[1])
    ax.set_ylim(bottom = ylim[0], top = ylim[1])

    ax.legend(fontsize = 12,loc=2)
    #plt.grid(True)
    
    plt.xlabel(r"$\rm M_{*,sph}$ / $M_{\odot}$",fontsize=16)
    plt.ylabel("$R_{e,sph}$ (kpc)",fontsize=16)
    
    plt.show()
    
    
#plot_sizemass_z0comparison()


########################