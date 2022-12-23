#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 03:02:36 2020
z
@author: dexter

The script to make the table for host galaxies selection and the spheroids parameters.
"""


import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAna

import numpy as np
import pickle

from astropy.cosmology import FlatLambdaCDM


M_sun = 4.53

#/home/dexter/result/stat/Into_selection_all
#/home/dexter/result/distance/
#vel_disp_list_all_mag.txt

# Name & RA (2000) & Dec (2000) & $Dist$  & $g$ & $i$ & scale  & Seeing & $\sigma_{velocity}$ & $morph$ & $Mag_{i-band}$ & \log $(M_*/M_\odot)_{total}$\\


def ScS_type_generate(file_name):
    table = SRead.read_list(file_name)
    ScStype= []
    for i in range(len(table)):
        if "Bulge" in table[i]:
            ScStype.append("S")
        elif "CoreBulge" in table[i]:
            ScStype.append("cS")
    return ScStype

def cal_mu_e(file_name):
    table = SRead.read_list(file_name)
    mu_e_list = []
    for i in range(len(table)):
        para_whole = SRead.grab_parameter_whole(file_name, ["Bulge","CoreBulge"])
        if "Bulge" in table[i]:
            mu_e = para_whole[i][0]
            mu_e_list.append(mu_e)
            
            #print(i)
            #print("mu_e, r_e, n_ser")
            #print(mu_e, para_whole[i][1], para_whole[i][2])
        elif "CoreBulge" in table[i]:
            mu_p = para_whole[i][0] 
            r_e=  para_whole[i][1] 
            n_ser =  para_whole[i][2] 
            r_b =  para_whole[i][3] 
            al =  para_whole[i][4] 
            ga =  para_whole[i][5] 
            #print(i)
            #print("mu_p, r_e, n_ser, r_b, al, ga")
            #print(mu_p, r_e, n_ser, r_b, al, ga)
            mu_e = SAna.AnalyticFunctions.mu_core_sersic_func(r_e, mu_p, r_e, n_ser, r_b, al, ga)
            mu_e_list.append(mu_e)
    mu_e_list = np.array(mu_e_list)
    return mu_e_list




##############################
# make table for the parent data
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

ScStype1 = ScS_type_generate("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
ScStype2 = ScS_type_generate("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
ScStype3 = ScS_type_generate("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

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

# the corrected mag g and i, Kcorrection+EXTINCTIOn(extinction already accounted for in the magnitude)
mag_g1_corr, mag_i1_corr = mag_g1-g1_kcorr, mag_i1-i1_kcorr
mag_g2_corr, mag_i2_corr = mag_g2-g2_kcorr, mag_i2-i2_kcorr
mag_g3_corr, mag_i3_corr = mag_g3-g3_kcorr, mag_i3-i3_kcorr

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

print('Check', mag_g1_corr[1], mag_i1_corr[1],  ML_select1_IP13_corr[1],
      corr_D1[1], Mag_i1_kcorr_cDis[1], E1_gal_cDis[1])


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

############# calculating total stellar mass, both SDSS mag and profiler mag
#Bin1
#SDSS M/L
M1_K_SE_SDSS = SPlot.MassCalculation(mag_i1, D1, 4.53,mag_g1_corr,mag_i1_corr)
#profiler M/L
M1_K_SE_prof = SPlot.MassCalculation(total_mag1, D1, 4.53,mag_g1_corr,mag_i1_corr)

E1_IP13_K_SE_SDSS = M1_K_SE_SDSS.cal_Mass(ML_select1_IP13_K)
E1_R15BC_K_SE_SDSS = M1_K_SE_SDSS.cal_Mass(ML_select1_R15BC_K)
E1_Z09_K_SE_SDSS = M1_K_SE_SDSS.cal_Mass(ML_select1_Z09_K)
E1_T11_K_SE_SDSS = M1_K_SE_SDSS.cal_Mass(ML_select1_T11_K)

E1_IP13_K_SE_prof = M1_K_SE_prof.cal_Mass(ML_select1_IP13_K)
E1_R15BC_K_SE_prof = M1_K_SE_prof.cal_Mass(ML_select1_R15BC_K)
E1_Z09_K_SE_prof = M1_K_SE_prof.cal_Mass(ML_select1_Z09_K)
E1_T11_K_SE_prof = M1_K_SE_prof.cal_Mass(ML_select1_T11_K)

#Bin2
#SDSS M/L
M2_K_SE_SDSS = SPlot.MassCalculation(mag_i2, D2, 4.53,mag_g2_corr,mag_i2_corr)
#profiler M/L
M2_K_SE_prof = SPlot.MassCalculation(total_mag2, D2, 4.53,mag_g2_corr,mag_i2_corr)

E2_IP13_K_SE_SDSS = M2_K_SE_SDSS.cal_Mass(ML_select2_IP13_K)
E2_R15BC_K_SE_SDSS = M2_K_SE_SDSS.cal_Mass(ML_select2_R15BC_K)
E2_Z09_K_SE_SDSS = M2_K_SE_SDSS.cal_Mass(ML_select2_Z09_K)
E2_T11_K_SE_SDSS = M2_K_SE_SDSS.cal_Mass(ML_select2_T11_K)

E2_IP13_K_SE_prof = M2_K_SE_prof.cal_Mass(ML_select2_IP13_K)
E2_R15BC_K_SE_prof = M2_K_SE_prof.cal_Mass(ML_select2_R15BC_K)
E2_Z09_K_SE_prof = M2_K_SE_prof.cal_Mass(ML_select2_Z09_K)
E2_T11_K_SE_prof = M2_K_SE_prof.cal_Mass(ML_select2_T11_K)

#Bin3
#SDSS M/L
M3_K_SE_SDSS = SPlot.MassCalculation(mag_i3, D3, 4.53,mag_g3_corr,mag_i3_corr)
#profiler M/L
M3_K_SE_prof = SPlot.MassCalculation(total_mag3, D3, 4.53,mag_g3_corr,mag_i3_corr)

E3_IP13_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_Z09_K)
E3_T11_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_T11_K)

E3_IP13_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_Z09_K)
E3_T11_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_T11_K)

##########################################
#K correction and extinction for sph mag
sph_mag1 = sph_mag1 - i1_EXT - i1_kcorr
sph_mag2 = sph_mag2 - i2_EXT - i2_kcorr
sph_mag3 = sph_mag3 - i3_EXT - i3_kcorr

sph_abs_mag1 = SPlot.MassCalculation(mag_g1, mag_i1,sph_mag1, D1, M_sun).cal_abs_mag(sph_mag1, D1)
sph_abs_mag2 = SPlot.MassCalculation(mag_g2, mag_i2,sph_mag2, D2, M_sun).cal_abs_mag(sph_mag2, D2)
sph_abs_mag3 = SPlot.MassCalculation(mag_g3, mag_i3,sph_mag3, D3, M_sun).cal_abs_mag(sph_mag3, D3)


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

################################
#calculate the mass error

mag_e = 0.3 #magnitude error
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


mu_e_1 = cal_mu_e("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
mu_e_2 = cal_mu_e("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
mu_e_3 = cal_mu_e("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")


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

#########################################
#########################################
#########################################

# making the table

name = name1
RA, DEC = RA_1, DEC_1
Dist, Dist_lerr, Dist_uerr = corr_D1, D1_lerr, D1_uerr
mag_g, mag_i = mag_g1_corr, mag_i1_corr
scale = scale1
Mag_i = Mag_i1_kcorr_cDis
seeing = seeing1
morph = morph1
morph_new = morph1_new
mass = E1_gal_cDis

Dist_final = Dist
Rmax = Rmax1
#Dist_final = SSort.str_zipping_generic( "$",
#                                        list(np.around(Dist,decimals=2)),
#                                       "\substack{+",
#                                       list(np.around(Dist_uerr,decimals=2)),
#                                       "\\\-",
#                                       list(np.around(Dist_lerr,decimals=2)),
#                                       "}$")


table = {"Name": name,
         "RA": np.around(RA, decimals=1),
         "DEC": np.around(DEC, decimals=1),
         "Dist": np.around(Dist_final,decimals=1),
         "mag_g": np.around(mag_g, decimals=1) ,
         "mag_i": np.around(mag_i, decimals=1) ,
         #"Scale": np.around(scale, decimals=2) ,
         "seeing": np.around(seeing, decimals=1),
         "morph (old)": morph,
         "morph (new)": morph_new,
         "Mag_i": np.around(Mag_i, decimals=1),
         "mass": np.around(mass/1e10, decimals=1),
         "R_max": np.around(Rmax1, decimals=1)}

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))

# save the table in dict form    
with open("Gal_table1_bin1V", 'wb') as f:
    pickle.dump(table, f)

# convert the dict to ascii file
SRead.convert_dict_ascii("Gal_table1_bin1V","Gal_table1_bin1V.txt")

##################################
#Bundle = "/home/dexter/result/Gal_bundle_equvi_bin4_cpt"

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")



name = name1
DDist = D1
ScStype = ScStype1
mu_e = mu_e_1
R_e = Re_1_kpc
n = n1
mag_sph =  sph_mag1
Mag_sph = sph_abs_mag1



Taylor_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E1_T11_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(mass_T11_uerr1,decimals=1)),
                                       "$")


Zibetti_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E1_Z09_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(mass_Z09_uerr1,decimals=1)),
                                       "$")

Roediger_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E1_R15BC_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(mass_R15BC_uerr1,decimals=1)),
                                       "$")

Into_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E1_IP13_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(mass_IP13_uerr1,decimals=1)),
                                       "$")




table = {"Name": name,
         "Dist": np.around(DDist,decimals=1),
         "type": ScStype,
         "mu_e": np.around(mu_e, decimals=1),
         "R_e": np.around(R_e, decimals=1),
         "n": np.around(n, decimals=1),
         "mag_sph": np.around(mag_sph, decimals=1),
         "Mag_sph": np.around(Mag_sph, decimals=1),
         "Taylor_mass": Taylor_mass,
         "Zibetti_mass":Zibetti_mass, 
         "Roediger_mass": Roediger_mass,
         "Into_mass": Into_mass} 

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table2_bin1V", 'wb') as f:
    pickle.dump(table, f)

SRead.convert_dict_ascii("Gal_table2_bin1V","Gal_table2_bin1V.txt")


#########################################

# making the table
# Bin2
name = name2
RA, DEC = RA_2, DEC_2
Dist, Dist_lerr, Dist_uerr = corr_D2, D2_lerr, D2_uerr
mag_g, mag_i = mag_g2_corr, mag_i2_corr
scale = scale2
Mag_i = Mag_i2_kcorr_cDis
seeing = seeing2
morph = morph2
morph_new = morph2_new
mass = E2_gal_cDis


Dist_final = Dist
Rmax = Rmax2
#Dist_final = SSort.str_zipping_generic( "$",
#                                        list(np.around(Dist,decimals=2)),
#                                       "\substack{+",
#                                       list(np.around(Dist_uerr,decimals=2)),
#                                       "\\\-",
#                                       list(np.around(Dist_lerr,decimals=2)),
#                                       "}$")


table = {"Name": name,
         "RA": np.around(RA, decimals=1),
         "DEC": np.around(DEC, decimals=1),
         "Dist": np.around(Dist_final,decimals=1),
         "mag_g": np.around(mag_g, decimals=1),
         "mag_i": np.around(mag_i, decimals=1),
         #"Scale": np.around(scale, decimals=2) ,
         "seeing": np.around(seeing, decimals=1),
         "morph (old)": morph,
         "morph (new)": morph_new,
         "Mag_i": np.around(Mag_i, decimals=1),
         "mass": np.around(mass/1e10, decimals=1),
         "R_max": np.around(Rmax, decimals=1)}

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))

# save the table in dict form    
with open("Gal_table1_bin2V", 'wb') as f:
    pickle.dump(table, f)

# convert the dict to ascii file
SRead.convert_dict_ascii("Gal_table1_bin2V","Gal_table1_bin2V.txt")

##################################
#Bundle = "/home/dexter/result/Gal_bundle_equvi_bin4_cpt"
name = name2
DDist = D2
ScStype = ScStype2
mu_e = mu_e_2
R_e = Re_2_kpc
n = n2
mag_sph =  sph_mag2
Mag_sph = sph_abs_mag2


Taylor_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E2_T11_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E2_T11_K+mass_T11_uerr2)-np.log10(E2_T11_K),decimals=1)),
                                       "$")


Zibetti_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E2_Z09_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E2_Z09_K+mass_Z09_uerr2)-np.log10(E2_Z09_K),decimals=1)),
                                       "$")

Roediger_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E2_R15BC_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E2_R15BC_K+mass_R15BC_uerr2)-np.log10(E2_R15BC_K),decimals=1)),
                                       "$")

Into_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E2_IP13_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E2_IP13_K+mass_IP13_uerr2)-np.log10(E2_IP13_K),decimals=1)),
                                       "$")


table = {"Name": name,
         "Dist": np.around(DDist,decimals=1),
         "type": ScStype,
         "mu_e": np.around(mu_e, decimals=1),
         "R_e": np.around(R_e, decimals=1),
         "n": np.around(n, decimals=1),
         "mag_sph": np.around(mag_sph, decimals=1),
         "Mag_sph": np.around(Mag_sph, decimals=1),
         "Taylor_mass": Taylor_mass,
         "Zibetti_mass":Zibetti_mass, 
         "Roediger_mass": Roediger_mass,
         "Into_mass": Into_mass} 

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table2_bin2V", 'wb') as f:
    pickle.dump(table, f)

SRead.convert_dict_ascii("Gal_table2_bin2V","Gal_table2_bin2V.txt")

#########################################

# making the table
# Bin3
name = name3
RA, DEC = RA_3, DEC_3
Dist, Dist_lerr, Dist_uerr = corr_D3, D3_lerr, D3_uerr
mag_g, mag_i = mag_g3_corr, mag_i3_corr
scale = scale3
Mag_i = Mag_i3_kcorr_cDis
seeing = seeing3
morph = morph3
morph_new = morph3_new
mass = E3_gal_cDis

Dist_final = Dist
Rmax = Rmax3
#Dist_final = SSort.str_zipping_generic( "$",
#                                        list(np.around(Dist,decimals=2)),
#                                       "\substack{+",
#                                      list(np.around(Dist_uerr,decimals=2)),
#                                       "\\\-",
#                                       list(np.around(Dist_lerr,decimals=2)),
#                                       "}$")


table = {"Name": name,
         "RA": np.around(RA, decimals=1),
         "DEC": np.around(DEC, decimals=1),
         "Dist": np.around(Dist_final,decimals=1),
         "mag_g": np.around(mag_g, decimals=1) ,
         "mag_i": np.around(mag_i, decimals=1) ,
         #"Scale": np.around(scale, decimals=2) ,
         "seeing": np.around(seeing, decimals=1),
         "morph (old)": morph,
         "morph (new)": morph_new,
         "Mag_i": np.around(Mag_i, decimals=1),
         "mass": np.around(mass/1e10, decimals=1),
         "R_max": np.around(Rmax, decimals=1)}         

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))

# save the table in dict form    
with open("Gal_table1_bin3V", 'wb') as f:
    pickle.dump(table, f)

# convert the dict to ascii file
SRead.convert_dict_ascii("Gal_table1_bin3V","Gal_table1_bin3V.txt")

##################################
#Bundle = "/home/dexter/result/Gal_bundle_equvi_bin4_cpt"
name = name3
DDist = D3
ScStype = ScStype3
mu_e = mu_e_3
R_e = Re_3_kpc
n = n3
mag_sph =  sph_mag3
Mag_sph = sph_abs_mag3


Taylor_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E3_T11_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E3_T11_K+mass_T11_uerr3)-np.log10(E3_T11_K),decimals=1)),
                                       "$")


Zibetti_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E3_Z09_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E3_Z09_K+mass_Z09_uerr3)-np.log10(E3_Z09_K),decimals=1)),
                                       "$")

Roediger_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E3_R15BC_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E3_R15BC_K+mass_R15BC_uerr3)-np.log10(E3_R15BC_K),decimals=1)),
                                       "$")

Into_mass = SSort.str_zipping_generic('$',list(np.around(np.log10(E3_IP13_K),decimals=1)), 
                                       "\pm",
                                       list(np.around(np.log10(E3_IP13_K+mass_IP13_uerr3)-np.log10(E3_IP13_K),decimals=1)),
                                       "$")


table = {"Name": name,
         "Dist": np.around(DDist,decimals=1),
         "type": ScStype,
         "mu_e": np.around(mu_e, decimals=1),
         "R_e": np.around(R_e, decimals=1),
         "n": np.around(n, decimals=1),
         "mag_sph": np.around(mag_sph, decimals=1),
         "Mag_sph": np.around(Mag_sph, decimals=1),
         "Taylor_mass": Taylor_mass,
         "Zibetti_mass":Zibetti_mass, 
         "Roediger_mass": Roediger_mass,
         "Into_mass": Into_mass} 

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table2_bin3V", 'wb') as f:
    pickle.dump(table, f)

SRead.convert_dict_ascii("Gal_table2_bin3V","Gal_table2_bin3V.txt")

#Mag_i = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")


#ML_select_4 = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
#M_4 = SPlot.MassCalculation(Mag_i, Dist, M_sun,mag_g,mag_i)
#mass_4 = M_4.cal_Mass(ML_select_4)


#   Name & $\mu_e$ & $R_e$ & $n$  & $mag_{Sph}$ & $Abs Mag_{Sph}$ & \log$(M_*/M_\odot)_{Taylor}$  & log$(M_*/M_\odot)_{Into}$ & log$(M_*/M_\odot)$ & log$(M_*/M_\odot)$  \\
