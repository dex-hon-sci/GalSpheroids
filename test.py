#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 05:57:25 2021

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import numpy as np
import pickle

from astropy.cosmology import FlatLambdaCDM


M_sun = 4.53

#/home/dexter/result/stat/Into_selection_all
#/home/dexter/result/distance/
#vel_disp_list_all_mag.txt

# Name & RA (2000) & Dec (2000) & $Dist$  & $g$ & $i$ & scale  & Seeing & $\sigma_{velocity}$ & $morph$ & $Mag_{i-band}$ & \log $(M_*/M_\odot)_{total}$\\


##############################The Mini Mini no Mi which The Mini Mini no Mi which turns the user smaller, while Wolf's Devil Fruit makes him bigger.turns the user smaller, while Wolf's Devil Fruit makes him bigger.
# make table for the parent data


D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4.txt",
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

# calculate the ellipticity of the Sersic dist2D fit
b_a_1 = D0_Bin1_table[:,34]
b_a_2 = D0_Bin2_table[:,34]
b_a_3 = D0_Bin3_table[:,34]

seeing1 = D0_Bin1_table[:,15]
seeing2 = D0_Bin2_table[:,15]
seeing3 =  D0_Bin3_table[:,15]

# calculate the Radius in equivalent axis 
#Sersic2D_50rad_1 = D0_Bin1_table[:,33]*np.sqrt(1-(1-(b_a_1)**2))
#Sersic2D_50rad_2 = D0_Bin2_table[:,33]*np.sqrt(1-(1-(b_a_2)**2))
#Sersic2D_50rad_3 = D0_Bin3_table[:,33]*np.sqrt(1-(1-(b_a_3)**2))

Sersic2D_50rad_1 = D0_Bin1_table[:,33]
Sersic2D_50rad_2 = D0_Bin2_table[:,33]
Sersic2D_50rad_3 = D0_Bin3_table[:,33]

#Get the name of the galaxies
name_D1 = D0_Bin1_table_n[:,0]
name_D2 = D0_Bin2_table_n[:,0]
name_D3 = D0_Bin3_table_n[:,0]

#Get the morphology of the galaxies
morph1 = D0_Bin1_table_n[:,17]
morph2 = D0_Bin2_table_n[:,17]
morph3 = D0_Bin3_table_n[:,17]

morph1_new = D0_Bin1_table_n[:,-2]
morph2_new = D0_Bin2_table_n[:,-2]
morph3_new = D0_Bin3_table_n[:,-2]

corr_mass1 = D0_Bin1_table[:,13]
corr_mass2 = D0_Bin2_table[:,13]
corr_mass3 = D0_Bin3_table[:,13]

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

sph_abs_mag1 = SPlot.MassCalculation(mag_g1, mag_i1,sph_mag1, D1, M_sun).cal_abs_mag()
sph_abs_mag2 = SPlot.MassCalculation(mag_g2, mag_i2,sph_mag2, D2, M_sun).cal_abs_mag()
sph_abs_mag3 = SPlot.MassCalculation(mag_g3, mag_i3,sph_mag3, D3, M_sun).cal_abs_mag()
print(D1)
print("sph_mag1",sph_mag1)
print("sph_abs_mag1",sph_mag1-25-5*np.log10(D1))
print("sph_abs_mag1",sph_abs_mag1)

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

Sersic2D_50rad_1_kpc = Sersic2D_50rad_1*scale1
Sersic2D_50rad_2_kpc = Sersic2D_50rad_2*scale2
Sersic2D_50rad_3_kpc = Sersic2D_50rad_3*scale3

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