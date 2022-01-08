#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 01:48:10 2021

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

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

##Define the selection volume
D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180))-np.cos(np.pi/2))

#V1,V2,V3=voll[2],voll[1],voll[0]

V1_V = voll[2]-voll[1]
V2_V = voll[1]-voll[0]
V3_V = voll[0]

V1,V2,V3 = V1_V,V2_V,V3_V
#%%
# Take away the VCC gaalxies from Bin3
D0_Bin3_table_noVirgo = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_noVirgo.txt")

D0_Bin3_table_noVirgo_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_noVirgo.txt",
    dtype = 'str')

K_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr_EXT_noVirgo.dat")

# The apparant magnitude of the galaxy in g- and i-band.
mag_g3, mag_i3 = K_table3[:,10], K_table3[:,9]
g3_EXT, i3_EXT = K_table3[:,23], K_table3[:,24]
g3_kcorr, i3_kcorr = K_table3[:,25], K_table3[:,26]

# the corrected mag g and i, Kcorrection+EXTINCTIOn
mag_g3_corr, mag_i3_corr = mag_g3-g3_kcorr, mag_i3-i3_kcorr

# The distance (final decision), lower limit and upper limit
D3, D3_lerr, D3_uerr = D0_Bin3_table_noVirgo[:,29], D0_Bin3_table_noVirgo[:,30], D0_Bin3_table_noVirgo[:,31]

#extended disk ellipicity
elle3 = D0_Bin3_table_noVirgo[:,-1] 

#Get the name of the galaxies
name_D3 = D0_Bin3_table_noVirgo_n[:,0]

#Get the morphology of the galaxies in RC3
morph3 = D0_Bin3_table_noVirgo_n[:,17]

# Get the new morphology given by us 
morph3_new = D0_Bin3_table_noVirgo_n[:,-2]


name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_noVirgo_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_noVirgo_cpt")

sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_noVirgo_cpt", ["Bulge","CoreBulge"])

Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_noVirgo_cpt", ["Bulge","CoreBulge"], 1) #get Re

Sersic_n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_noVirgo_cpt", ["Bulge","CoreBulge"], 2) #get n

mu_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_noVirgo_cpt", ["Bulge","CoreBulge"], 0) #get n

core_sersic_mu_p3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_noVirgo_cpt", ["CoreBulge"], 0)

#Get Re major axis
Re_3_major = SRead.grab_parameter("F_Gal_bundle_major_Bin3V_noVirgo_cpt", ["Bulge","CoreBulge"], 1) #get Re

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale3 = D3* ars
scale3_lerr, scale3_uerr =  (D3-D3_lerr)*ars, (D3+D3_uerr)*ars

Re_3_kpc = Re_3* scale3
Re_3_kpc_major = Re_3_major* scale3
Re_3_kpc_lerr, Re_3_kpc_uerr = abs(Re_3* scale3_lerr - Re_3_kpc), abs(Re_3* scale3_uerr - Re_3_kpc)

Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]

# spheroid mag correction
sph_mag3 = sph_mag3 - i3_EXT - i3_kcorr
################################
#Calculate mass 

# "F_Gal_bundle_equvi_Bin3V_noVirgo_cpt"


ML_select3_IP13_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Into13_MassRatio
ML_select3_R15BC_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Roediger15BC03_MassRatio
ML_select3_Z09_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Zibetti09_MassRatio
ML_select3_T11_K = SPlot.MLRelationIband(mag_g3_corr,mag_i3_corr).Taylor11_MassRatio

M3_K = SPlot.MassCalculation(sph_mag3, D3, 4.53, mag_g3_corr,mag_i3_corr)

E3_IP13_K = M3_K.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K = M3_K.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K = M3_K.cal_Mass(ML_select3_Z09_K)
E3_T11_K = M3_K.cal_Mass(ML_select3_T11_K)

################################
#calculate the mass error

MLR3 = ML_select3_T11_K
MLR_e3 = 10**0.1

mag_e = 0.3 #magnitude error

mass_uerr3 = np.sqrt(((mag_e/2.5)**2)+((2*D3_uerr/(D3*np.log(10)))**2)+((MLR_e3/(MLR3*np.log(10)))**2))
mass_err3 = mass_uerr3
mass3 = E3_IP13_K
    
# perform the size-mass cut
Bcut3_Barro = SPlot.SelectionCut(mass3, D3).Barro13_cut()
S3_Barro = SSort.selection_generic(mass3, Re_3_kpc, Bcut3_Barro)    
n3_Barro = np.size(S3_Barro["bag_y"]) 

Bcut3_vdWel = SPlot.SelectionCut(mass3, D3).vdWel14_cut()
S3_vdWel = SSort.selection_generic(mass3, Re_3_kpc_major, Bcut3_vdWel)    
n3_vdWel = np.size(S3_vdWel["bag_y"]) 

Bcut3_vDokkum = SPlot.SelectionCut(mass3, D3).vDokkum15_cut()
S3_vDokkum = SSort.selection_generic(mass3, Re_3_kpc, Bcut3_vDokkum)    
n3_vDokkum = np.size(S3_vDokkum["bag_y"]) 

Bcut3_Dam = SPlot.SelectionCut(mass3, D3).Damjanov14_cut()
S3_Dam = SSort.selection_generic(mass3, Re_3_kpc, Bcut3_Dam)
n3_Dam = np.size(S3_Dam["bag_y"]) 
    
    

print("N_Barro", n3_Barro, "n",n3_Barro/V3)
print("N_vdWel", n3_vdWel, "n",n3_vdWel/V3)
print("N_vDokkum",n3_vDokkum, "n",n3_vDokkum/V3)
print("N_Dam",n3_Dam, "n",n3_Dam/V3)
