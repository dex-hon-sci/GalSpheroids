#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 18:22:13 2021

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

# Read all galaxies

All_gal = SRead.read_table("/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt")

all_gal_RA, all_gal_DEC = All_gal[:,1], All_gal[:,2]

# Read the host galaxy candidate 

Host_gal = SRead.read_table("/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

host_gal_RA, host_gal_DEC = Host_gal[:,1], Host_gal[:,2]

# Read the Virgo CLuster Galaxy

VCC_gal = SRead.read_table("/home/dexter/result/stat/completeness/VCC_cross_list.txt")

VCC_gal_RA, VCC_gal_DEC = VCC_gal[:,2], VCC_gal[:,3]

# plot the RA DEC diagram

def plot_RA_DEC_VCC():
    fig = plt.figure(figsize=(6.4, 5.8))

    plt.plot(all_gal_RA, all_gal_DEC,'o', color = "#8e8e8e" ,ms = 1.0, alpha = 0.6, label= r"$\rm All$")
    plt.plot(host_gal_RA, host_gal_DEC,'o', color = "b",  ms =10, label = r"$\rm Our~galaxies$")
    plt.plot(VCC_gal_RA, VCC_gal_DEC,'ro',ms=10, label= r"$\rm Our~Virgo~cluster~galaxies$")
    plt.xlabel(r"$R.A.~(\rm deg)$", fontsize=16)
    plt.ylabel(r"$Dec.~(\rm deg)$", fontsize=16)    
    plt.legend(numpoints = 1, loc=2,fontsize=12)
    
    plt.xlim(135,216)
    plt.ylim(-5,74)
    
    plt.tight_layout()

    return fig

plot_RA_DEC_VCC()

# Load the size-mass stuff
mass0 = np.linspace(2e8,0.5e13,2000)
Dist0 = np.linspace(0,120,2000)

D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_3.txt")

D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_3.txt",
    dtype = 'str')

# The apparant magnitude of the galaxy in g- and i-band.
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

# The distance (final decision), lower limit and upper limit
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

#extended disk ellipicity
elle3 = D0_Bin3_table[:,-1] 

#Get the name of the galaxies
name_D3 = D0_Bin3_table_n[:,0]

#Get the morphology of the galaxies in RC3
morph3 = D0_Bin3_table_n[:,17]

# Get the new morphology given by us 
morph3_new = D0_Bin3_table_n[:,-2]

############### reading result files###############
name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])

Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

Sersic_n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 2) #get n

mu_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 0) #get n

core_sersic_mu_p3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["CoreBulge"], 0)

#Get Re major axis
Re_3_major = SRead.grab_parameter("F_Gal_bundle_major_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale3 = D3* ars
scale3_lerr, scale3_uerr =  (D3-D3_lerr)*ars, (D3+D3_uerr)*ars

Re_3_kpc = Re_3* scale3
Re_3_kpc_major = Re_3_major* scale3
Re_3_kpc_lerr, Re_3_kpc_uerr = abs(Re_3* scale3_lerr - Re_3_kpc), abs(Re_3* scale3_uerr - Re_3_kpc)

Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]

################################
#Calculate mass with K-correction

# "F_Gal_bundle_equvi_Bin3V_noVirgo_cpt"

K_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr.dat")
K_table3_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr.dat", 
    dtype='str')

K_name3 = K_table3_n[:,4]

mag_g3_kcorr, mag_i3_kcorr = K_table3[:,19], K_table3[:,18]

ML_select3_IP13_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Into13_MassRatio
ML_select3_R15BC_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Roediger15BC03_MassRatio
ML_select3_Z09_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Zibetti09_MassRatio
ML_select3_T11_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Taylor11_MassRatio

M3_K = SPlot.MassCalculation(sph_mag3, D3, 4.53, mag_g3_kcorr,mag_i3_kcorr)

E3_IP13_K = M3_K.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K = M3_K.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K = M3_K.cal_Mass(ML_select3_Z09_K)
E3_T11_K = M3_K.cal_Mass(ML_select3_T11_K)

################################
#calculate the mass error

MLR3 = ML_select3_R15BC_K
MLR_e3 = 10**0.1

mag_e = 0.3 #magnitude error

mass_uerr3 = np.sqrt(((mag_e/2.5)**2)+((2*D3_uerr/(D3*np.log(10)))**2)+((MLR_e3/(MLR3*np.log(10)))**2))

mass_err3 = mass_uerr3

