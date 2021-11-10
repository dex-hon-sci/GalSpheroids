#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 17:47:30 2021

@author: dexter
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAnalysis

import matplotlib.pyplot as plt
import numpy as np

## style setting
import matplotlib.style
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

plt.style.use('classic')

matplotlib.rcParams['axes.unicode_minus'] = True

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0

markers = ["E0","E1","E2","E3","EAS","EABS","EBS","SA0","SAB0", "SB0", "SAa", 
           "SABa", "SBa", "SAb","SABb", "SBb", "SAc","SABc", "SBb","SAc","SABc",
           "SBc","SAd","SABd","SBd","SAm","SABm","SBm"]

markers_group = ["E","ES","S0","Sa","Sb","Sc","Sd","Sm"]



#read the data
mass0 = np.linspace(2e8,0.5e13,2000)
Dist0 = np.linspace(0,120,2000)

D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_3.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_3.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_3.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_3.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_3.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_3.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_3.txt",
    dtype = 'str')

D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_3.txt")

# The apparant magnitude of the galaxy in g- and i-band.
mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

# The distance (final decision), lower limit and upper limit
D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

#extended disk ellipicity
elle1 = D0_Bin1_table[:,-1] 
elle2 = D0_Bin2_table[:,-1] 
elle3 = D0_Bin3_table[:,-1] 

# calculate the ellipticity of the Sersic 2D fit
b_a_1 = D0_Bin1_table[:,34]
b_a_2 = D0_Bin2_table[:,34]
b_a_3 = D0_Bin3_table[:,34]

# Calculate the Radius in equivalent axis given by the SDSS ATLAS catalogue
# Accuracy questionable
Sersic2D_50rad_1 = D0_Bin1_table[:,33]*np.sqrt(1-(1-(b_a_1)**2))
Sersic2D_50rad_2 = D0_Bin2_table[:,33]*np.sqrt(1-(1-(b_a_2)**2))
Sersic2D_50rad_3 = D0_Bin3_table[:,33]*np.sqrt(1-(1-(b_a_3)**2))

Sersic2D_50rad_1 = D0_Bin1_table[:,33]
Sersic2D_50rad_2 = D0_Bin2_table[:,33]
Sersic2D_50rad_3 = D0_Bin3_table[:,33]

#Get the name of the galaxies
name_D1 = D0_Bin1_table_n[:,0]
name_D2 = D0_Bin2_table_n[:,0]
name_D3 = D0_Bin3_table_n[:,0]

#Get the morphology of the galaxies in RC3
morph1 = D0_Bin1_table_n[:,17]
morph2 = D0_Bin2_table_n[:,17]
morph3 = D0_Bin3_table_n[:,17]

# Get the new morphology given by us 
morph1_new = D0_Bin1_table_n[:,-2]
morph2_new = D0_Bin2_table_n[:,-2]
morph3_new = D0_Bin3_table_n[:,-2]

############### reading result files###############
name1 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
name2 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])

Re_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

Sersic_n1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 2) #get n
Sersic_n2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 2) #get n
Sersic_n3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 2) #get n

mu_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 0) #get n
mu_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 0) #get n
mu_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 0) #get n

core_sersic_mu_p1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["CoreBulge"], 0)
core_sersic_mu_p2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["CoreBulge"], 0)
core_sersic_mu_p3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["CoreBulge"], 0)

#Get Re major axis
Re_1_major = SRead.grab_parameter("F_Gal_bundle_major_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2_major = SRead.grab_parameter("F_Gal_bundle_major_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3_major = SRead.grab_parameter("F_Gal_bundle_major_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale

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

################################
#Calculate mass with K-correction
K_table1 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1V_Kcorr.dat")
K_table2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2V_Kcorr.dat")
K_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr.dat")

K_table1_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1V_Kcorr.dat", 
    dtype='str')
K_table2_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2V_Kcorr.dat",
    dtype='str')
K_table3_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr.dat", 
    dtype='str')

K_name1, K_name2, K_name3 = K_table1_n[:,4], K_table2_n[:,4], K_table3_n[:,4]

mag_g1_kcorr, mag_i1_kcorr = K_table1[:,19], K_table1[:,18]
mag_g2_kcorr, mag_i2_kcorr = K_table2[:,19], K_table2[:,18]
mag_g3_kcorr, mag_i3_kcorr = K_table3[:,19], K_table3[:,18]

ML_select1_IP13_K = SPlot.MLRelationIband(mag_g1_kcorr,mag_i1_kcorr).Into13_MassRatio
ML_select1_R15BC_K = SPlot.MLRelationIband(mag_g1_kcorr,mag_i1_kcorr).Roediger15BC03_MassRatio
ML_select1_Z09_K = SPlot.MLRelationIband(mag_g1_kcorr,mag_i1_kcorr).Zibetti09_MassRatio
ML_select1_T11_K = SPlot.MLRelationIband(mag_g1_kcorr,mag_i1_kcorr).Taylor11_MassRatio

M1_K = SPlot.MassCalculation(sph_mag1, D1, 4.53,mag_g1_kcorr,mag_i1_kcorr)

E1_IP13_K = M1_K.cal_Mass(ML_select1_IP13_K)
E1_R15BC_K = M1_K.cal_Mass(ML_select1_R15BC_K)
E1_Z09_K = M1_K.cal_Mass(ML_select1_Z09_K)
E1_T11_K = M1_K.cal_Mass(ML_select1_T11_K)


ML_select2_IP13_K = SPlot.MLRelationIband(mag_g2_kcorr,mag_i2_kcorr).Into13_MassRatio
ML_select2_R15BC_K = SPlot.MLRelationIband(mag_g2_kcorr,mag_i2_kcorr).Roediger15BC03_MassRatio
ML_select2_Z09_K = SPlot.MLRelationIband(mag_g2_kcorr,mag_i2_kcorr).Zibetti09_MassRatio
ML_select2_T11_K = SPlot.MLRelationIband(mag_g2_kcorr,mag_i2_kcorr).Taylor11_MassRatio

M2_K = SPlot.MassCalculation(sph_mag2, D2, 4.53,mag_g2_kcorr,mag_i2_kcorr)

E2_IP13_K = M2_K.cal_Mass(ML_select2_IP13_K)
E2_R15BC_K = M2_K.cal_Mass(ML_select2_R15BC_K)
E2_Z09_K = M2_K.cal_Mass(ML_select2_Z09_K)
E2_T11_K = M2_K.cal_Mass(ML_select2_T11_K)


ML_select3_IP13_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Into13_MassRatio
ML_select3_R15BC_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Roediger15BC03_MassRatio
ML_select3_Z09_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Zibetti09_MassRatio
ML_select3_T11_K = SPlot.MLRelationIband(mag_g3_kcorr,mag_i3_kcorr).Taylor11_MassRatio

M3_K = SPlot.MassCalculation(sph_mag3, D3, 4.53, mag_g3_kcorr,mag_i3_kcorr)

E3_IP13_K = M3_K.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K = M3_K.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K = M3_K.cal_Mass(ML_select3_Z09_K)
E3_T11_K = M3_K.cal_Mass(ML_select3_T11_K)

#####################

#stich the Bins together
E_R15BC_K = np.concatenate((E1_R15BC_K,E2_R15BC_K, E3_R15BC_K))
E_T11_K = np.concatenate((E1_T11_K,E2_T11_K, E3_T11_K))

################################
#calculate the mass error

MLR1 = ML_select1_R15BC_K
MLR_e1 = 10**0.1    

MLR2 = ML_select2_R15BC_K
MLR_e2 = 10**0.1

MLR3 = ML_select3_R15BC_K
MLR_e3 = 10**0.1

mag_e = 0.3 #magnitude error

mass_uerr1 = np.sqrt(((mag_e/2.5)**2)+((2*D1_uerr/(D1*np.log(10)))**2)+((MLR_e1/(MLR1*np.log(10)))**2))
mass_uerr2 = np.sqrt(((mag_e/2.5)**2)+((2*D2_uerr/(D2*np.log(10)))**2)+((MLR_e2/(MLR2*np.log(10)))**2))
mass_uerr3 = np.sqrt(((mag_e/2.5)**2)+((2*D3_uerr/(D3*np.log(10)))**2)+((MLR_e3/(MLR3*np.log(10)))**2))

mass_err1 = mass_uerr1
mass_err2 = mass_uerr2
mass_err3 = mass_uerr3

##########################################
#Host Galaxy info from 1Sersic
name1_1Sersic = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin1V_cpt")
name2_1Sersic = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin2V_cpt")
name3_1Sersic = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin3V_cpt")

total_mag1_1Sersic = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin1V_cpt")
total_mag2_1Sersic = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin2V_cpt")
total_mag3_1Sersic = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_Bin3V_cpt")

sph_mag1_1Sersic = SRead.grab_mag("F_Gal_bundle_1Sersic_equvi_Bin1V_cpt", ["Bulge"])
sph_mag2_1Sersic = SRead.grab_mag("F_Gal_bundle_1Sersic_equvi_Bin2V_cpt", ["Bulge"])
sph_mag3_1Sersic = SRead.grab_mag("F_Gal_bundle_1Sersic_equvi_Bin3V_cpt", ["Bulge"])

total_mag1_1Sersic = SRead.grab_total_mag("F_Gal_bundle_1Sersic_equvi_Bin1V_cpt")
total_mag2_1Sersic = SRead.grab_total_mag("F_Gal_bundle_1Sersic_equvi_Bin2V_cpt")
total_mag3_1Sersic = SRead.grab_total_mag("F_Gal_bundle_1Sersic_equvi_Bin3V_cpt")

Re_1_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin1V_cpt", ["Bulge"], 1) #get Re
Re_2_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin2V_cpt", ["Bulge"], 1) #get Re
Re_3_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin3V_cpt", ["Bulge"], 1) #get Re

Sersic_n1_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin1V_cpt", ["Bulge"], 2) #get n
Sersic_n2_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin2V_cpt", ["Bulge"], 2) #get n
Sersic_n3_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin3V_cpt", ["Bulge"], 2) #get n

Re_1_kpc_1Sersic = Re_1_1Sersic* scale1
Re_2_kpc_1Sersic = Re_2_1Sersic* scale2
Re_3_kpc_1Sersic = Re_3_1Sersic* scale3

Re_1_kpc_lerr_1Sersic, Re_1_kpc_uerr_1Sersic = abs(Re_1_kpc_1Sersic - Re_1_1Sersic* scale1_lerr) , abs(Re_1_1Sersic* scale1_uerr - Re_1_kpc_1Sersic)
Re_2_kpc_lerr_1Sersic, Re_2_kpc_uerr_1Sersic = abs(Re_2_1Sersic* scale2_lerr - Re_2_kpc) , abs(Re_2_1Sersic* scale2_uerr - Re_2_kpc_1Sersic)
Re_3_kpc_lerr_1Sersic, Re_3_kpc_uerr_1Sersic = abs(Re_3_1Sersic* scale3_lerr - Re_3_kpc), abs(Re_3_1Sersic* scale3_uerr - Re_3_kpc_1Sersic)

Re_1_kpc_err_1Sersic =[Re_1_kpc_lerr_1Sersic, Re_1_kpc_uerr_1Sersic]
Re_2_kpc_err_1Sersic =[Re_2_kpc_lerr_1Sersic, Re_2_kpc_uerr_1Sersic]
Re_3_kpc_err_1Sersic =[Re_3_kpc_lerr_1Sersic, Re_3_kpc_uerr_1Sersic]


# Galaxy masses (from 1Sersic)
M1_K_1Sersic = SPlot.MassCalculation(total_mag1_1Sersic, D1, 4.53,mag_g1_kcorr,mag_i1_kcorr)

E1_IP13_K_1Sersic = M1_K_1Sersic.cal_Mass(ML_select1_IP13_K)
E1_R15BC_K_1Sersic = M1_K_1Sersic.cal_Mass(ML_select1_R15BC_K)
E1_Z09_K_1Sersic = M1_K_1Sersic.cal_Mass(ML_select1_Z09_K)
E1_T11_K_1Sersic = M1_K_1Sersic.cal_Mass(ML_select1_T11_K)

M2_K_1Sersic = SPlot.MassCalculation(total_mag2_1Sersic, D2, 4.53,mag_g2_kcorr,mag_i2_kcorr)

E2_IP13_K_1Sersic = M2_K_1Sersic.cal_Mass(ML_select2_IP13_K)
E2_R15BC_K_1Sersic = M2_K_1Sersic.cal_Mass(ML_select2_R15BC_K)
E2_Z09_K_1Sersic = M2_K_1Sersic.cal_Mass(ML_select2_Z09_K)
E2_T11_K_1Sersic = M2_K_1Sersic.cal_Mass(ML_select2_T11_K)

M3_K_1Sersic = SPlot.MassCalculation(total_mag3_1Sersic, D3, 4.53,mag_g3_kcorr,mag_i3_kcorr)

E3_IP13_K_1Sersic = M3_K_1Sersic.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K_1Sersic = M3_K_1Sersic.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K_1Sersic = M3_K_1Sersic.cal_Mass(ML_select3_Z09_K)
E3_T11_K_1Sersic = M3_K_1Sersic.cal_Mass(ML_select3_T11_K)

######################################
# The final size-mass input

mass1, mass2, mass3 = E1_R15BC_K, E2_R15BC_K, E3_R15BC_K
mass1_gal, mass2_gal, mass3_gal = E1_R15BC_K_1Sersic, E2_R15BC_K_1Sersic, E3_R15BC_K_1Sersic

R1,R2,R3 = Re_1_kpc, Re_2_kpc, Re_3_kpc
R1_gal, R2_gal, R3_gal = Re_1_kpc_1Sersic, Re_2_kpc_1Sersic, Re_3_kpc_1Sersic


######################################
# compare Re and Rmax

Re_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_V_cpt", ["Bulge"], 1) #get Re
Re_1Sersic_Bin1 = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin1V_cpt", ["Bulge"], 1)
Re_1Sersic_Bin2 = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin2V_cpt", ["Bulge"], 1)
Re_1Sersic_Bin3 = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_Bin3V_cpt", ["Bulge"], 1)

name = SRead.grab_name("F_Gal_bundle_1Sersic_equvi_V_cpt") #get Re
D, D_lerr, D_uerr = D0_all_table[:,29], D0_all_table[:,30], D0_all_table[:,31]
scale = D* ars

total_mag_1Sersic = SRead.grab_total_mag("F_Gal_bundle_1Sersic_equvi_V_cpt")
total_mag= SRead.grab_total_mag("F_Gal_bundle_equvi_V_cpt")


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

l,q =0,0

print("============name/Re/Rmax/mag1Sersic/magMulti===========")
for i in range(len(Rmax)):
    if Re_1Sersic[i]/Rmax[i] < 1.0:
        l = l+1
    if Re_1Sersic[i]/Rmax[i] < 0.:
        q= q+1 
        
    print(name[i],Re_1Sersic[i],Rmax[i], total_mag_1Sersic[i], total_mag[i])
print("l", l , "q",q)

#########plot transition##############
import matplotlib.patches as mpatches

def count_gal(x,y,x1,x2,y1,y2):
    # quick function to count how many sph live within x1,x2,y1,y2 boundary
    # x ,y are 1D arrays
    count =0
    for i in range(len(x)):
        if x[i] < x2 and x[i] >x1 and y[i] < y2 and y[i] > y1:
            count = count +1 
        else:
            pass
    print('count',count)

#count_gal(E_R15BC_K,Re_kpc_mine, 1e10,2.5e11,0.4,6.0)

def add_arrow(A,x0,y0,x1,y1):
    N = len([x0])
    print('Num',N)
    if N > 0:
        i=0
        for i in range(N):
            print(i)
            #A.arrow(x0[i], y0[i], x1[i], y1[i],
            #        head_width=0.05, head_length=0.1, fc='k', ec='k')
            arrow = mpatches.FancyArrowPatch((x0[i],y0[i]), ((x1[i]),(y1[i])),
                                         lw = 2,
                                         mutation_scale=30, arrowstyle="-")
            A.add_patch(arrow)
    
def add_lines(A,x0,y0,x1,y1):
    N = len([x0])

    i=0
    for i in range(N):
        A.plot([x0[i],y0[i]], [x1[i],y1[i]],'--',lw=3)


def plot_sizemass_trans(A, mass_old, R_old, mass_new, R_new,
                        label_old = None, label_new=None):

    #Bin1
    A.scatter(mass_new, R_new,marker='o',c='#a5200b',label=r"$\rm Hidden~c,Sph$", 
              s =360, alpha=0.8)
    
    #A.scatter(E1_R15BC_K_SE_prof, Re_1_kpc,marker='x',c='#a5200b',label='Bin1', 
     #         s =70, alpha=0.7)
    
    A.scatter(mass_old, R_old,marker='P',c='g',label=r"$\rm Host~galaxies$", 
              s =360, alpha=0.8)
    
    
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox


def add_morph_marker(A,x,y):
    E0_img = mpimg.imread('./morph/E0.png')    

    imagebox = OffsetImage(E0_img, zoom=0.1)
                           

    N = len(x)
    for i in range(N):
        ab = AnnotationBbox(imagebox, (x[i], y[i]), 
                            bboxprops =dict(edgecolor='None'))

        A.add_artist(ab)

xlim_mo = [4e9,2e12]
ylim_mo = [0.18,600]
#plot morphology based selection
text_location=[0.8e10,24]
delta_text = 8


def plot_sizemass_trans_2plots(xo_0=None,yo_0=None,xn_0=None,yn_0=None, name0=None,
                               xo_1=None,yo_1=None,xn_1=None,yn_1=None, name1=None):
    fig = plt.figure(figsize=(8.2, 11.6))
    gs = gridspec.GridSpec(ncols=1, nrows=2,
                               hspace=0, wspace=0.0) 

    axt0 = plt.subplot(gs[0])
    #add_lines(axt0,xo_0,yo_0,xn_0,yn_0)

    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      "", 
                                                      alpha0=0,AX=axt0)

    med0_xo, med0_yo = np.median(xo_0),np.median(yo_0)
    med0_xn, med0_yn = np.median(xn_0),np.median(yn_0)
    
    arrow0 = mpatches.FancyArrowPatch((med0_xn,med0_yn), (med0_xo,med0_yo),
                                 mutation_scale=100,alpha=0.5,color="b",ec='k')
    
    # Find the outliers
    outlier_S0 = SSort.selection_generic(E_host_Bin_S0, R_host_Bin_S0, np.repeat(30.0,len(R_host_Bin_S0)),
                                     direction="high", axis="y")
    outlier_S = SSort.selection_generic(E_host_Bin_S, R_host_Bin_S, np.repeat(30.0,len(R_host_Bin_S)),
                                     direction="high", axis="y")
    
    
    axt0.add_patch(arrow0)
    
    plot_sizemass_trans(axt0,xo_0,yo_0,xn_0,yn_0)
    
    # circle all the outliers
    axt0.scatter(outlier_S0["bag_x"],outlier_S0["bag_y"],
                 facecolors='none', edgecolors='r', s = 650)

    #SPlot.ShowcaseIndi.show_name(xo_0,yo_0, name0, A=axt0,size=16)

    axt0.text(text_location[0],text_location[1],r"$\rm S0~hosts$",fontsize=32,color="k")    
    
    axt0.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt0.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt0.set_xscale( 'log' )
    axt0.set_yscale( 'log' )
    
    axt0.set_ylabel(r"$R_\mathrm{e,equ}~\rm(kpc)$", fontsize=16)
    axt0.set_xlabel(r"$M_{*} / \rm M_{\odot}(RC15)$", fontsize=16)
    
    
    #axt0.legend(loc=4)
    #axt0.grid(True)

    axt1 = plt.subplot(gs[1],sharey=axt0)


    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt1)

    #SPlot.ShowcaseIndi.show_name(xo_1,yo_1, name1, A=axt1,size=16)
    
    med1_xo, med1_yo = np.median(xo_1),np.median(yo_1)
    med1_xn, med1_yn = np.median(xn_1),np.median(yn_1)
    
    arrow1 = mpatches.FancyArrowPatch((med1_xn,med1_yn), (med1_xo,med1_yo),
                                 mutation_scale=100,alpha=0.5,color="b",ec='k')
    axt1.add_patch(arrow1)
   # add_arrow(axt1,med_xo,med_yo,med_xn,med_yn)

    plot_sizemass_trans(axt1,xo_1,yo_1,xn_1,yn_1)
    
    axt1.scatter(outlier_S["bag_x"],outlier_S["bag_y"],
                 facecolors='none', edgecolors='r', s = 650)

    axt1.text(text_location[0],text_location[1],r"$\rm S~hosts$",fontsize=32,color="k") 
    
    axt1.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt1.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt1.set_xscale( 'log' )
    axt1.set_yscale( 'log' )
    
    axt1.set_ylabel(r"$R_\mathrm{e,equ}~\rm(kpc)$", fontsize=16)    
    axt1.set_xlabel(r"$M_{*} / \rm M_{\odot}~(RC15)$", fontsize=16)
    axt1.legend(loc=4)
    #axt1.grid(True)
    

    plt.show()


# Apply Barro cut, select the subsample
Bcut1 = SPlot.SelectionCut(mass1, D1).Barro13_cut()
Bcut2 = SPlot.SelectionCut(mass2, D2).Barro13_cut()
Bcut3 = SPlot.SelectionCut(mass3, D3).Barro13_cut()

S1 = SSort.selection_generic(mass1, Re_1_kpc, Bcut1)
S2 = SSort.selection_generic(mass2, Re_2_kpc, Bcut2)
S3 = SSort.selection_generic(mass3, Re_3_kpc, Bcut3)

#select the morphology of the host that passes the selection, seperated by Bins
morph1_sub = SSort.cherry_pick(S1['index'], morph1_new)
morph2_sub = SSort.cherry_pick(S2['index'], morph2_new)
morph3_sub = SSort.cherry_pick(S3['index'], morph3_new)

# Seperate by morphology
# Bin1 size and mass
# Find the indices
index_Bin1_E = SSort.morph_str_selection(S1['index'], morph1_sub)["E"]
index_Bin1_S0 = SSort.morph_str_selection(S1['index'], morph1_sub)["S0"]
index_Bin1_S = SSort.morph_str_selection(S1['index'], morph1_sub)["S"]

name_Bin1_E = SSort.cherry_pick(index_Bin1_E, name1)
name_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, name1)
name_Bin1_S = SSort.cherry_pick(index_Bin1_S, name1)

Rmax1_E = SSort.cherry_pick(index_Bin1_E, Rmax1)
Rmax1_S0 = SSort.cherry_pick(index_Bin1_S0, Rmax1)
Rmax1_S = SSort.cherry_pick(index_Bin1_S, Rmax1)

morph_Bin1_E = SSort.cherry_pick(index_Bin1_E, morph1_new)
morph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, morph1_new)
morph_Bin1_S = SSort.cherry_pick(index_Bin1_S, morph1_new)

E_host_Bin1_E = SSort.cherry_pick(index_Bin1_E, mass1_gal)
R_host_Bin1_E = SSort.cherry_pick(index_Bin1_E, R1_gal)
E_sph_Bin1_E = SSort.cherry_pick(index_Bin1_E, mass1)
R_sph_Bin1_E = SSort.cherry_pick(index_Bin1_E, Re_1_kpc)

E_host_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, mass1_gal)
R_host_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, R1_gal)
E_sph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, mass1)
R_sph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, R1)

E_host_Bin1_S = SSort.cherry_pick(index_Bin1_S, mass1_gal)
R_host_Bin1_S = SSort.cherry_pick(index_Bin1_S, R1_gal)
E_sph_Bin1_S = SSort.cherry_pick(index_Bin1_S, mass1)
R_sph_Bin1_S = SSort.cherry_pick(index_Bin1_S, R1)
 
#Bin2 size and mass
index_Bin2_E = SSort.morph_str_selection(S2['index'], morph2_sub)["E"]
index_Bin2_S0 = SSort.morph_str_selection(S2['index'], morph2_sub)["S0"]
index_Bin2_S = SSort.morph_str_selection(S2['index'], morph2_sub)["S"]

name_Bin2_E = SSort.cherry_pick(index_Bin2_E, name2)
name_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, name2)
name_Bin2_S = SSort.cherry_pick(index_Bin2_S, name2)

Rmax2_E = SSort.cherry_pick(index_Bin2_E, Rmax2)
Rmax2_S0 = SSort.cherry_pick(index_Bin2_S0, Rmax2)
Rmax2_S = SSort.cherry_pick(index_Bin2_S, Rmax2)

morph_Bin2_E = SSort.cherry_pick(index_Bin2_E, morph2_new)
morph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, morph2_new)
morph_Bin2_S = SSort.cherry_pick(index_Bin2_S, morph2_new)

E_host_Bin2_E = SSort.cherry_pick(index_Bin2_E, mass2_gal)
R_host_Bin2_E = SSort.cherry_pick(index_Bin2_E, R2_gal)
E_sph_Bin2_E = SSort.cherry_pick(index_Bin2_E, mass2)
R_sph_Bin2_E = SSort.cherry_pick(index_Bin2_E, R2)

E_host_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, mass2_gal)
R_host_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, R2_gal)
E_sph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, mass2)
R_sph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, R2)

E_host_Bin2_S = SSort.cherry_pick(index_Bin2_S, mass2_gal)
R_host_Bin2_S = SSort.cherry_pick(index_Bin2_S, R2_gal)
E_sph_Bin2_S =  SSort.cherry_pick(index_Bin2_S, mass2)
R_sph_Bin2_S = SSort.cherry_pick(index_Bin2_S, R2)

# Bin3 size and mass
index_Bin3_E = SSort.morph_str_selection(S3['index'], morph3_sub)["E"]
index_Bin3_S0 = SSort.morph_str_selection(S3['index'], morph3_sub)["S0"]
index_Bin3_S = SSort.morph_str_selection(S3['index'], morph3_sub)["S"]

name_Bin3_E = SSort.cherry_pick(index_Bin3_E, name3)
name_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, name3)
name_Bin3_S = SSort.cherry_pick(index_Bin3_S, name3)

Rmax3_E = SSort.cherry_pick(index_Bin3_E, Rmax3)
Rmax3_S0 = SSort.cherry_pick(index_Bin3_S0, Rmax3)
Rmax3_S = SSort.cherry_pick(index_Bin3_S, Rmax3)

morph_Bin3_E = SSort.cherry_pick(index_Bin3_E, morph3_new)
morph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, morph3_new)
morph_Bin3_S = SSort.cherry_pick(index_Bin3_S, morph3_new)

E_host_Bin3_E = SSort.cherry_pick(index_Bin3_E, mass3_gal)
R_host_Bin3_E = SSort.cherry_pick(index_Bin3_E, R3_gal)
E_sph_Bin3_E = SSort.cherry_pick(index_Bin3_E, mass3)
R_sph_Bin3_E = SSort.cherry_pick(index_Bin3_E, R3)

E_host_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, mass3_gal)
R_host_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, R3_gal)
E_sph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, mass3)
R_sph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, R3)

E_host_Bin3_S = SSort.cherry_pick(index_Bin3_S, mass3_gal)
R_host_Bin3_S = SSort.cherry_pick(index_Bin3_S, R3_gal)
E_sph_Bin3_S =  SSort.cherry_pick(index_Bin3_S, mass3)
R_sph_Bin3_S = SSort.cherry_pick(index_Bin3_S, R3)


# combine Bin1-3
index_Bin_E = index_Bin1_E+index_Bin2_E+index_Bin3_E
index_Bin_S0 = index_Bin1_S0+index_Bin2_S0+index_Bin3_S0
index_Bin_S = index_Bin1_S+index_Bin2_S+index_Bin3_S

name_E = name_Bin1_E + name_Bin2_E + name_Bin3_E
name_SO = name_Bin1_S0 + name_Bin2_S0 + name_Bin3_S0
name_S = name_Bin1_S + name_Bin2_S + name_Bin3_S

Rmax_E = Rmax1_E + Rmax2_E + Rmax3_E 
Rmax_S0 = Rmax1_S0 + Rmax2_S0 + Rmax3_S0
Rmax_S = Rmax1_S + Rmax2_S + Rmax3_S 

E_host_Bin_E = E_host_Bin1_E + E_host_Bin2_E + E_host_Bin3_E
E_host_Bin_S0 = E_host_Bin1_S0 + E_host_Bin2_S0 + E_host_Bin3_S0
E_host_Bin_S = E_host_Bin1_S + E_host_Bin2_S + E_host_Bin3_S

R_host_Bin_E = R_host_Bin1_E + R_host_Bin2_E + R_host_Bin3_E
R_host_Bin_S0 = R_host_Bin1_S0 + R_host_Bin2_S0 + R_host_Bin3_S0
R_host_Bin_S = R_host_Bin1_S + R_host_Bin2_S + R_host_Bin3_S

E_sph_Bin_E = E_sph_Bin1_E + E_sph_Bin2_E + E_sph_Bin3_E
E_sph_Bin_S0 = E_sph_Bin1_S0 + E_sph_Bin2_S0 + E_sph_Bin3_S0
E_sph_Bin_S = E_sph_Bin1_S + E_sph_Bin2_S + E_sph_Bin3_S

R_sph_Bin_E = R_sph_Bin1_E + R_sph_Bin2_E + R_sph_Bin3_E
R_sph_Bin_S0 = R_sph_Bin1_S0 + R_sph_Bin2_S0 + R_sph_Bin3_S0
R_sph_Bin_S = R_sph_Bin1_S + R_sph_Bin2_S + R_sph_Bin3_S

morph_Bin_E = morph_Bin1_E + morph_Bin2_E + morph_Bin3_E
morph_Bin_S0 = morph_Bin1_S0 + morph_Bin2_S0 + morph_Bin3_S0
morph_Bin_S = morph_Bin1_S + morph_Bin2_S + morph_Bin3_S

# plot the comaprison
plot_sizemass_trans_2plots(xo_0 = E_host_Bin_S0, yo_0 = R_host_Bin_S0, xn_0 = E_sph_Bin_S0, yn_0 = R_sph_Bin_S0, name0 = morph_Bin_S0, 
                           xo_1 = E_host_Bin_S, yo_1 = R_host_Bin_S, xn_1 = E_sph_Bin_S, yn_1 = R_sph_Bin_S, name1 = morph_Bin_S)
                           
outlier_S0 = SSort.selection_generic(E_host_Bin_S0, R_host_Bin_S0, np.repeat(30.0,len(R_host_Bin_S0)),
                                     direction="high", axis="y")

outlier_S = SSort.selection_generic(E_host_Bin_S, R_host_Bin_S, np.repeat(30.0,len(R_host_Bin_S)),
                                     direction="high", axis="y")

print("S0", len(E_sph_Bin1_S0)+len(E_sph_Bin2_S0)+len(E_sph_Bin3_S0))
print("S", len(E_sph_Bin1_S)+len(E_sph_Bin2_S)+len(E_sph_Bin3_S))


print("outlier_S0",outlier_S0)
print("outlier_S",outlier_S)


print(name_SO[7])
print(name_S[0])
print(name_S[10])

#reserve
#
#fig, ax = plt.subplots()        
##SPlot.SelectionCut(mass0,Dist0).plot_cut()
#SPlot.ShowcaseIndi.Mass_Re_plot(Benzanson2015_mass,Benzanson2015_Re, colour = '#4277cf',
 #                               legend='Benzanson et al. 2015',ms=8,alpha0 = 0.2, lw=3, marker = "^")

#SPlot.ShowcaseIndi.Mass_Re_plot(Zahid2015_mass,Zahid2015_Re, colour = '#c17dd5',
#                                legend='Zahid et al. 2015',ms=8,alpha0 = 0.2, lw=3, marker = "^")

#plot_dexter_sample_all()
#
#plt.ylim(ylim[0],ylim[1])
#plt.xlim(xlim[0],xlim[1])
#plt.show()
#
#
#fig, ax = plt.subplots()        
#plt.plot(E1_R15BC,mag_g1-mag_i1,'o')
#plt.plot(E2_R15BC,mag_g2-mag_i2,'o')
#plt.plot(E3_R15BC,mag_g3-mag_i3,'o')
#
#plt.show()


# plot mag_mag 1Sersic vs multi

# plot Re_Re

#plot the size mass relation transision

ReRmax = SRead.read_table("/home/dexter/result/stat/completeness/gal_ReRmaxMag2.txt")
ReRmax_n = SRead.read_table("/home/dexter/result/stat/completeness/gal_ReRmaxMag2.txt",dtype="str")

mag_one_sersic = ReRmax[:,3]
mag_multi_cpt= ReRmax[:,4]
name_ReRmax = ReRmax_n[:,0]

SPlot.ShowcaseCompare2.plot_compare_generic(mag_one_sersic, 
                                            mag_multi_cpt,para_name="mag", 
                                            name=name_ReRmax,label=["{1-S{\'e}rsic}","{multi-cpt}"])

# Find the insider in Re/Rmax
ReRmax_insiders = SSort.selection_generic(
    np.repeat(0,len(Re_1Sersic)), Re_1Sersic/Rmax, np.repeat(
        3.0,len(Re_1Sersic)), direction="low", axis="y")

# Find the outliers in Re/Rmax
ReRmax_outliers = SSort.selection_generic(
    np.repeat(0,len(Re_1Sersic)), Re_1Sersic/Rmax, np.repeat(
        3.0,len(Re_1Sersic)), direction="high", axis="y")

print((ReRmax_outliers['bag_y']))
print(len(ReRmax_insiders['bag_y']))

# Seperate Re and Rmax by morphology
Re1_spc = SSort.morph_str_selection(Re_1Sersic_Bin1,morph1_new)
Re2_spc = SSort.morph_str_selection(Re_1Sersic_Bin2,morph2_new)
Re3_spc = SSort.morph_str_selection(Re_1Sersic_Bin3,morph3_new)

Rmax1_spc = SSort.morph_str_selection(Rmax1,morph1_new)
Rmax2_spc = SSort.morph_str_selection(Rmax2,morph2_new)
Rmax3_spc = SSort.morph_str_selection(Rmax3,morph3_new)

name1_spc = SSort.morph_str_selection(name1,morph1_new)
name2_spc = SSort.morph_str_selection(name2,morph2_new)
name3_spc = SSort.morph_str_selection(name3,morph3_new)


name_E = name1_spc["E"]+name2_spc["E"]+name3_spc["E"]
name_S0 = name1_spc["S0"]+name2_spc["S0"]+name3_spc["S0"]
name_S = name1_spc["S"]+name2_spc["S"]+name3_spc["S"]

Re_E = np.concatenate((Re1_spc["E"],Re2_spc["E"],Re3_spc["E"]))
Re_S0 = np.concatenate((Re1_spc["S0"], Re2_spc["S0"], Re3_spc["S0"]))
Re_S = np.concatenate((Re1_spc["S"], Re2_spc["S"], Re3_spc["S"]))

Rmax_E = np.concatenate((Rmax1_spc["E"],Rmax2_spc["E"],Rmax3_spc["E"]))
Rmax_S0 = np.concatenate((Rmax1_spc["S0"], Rmax2_spc["S0"], Rmax3_spc["S0"]))
Rmax_S = np.concatenate((Rmax1_spc["S"], Rmax2_spc["S"], Rmax3_spc["S"]))


print(np.array(Re_E)/np.array(Rmax_E))
print(np.array(Re_S0)/np.array(Rmax_S0))
print(np.array(Re_S)/np.array(Rmax_S))

def plot_ReRmax_hist():
    fig = plt.figure(figsize=(6.4, 4.8))

    # Divided by morphology
    # 20,20,160
    
    space = np.linspace(0.0,50.0,num=400)
    
    plt.hist(np.array(Re_S)/np.array(Rmax_S),space, 
                                facecolor='blue',
                                alpha=0.6, histtype ='stepfilled',label="S")
 
    plt.hist(np.array(Re_S0)/np.array(Rmax_S0),space, 
                                facecolor='orange',
                                alpha=0.6, histtype ='stepfilled',label="S0")
    plt.hist(np.array(Re_E)/np.array(Rmax_E),space, 
                                facecolor='red',
                                alpha=0.6, histtype ='stepfilled',label="E")
    
    
    #plt.hist(np.array(R_sph_Bin_S)/np.array(Rmax_S),650, 
    #                            facecolor='orange',
    #                            alpha=0.7, histtype ='stepfilled',label="Our galaxies")
    #plt.hist(ReRmax_outliers['bag_y'],662, 
    #                            facecolor='red',
    #                            alpha=0.7, histtype ='stepfilled',label="Outliers")    

    #print(Re_1Sersic/Rmax)
    plt.xlabel(r"$R_\mathrm{e} / \rm R_{max}$", fontsize=16)
    plt.ylabel(r"$\rm Number~of~galaxies$", fontsize=16)    
    
    plt.xlim(0,3)
    plt.ylim(0.5,30)
    
    plt.yscale( 'log' )  
    #plt.xscale( 'log' )  
    plt.legend()

    
    plt.tight_layout()


plot_ReRmax_hist()
