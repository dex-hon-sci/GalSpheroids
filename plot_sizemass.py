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
##
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

###

D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V.txt",
    dtype = 'str')

D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

# calculate the ellipticity of the Sersic 2D fit
b_a_1 = D0_Bin1_table[:,34]
b_a_2 = D0_Bin2_table[:,34]
b_a_3 = D0_Bin3_table[:,34]

# calculate the Radius in equivalent axis 
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

#Get the morphology of the galaxies
morph1 = D0_Bin1_table_n[:,17]
morph2 = D0_Bin2_table_n[:,17]
morph3 = D0_Bin3_table_n[:,17]

morph1_new = D0_Bin1_table_n[:,-1]
morph2_new = D0_Bin2_table_n[:,-1]
morph3_new = D0_Bin3_table_n[:,-1]



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

print("core_sersic_mu_p1",core_sersic_mu_p1)
print("core_sersic_mu_p2",core_sersic_mu_p2)
print("core_sersic_mu_p3",core_sersic_mu_p3)

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


Sersic2D_50rad_1_kpc = Sersic2D_50rad_1*scale1
Sersic2D_50rad_2_kpc = Sersic2D_50rad_2*scale2
Sersic2D_50rad_3_kpc = Sersic2D_50rad_3*scale3

for i in range(len(name1)):
    print(name1[i],morph1[i], morph1_new[i])
    
for i in range(len(name1)):
    print(name1[i],Sersic2D_50rad_1_kpc[i],Sersic2D_50rad_1[i],scale1[i])



############# calculating spheroid mass

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

############# calculating spheroid absoulte magnitude
Abs_sph_mag1 = M1.cal_abs_mag()
Abs_sph_mag2 = M2.cal_abs_mag()
Abs_sph_mag3 = M3.cal_abs_mag()

print("Abs_sph_mag1",Abs_sph_mag1)

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
E_R15BC_K = np.concatenate((E1_R15BC_K,E2_R15BC_K, E3_R15BC_K))
E_T11_K = np.concatenate((E1_T11_K,E2_T11_K, E3_T11_K))

print("median",np.median(E_R15BC_K)-np.median(E_T11_K))
print(np.median(E_R15BC_K-E_T11_K))
########################
def plot_mass_compare_Kcor():
    SPlot.ShowcaseCompare2.plot_compare_generic(E1_R15BC_K, E1_R15BC, 
                                            sub=False , div=True,
                                            para_name="Spheroid Mass", colour="blue", 
                                            name=K_name1, label=["Kcorr",""])

    SPlot.ShowcaseCompare2.plot_compare_generic(E2_R15BC_K, E2_R15BC, 
                                            sub=False , div=True,
                                            para_name="Spheroid Mass", colour="blue", 
                                            name=K_name2, label=["Kcorr",""])
    SPlot.ShowcaseCompare2.plot_compare_generic(E3_R15BC_K, E3_R15BC, 
                                                sub=False , div=True,
                                                para_name="Spheroid Mass", colour="blue", 
                                                name=K_name3, label=["Kcorr",""])
    plt.show()

############# calculating total stellar mass, both SDSS mag and profiler mag
#Bin1
#SDSS M/L
M1_K_SE_SDSS = SPlot.MassCalculation(mag_i1, D1, 4.53,mag_g1_kcorr,mag_i1_kcorr)
#profiler M/L
M1_K_SE_prof = SPlot.MassCalculation(total_mag1, D1, 4.53,mag_g1_kcorr,mag_i1_kcorr)

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
M2_K_SE_SDSS = SPlot.MassCalculation(mag_i2, D2, 4.53,mag_g2_kcorr,mag_i2_kcorr)
#profiler M/L
M2_K_SE_prof = SPlot.MassCalculation(total_mag2, D2, 4.53,mag_g2_kcorr,mag_i2_kcorr)

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
M3_K_SE_SDSS = SPlot.MassCalculation(mag_i3, D3, 4.53,mag_g3_kcorr,mag_i3_kcorr)
#profiler M/L
M3_K_SE_prof = SPlot.MassCalculation(total_mag3, D3, 4.53,mag_g3_kcorr,mag_i3_kcorr)

E3_IP13_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_Z09_K)
E3_T11_K_SE_SDSS = M3_K_SE_SDSS.cal_Mass(ML_select3_T11_K)

E3_IP13_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_IP13_K)
E3_R15BC_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_R15BC_K)
E3_Z09_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_Z09_K)
E3_T11_K_SE_prof = M3_K_SE_prof.cal_Mass(ML_select3_T11_K)


print ("T11 average M/L ratio", np.average(ML_select1_T11_K),
       np.average(ML_select2_T11_K),np.average(ML_select3_T11_K))
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
# Define ploting in RDJ15
def plot_dexter_sample_Bin():
    
    
    SPlot.ShowcaseIndi.Mass_Re_plot(E1_R15BC_K, Re_1_kpc, yerr = Re_1_kpc_err, 
                                xerr = mass_err1*E1_R15BC_K,
                                colour = '#a5200b',
                                name=name1,legend='Bin1',ms=8,alpha0 = 0.2, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E2_R15BC_K, Re_2_kpc, yerr = Re_2_kpc_err, 
                                xerr = mass_err2*E2_R15BC_K,
                                colour = '#0b5786',
                                name=name2,legend='Bin2',ms=8,alpha0 = 0.2, lw=3)
    SPlot.ShowcaseIndi.Mass_Re_plot(E3_R15BC_K, Re_3_kpc, yerr = Re_3_kpc_err, 
                                xerr = mass_err3*E3_R15BC_K,
                                colour='#2a3236',
                                name=name1,legend='Bin3',ms=8,alpha0 = 0.2, lw=3)
    
def plot_dexter_sample_all():
    
    SPlot.ShowcaseIndi.Mass_Re_plot(E1_R15BC_K, Re_1_kpc, yerr = Re_1_kpc_err,
                                xerr = mass_err1*E1_R15BC_K,
                                colour = '#a5200b',
                                name=name1,legend='This work',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E2_R15BC_K, Re_2_kpc, yerr = Re_2_kpc_err,
                                xerr = mass_err2*E2_R15BC_K,
                                colour = '#a5200b',
                                name=name2,legend='',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E3_R15BC_K, Re_3_kpc, yerr = Re_3_kpc_err,
                                xerr = mass_err3*E3_R15BC_K,
                                colour='#a5200b',
                                name=name1,legend='',ms=8,alpha0 = 0.4,lw=3 )
    
    
txsep1 = 40#10*5.5
txsep2 = 65#10*6.0
txsep3 = 80#10*6.5

x,y = 8.6e8,100

circle_size = 310

def plot_Barro_cut_all(AX):
    
    Bcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).Barro13_cut()
    Bcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).Barro13_cut()
    Bcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).Barro13_cut()

    S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc, Bcut1)
    S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc, Bcut2)
    S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc, Bcut3)
    
    
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

    vdWcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).vdWel14_cut()
    vdWcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).vdWel14_cut()
    vdWcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).vdWel14_cut()

    S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc_major, vdWcut1)
    S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc_major, vdWcut2)
    S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc_major, vdWcut3)
    
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
    
    vDcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).vDokkum15_cut()
    vDcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).vDokkum15_cut()
    vDcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).vDokkum15_cut()

    S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc, vDcut1)
    S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc, vDcut2)
    S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc, vDcut3)
    
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
    

    Damcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).Damjanov14_cut()
    Damcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).Damjanov14_cut()
    Damcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).Damjanov14_cut()

    S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc, Damcut1)
    S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc, Damcut2)
    S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc, Damcut3)

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
    

    Gcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).Graham15_broad_cut()
    Gcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).Graham15_broad_cut()
    Gcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).Graham15_broad_cut()

    S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc, Gcut1)
    S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc, Gcut2)
    S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc, Gcut3)

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

def plot_sizemass_SDSS_mine():
    fig1, ax1 = plt.subplots()        


    SPlot.ShowcaseIndi.Mass_Re_plot(nsa_mass_R15BC, nsa_Re, marker='x', colour = 'grey',
                                ms = 5, legend='SDSS galaxies',alpha0 = 0.2)

    #SPlot.SelectionCut(mass0,Dist0).plot_cut()
    plot_dexter_sample_Bin()
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
    A.scatter(E1_R15BC, Re_1_kpc,marker='o',c='#a5200b',label='Bin1', 
              s =70, alpha=0.7)
    A.errorbar(E1_R15BC, Re_1_kpc, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*E1_R15BC, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=0.65, marker='o')
    #Bin2
    A.scatter(E2_R15BC, Re_2_kpc,marker='o',c='#0b5786',label='Bin2', 
              s =70,alpha=0.71)
    A.errorbar(E2_R15BC, Re_2_kpc, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*E2_R15BC,ls='none',linewidth=4,
                  color = '#0b5786',
                  ecolor='#0b5786',capsize=0, alpha=0.65, marker='o')
    #Bin3
    A.scatter(E3_R15BC, Re_3_kpc,marker='o',c='#2a3236',label='Bin3', 
              s =70, alpha=0.7)
    A.errorbar(E3_R15BC, Re_3_kpc, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*E3_R15BC,ls='none',linewidth=4,
                  color = '#2a3236',
                  ecolor='#2a3236',capsize=0,
                  alpha=0.65, marker='o')   
    
    A.set_ylim(ylim[0],ylim[1])
    A.set_xlim(xlim[0],xlim[1])
    
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )    

def plot_dexter_sample_all2(A,alpha=0.65):
    #Bin1
    A.scatter(E1_R15BC_K, Re_1_kpc,marker='o',c='#a5200b',label='This work', 
              s =70, alpha=0.1)
    A.errorbar(E1_R15BC_K, Re_1_kpc, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*E1_R15BC_K, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=alpha, marker='o')
    #Bin2
    A.scatter(E2_R15BC_K, Re_2_kpc,marker='o',c='#a5200b',label='', 
              s =70,alpha=0.71)
    A.errorbar(E2_R15BC_K, Re_2_kpc, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*E2_R15BC_K,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0, alpha=alpha, marker='o')
    #Bin3
    A.scatter(E3_R15BC_K, Re_3_kpc,marker='o',c='#a5200b',label='', 
              s =70, alpha=0.7)
    A.errorbar(E3_R15BC_K, Re_3_kpc, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*E3_R15BC_K,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0,
                  alpha=alpha, marker='o')   

    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale( 'log' )
    A.set_yscale( 'log' )
    
    
    
def plot_dexter_sample_all2_T11(A,scale='log',alpha=0.65):
    #Bin1
    A.scatter(E1_T11, Re_1_kpc,marker='o',c='#a5200b',label='This work', 
              s =70, alpha=0.7)
    A.errorbar(E1_T11, Re_1_kpc, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*E1_T11, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=alpha, marker='o')
    #Bin2
    A.scatter(E2_T11, Re_2_kpc,marker='o',c='#a5200b',label='', 
              s =70,alpha=0.71)
    A.errorbar(E2_T11, Re_2_kpc, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*E2_T11,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0, alpha=alpha, marker='o')
    #Bin3
    A.scatter(E3_T11, Re_3_kpc,marker='o',c='#a5200b',label='', 
              s =70, alpha=0.7)
    A.errorbar(E3_T11, Re_3_kpc, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*E3_T11,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0,
                  alpha=alpha, marker='o')   

    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale(scale)
    A.set_yscale(scale)


def plot_dexter_sample_all2_major(A,scale='log',alpha=0.65):
    #Bin1
    A.scatter(E1_R15BC_K, Re_1_kpc_major, marker='o',c='#a5200b',label='This work', 
              s =70, alpha=0.7)
    A.errorbar(E1_R15BC_K, Re_1_kpc_major, yerr = Re_1_kpc_err, 
                  xerr = mass_err1*E1_R15BC, ls='none',linewidth=4, 
                  color ='#a5200b',
                  ecolor='#a5200b', capsize=0, alpha=alpha, marker='o')
    #Bin2
    A.scatter(E2_R15BC_K, Re_2_kpc_major,marker='o',c='#a5200b',label='', 
              s =70,alpha=0.71)
    A.errorbar(E2_R15BC_K, Re_2_kpc_major, yerr = Re_2_kpc_err, 
                  xerr = mass_err2*E2_R15BC,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0, alpha=alpha, marker='o')
    #Bin3
    A.scatter(E3_R15BC_K, Re_3_kpc_major,marker='o',c='#a5200b',label='', 
              s =70, alpha=0.7)
    A.errorbar(E3_R15BC_K, Re_3_kpc_major, yerr = Re_3_kpc_err, 
                  xerr = mass_err3*E3_R15BC,ls='none',linewidth=4,
                  color = '#a5200b',
                  ecolor='#a5200b',capsize=0,
                  alpha=alpha, marker='o')   

    
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
    
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=2, nrows=3,
                               hspace=0, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])
    plot_dexter_sample_Bin2(axs0)
    axs0.legend(loc=2)
    axs0.set_ylabel("$R_{e,sph}$ (kpc)",fontsize=16)
    #axs0.grid(True)

    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharey=axs0)
    plot_sizemass_DamCut_mine(axs1)
    axs1.legend(loc=4)
    plt.setp(axs1.get_yticklabels(), visible=False)
    #axs1.set_yticks([])
    #axs1.grid(True)
    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2])
    plot_sizemass_BarroCut_mine(axs2)
    axs2.legend(loc=4)
    axs2.set_ylabel("$R_{e,sph}$ (kpc)",fontsize=16)
    #axs2.grid(True)
    
    #plot Panel (4)
    axs3 = plt.subplot(gs[3],sharey=axs2)
    plot_sizemass_vDokkumCut_mine(axs3)
    axs3.legend(loc=4)
    #axs3.set_yticks([])
    plt.setp(axs3.get_yticklabels(), visible=False)
    #axs3.grid(True)

    #plot Panel (5)
    axs4 = plt.subplot(gs[4])
    plot_sizemass_vdWelCut_mine(axs4)
    axs4.legend(loc=4)
    axs4.set_ylabel("$ R_{e,sph}$ (kpc)", fontsize=16)
    axs4.set_xlabel(r"$\rm M_{*,sph} / M_{\odot}$", fontsize=16)
    #axs4.grid(True)
   
    #plot Panel (6)
    axs5 = plt.subplot(gs[5],sharey=axs4)
    plot_sizemass_GrahamCut_mine(axs5)
    axs5.legend(loc=4)
    plt.setp(axs5.get_yticklabels(), visible=False)
    #axs5.set_yticks([])
    axs5.set_xlabel(r"$\rm M_{*,sph} / M_{\odot}$", fontsize=16)
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
Abs_sph_mag_mine = np.concatenate((Abs_sph_mag1, Abs_sph_mag2, Abs_sph_mag3))
mu_mine = np.concatenate((mu_1, mu_2, mu_3))

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
    
    ax.legend(loc=2)
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
    
    ax.legend(loc=2)
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

print(len(E_T11),len(Re_kpc))

##fitting
popt,pcov = curve_fit(s, E_T11, Re_kpc)
popt2,pcov2 = curve_fit(s,E_T11_mine, Re_kpc_mine)
popt3,pcov3 = curve_fit(s,Savorgnan_mass_T11, Savorgnan_size_eq_kpc)
popt4,pcov4 = curve_fit(s,Davis_mass_T11, Davis_size_eq_kpc)
popt5,pcov5 = curve_fit(s,Sahu_mass_T11, Sahu_size_eq_kpc)

popt_g,pcov_g = curve_fit(ss, E_T11, Re_kpc)

print("graham_equ_para", *popt_g )
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

    ax.legend(loc=2)
    #plt.grid(True)
    
    plt.xlabel(r"$\rm M_{*,sph}$ / $M_{\odot}$",fontsize=16)
    plt.ylabel("$R_{e,sph}$ (kpc)",fontsize=16)
    
    plt.show()
    
    
plot_sizemass_z0comparison()


########################


import matplotlib.patches as mpatches

def add_arrow(A,x0,y0,x1,y1):
    N = len(x0)

    i=0
    for i in range(N):
        #A.arrow(x0[i], y0[i], x1[i], y1[i],
        #        head_width=0.05, head_length=0.1, fc='k', ec='k')
        arrow = mpatches.FancyArrowPatch((x0[i],y0[i]), ((x1[i]),(y1[i])),
                                         lw = 2,
                                         mutation_scale=30, arrowstyle="-|>")
        A.add_patch(arrow)
    


def plot_sizemass_trans(A, mass_old, R_old, mass_new, R_new,
                        label_old = None, label_new=None):

    #Bin1
    A.scatter(mass_new, R_new,marker='o',c='#a5200b',label=label_new, 
              s =120, alpha=0.7)
    
    #A.scatter(E1_R15BC_K_SE_prof, Re_1_kpc,marker='x',c='#a5200b',label='Bin1', 
     #         s =70, alpha=0.7)
    
    A.scatter(mass_old, R_old,marker='x',c='g',label=label_old, 
              s =70, alpha=0.0)
    
    
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


#########plot transition##############

xlim_mo = [4.8e9,1.24e12]
ylim_mo = [0.28,44]
#plot morphology based selection
text_location=[1.61e10,10]
delta_text = 3


def plot_sizemass_trans_3plots(xo_0=None,yo_0=None,xn_0=None,yn_0=None, name0=None,
                               xo_1=None,yo_1=None,xn_1=None,yn_1=None, name1=None,
                               xo_2=None,yo_2=None,xn_2=None,yn_2=None, name2=None,
                               xo_3=None,yo_3=None,xn_3=None,yn_3=None, name3=None,
                               xo_4=None,yo_4=None,xn_4=None,yn_4=None, name4=None,
                               xo_5=None,yo_5=None,xn_5=None,yn_5=None, name5=None,
                               xo_6=None,yo_6=None,xn_6=None,yn_6=None, name6=None,
                               xo_7=None,yo_7=None,xn_7=None,yn_7=None, name7=None,
                               xo_8=None,yo_8=None,xn_8=None,yn_8=None, name8=None):
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=3, nrows=3,
                               hspace=0, wspace=0.0) 

    axt0 = plt.subplot(gs[0])

    plot_sizemass_trans(axt0,xo_0,yo_0,xn_0,yn_0)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      "", 
                                                      alpha0=0,AX=axt0)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      "", 
                                                      alpha0=0,AX=axt0)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      "", 
                                                      alpha0=0,AX=axt0)             
    
    add_arrow(axt0,xo_0,yo_0,xn_0,yn_0)
    
    SPlot.ShowcaseIndi.show_name(xo_0,yo_0, name0, A=axt0,size=16)

    axt0.text(text_location[0],text_location[1],r"$\rm Bin 1$",fontsize=22,color="#a5200b")    
    axt0.text(text_location[0],text_location[1]-delta_text,r"$\rm E -> E$",fontsize=22,color="k")    
    
    axt0.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt0.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt0.set_xscale( 'log' )
    axt0.set_yscale( 'log' )
    
    axt0.set_ylabel(r"$\rm R_{e}$ (kpc)", fontsize=16)
    axt0.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)
    
    
    #axt0.legend(loc=4)
    #axt0.grid(True)

    axt1 = plt.subplot(gs[1],sharey=axt0)

    plot_sizemass_trans(axt1,xo_1,yo_1,xn_1,yn_1)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt1)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt1)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt1)
    add_arrow(axt1,xo_1,yo_1,xn_1,yn_1)

    SPlot.ShowcaseIndi.show_name(xo_1,yo_1, name1, A=axt1,size=16)
    
    axt1.text(text_location[0],text_location[1],r"$\rm Bin 1$",fontsize=22,color="#a5200b") 
    axt1.text(text_location[0],text_location[1]-delta_text,r"$\rm S0 -> E$",fontsize=22,color="k")    
    
    axt1.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt1.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt1.set_xscale( 'log' )
    axt1.set_yscale( 'log' )
    
    axt1.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)
    #axt1.legend(loc=4)
    #axt1.grid(True)
    plt.setp(axt1.get_yticklabels(), visible=False)
    
    axt2 = plt.subplot(gs[2],sharey=axt0)

    plot_sizemass_trans(axt2,xo_2,yo_2,xn_2,yn_2)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt2)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt2)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt2)
    add_arrow(axt2,xo_2,yo_2,xn_2,yn_2)
       
    #add_morph_marker(axt2,E3_R15BC,Re_3_kpc)
    #add_arrow(axt2,E3_R15BC_K_SE_SDSS,Sersic2D_50rad_3_kpc,E3_R15BC,Re_3_kpc)
    SPlot.ShowcaseIndi.show_name(xo_2,yo_2, name2, A=axt2,size=16)
    
    axt2.text(text_location[0],text_location[1],r"$\rm Bin 1$",fontsize=22,color="#a5200b") 
    axt2.text(text_location[0],text_location[1]-delta_text,r"$\rm S -> E$",fontsize=22,color="k")    
    
    axt2.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt2.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt2.set_xscale( 'log' )
    axt2.set_yscale( 'log' )
        
    axt2.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)
    #axt2.legend(loc=4)
    #axt2.grid(True)
    plt.setp(axt2.get_yticklabels(), visible=False)


    axt3 = plt.subplot(gs[3])
    plot_sizemass_trans(axt3,xo_3,yo_3,xn_3,yn_3)
    SPlot.ShowcaseIndi.show_name(xo_3,yo_3, name3, A=axt3,size=16)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt3)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt3)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt3)    
    add_arrow(axt3,xo_3,yo_3,xn_3,yn_3)

    axt3.text(text_location[0],text_location[1],r"$\rm Bin 2$",fontsize=22,color="#0b5786") 
    axt3.text(text_location[0],text_location[1]-delta_text,r"$\rm E -> E$",fontsize=22,color="k")    
   
    axt3.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt3.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt3.set_xscale( 'log' )
    axt3.set_yscale( 'log' )
        
    axt3.set_ylabel(r"$\rm R_{e}$ (kpc)", fontsize=16)
    #axt3.legend(loc=4)
    #axt3.grid(True)
    plt.setp(axt3.get_xticklabels(), visible=False)
    
    axt4 = plt.subplot(gs[4])
    
    plot_sizemass_trans(axt4,xo_4,yo_4,xn_4,yn_4)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt4)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt4)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt4)    
    add_arrow(axt4,xo_4,yo_4,xn_4,yn_4)
    SPlot.ShowcaseIndi.show_name(xo_4,yo_4, name4, A=axt4,size=16)

    axt4.text(text_location[0],text_location[1],
              r"$\rm Bin 2$",fontsize=22,color="#0b5786") 
    axt4.text(text_location[0],text_location[1]-delta_text,
              r"$\rm S0 -> E$",fontsize=22,color="k")    

    axt4.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt4.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt4.set_xscale( 'log' )
    axt4.set_yscale( 'log' )
        
    #axt4.legend(loc=4)
    #axt4.grid(True)
    plt.setp(axt4.get_yticklabels(), visible=False)

    axt5 = plt.subplot(gs[5])
    
    plot_sizemass_trans(axt5,xo_5,yo_5,xn_5,yn_5)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt5)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt5)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt5)
    add_arrow(axt5,xo_5,yo_5,xn_5,yn_5)
    SPlot.ShowcaseIndi.show_name(xo_5,yo_5, name5, A=axt5,size=16)

    axt5.text(text_location[0],text_location[1],r"$\rm Bin 2$",fontsize=22,color="#0b5786") 
    axt5.text(text_location[0],text_location[1]-delta_text,r"$\rm S -> E$",fontsize=22,color="k")    
    
    axt5.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt5.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt5.set_xscale( 'log' )
    axt5.set_yscale( 'log' )
        
    #axt5.legend(loc=4)
    #axt5.grid(True)
    plt.setp(axt5.get_yticklabels(), visible=False)

    axt6 = plt.subplot(gs[6])
    
    plot_sizemass_trans(axt6,xo_6,yo_6,xn_6,yn_6,label_new="Spheroids")
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      "Barro et al. 2013", 
                                                      alpha0=0,AX=axt6)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      "van der Wel et al. 2014", 
                                                      alpha0=0,AX=axt6)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      "van Dokkum et al. 2015", 
                                                      alpha0=0,AX=axt6)      
    add_arrow(axt6,xo_6,yo_6,xn_6,yn_6)
    SPlot.ShowcaseIndi.show_name(xo_6,yo_6, name6, A=axt6,size=16)

    axt6.text(text_location[0],text_location[1],r"$\rm Bin 3$",fontsize=22,color="#2a3236") 
    axt6.text(text_location[0],text_location[1]-delta_text,r"$\rm E -> E$",fontsize=22,color="k")    

    axt6.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt6.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt6.set_xscale( 'log' )
    axt6.set_yscale( 'log' )
        
    axt6.legend(loc=3)
    #axt6.grid(True)
    axt6.set_ylabel(r"$\rm R_{e}$ (kpc)", fontsize=16)
    axt6.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)    


    axt7 = plt.subplot(gs[7])
      
    plot_sizemass_trans(axt7,xo_7,yo_7,xn_7,yn_7)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt7)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt7)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt7)     
    #add_morph_marker(axt2,E3_R15BC,Re_3_kpc)
    add_arrow(axt7,xo_7,yo_7,xn_7,yn_7)
    SPlot.ShowcaseIndi.show_name(xo_7,yo_7, name7, A=axt7,size=16)

    axt7.text(text_location[0],text_location[1],r"$\rm Bin 3$",fontsize=22,color="#2a3236") 
    axt7.text(text_location[0],text_location[1]-delta_text,r"$\rm S0 -> E$",fontsize=22,color="k")    

    
    axt7.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt7.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt7.set_xscale( 'log' )
    axt7.set_yscale( 'log' )
        
    axt7.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)
    #axt7.legend(loc=4)
    #axt7.grid(True)
    plt.setp(axt7.get_yticklabels(), visible=False)

    axt8 = plt.subplot(gs[8])
    plot_sizemass_trans(axt8,xo_8,yo_8,xn_8,yn_8)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      None, 
                                                      alpha0=0,AX=axt8)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", 
                                                      None, 
                                                      alpha0=0,AX=axt8)
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", 
                                                      None, 
                                                      alpha0=0,AX=axt8)    
    #add_morph_marker(axt2,E3_R15BC,Re_3_kpc)
    add_arrow(axt8,xo_8,yo_8,xn_8,yn_8)
    SPlot.ShowcaseIndi.show_name(xo_8,yo_8, name8, A=axt8,size=18)

    axt8.text(text_location[0],text_location[1],r"$\rm Bin 3$",fontsize=22,color="#2a3236") 
    axt8.text(text_location[0],text_location[1]-delta_text,r"$\rm S -> E$",fontsize=22,color="k")    

    
    axt8.set_xlim(left = xlim_mo[0], right = xlim_mo[1])
    axt8.set_ylim(bottom = ylim_mo[0], top = ylim_mo[1])
       
    axt8.set_xscale( 'log' )
    axt8.set_yscale( 'log' )
        
    axt8.set_xlabel(r"$\rm M_{*} / M_{\odot}$", fontsize=16)
    #axt8.legend(loc=4)
    #axt8.grid(True)
    plt.setp(axt8.get_yticklabels(), visible=False)
    plt.show()


# Apply Barro cut, select the subsample
Bcut1 = SPlot.SelectionCut(E1_R15BC_K, D1).Barro13_cut()
Bcut2 = SPlot.SelectionCut(E2_R15BC_K, D2).Barro13_cut()
Bcut3 = SPlot.SelectionCut(E3_R15BC_K, D3).Barro13_cut()

S1 = SSort.selection_generic(E1_R15BC_K, Re_1_kpc, Bcut1)
S2 = SSort.selection_generic(E2_R15BC_K, Re_2_kpc, Bcut2)
S3 = SSort.selection_generic(E3_R15BC_K, Re_3_kpc, Bcut3)

#select the morphology of the host that passes the selection
morph1_sub = SSort.cherry_pick(S1['index'], morph1_new)
morph2_sub = SSort.cherry_pick(S2['index'], morph2_new)
morph3_sub = SSort.cherry_pick(S3['index'], morph3_new)

# Bin1 size and mass

index_Bin1_E = SSort.morph_str_selection(S1['index'], morph1_sub)["E"]
index_Bin1_S0 = SSort.morph_str_selection(S1['index'], morph1_sub)["S0"]
index_Bin1_S = SSort.morph_str_selection(S1['index'], morph1_sub)["S"]

morph_Bin1_E = SSort.cherry_pick(index_Bin1_E, morph1_new)
morph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, morph1_new)
morph_Bin1_S = SSort.cherry_pick(index_Bin1_S, morph1_new)

E_host_Bin1_E = SSort.cherry_pick(index_Bin1_E, E1_R15BC_K_SE_SDSS)
R_host_Bin1_E = SSort.cherry_pick(index_Bin1_E, Sersic2D_50rad_1_kpc)
E_sph_Bin1_E = SSort.cherry_pick(index_Bin1_E, E1_R15BC)
R_sph_Bin1_E = SSort.cherry_pick(index_Bin1_E, Re_1_kpc)

E_host_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, E1_R15BC_K_SE_SDSS)
R_host_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, Sersic2D_50rad_1_kpc)
E_sph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, E1_R15BC)
R_sph_Bin1_S0 = SSort.cherry_pick(index_Bin1_S0, Re_1_kpc)

E_host_Bin1_S = SSort.cherry_pick(index_Bin1_S, E1_R15BC_K_SE_SDSS)
R_host_Bin1_S = SSort.cherry_pick(index_Bin1_S, Sersic2D_50rad_1_kpc)
E_sph_Bin1_S = SSort.cherry_pick(index_Bin1_S, E1_R15BC)
R_sph_Bin1_S = SSort.cherry_pick(index_Bin1_S, Re_1_kpc)
 
#Bin2 size and mass
index_Bin2_E = SSort.morph_str_selection(S2['index'], morph2_sub)["E"]
index_Bin2_S0 = SSort.morph_str_selection(S2['index'], morph2_sub)["S0"]
index_Bin2_S = SSort.morph_str_selection(S2['index'], morph2_sub)["S"]

morph_Bin2_E = SSort.cherry_pick(index_Bin2_E, morph2_new)
morph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, morph2_new)
morph_Bin2_S = SSort.cherry_pick(index_Bin2_S, morph2_new)

E_host_Bin2_E = SSort.cherry_pick(index_Bin2_E, E2_R15BC_K_SE_SDSS)
R_host_Bin2_E = SSort.cherry_pick(index_Bin2_E, Sersic2D_50rad_2_kpc)
E_sph_Bin2_E = SSort.cherry_pick(index_Bin2_E, E2_R15BC)
R_sph_Bin2_E = SSort.cherry_pick(index_Bin2_E, Re_2_kpc)

E_host_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, E2_R15BC_K_SE_SDSS)
R_host_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, Sersic2D_50rad_2_kpc)
E_sph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, E2_R15BC)
R_sph_Bin2_S0 = SSort.cherry_pick(index_Bin2_S0, Re_2_kpc)

E_host_Bin2_S = SSort.cherry_pick(index_Bin2_S, E2_R15BC_K_SE_SDSS)
R_host_Bin2_S = SSort.cherry_pick(index_Bin2_S, Sersic2D_50rad_2_kpc)
E_sph_Bin2_S =  SSort.cherry_pick(index_Bin2_S, E2_R15BC)
R_sph_Bin2_S = SSort.cherry_pick(index_Bin2_S, Re_2_kpc)

# Bin3 size and mass
index_Bin3_E = SSort.morph_str_selection(S3['index'], morph3_sub)["E"]
index_Bin3_S0 = SSort.morph_str_selection(S3['index'], morph3_sub)["S0"]
index_Bin3_S = SSort.morph_str_selection(S3['index'], morph3_sub)["S"]

morph_Bin3_E = SSort.cherry_pick(index_Bin3_E, morph3_new)
morph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, morph3_new)
morph_Bin3_S = SSort.cherry_pick(index_Bin3_S, morph3_new)

E_host_Bin3_E = SSort.cherry_pick(index_Bin3_E, E3_R15BC_K_SE_SDSS)
R_host_Bin3_E = SSort.cherry_pick(index_Bin3_E, Sersic2D_50rad_3_kpc)
E_sph_Bin3_E = SSort.cherry_pick(index_Bin3_E, E3_R15BC)
R_sph_Bin3_E = SSort.cherry_pick(index_Bin3_E, Re_3_kpc)

E_host_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, E3_R15BC_K_SE_SDSS)
R_host_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, Sersic2D_50rad_3_kpc)
E_sph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, E3_R15BC)
R_sph_Bin3_S0 = SSort.cherry_pick(index_Bin3_S0, Re_3_kpc)

E_host_Bin3_S = SSort.cherry_pick(index_Bin3_S, E3_R15BC_K_SE_SDSS)
R_host_Bin3_S = SSort.cherry_pick(index_Bin3_S, Sersic2D_50rad_3_kpc)
E_sph_Bin3_S =  SSort.cherry_pick(index_Bin3_S, E3_R15BC)
R_sph_Bin3_S = SSort.cherry_pick(index_Bin3_S, Re_3_kpc)

# plot the comaprison
plot_sizemass_trans_3plots(xo_0 = E_host_Bin1_E, yo_0 = R_host_Bin1_E, xn_0 = E_sph_Bin1_E, yn_0 = R_sph_Bin1_E, name0 = morph_Bin1_E, 
                           xo_1 = E_host_Bin1_S0, yo_1 = R_host_Bin1_S0, xn_1 = E_sph_Bin1_S0, yn_1 = R_sph_Bin1_S0, name1 = morph_Bin1_S0,
                           xo_2 = E_host_Bin1_S, yo_2= R_host_Bin1_S, xn_2 = E_sph_Bin1_S, yn_2 = R_sph_Bin1_S, name2 = morph_Bin1_S,
                           xo_3 = E_host_Bin2_E, yo_3 = R_host_Bin2_E, xn_3 = E_sph_Bin2_E, yn_3 = R_sph_Bin2_E, name3 = morph_Bin2_E,
                           xo_4 = E_host_Bin2_S0, yo_4 = R_host_Bin2_S0, xn_4 = E_sph_Bin2_S0, yn_4 = R_sph_Bin2_S0, name4 = morph_Bin2_S0,
                           xo_5 = E_host_Bin2_S, yo_5 = R_host_Bin2_S, xn_5 = E_sph_Bin2_S, yn_5 = R_sph_Bin2_S, name5 = morph_Bin2_S,
                           xo_6 = E_host_Bin3_E, yo_6 = R_host_Bin3_E, xn_6 = E_sph_Bin3_E, yn_6 = R_sph_Bin3_E, name6 = morph_Bin3_E,
                           xo_7 = E_host_Bin3_S0, yo_7 = R_host_Bin3_S0, xn_7 = E_sph_Bin3_S0, yn_7 = R_sph_Bin3_S0, name7 = morph_Bin3_S0,
                           xo_8 = E_host_Bin3_S, yo_8 = R_host_Bin3_S, xn_8 = E_sph_Bin3_S, yn_8 = R_sph_Bin3_S, name8 = morph_Bin3_S)

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

