#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 15:33:40 2021

@author: dexter
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:14:07 2021

@author: dexter

This is a script dedicated to
visualise the spheroid parameters.

This script produce the following plots:
0-1) stacked radial profile, as well as listing the Sersic mu0 (Done) 
0-2) Plot the average ellipticity of each galaxy
1) mu_0 - n plots (Done) 
2) mu_0 - Re plots (Done)
3) size-mass plot with curved fit
4) Seperation between Core and Sersic (Done)
5) Seperation between ETG-LTG (Done)
6) Seperation between E,S0, S (Done)
7) Tree arrangement, E,S0,S then core and Sersic
8) Tree arrangement nuclear cpt and broken disk
    
    
A) Pre-processing
    - input 
    - calculation
B) Outlier treatement
    - Find them
    - manual deletion
C) Separation
    - Morphology separation
D) Plotting
    
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
from astropy.cosmology import FlatLambdaCDM

from scipy.optimize import curve_fit
from scipy import stats
import pickle
import cmasher as cmr

"""
A) Pre-processing
    - input 
        1 read galaxy bundle, geom file, 
    - calculation
"""

plt.style.use('classic')
#mpl.rcParams['text.usetex'] = True

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale
## Read in the data, base data#################################################
# Read the name of the ISOFIT output from a list
outlist = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')

# Read the geometry file for the sample
geom_file = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat")
geom_file_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')

bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt"
bundle = SRead.read_list(bundle_name)

bundle_BD_name = "/home/dexter/SphProject/F_Gal_bundle_BD_equvi_V_cpt"
bundle_BD = SRead.read_list(bundle_BD_name)


bundle_1Sersic_name = "/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_V_cpt"
bundle_1Sersic = SRead.read_list(bundle_BD_name)


name = geom_file_n[:,0]
mu0 = geom_file[:,6]

# Read the data from the galaxy bundle
sph_mag = SRead.grab_mag("F_Gal_bundle_equvi_V_cpt", ["Bulge","CoreBulge"])
Re = SRead.grab_parameter("F_Gal_bundle_equvi_V_cpt", ["Bulge","CoreBulge"], 1) 
Sersic_n = SRead.grab_parameter("F_Gal_bundle_equvi_V_cpt", ["Bulge","CoreBulge"], 2) 
mu_e = SRead.grab_parameter("F_Gal_bundle_equvi_V_cpt", ["Bulge","CoreBulge"], 0) 

core_sersic_mu_p = SRead.grab_parameter("F_Gal_bundle_equvi_V_cpt", ["CoreBulge"], 0)

total_mag = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt")

# Read the data from the galaxy bundle from B+D model
sph_mag_BD = SRead.grab_mag("F_Gal_bundle_BD_equvi_V_cpt", ["Bulge","CoreBulge"])
Re_BD = SRead.grab_parameter("F_Gal_bundle_BD_equvi_V_cpt", ["Bulge","CoreBulge"], 1) 
Sersic_n_BD = SRead.grab_parameter("F_Gal_bundle_BD_equvi_V_cpt", ["Bulge","CoreBulge"], 2) 
mu_e_BD = SRead.grab_parameter("F_Gal_bundle_BD_equvi_V_cpt", ["Bulge","CoreBulge"], 0) 

total_mag_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_BD_equvi_V_cpt")

# Read the data from the galaxy bundle from 1Sersic model
sph_mag_1Sersic = SRead.grab_mag("F_Gal_bundle_1Sersic_equvi_V_cpt", ["Bulge","CoreBulge"])
Re_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_V_cpt", ["Bulge","CoreBulge"], 1) 
Sersic_n_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_V_cpt", ["Bulge","CoreBulge"], 2) 
mu_e_1Sersic = SRead.grab_parameter("F_Gal_bundle_1Sersic_equvi_V_cpt", ["Bulge","CoreBulge"], 0) 

total_mag_1Sersic = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_1Sersic_equvi_V_cpt")

#%%
# Read cloud 
# read NSA-Sloan catalog
nsa = SRead.read_table('/home/dexter/result/stat/completeness/nsa_sizemass.dat')

nsa_z = nsa[:,0]

nsa_mass_R15BC = nsa[:,4]

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)
nsa_d = cosmo.comoving_distance(nsa_z).value

nsa_scale = nsa_d * ars

nsa_Re = nsa[:,1]*nsa_scale

#read van der wel 2014 z=1.45 points

vdWel_highz = SRead.read_table('/home/dexter/result/stat/completeness/van_der_Wel2014_z1_25Q.dat')
vdWel_highz_2 = SRead.read_table('/home/dexter/result/stat/completeness/van_der_Wel2014_z2_25Q.dat')

vdWel_highz_mass = vdWel_highz[:,0]
vdWel_highz_size = vdWel_highz[:,1]

vdWel_highz_mass2 = vdWel_highz_2[:,0]
vdWel_highz_size2 = vdWel_highz_2[:,1]
#%%
# Read Sahu, Davis, and Savorgnan

#read data
Savorgnan_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Savorgnan_para2.dat")
Savorgnan_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Savorgnan_para2.dat",dtype='str')

Davis_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Davis_para2.dat")
Davis_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Davis_para2.dat",dtype='str')

Sahu_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/Sahu_para2.dat")
Sahu_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Sahu_para2.dat",dtype='str')

ALL3_data = SRead.read_table(
    "/home/dexter/result/stat/completeness/BHsample_updated_SPITZER.dat")
ALL3_data_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/BHsample_updated_SPITZER.dat",dtype='str')

Savorgnan_name = Savorgnan_data_n[:,0]
Savorgnan_size_eq_kpc = Savorgnan_data[:,26]

Savorgnan_morph = Savorgnan_data_n[:,4]

Savorgnan_BH_mass = 10**Savorgnan_data[:,19]

Savorgnan_Sersic_n = Savorgnan_data[:,21]
Savorgnan_b = Savorgnan_data[:,22]
Savorgnan_mu_e = Savorgnan_data[:,23]

Savorgnan_mu_0 = Savorgnan_mu_e-2.5*Savorgnan_b/np.log(10)

Savorgnan_dist = Savorgnan_data[:,17]

Savorgnan_sph_mag = Savorgnan_data[:,13]
Savorgnan_gal_Mag = Savorgnan_data[:,15]

Savorgnan_Mag = Savorgnan_sph_mag-25-5*np.log10(Savorgnan_dist)
Savorgnan_gal_Mag = Savorgnan_gal_Mag-25-5*np.log10(Savorgnan_dist)

Savorgnan_MLR = np.repeat(0.6, len(Savorgnan_Mag))
Savorgnan_MLR_IP13 = Savorgnan_data[:,28] # Nandini's corrected IP13

Savorgnan_Mo = Savorgnan_data[:,10]
Savorgnan_scale = Savorgnan_data[:,25]*(10**3) #in pc/arcsec

#calculate spheroid mass and galaxy mass
Savorgnan_mass_36 = Savorgnan_MLR*(10**((Savorgnan_Mo-Savorgnan_Mag)/2.5))
Savorgnan_gal_mass =Savorgnan_MLR*(10**((Savorgnan_Mo-Savorgnan_gal_Mag)/2.5))

Savorgnan_mass_IP13 = 10**Savorgnan_data[:,29]
Savorgnan_gal_mass_IP13 = 10**Savorgnan_data[:,31]

#Savorgnan_mass_T11 = 10**(0.88 * Savorgnan_mass_36+1.02)
#Savorgnan_mass_T11 = 10**Savorgnan_mass_36

Davis_name = Davis_data_n[:,0]
Davis_size_eq_kpc = Davis_data[:,26]
Davis_morph = Davis_data_n[:,4]

Davis_BH_mass = 10**Davis_data[:,19]

Davis_Sersic_n = Davis_data[:,21]
Davis_b = Davis_data[:,22]
Davis_mu_e = Davis_data[:,23]

Davis_mu_0 = Davis_mu_e-2.5*Davis_b/np.log(10)

Davis_dist = Davis_data[:,17]

Davis_sph_mag = Davis_data[:,13]
Davis_gal_mag = Davis_data[:,15]

Davis_Mag = Davis_sph_mag-25-5*np.log10(Davis_dist)  #V
Davis_gal_Mag = Davis_gal_mag-25-5*np.log10(Davis_dist)  #V

Davis_MLR = np.repeat(0.453, len(Davis_Mag))
Davis_MLR_IP13 = Davis_data[:,28]

Davis_Mo = Davis_data[:,10]
Davis_scale = Davis_data[:,25]*(10**3) #in pc/arcsec

Davis_mass_36 = Davis_MLR*(10**((Davis_Mo-Davis_Mag)/2.5))
Davis_gal_mass = Davis_MLR*(10**((Davis_Mo-Davis_Mag)/2.5))

Davis_mass_IP13 = 10**Davis_data[:,29]
Davis_gal_mass_IP13 = 10**Davis_data[:,31]

#Davis_mass_T11 = 10**(0.88 * Davis_mass_36+1.02)
#Davis_mass_T11 = 10**Davis_mass_36

Sahu_name = Sahu_data_n[:,0]
Sahu_size_eq_kpc = Sahu_data[:,26]
Sahu_morph = Sahu_data_n[:,4]

Sahu_BH_mass = 10**Sahu_data[:,19]

Sahu_Sersic_n = Sahu_data[:,21]
Sahu_b = Sahu_data[:,22]
Sahu_mu_e = Sahu_data[:,23]

Sahu_mu_0 = Sahu_mu_e-2.5*Sahu_b/np.log(10)

Sahu_dist = Sahu_data[:,17]

Sahu_sph_mag = Sahu_data[:,13]
Sahu_gal_mag = Sahu_data[:,15]

Sahu_Mag = Sahu_sph_mag-25-5*np.log10(Sahu_dist)  #V
Sahu_gal_Mag = Sahu_gal_mag-25-5*np.log10(Sahu_dist)

Sahu_MLR = np.repeat(0.6, len(Sahu_Mag))
Sahu_MLR_IP13 = Sahu_data[:,28]

Sahu_Mo = Sahu_data[:,10]
Sahu_scale = Sahu_data[:,25]*(10**3) #in pc/arcsec

Sahu_mass_36 = Sahu_MLR*(10**((Sahu_Mo-Sahu_Mag)/2.5))
Sahu_gal_mass = Sahu_MLR*(10**((Sahu_Mo-Sahu_gal_Mag)/2.5))

Sahu_mass_IP13 = 10**Sahu_data[:,29]
Sahu_gal_mass_IP13 = 10**Sahu_data[:,31]

#Sahu_mass_T11 = 10**(0.88 * Sahu_mass_36+1.02)
#Sahu_mass_T11 = 10**Sahu_mass_36


# read ALL3 in one go
ALL3_name = ALL3_data_n[:,0] #V
ALL3_size_eq_kpc = ALL3_data[:,26] #V
ALL3_morph = ALL3_data_n[:,4]

ALL3_BH_mass = 10**ALL3_data[:,19] #V

ALL3_Sersic_n = ALL3_data[:,21] #V
ALL3_b = ALL3_data[:,22] #V
ALL3_mu_e = ALL3_data[:,23] #V

ALL3_mu_0 = ALL3_mu_e-2.5*ALL3_b/np.log(10)

ALL3_dist = ALL3_data[:,17] #V

ALL3_sph_mag = ALL3_data[:,13]  #V
ALL3_gal_mag = ALL3_data[:,15]  #V

ALL3_Mag = ALL3_sph_mag-25-5*np.log10(ALL3_dist)  #V
ALL3_gal_Mag = ALL3_gal_mag-25-5*np.log10(ALL3_dist) #V

ALL3_MLR = np.concatenate((Savorgnan_MLR,Sahu_MLR,Davis_MLR))#V
ALL3_MLR_IP13 = 10**ALL3_data[:,28] #V

ALL3_Mo = ALL3_data[:,10] #V
ALL3_scale = ALL3_data[:,25]*(10**3) #in pc/arcsec

ALL3_mass = ALL3_MLR*(10**((ALL3_Mo-ALL3_Mag)/2.5))
ALL3_gal_mass = ALL3_MLR*(10**((ALL3_Mo-ALL3_gal_Mag)/2.5))

ALL3_mass_IP13 = 10**ALL3_data[:,29] #V
ALL3_gal_mass_IP13 = 10**ALL3_data[:,31] #V

def get_bn_iterate(array,target = 0.5):
    storage = []
    for i in range(len(array)):
        storage.append(SAna.get_bn_Nandini(array[i],z=target))
        print(i)

    return np.array(storage)

def get_mu0(mu_e,Re,Sersic_n):
    """
    For only Sersic bulges

    Parameters
    ----------
    mu_e : TYPE
        DESCRIPTION.
    Re : TYPE
        DESCRIPTION.
    Sersic_n : TYPE
        DESCRIPTION.

    Returns
    -------
    mu0 : TYPE
        DESCRIPTION.

    """
    mu0 = []
    for i in range(len(mu_e_BD)):
        mu0_i = SAna.AnalyticFunctions.mu_sersic_func(0, mu_e[i],Re[i],Sersic_n[i])
        mu0.append(mu0_i)
    mu0 = np.array(mu0)
    return mu0
    
# get the mu0
mu0_BD = get_mu0(mu_e_BD,Re_BD,Sersic_n_BD)
mu0_1Sersic = get_mu0(mu_e_1Sersic,Re_1Sersic,Sersic_n_1Sersic)

# zip the basic parameters

#%% Read in early decomposition work info
Baldry2009_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Barway_BD_decom.dat")
Baldry2009_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Barway_BD_decom.dat", dtype="str")

Baldryname = Baldry2009_table_n[:,0]

Baldry_Sersic_n = Baldry2009_table[:,5]
Baldry_mu0 = Baldry2009_table[:,3]
Baldry_R_kpc = Baldry2009_table[:,4]

Baldry_total_Mag = Baldry2009_table[:,11] #Tband
Baldry_Mag_B_T = Baldry2009_table[:,10]

offset_Baldry = 2.6
# total sph abs mag
Baldry_sph_Mag = Baldry_total_Mag-2.5*np.log10(Baldry_Mag_B_T) + offset_Baldry
# Khosroshahi 2000
Khosroshahi_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Khosroshahi_BD_decom.dat")
Khosroshahi_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Khosroshahi_BD_decom.dat", dtype="str")

Khosroshahi_mag = Khosroshahi_table[:,2]

Khosroshahi_z = Khosroshahi_table[:,3]

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)
Khosroshahi_Dist = cosmo.comoving_distance(Khosroshahi_z).value

Khosroshahi_Abs_total_mag = Khosroshahi_mag-25-5*np.log10(Khosroshahi_Dist) 

Khosroshahi_Sersic_n = Khosroshahi_table[:,6]
Khosroshahi_R_kpc = Khosroshahi_table[:,8]
Khosroshahi_mu0 = Khosroshahi_table[:,10]

Khosroshahi_B_D = Khosroshahi_table[:,12]
Khosroshahi_B_T = ((Khosroshahi_B_D**-1)+1)**-1

offset_Khosroshahi = -3.6
# total sph abs mag
Khosroshahi_sph_mag = Khosroshahi_Abs_total_mag-2.5*np.log10(Khosroshahi_B_T) +\
    offset_Khosroshahi
    
#Laurikainen 2010
Laurikainen_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Laurikainen2010_BD_decom.dat")

Laurikainen_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/Laurikainen2010_BD_decom.dat", 
    dtype="str")

Laurikainen_name = Laurikainen_table_n[:,0]

Laurikainen_Sersic_n = Laurikainen_table[:,2]
Laurikainen_R_kpc = Laurikainen_table[:,3]
Laurikainen_mu0 = Laurikainen_table[:,6]

Laurikainen_Abs_total_mag = Laurikainen_table[:,4]
Laurikainen_B_T = Laurikainen_table[:,8]

#Laurikainen_mu0 = Laurikainen_mu0-2.5*np.log10(Laurikainen_B_T)


offset_Laurikainen = 1.6
Laurikainen_sph_mag = Laurikainen_Abs_total_mag-2.5*np.log10(Laurikainen_B_T) +\
    offset_Laurikainen
    
#%%
# Get the distance from the parent sample completeness
D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt")
D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt",
    dtype = "str")

K_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin_all_Kcorr_byName.dat")

mag_g, mag_i = D0_all_table[:,11], D0_all_table[:,10] #*Galaxy colours
mag_g_kcorr, mag_i_kcorr = K_table[:,19], K_table[:,18] #Galaxy colours with K correction

D, D_lerr, D_uerr = D0_all_table[:,29], D0_all_table[:,30], D0_all_table[:,31]

morph = D0_all_table_n[:,-2] 

elle = geom_file[:,7] #extended disk ellipicity

# Read the data from the galaxy BD bundle for dust correction later
sph_mag_i_BD = SRead.grab_mag(bundle_BD_name, ["Bulge"])
disc_mag_i_BD = SRead.grab_mag(bundle_BD_name, ["Disk"])

Abs_sph_mag_i_BD = sph_mag_i_BD-25-5*np.log10(D) 
Abs_disc_mag_i_BD = disc_mag_i_BD-25-5*np.log10(D)

############################################################################################
# Calculate the absoulte magnitude and the stellar mass
ML_select_T11 = SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Taylor11_MassRatio
ML_select_RC15 = SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Roediger15BC03_MassRatio
ML_select_Z09 = SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Zibetti09_MassRatio
ML_select_IP13 =  SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Into13_MassRatio

M = SPlot.MassCalculation(sph_mag, D, 4.53,mag_g_kcorr,mag_i_kcorr)
Abs_sph_mag = M.cal_abs_mag(sph_mag, D)
Abs_total_mag = total_mag-25-5*np.log10(D) 

M_BD = SPlot.MassCalculation(sph_mag_BD, D, 4.53,mag_g_kcorr,mag_i_kcorr)
Abs_sph_mag_BD = M.cal_abs_mag(sph_mag_BD, D)
Abs_total_mag_BD = total_mag_BD-25-5*np.log10(D) 

M_1Sersic = SPlot.MassCalculation(sph_mag_1Sersic, D, 4.53,mag_g_kcorr,mag_i_kcorr)
Abs_sph_mag_1Sersic = M.cal_abs_mag(sph_mag_1Sersic, D)
Abs_total_mag_1Sersic = total_mag_1Sersic-25-5*np.log10(D) 

#establish B/T ratio
B_T_ratio_i_old = 10**((sph_mag-total_mag)/(-2.5)) #old B/T ratio
D_B_ratio_i_old = ((B_T_ratio_i_old)**-1)  - np.repeat(1,len(B_T_ratio_i_old))  #old D/B ratio
B_D_ratio_i_old = D_B_ratio_i_old**-1

#print('D_B_ratio_i_old',D_B_ratio_i_old)

# Calculate the dust corrected version of abs mag of sph for ALL sample 
# E and S0 does not require that but I calculate them nonethesless
Abs_sph_mag_dustCorr = M.dust_correction_Driver08(Abs_sph_mag,elle)

# Calculate the dust corrected version of i-band galaxy total mag for ALL sample 
#sph_mag_i_BD_dustcorr = M.dust_correction_Driver08(Abs_sph_mag_i_BD,elle,struc = "Bulge", band = "i")
#disc_mag_i_BD_dustcorr = M.dust_correction_Driver08(Abs_disc_mag_i_BD,elle,struc = "Disk", band = "i")

#****
disc_mag = sph_mag-2.5*np.log10(D_B_ratio_i_old) #disc_mag based on old D/B ratio
disc_mag2 = Abs_sph_mag_dustCorr-2.5*np.log10(D_B_ratio_i_old)
#***

Abs_disc_mag = disc_mag-25-5*np.log10(D)  
Abs_disc_mag_dustCorr = M.dust_correction_Driver08(Abs_disc_mag,elle,struc="Disk",band="i")

# assume B/T ratio to be the same in g-band, apply the same dust correction on g-band
B_D_ratio_i = 10**((Abs_sph_mag_dustCorr-Abs_disc_mag_dustCorr)/(-2.5))
B_T_ratio_i = (1+B_D_ratio_i**-1)**-1
#B_D_ratio_i = B_D_ratio_i_old

#The total magnitude of the dust corrected, total i-band galaxy magnitude
#mag_i_dustcorr = -2.5*np.log10(10**(sph_mag_i_BD_dustcorr/-2.5) + 10**(disc_mag_i_BD_dustcorr/-2.5))
mag_i_dustcorr = -2.5*np.log10(10**(Abs_sph_mag_dustCorr/-2.5) + 10**(Abs_disc_mag_dustCorr/-2.5))


#Calculate the dust corrected version of g-band galaxy total mag for ALL sample 
L_gal_g = 10**(mag_g_kcorr/-2.5)
L_sph_g = L_gal_g * B_T_ratio_i_old
L_disc_g = L_gal_g * (np.repeat(1,len(B_T_ratio_i_old))-B_T_ratio_i_old)

sph_mag_g_BD = -2.5*np.log10(L_sph_g)
disc_mag_g_BD = -2.5*np.log10(L_disc_g)

sph_mag_g_BD_dustcorr = M.dust_correction_Driver08(sph_mag_g_BD,elle,struc = "Bulge", band = "g")
disc_mag_g_BD_dustcorr = M.dust_correction_Driver08(disc_mag_g_BD,elle,struc = "Disk", band = "g")

#The total magnitude of the dust corrected, total g-band galaxy magnitude
mag_g_dustcorr = -2.5*np.log10(10**(sph_mag_g_BD_dustcorr/-2.5) + 10**(disc_mag_g_BD_dustcorr/-2.5))

#print(mag_g_dustcorr,mag_i_dustcorr)

K_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_BinV_Kcorr_EXT.dat")

K_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_BinV_Kcorr_EXT.dat", 
    dtype='str')

# the corrected mag g and i, Kcorrection+EXTINCTIOn
mag_g, mag_i = K_table[:,10], K_table[:,9]
g_EXT, i_EXT = K_table[:,23], K_table[:,24]
g_kcorr, i_kcorr = K_table[:,25], K_table[:,26]

mag_g_corr, mag_i_corr = mag_g-g_kcorr, mag_i-i_kcorr


# Calculate the absoulte magnitude and the stellar mass
ML_select_T11_dustcorr = SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Taylor11_MassRatio
ML_select_RC15_dustcorr = SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Roediger15BC03_MassRatio
ML_select_Z09_dustcorr = SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Zibetti09_MassRatio
ML_select_IP13_dustcorr =  SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Into13_MassRatio

# Calculate the error of the Re and stellar mass
MLR = ML_select_RC15_dustcorr
MLR_dust = ML_select_RC15_dustcorr
MLR_e = 10**0.1
mag_e = 0.3 #magnitude errorNote that the y-axis of (1) is the same range of the x-axis of (2)

# Calculate the stellar mass, without dust correction
E_T11_K = M.cal_Mass(MLR)
E_BD = M_BD.cal_Mass(MLR)
E_1Sersic = M_1Sersic.cal_Mass(MLR)

# Calculate the dust corrected version of stellar mass for ALL sample 
E_T11_K_dustCorr = MLR_dust*(10**((4.53-Abs_sph_mag_dustCorr)/2.5))

#print(MLR_dust,MLR)
E_T11_K_dustCorr_old = MLR*(10**((4.53-Abs_sph_mag_dustCorr)/2.5))

mass_uerr = np.sqrt(((mag_e/2.5)**2)+((2*D_uerr/(D*np.log(10)))**2)+((MLR_e/(MLR*np.log(10)))**2))
mass_err = mass_uerr

# calculate the scale
scale = D* ars

scale_other = scale * 10**3

scale_lerr, scale_uerr = (D-D_lerr)*ars, (D+D_uerr)*ars

# Calculate Re with multi-cpt, B+D, and 1-Sersic
Re_kpc = Re* scale
Re_kpc_BD = Re_BD* scale
Re_kpc_1Sersic = Re_1Sersic* scale

Re_kpc_lerr, Re_kpc_uerr = abs(Re_kpc - Re* scale_lerr) , abs(Re* scale_uerr - Re_kpc)
Re_kpc_err =[Re_kpc_lerr, Re_kpc_uerr]

print(Sersic_n[26])
#print(name[5],mu_e[5],Re_kpc[5],Sersic_n[5],sph_mag[5],Abs_total_mag[5],np.log10(E_T11_K[5]))

#print(name[42],mu_e[42],Re_kpc[42],Sersic_n[42],sph_mag[42],Abs_total_mag[42],np.log10(E_T11_K[42]))

#print(name[54],mu_e[54],Re_kpc[54],Sersic_n[54],sph_mag[54],Abs_total_mag[54],np.log10(E_T11_K[54]))

#%%
# Read in the data from the others
core = SRead.grab_parameter_whole(bundle_name, ["CoreBulge"])
#######Seperation 1###########################################################
#core vs not core seperation here
S = SSort.seperator_label_generic(bundle, ["Bulge","CoreBulge"])
S_bar = SSort.seperator_label_generic(bundle, ["PrimBar"])

# Extract the info from Bulge only
sph_bulge_mag = list(SRead.grab_mag(S[0]["Bulge"], ["Bulge"]))
Re_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 1) )
Sersic_n_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 2))
mu_e_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 0))
total_mag_bulge = list(SRead.grab_total_mag(S[0]["Bulge"]))

# cherry pick info based on the index list
name_bulge = SSort.cherry_pick(S[0]["Bulge_index"], name)
mu0_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mu0)

mag_g_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_g)
mag_i_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_i)
mag_g_kcorr_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_g_kcorr)
mag_i_kcorr_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_i_kcorr)

D_bulge = SSort.cherry_pick(S[0]["Bulge_index"], D)
D_lerr_bulge =  SSort.cherry_pick(S[0]["Bulge_index"], D_lerr)
D_uerr_bulge =  SSort.cherry_pick(S[0]["Bulge_index"], D_uerr)

Abs_sph_mag_bulge = SSort.cherry_pick(S[0]["Bulge_index"], Abs_sph_mag)
E_T11_K_bulge = SSort.cherry_pick(S[0]["Bulge_index"], E_T11_K)

MLR_bulge = SSort.cherry_pick(S[0]["Bulge_index"], MLR)
mass_uerr_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mass_uerr)
mass_err_bulge = np.array(mass_uerr_bulge)

scale_bulge = np.array(D_bulge)* ars
scale_lerr_bulge, scale_uerr_bulge = (np.array(D_bulge)-np.array(D_lerr_bulge))*ars, (np.array(D_bulge)+np.array(D_uerr_bulge))*ars
Re_kpc_bulge = Re_bulge* scale_bulge

Re_kpc_lerr_bulge, Re_kpc_uerr_bulge = abs(Re_kpc_bulge - Re_bulge* scale_lerr_bulge) , abs(Re_bulge* scale_uerr_bulge - Re_kpc_bulge)
Re_kpc_err_bulge =[Re_kpc_lerr_bulge, Re_kpc_uerr_bulge]

Abs_total_mag_bulge = SSort.cherry_pick(S[0]["Bulge_index"], Abs_total_mag)

# Extract the info from Core Bulge only
sph_corebulge_mag = list(SRead.grab_mag(S[1]["CoreBulge"], ["CoreBulge"]))
Re_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 1))
Sersic_n_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 2)) 
mu_e_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 0))
total_mag_corebulge = list(SRead.grab_total_mag(S[1]["CoreBulge"]))

# cherry pick info based on the index list (output list)
name_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], name)
mu0_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mu0)

#print("coreBulge",name_corebulge[12],name_corebulge[14])

mag_g_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_g)
mag_i_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_i)
mag_g_kcorr_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_g_kcorr)
mag_i_kcorr_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_i_kcorr)

D_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], D)
D_lerr_corebulge =  SSort.cherry_pick(S[1]["CoreBulge_index"], D_lerr)
D_uerr_corebulge =  SSort.cherry_pick(S[1]["CoreBulge_index"], D_uerr)

Abs_sph_mag_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], Abs_sph_mag)
E_T11_K_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], E_T11_K)

MLR_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], MLR)
mass_uerr_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mass_uerr)
mass_err_corebulge = np.array(mass_uerr_corebulge)

elle_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], elle)

scale_corebulge = np.array(D_corebulge)* ars
scale_lerr_corebulge, scale_uerr_corebulge = (np.array(D_corebulge)-np.array(D_lerr_corebulge))*ars, (np.array(D_corebulge)+np.array(D_uerr_corebulge))*ars
Re_kpc_corebulge = Re_corebulge* scale_corebulge

Abs_total_mag_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], Abs_total_mag)

Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge = abs(Re_kpc_corebulge - Re_corebulge* scale_lerr_corebulge) , abs(Re_corebulge* scale_uerr_corebulge - Re_kpc_corebulge)
Re_kpc_err_corebulge =[Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge]

# put the two array together  
n_combine = np.array([Sersic_n_corebulge,Sersic_n_bulge])
mu0_combine = np.array([mu0_corebulge,mu0_bulge])
Mag_combine = np.array([Abs_sph_mag_corebulge,Abs_sph_mag_bulge])

Abs_total_mag_combine = np.array([Abs_total_mag_corebulge,Abs_total_mag_bulge])

#K correction and extinction for sph mag
sph_mag = sph_mag - i_EXT - i_kcorr
#%%
def get_bn_Re_mu_1090():
    # calculate b_z,n with z = 0.1, 0.5, and 0.9 
    bn_10 = get_bn_iterate(Sersic_n,target=0.05)
    bn_50 = get_bn_iterate(Sersic_n,target=0.5)
    bn_90 = get_bn_iterate(Sersic_n,target=0.95)
    
    #print('bn_90',bn_90)
    
    # calculate R_z with z = 0.1 and 0.9
    Re_kpc_10 = Re_kpc*(bn_10/bn_50)**Sersic_n
    Re_kpc_90 = Re_kpc*(bn_90/bn_50)**Sersic_n
    
    # calculate the mu_z with z = 0.1, 0.9
    mu_10 = mu0 + (2.5*bn_10)/np.log(10)
    mu_50 = mu0 + (2.5*bn_50)/np.log(10)
    mu_90 = mu0 + (2.5*bn_90)/np.log(10)
    
    Re_kpc_10_combine_morph = seperate_morph_simple(Re_kpc_10,index, morph) 
    Re_kpc_90_combine_morph = seperate_morph_simple(Re_kpc_90,index, morph) 

    mu0_10_combine_morph = seperate_morph_simple(mu_10,index, morph) 
    mu0_50_combine_morph = seperate_morph_simple(mu_50,index, morph)
    mu0_90_combine_morph = seperate_morph_simple(mu_90,index, morph) 
    
    bn_dict = {"bn_10":bn_10,
               "bn_50":bn_50,
               "bn_90":bn_90,
               "Re_kpc_10":Re_kpc_10,
               "Re_kpc_90":Re_kpc_90,
               "mu_10": mu_10,
               "mu_50": mu_50,
               "mu_90": mu_90,
               "Re_kpc_10_combine_morph": Re_kpc_10_combine_morph, 
               "Re_kpc_90_combine_morph": Re_kpc_90_combine_morph,
               "mu0_10_combine_morph": mu0_10_combine_morph,
               "mu0_50_combine_morph": mu0_50_combine_morph,
               "mu0_90_combine_morph": mu0_90_combine_morph}
    
    with open('bn_dict.pkl', 'wb') as f:
        pickle.dump(bn_dict, f)
        
    return (Re_kpc_10_combine_morph,Re_kpc_90_combine_morph, 
            mu0_10_combine_morph, mu0_50_combine_morph, mu0_90_combine_morph)


def get_bn_Re_mu_1090_BIG(Sersic_n,Re_kpc,mu0,index,morph):
    # calculate b_z,n with z = 0.1, 0.5, and 0.9 
    bn_05 = get_bn_iterate(Sersic_n,target=0.05)
    bn_10 = get_bn_iterate(Sersic_n,target=0.10)
    bn_20 = get_bn_iterate(Sersic_n,target=0.20)
    bn_30 = get_bn_iterate(Sersic_n,target=0.30)
    bn_40 = get_bn_iterate(Sersic_n,target=0.40)
    bn_50 = get_bn_iterate(Sersic_n,target=0.50)
    bn_60 = get_bn_iterate(Sersic_n,target=0.60)
    bn_70 = get_bn_iterate(Sersic_n,target=0.70)
    bn_80 = get_bn_iterate(Sersic_n,target=0.80)
    bn_90 = get_bn_iterate(Sersic_n,target=0.90)
    bn_95 = get_bn_iterate(Sersic_n,target=0.95)
    bn_100 = get_bn_iterate(Sersic_n,target=1.0)
    
    #print('bn_90',bn_90)
    
    # calculate R_z with z = 0.1 and 0.9
    Re_kpc_05 = Re_kpc*(bn_05/bn_50)**Sersic_n
    Re_kpc_10 = Re_kpc*(bn_10/bn_50)**Sersic_n
    Re_kpc_20 = Re_kpc*(bn_20/bn_50)**Sersic_n
    Re_kpc_30 = Re_kpc*(bn_30/bn_50)**Sersic_n
    Re_kpc_40 = Re_kpc*(bn_40/bn_50)**Sersic_n
    Re_kpc_50 = Re_kpc*(bn_50/bn_50)**Sersic_n
    Re_kpc_60 = Re_kpc*(bn_60/bn_50)**Sersic_n
    Re_kpc_70 = Re_kpc*(bn_70/bn_50)**Sersic_n
    Re_kpc_80 = Re_kpc*(bn_80/bn_50)**Sersic_n
    Re_kpc_90 = Re_kpc*(bn_90/bn_50)**Sersic_n
    Re_kpc_95 = Re_kpc*(bn_95/bn_50)**Sersic_n
    Re_kpc_100 = Re_kpc*(bn_100/bn_50)**Sersic_n

    # calculate the mu_z with z = 0.1, 0.9
    mu_05 = mu0 + (2.5*bn_05)/np.log(10)
    mu_10 = mu0 + (2.5*bn_10)/np.log(10)
    mu_20 = mu0 + (2.5*bn_20)/np.log(10)
    mu_30 = mu0 + (2.5*bn_30)/np.log(10)
    mu_40 = mu0 + (2.5*bn_40)/np.log(10)
    mu_50 = mu0 + (2.5*bn_50)/np.log(10)
    mu_60 = mu0 + (2.5*bn_60)/np.log(10)
    mu_70 = mu0 + (2.5*bn_70)/np.log(10)
    mu_80 = mu0 + (2.5*bn_80)/np.log(10)
    mu_90 = mu0 + (2.5*bn_90)/np.log(10)
    mu_95 = mu0 + (2.5*bn_95)/np.log(10)
    mu_100 = mu0 + (2.5*bn_100)/np.log(10)
    
    #Re_kpc_05_combine_morph = seperate_morph_simple(Re_kpc_05,index, morph) 
    #Re_kpc_95_combine_morph = seperate_morph_simple(Re_kpc_95,index, morph) 

    #mu0_05_combine_morph = seperate_morph_simple(mu_05,index, morph) 
    #mu0_50_combine_morph = seperate_morph_simple(mu_50,index, morph)
    #mu0_95_combine_morph = seperate_morph_simple(mu_95,index, morph) 
    
    bn_dict = {"bn_05":bn_05, "bn_10":bn_10, "bn_20":bn_20, "bn_30":bn_30,
               "bn_40":bn_40, "bn_50":bn_50, "bn_60":bn_60, "bn_70":bn_70,
               "bn_80":bn_80, "bn_90":bn_90, "bn_95":bn_95, "bn_100":bn_100,
               "Re_kpc_05":Re_kpc_05, "Re_kpc_10":Re_kpc_10, "Re_kpc_20":Re_kpc_20,
               "Re_kpc_30":Re_kpc_30, "Re_kpc_40":Re_kpc_40, "Re_kpc_50":Re_kpc_50,
               "Re_kpc_60":Re_kpc_60, "Re_kpc_70":Re_kpc_70, "Re_kpc_80":Re_kpc_80,
               "Re_kpc_90":Re_kpc_90, "Re_kpc_05":Re_kpc_95, "Re_kpc_05":Re_kpc_100,
               "mu_05": mu_05, "mu_10": mu_10,"mu_20": mu_20, "mu_30": mu_30, 
               "mu_40": mu_40, "mu_50": mu_50, "mu_60": mu_60, "mu_70": mu_70, 
               "mu_80": mu_80, "mu_90": mu_90,  "mu_95": mu_95, "mu_100": mu_100 
               }
    
    with open('bn_dict_BIG.pkl', 'wb') as f:
        pickle.dump(bn_dict, f)
        
    return bn_dict

def plot_stack_surface_brightness_profile(r):
    bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt"
    bundle=SRead.read_list(bundle_name)
    name = SRead.grab_name(bundle_name)
    
    fig = plt.figure(figsize=(6.4, 4.8))
    new_mu_list = []
    
    #plot individual curve
    for i in range(len(bundle)):        
        for j in range(len(bundle[i])):
            identifier = bundle[i][j]
            if identifier == "Bulge": #if the spheroid is a Sersic bulge
                # Load the Sersic fits parameters
                para = bundle[i][j+1]
                
                #plot the radial SB profile using the Sersic parameters
                line = SAna.AnalyticFunctions.mu_sersic_func(r,*list(para))
                
                plt.plot(r,line,"r--",lw =3)
                
                #get b and mu_e
                mu_e, Re ,n = para[0], para[1], para[2]
                print(mu_e, Re ,n)
                
                b = SAna.get_bn(n) # get_bn is the exact function used in profiler
                # obtain mu_0 by equ7 from Graham2005
                new_mu = mu_e-2.5*b/np.log(10)
              #  print(name[i],new_mu,Re[i],n,sph_mag[i]) # list mu0
                new_mu_list.append(new_mu) # save the new mu_0
                
                #list name, first element of from the radial SB, and the mu_0 
                #from equ 7
                print(name[i],line[0])#, new_mu
                pass
            elif identifier == "CoreBulge": #if the spheroid is a Core Sersic bulge
                para = bundle[i][j+1]
                
                # calculate the mu_e of the core Sersic funtion
                Re, n = para[1], para[2]
                mu_e = SAna.AnalyticFunctions.simple_mu_core_sersic(Re,*list(para))
                print(mu_e, Re ,n)

                # put in the core Sersic fit parameter mu_e, r_e , n_ser
                line = SAna.AnalyticFunctions.mu_sersic_func(r,mu_e,para[1],para[2])
                plt.plot(r,line,"k-",lw= 3)
                new_mu_list.append(line[0])

                #print(name[i],line[0], Re[i],n,sph_mag[i]) #list mu0
                print(name[i],line[0])#, new_mu

                pass
    plt.gca().invert_yaxis()
    #plt.ylim(np.log10(24),np.log(10))
    plt.ylim(24,12)
    plt.xlim(0,120)
    
    plt.tight_layout()
    plt.show()
    
    return fig,new_mu_list

def list_mu0_extrapolate(new_mu_list):
    """
    List the name and mu0
    mu0 is an extrapolation to the centre

    """
    for i in range(len(name)):
        print(name[i],new_mu_list[i])


#%% Read other works
#size_mass_Graham_equ(Mass,ML,n,A,B)
#print("fitting", curve_fit(SAna.AnalyticFunctions.size_mass_Graham_equ,(E_T11_K,MLR),Re_kpc))

#mass_gen = np.linspace()
#MLR_gen = np.linspace()
#Sersic_n_gen = np.

mass0 = np.linspace(8.5,12.5,15)

ZOE = 0.53*(((10**mass0)/(3e10))**-0.2)*(0.5+0.5*((10**mass0)/(3e10))**8)**0.119
# Shen 2003 size-mass relatio (z-band)
#Shen2003 = SAna.AnalyticFunctions.size_mass_powerlaw(10**mass0,0.347e-5,0.56)

mass_shen03 = np.linspace(10,12.5,15)

Shen2003 = (0.347e-5)*(10**mass_shen03)**0.56 # equ 17 for Fig11 z-band

# Lange 2015 ETG by morphology
mass_Lange15 = np.linspace(9.4,12.5,15)
Lange2015_ETG_morph = (2.44e-5)*(10**mass_Lange15)**0.49 #morphology cut i-band Lange2015


mass1, mass2 = np.linspace(10.4,12.5,7), np.linspace(9,10.4,7)
Saracco2018_high = 2.7e-7*(10**mass1)**0.64
Saracco2018_low = 26*(10**mass2)**-0.13

# Lange 2016 size-mass relation from table 1
#Lange2016_E = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,7.4e-5,0.44)
#Lange2016_E_M1e10 = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,1.382,0.643)
#Lange2016_ETG_bulge = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,1.836,0.267)

mass_Lange16 = np.linspace(9,12.5,15)

Lange2016_E = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Lange16)/1e10,2.114,0.329)
Lange2016_E_M1e10 = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Lange16)/1e10,1.382,0.643)
Lange2016_ETG_bulge = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Lange16)/1e10,1.836,0.267)
Lange2016_final_bulge = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Lange16)/1e10,2.063,0.263)


mass_Nedkova21 = np.linspace(10.3,11.6,15)
# 0.2< z < 0.
Nedkova2021_lowz = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Nedkova21)/5e10,10**0.61,0.68)
Nedkova2021_highz1 = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Nedkova21)/5e10,10**0.45,0.64)
Nedkova2021_highz2 = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Nedkova21)/5e10,10**0.28,0.63)
Nedkova2021_highz3 = SAna.AnalyticFunctions.size_mass_powerlaw((10**mass_Nedkova21)/5e10,10**0.18,0.61)

mass_Baldry21 = np.linspace(10,11.5,15)
Baldry2020 = 10**((1/1.5)*mass_Baldry21-(9.85/1.5))


#Dimauro2018
mass_Dimauro18 = np.linspace(10.3,12.5,15)
Dimauro2018_bulge_lowz1 = (10**0.43)*((10**mass_Dimauro18)/5e10)**0.38

#Méndez-Abreu 2021
mass_Mendez_Abreu21 = np.linspace(10.5,12,15)
Mendez_Abreu2021_bulge_highmass = (10**(-4.22+0.71*mass_Mendez_Abreu21))*(1/1e3)
Mendez_Abreu2021_bulge_lowmass = 10**(-0.39+0.34*mass_Mendez_Abreu21)*(1/1e3)

## End in reading base data#################################################
#%% 
"""
B) Outlier treatment
    
"""
# investigation
# M_i > -20
# weird core-Sersic gal, NGC4382, NGC4636
# Core+ disk nothing special??
Q = SSort.selection_generic(name_corebulge, Abs_sph_mag_corebulge, np.repeat(-20,len(name_corebulge)),
                            direction="high")

#print(len(sph_corebulge_mag),len(Abs_sph_mag_corebulge))
Q2 = SSort.selection_generic(sph_corebulge_mag, Abs_sph_mag_corebulge, np.repeat(-20,len(name_corebulge)),
                            direction="high")
print(Q,Q2) #show the core spheroid magnitude


Q33 = SSort.selection_generic(mu0_bulge, name_bulge, np.repeat(5,len(mu0_bulge)),
                            direction="low",axis='x')

print("Q33",Q33)

print(D_corebulge[12],D_corebulge[14]) # show the distance

print(SPlot.MassCalculation(10.006806, D_corebulge[12], 4.53,np.zeros(1),np.zeros(1)).cal_abs_mag(10.006806, D_corebulge[12]))
print(SPlot.MassCalculation(10.05999, D_corebulge[14], 4.53,np.zeros(1),np.zeros(1)).cal_abs_mag(10.05999, D_corebulge[14]))

ML_T11 = SPlot.MLRelationIband(np.array([mag_g_corebulge[12]]),np.array([mag_i_corebulge[12]])).Taylor11_MassRatio
ML_RC15 = SPlot.MLRelationIband(np.array([mag_g_corebulge[12]]),np.array([mag_i_corebulge[12]])).Roediger15BC03_MassRatio
ML_IP13 = SPlot.MLRelationIband(np.array([mag_g_corebulge[12]]),np.array([mag_i_corebulge[12]])).Into13_MassRatio


#print(np.log10(SPlot.MassCalculation(10.006806, 17.881, 4.53,
#                            mag_g_corebulge[12],mag_i_corebulge[12]).cal_Mass(ML_T11)))
#print(np.log10(SPlot.MassCalculation(10.006806, 17.881, 4.53,
#                            mag_g_corebulge[12],mag_i_corebulge[12]).cal_Mass(ML_RC15)))
#print(np.log10(SPlot.MassCalculation(10.006806, D_corebulge[12], 4.53,
#                            mag_g_corebulge[12],mag_i_corebulge[12]).cal_Mass(ML_IP13)))
# Find 3.2,4.6, 4e10 5.2e10 S0

# weird S gal, high bugle size and mass NGC3270, paticular high size. nuclear cpt
#Q3 = SSort.zone_in_2D(np.array([list(E_T11_K_S),list(Re_kpc_S)]),[4e10,3.7],[5.2e10,4.6])

#print(Q3)
#print(name[Q3["index"]])

#weird dim S0, NGC3665 nuclear **Dusty disk
#Q4 = SSort.zone_in_2D([Sersic_n_S0,Abs_sph_mag_S0],[0,-18],[1,-16])

#print(Q4,name[Q4["index"]])

print("outlier index", SSort.outlier_detect(mu0, 20,direction=">"))
print("outlier index", SSort.outlier_detect(np.array(mu0_corebulge), 20,direction=">"))


#%% Delete the outliers
# group the data by intervals
C1_n = SSort.group_to_bin1D(Abs_total_mag, Sersic_n, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])

C2_mu = SSort.group_to_bin1D(Abs_total_mag, mu0, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])
    
C3_Re = SSort.group_to_bin1D(Abs_total_mag, Re, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])    


# Remove outliers via sigma clipping
def sigma_clip_loop(C_x,C_y):
    out_Cx, out_Cy = [], []
    for i in range(len(C_x)):
        out_Cx_i, out_Cy_i = SSort.sigma_clip2D(C_x[i],C_y[i],sigma=2.0)
        
        out_Cx.append(out_Cx_i)
        out_Cy.append(out_Cy_i)
    out_Cx, out_Cy = np.concatenate(out_Cx), np.concatenate(out_Cy)
    
    return out_Cx, out_Cy

out_Cxn, out_Cyn = sigma_clip_loop(np.array(C1_n['y']),np.array(C1_n['x']))
out_Cxmu, out_Cymu = sigma_clip_loop(np.array(C2_mu['y']),np.array(C2_mu['x']))
out_CxRe, out_CyRe = sigma_clip_loop(np.array(C3_Re['y']),np.array(C3_Re['x']))

#print("Length",len(out_Cxn),len(out_Cyn),len(out_Cxmu), 
#      len(out_Cymu), len(out_CxRe), len(out_CyRe)) 

#Remove outliers manually
def outlier_removal_manual(morph, sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag, scale, outlier_list):
    morph = np.delete(morph, outlier_list)
    sph_mag = np.delete(sph_mag, outlier_list)
    Re = np.delete(Re, outlier_list)
    Sersic_n = np.delete(Sersic_n, outlier_list)
    mu_e = np.delete(mu_e, outlier_list)
    total_mag = np.delete(total_mag, outlier_list)
    name = np.delete(name, outlier_list)
    mu0 = np.delete(mu0, outlier_list)
    mag_g_kcorr = np.delete(mag_g_kcorr, outlier_list)
    mag_i_kcorr = np.delete(mag_i_kcorr, outlier_list)
    D = np.delete(D, outlier_list)              
    D_lerr =np.delete(D_lerr, outlier_list)   
    D_uerr =np.delete(D_uerr, outlier_list)    
    Abs_sph_mag = np.delete(Abs_sph_mag, outlier_list)
    elle =np.delete(elle, outlier_list)   
    E_T11_K =np.delete(E_T11_K, outlier_list)   
    MLR = np.delete(MLR, outlier_list)   
    mass_uerr = np.delete(mass_uerr, outlier_list)   
    Abs_total_mag = np.delete(Abs_total_mag, outlier_list)
    #Re_kpc = np.delete(Re_kpc, [43,73])
    scale = np.delete(scale, outlier_list)   
    
    return (morph, sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag, scale)
     

morph, sph_mag, Re, Sersic_n, mu_e,total_mag, \
    name, mu0, mag_g_kcorr, mag_i_kcorr, D,D_lerr, \
        D_uerr, Abs_sph_mag, elle, E_T11_K,\
                     MLR, mass_uerr,Abs_total_mag, scale = \
                         outlier_removal_manual(\
                         morph, sph_mag, Re, Sersic_n, mu_e,total_mag, \
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D, \
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K, 
                     MLR, mass_uerr,Abs_total_mag, scale, [43,73])                             

morph_corebulge, sph_corebulge_mag, Re_corebulge, Sersic_n_corebulge, \
mu_e_corebulge,total_mag_corebulge, \
    name_corebulge, mu0_corebulge, mag_g_kcorr_corebulge, \
    mag_i_kcorr_corebulge, D_corebulge,D_lerr_corebulge, \
        D_uerr_corebulge, Abs_sph_mag_corebulge, elle_corebulge, E_T11_K_corebulge,\
                     MLR_corebulge, mass_uerr_corebulge, \
                         Abs_total_mag_corebulge, scale_corebulge = \
                         outlier_removal_manual(\
                         morph, sph_corebulge_mag, Re_corebulge, Sersic_n_corebulge, 
                         mu_e_corebulge,total_mag_corebulge, \
                     name_corebulge, mu0_corebulge, mag_g_kcorr_corebulge, \
                     mag_i_kcorr_corebulge, D_corebulge, \
                     D_lerr_corebulge, D_uerr_corebulge, Abs_sph_mag_corebulge, \
                     elle_corebulge, E_T11_K_corebulge, \
                     MLR_corebulge, mass_uerr_corebulge,Abs_total_mag_corebulge, \
                     scale_corebulge, [8,14])  


scale_other = scale * 10**3

#recalculate some qunatities
Re_kpc = Re* scale
#Re_kpc_BD = Re_BD* scale
#Re_kpc_1Sersic = Re_1Sersic* scale

#scale_corebulge = np.array(D_corebulge)* ars
#scale_lerr_corebulge, scale_uerr_corebulge = (np.array(D_corebulge)-np.array(D_lerr_corebulge))*ars, (np.array(D_corebulge)+np.array(D_uerr_corebulge))*ars
Re_kpc_corebulge = Re_corebulge* scale_corebulge


#Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge = abs(Re_kpc_corebulge - Re_corebulge* scale_lerr_corebulge) , abs(Re_corebulge* scale_uerr_corebulge - Re_kpc_corebulge)
#Re_kpc_err_corebulge =[Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge]

#Re_kpc_lerr, Re_kpc_uerr = abs(Re_kpc - Re* scale_lerr) , abs(Re* scale_uerr - Re_kpc)
#Re_kpc_err =[Re_kpc_lerr, Re_kpc_uerr]

index= np.arange(0,len(Re),1)
print(len(Abs_total_mag),len(Sersic_n),len(mu0),len(Re))

#%%
"""
C) Separation
"""
####Start Seperation 2 ########################################################
index = np.arange(len(morph))
index_SG16 = np.arange(len(Savorgnan_morph))
index_D19 = np.arange(len(Davis_morph))
index_S19 = np.arange(len(Sahu_morph))
index_A3 = np.arange(len(ALL3_morph))

def seperate_morph_simple(A,index, morph):
    morph_dict = SSort.morph_str_selection(index, morph)
    
    A_E = SSort.cherry_pick(morph_dict['E'], A)
    A_S0 = SSort.cherry_pick(morph_dict['S0'], A)
    A_S = SSort.cherry_pick(morph_dict['S'], A)
    
    A_combine_morph = np.array([A_E,A_S0,A_S],dtype=object)
    return A_combine_morph


def seperation_morph(index, morph, sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag):

    #Seperate the sample based on morphlogy (ours)
    morph_dict = SSort.morph_str_selection(index, morph)
        
    sph_mag_E = SSort.cherry_pick(morph_dict['E'], sph_mag)
    sph_mag_S0 = SSort.cherry_pick(morph_dict['S0'], sph_mag)
    sph_mag_S = SSort.cherry_pick(morph_dict['S'], sph_mag)
    
    Re_E = SSort.cherry_pick(morph_dict['E'], Re)
    Re_S0 = SSort.cherry_pick(morph_dict['S0'], Re)
    Re_S = SSort.cherry_pick(morph_dict['S'], Re)
    
    Sersic_n_E =  SSort.cherry_pick(morph_dict['E'], Sersic_n)
    Sersic_n_S0 =  SSort.cherry_pick(morph_dict['S0'], Sersic_n)
    Sersic_n_S =  SSort.cherry_pick(morph_dict['S'], Sersic_n)
    
    mu_e_E = SSort.cherry_pick(morph_dict['E'], mu_e)
    mu_e_S0 = SSort.cherry_pick(morph_dict['S0'], mu_e)
    mu_e_S = SSort.cherry_pick(morph_dict['S'], mu_e)

    total_mag_E = SSort.cherry_pick(morph_dict['E'], total_mag)
    total_mag_S0 = SSort.cherry_pick(morph_dict['S0'], total_mag)
    total_mag_S = SSort.cherry_pick(morph_dict['S'], total_mag)

    # cherry pick info based on the index list
    name_E = SSort.cherry_pick(morph_dict['E'], name)
    name_S0 = SSort.cherry_pick(morph_dict['S0'], name)
    name_S = SSort.cherry_pick(morph_dict['S'], name)

    mu0_E = SSort.cherry_pick(morph_dict['E'], mu0)
    mu0_S0 = SSort.cherry_pick(morph_dict['S0'], mu0)
    mu0_S = SSort.cherry_pick(morph_dict['S'], mu0)

    mag_g_kcorr_E = SSort.cherry_pick(morph_dict['E'], mag_g_kcorr)
    mag_g_kcorr_S0 = SSort.cherry_pick(morph_dict['S0'], mag_g_kcorr)
    mag_g_kcorr_S = SSort.cherry_pick(morph_dict['S'], mag_g_kcorr)
    mag_i_kcorr_E = SSort.cherry_pick(morph_dict['E'], mag_i_kcorr)
    mag_i_kcorr_S0 = SSort.cherry_pick(morph_dict['S0'], mag_i_kcorr)
    mag_i_kcorr_S = SSort.cherry_pick(morph_dict['S'], mag_i_kcorr)
    
    D_E = SSort.cherry_pick(morph_dict['E'], D)
    D_S0 = SSort.cherry_pick(morph_dict['S0'], D)
    D_S = SSort.cherry_pick(morph_dict['S'], D)

    D_lerr_E =  SSort.cherry_pick(morph_dict['E'], D_lerr)
    D_lerr_S0 =  SSort.cherry_pick(morph_dict['S0'], D_lerr)
    D_lerr_S =  SSort.cherry_pick(morph_dict['S'], D_lerr)
    D_uerr_E =  SSort.cherry_pick(morph_dict['E'], D_uerr)
    D_uerr_S0 =  SSort.cherry_pick(morph_dict['S0'], D_uerr)
    D_uerr_S =  SSort.cherry_pick(morph_dict['S'], D_uerr)
    
    Abs_sph_mag_E = SSort.cherry_pick(morph_dict['E'], Abs_sph_mag)
    Abs_sph_mag_S0 = SSort.cherry_pick(morph_dict['S0'], Abs_sph_mag)
    Abs_sph_mag_S = SSort.cherry_pick(morph_dict['S'], Abs_sph_mag)
    
    elle_S = SSort.cherry_pick(morph_dict['S'],elle)
    #Abs_sph_mag_dustCorr_S = SSort.cherry_pick(morph_dict['S'],Abs_sph_mag_dustCorr)
    
    # my mass
    E_T11_K_E = SSort.cherry_pick(morph_dict['E'], E_T11_K)
    E_T11_K_S0 = SSort.cherry_pick(morph_dict['S0'], E_T11_K)
    E_T11_K_S = SSort.cherry_pick(morph_dict['S'], E_T11_K)

    #E_T11_K_dustcorr_E = SSort.cherry_pick(morph_dict['E'], E_T11_K_dustCorr)
    #E_T11_K_dustcorr_S0 = SSort.cherry_pick(morph_dict['S0'], E_T11_K_dustCorr)
    #E_T11_K_dustcorr_S = SSort.cherry_pick(morph_dict['S'], E_T11_K_dustCorr)

    #E_T11_K_dustcorr_E_old = SSort.cherry_pick(morph_dict['E'], E_T11_K_dustCorr_old)   
    #E_T11_K_dustcorr_S0_old = SSort.cherry_pick(morph_dict['S0'], E_T11_K_dustCorr_old)
    #E_T11_K_dustcorr_S_old = SSort.cherry_pick(morph_dict['S'], E_T11_K_dustCorr_old)

    MLR_E = SSort.cherry_pick(morph_dict['E'], MLR)
    MLR_S0 = SSort.cherry_pick(morph_dict['S0'], MLR)
    MLR_S = SSort.cherry_pick(morph_dict['S'], MLR)
    
    mass_uerr_E = SSort.cherry_pick(morph_dict['E'], mass_uerr)
    mass_uerr_S0 = SSort.cherry_pick(morph_dict['S0'], mass_uerr)
    mass_uerr_S = SSort.cherry_pick(morph_dict['S'], mass_uerr)
    
    mass_err_E = np.array(mass_uerr_E)
    mass_err_S0 = np.array(mass_uerr_S0)
    mass_err_S = np.array(mass_uerr_S)
    
    scale_E = np.array(D_E)*ars
    scale_S0 = np.array(D_S0)*ars
    scale_S = np.array(D_S)*ars
    
    scale_lerr_E, scale_uerr_E = (np.array(D_E)-np.array(D_lerr_E))*ars, (np.array(D_E)+np.array(D_uerr_E))*ars
    scale_lerr_S0, scale_uerr_S0 = (np.array(D_S0)-np.array(D_lerr_S0))*ars, (np.array(D_S0)+np.array(D_uerr_S0))*ars
    scale_lerr_S, scale_uerr_S = (np.array(D_S)-np.array(D_lerr_S))*ars, (np.array(D_S)+np.array(D_uerr_S))*ars
    
    Re_kpc_E = Re_E* scale_E
    Re_kpc_S0 = Re_S0* scale_S0
    Re_kpc_S = Re_S* scale_S
        
    Re_kpc_lerr_E, Re_kpc_uerr_E = abs(Re_kpc_E - Re_E* scale_lerr_E) , abs(Re_E* scale_uerr_E - Re_kpc_E)
    Re_kpc_lerr_S0, Re_kpc_uerr_S0 = abs(Re_kpc_S0 - Re_S0* scale_lerr_S0) , abs(Re_S0* scale_uerr_S0 - Re_kpc_S0)
    Re_kpc_lerr_S, Re_kpc_uerr_S = abs(Re_kpc_S - Re_S* scale_lerr_S) , abs(Re_S* scale_uerr_S - Re_kpc_S)
        
    Re_kpc_err_E =[Re_kpc_lerr_E, Re_kpc_uerr_E]
    Re_kpc_err_S0 =[Re_kpc_lerr_S0, Re_kpc_uerr_S0]
    Re_kpc_err_S =[Re_kpc_lerr_S, Re_kpc_uerr_S]

    Abs_total_mag_E = SSort.cherry_pick(morph_dict['E'], Abs_total_mag)
    Abs_total_mag_S0 = SSort.cherry_pick(morph_dict['S0'], Abs_total_mag)
    Abs_total_mag_S = SSort.cherry_pick(morph_dict['S'], Abs_total_mag)

    # put the Three array together  
    n_combine_morph = np.array([Sersic_n_E,Sersic_n_S0,Sersic_n_S])
    mu0_combine_morph = np.array([mu0_E,mu0_S0,mu0_S])
    Mag_combine_morph = np.array([Abs_sph_mag_E,Abs_sph_mag_S0, Abs_sph_mag_S])
    
    E_combine_morph = np.array([E_T11_K_E,E_T11_K_S0,E_T11_K_S])
    Re_kpc_combine_morph = np.array([Re_kpc_E,Re_kpc_S0,Re_kpc_S])
    Re_kpc_err_combine = np.array([Re_kpc_err_E, Re_kpc_err_S0, Re_kpc_err_S])
    mass_err_combine = np.array([mass_err_E,mass_err_S0,mass_err_S])
    
    Abs_total_mag_combine_morph = np.array([Abs_total_mag_E,Abs_total_mag_S0,Abs_total_mag_S])

    # put only the ETG array together  
    #n_combine_morph = np.array([Sersic_n_E,Sersic_n_S0])
    #mu0_combine_morph = np.array([mu0_E,mu0_S0])
    #Mag_combine_morph = np.array([Abs_sph_mag_E,Abs_sph_mag_S0])
    #
    #E_combine_morph = np.array([E_T11_K_E,E_T11_K_S0])
    #Re_kpc_combine_morph = np.array([Re_kpc_E,Re_kpc_S0])
    #Re_kpc_err_combine = np.array([Re_kpc_err_E, Re_kpc_err_S0])
    #mass_err_combine = np.array([mass_err_E,mass_err_S0])
    ##
    #Abs_total_mag_combine_morph = np.array([Abs_total_mag_E,Abs_total_mag_S0,Abs_total_mag_S])    


    return (n_combine_morph, mu0_combine_morph, Mag_combine_morph, 
            E_combine_morph, Re_kpc_combine_morph, Re_kpc_err_combine,
            mass_err_combine, Abs_total_mag_combine_morph)

#%%
"""
D) plotting
"""        
def plot_n_mu0_Mag_3plot(n,mu,Re,Mag, 
                         Sersic_n_corebulge,
                         mu0_corebulge, Re_kpc_corebulge, 
                         Abs_sph_mag_corebulge,
                         label=[],fit_instruc=0):
    """
    Plot a similar plot as in Graham 2019 figure 1 with my sample

    Parameters
    ----------
    n : numpy array
        Sersic indices
    mu0 : numpy array
        The central surface brightness.
    Mag : numpy array
        The absoulte magnitude.
    label : list, optional
        optional list of labels. The default is [].
    fit_instruc: float
        Fitting instruction. The default is the 0 element of the input array.
    """
    #Plotting
    fig = plt.figure(figsize=(6.4*3, 5.2))
    gs = gridspec.GridSpec(ncols=3, nrows=1,
                               hspace=0, wspace=0.0) 
    # Define colour and markers
    markers = ["o","X", "s","+","*"]
    colour = ["k","#d20d0d","#ebb800","#0f920f","#0f4592"]
    
    # The Mag vs Sersic n plot
    axt0 = plt.subplot(gs[0])   
    
    # Check the dimension of the input array 
    if n.shape[0] != len(label):
        axt0.plot(n,Mag,'ko',ms=10)
    elif n.shape[0] == len(label):
        for i in range(len(n)):
            axt0.plot(n[i],Mag[i], linestyle="None",
                      marker=markers[i], color = colour[i],
                      ms=10,label=label[i])
    
    # plot the scatter point     
    axt0.scatter(Sersic_n_corebulge,Abs_sph_mag_corebulge,
                 facecolors='none', edgecolors='r', s = 500,
                 label=r"$\rm core-S\'ersic$")

    # Calculate the Pearson correlation for the sample and its subgroups
    r_n = stats.pearsonr(np.log10(Sersic_n),Abs_sph_mag)
    r_E_n = stats.pearsonr(np.log10(n_combine_morph[0]),
                           Abs_total_mag_combine_morph[0])
    r_disc_n = stats.pearsonr(np.log10(np.concatenate((n_combine_morph[1],
                                                       n_combine_morph[2]))),\
                     np.concatenate((Abs_total_mag_combine_morph[1],
                                     Abs_total_mag_combine_morph[2])))
        
    #Calculate the Spearsman correlation for the sample and its subgroups
    rs_n = stats.spearmanr(np.log10(Sersic_n),Abs_sph_mag)
    rs_E_n = stats.spearmanr(np.log10(n_combine_morph[0]),
                           Abs_total_mag_combine_morph[0])
    rs_disc_n = stats.spearmanr(np.log10(np.concatenate((n_combine_morph[1],
                                                       n_combine_morph[2]))),\
                     np.concatenate((Abs_total_mag_combine_morph[1],
                                     Abs_total_mag_combine_morph[2])))
        
    axt0.text(4,-18.5-1.4, r'$\bf r_p(All)={:.2f}$'.format(r_n[0]),fontsize=18)
    axt0.text(4,-18.5-0.7, r'$r_p(E+ES)={:.2f}$'.format(r_E_n[0]),
              fontsize=18)
    axt0.text(4,-18.5, r"$r_p(S0+S)={:.2f}$".format(r_disc_n[0]),
              fontsize=18)
    
    axt0.text(4,-16-1.4, r'$\bf r_s(All)={:.2f}$'.format(rs_n[0]),fontsize=18)
    axt0.text(4,-16-0.7, r'$r_s(E+ES)={:.2f}$'.format(rs_E_n[0]),
              fontsize=18)
    axt0.text(4,-16, r"$r_s(S0+S)={:.2f}$".format(rs_disc_n[0]),
              fontsize=18)
    
    #axt0.legend(loc="lower right")
    axt0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n_\mathrm{Sph} $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
    
    axt0.legend(loc="upper left")
    axt0.tick_params(axis='x',labelsize=15)
    axt0.tick_params(axis='y',labelsize=15)
    # The Mag vs mu_0 plot
    axt1 = plt.subplot(gs[1],sharey=axt0)

    # Check the dimension of the input array  
    if mu.shape[0] != len(label):
        axt1.plot(mu0,Mag,'ko',ms=10)
    elif mu.shape[0] == len(label):
        for i in range(len(mu)):
            axt1.plot(mu[i],Mag[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i]) 
            
    # plot the scatter point     
    axt1.scatter(mu0_corebulge,Abs_sph_mag_corebulge,
                 facecolors='none', edgecolors='r', s = 500)
    
    # 3-sigma
    std_mu0 = np.std(mu0)
    print('median',np.median(mu0),'3sigma',3*std_mu0)
    
    # Calculate the Pearson correlation for the sample and its subgroups
    r_mu = stats.pearsonr(mu0,Abs_sph_mag)
    r_E_mu = stats.pearsonr(mu0_combine_morph[0],Abs_total_mag_combine_morph[0])
    r_disc_mu = stats.pearsonr(np.concatenate((mu0_combine_morph[1],mu0_combine_morph[2])),\
                     np.concatenate((Abs_total_mag_combine_morph[1],Abs_total_mag_combine_morph[2])))
 
    # Calculate the Spearsman correlation for the sample and its subgroups
    rs_mu = stats.spearmanr(mu0,Abs_sph_mag)
    rs_E_mu = stats.spearmanr(mu0_combine_morph[0],Abs_total_mag_combine_morph[0])
    rs_disc_mu = stats.spearmanr(np.concatenate((mu0_combine_morph[1],mu0_combine_morph[2])),\
                     np.concatenate((Abs_total_mag_combine_morph[1],Abs_total_mag_combine_morph[2])))
      
    # show the text on screen    
    axt1.text(8,-18.5-1.4, r"$\bf r_p(All)={:.2f}$".format(r_mu[0]),fontsize=18)
    axt1.text(8,-18.5-0.7, r"$r_p(E+ES)={:.2f}$".format(r_E_mu[0]),fontsize=18)
    axt1.text(8,-18.5, r"$r_p(S0+S)={:.2f}$".format(r_disc_mu[0]),fontsize=18)
        
    axt1.text(8,-16-1.4, r"$\bf r_s(All)={:.2f}$".format(rs_mu[0]),fontsize=18)
    axt1.text(8,-16-0.7, r"$r_s(E+ES)={:.2f}$".format(rs_E_mu[0]),fontsize=18)
    axt1.text(8,-16, r"$r_s(S0+S)={:.2f}$".format(rs_disc_mu[0]),fontsize=18)

    axt1.set_xlabel(r"$ \mu_\mathrm{0,Sph}~(\rm mag~arcsec^{-2})$", fontsize=22)
    axt1.set_xlim(19.2,-0.5)
    axt1.tick_params(axis='x',labelsize=15)

    # The Mag vs R_e
    axt2 = plt.subplot(gs[2],sharey=axt0)

    # Check the dimension of the input array  
    if Re.shape[0] != len(label):
        axt2.plot(Re,Mag,'ko',ms=10)
    elif Re.shape[0] == len(label):
        for i in range(len(Re)):
            axt2.plot(Re[i],Mag[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i])   
    axt2.scatter(Re_kpc_corebulge,Abs_sph_mag_corebulge,
                 facecolors='none', edgecolors='r', s = 500, 
                 label=r'$\rm core-S\'ersic$')
    
    # Calculate the Pearson correlation for the sample and its subgroups
    r_R = stats.pearsonr(np.log10(Re_kpc),Abs_sph_mag)
    r_E_R = stats.pearsonr(
        np.log10(Re_kpc_combine_morph[0]),Abs_total_mag_combine_morph[0])
    r_disc_R = stats.pearsonr(
        np.log10(np.concatenate((Re_kpc_combine_morph[1],
                                 Re_kpc_combine_morph[2]))),\
                     np.concatenate((Abs_total_mag_combine_morph[1],
                                     Abs_total_mag_combine_morph[2])))
    # Calculate the Spearsman correlation for the sample and its subgroups
    rs_R = stats.spearmanr(np.log10(Re_kpc),Abs_sph_mag)
    rs_E_R = stats.spearmanr(
        np.log10(Re_kpc_combine_morph[0]),Abs_total_mag_combine_morph[0])
    rs_disc_R = stats.spearmanr(
        np.log10(np.concatenate((Re_kpc_combine_morph[1],
                                 Re_kpc_combine_morph[2]))),\
                     np.concatenate((Abs_total_mag_combine_morph[1],
                                     Abs_total_mag_combine_morph[2])))
        
    #show text on screen
    axt2.text(7,-18.5-1.4, r"$\bf r_p(All)={:.2f}$".format(r_R[0]),fontsize=18)
    axt2.text(7,-18.5-0.7, r"$r_p(E+ES)={:.2f}$".format(r_E_R[0]),
              fontsize=18)
    axt2.text(7,-18.5, r"$r_p(S0+S)={:.2f}$".format(r_disc_R[0]),
              fontsize=18)
        
    axt2.text(7,-16-1.4, r"$\bf r_s(All)={:.2f}$".format(rs_R[0]),fontsize=18)
    axt2.text(7,-16-0.7, r"$r_s(E+ES)={:.2f}$".format(rs_E_R[0]),
              fontsize=18)
    axt2.text(7,-16, r"$r_s(S0+S)={:.2f}$".format(rs_disc_R[0]),
              fontsize=18)
  
    axt2.set_xlabel(r"$ R_\mathrm{e,Sph}~(\rm kpc)$", fontsize=22)
    axt2.set_xscale('log')
    axt2.set_xlim(0.09,120)
 
    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.setp(axt2.get_yticklabels(), visible=False)
    
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
    plt.tight_layout()

    axt2.tick_params(axis='x',labelsize=15)

    plt.show()
    return (fig)#,*popt_n,*popt_mu)
#%%
def plot_n_mu0_Mag_3Lplot(n,mu,Re, 
                         label=[],fit_instruc=0):
    """
    Plot a similar plot as in Graham 2019 figure 1 with my sample

    Parameters
    ----------
    n : numpy array
        Sersic indices
    mu0 : numpy array
        The central surface brightness.
    Mag : numpy array
        The absoulte magnitude.
    label : list, optional
        optional list of labels. The default is [].
    fit_instruc: float
        Fitting instruction. The default is the 0 element of the input array.
    """
    #Plotting
    fig = plt.figure(figsize=(6.4*2, 6.4*2))
    gs = gridspec.GridSpec(ncols=2, nrows=2,
                               hspace=0.2, wspace=0.0) 
    # Define colour and markers
    markers = ["o","X", "s","+","*"]
    colour = ["k","#d20d0d","#ebb800","#0f920f","#0f4592"]
    
    # The Mag vs Sersic n plot
    axt0 = plt.subplot(gs[0,0])   
    
    # Check the dimension of the input array 
    if n.shape[0] != len(label):
        axt0.plot(np.log10(n),np.log10(Re),'ko',ms=10)
    elif n.shape[0] == len(label):
        for i in range(len(n)):
            axt0.plot(np.log10(n[i]),np.log10(Re[i]), linestyle="None",
                      marker=markers[i], color = colour[i],
                      ms=10,label=label[i])
    
    # Calculate the Pearson correlation for the sample and its subgroups
    r_n = stats.pearsonr(np.log10(Sersic_n),
                         np.log10(Re_kpc))
    r_E_n = stats.pearsonr(np.log10(n[0]),np.log10(Re[0]))
    r_disc_n = stats.pearsonr(np.log10(np.concatenate((n[1],n[2]))),\
                     np.log10(np.concatenate((Re[1],Re[2]))))
        
    #Calculate the Spearsman correlation for the sample and its subgroups
    rs_n = stats.spearmanr(np.log10(Sersic_n),np.log10(Re_kpc))
    rs_E_n = stats.spearmanr(np.log10(n[0]),np.log10(Re[0]))
    rs_disc_n = stats.spearmanr(np.log10(np.concatenate((n[1], n[2]))),\
                     np.log10(np.concatenate((Re[1],Re[2]))))
        
    axt0.text(-0.3,np.log10(70), r'$\bf r_p(All)={:.2f}$'.format(r_n[0]),fontsize=18)
    axt0.text(-0.3,np.log10(70-22), r'$r_p(E+ES)={:.2f}$'.format(r_E_n[0]),
              fontsize=18)
    axt0.text(-0.3,np.log10(70-37), r"$r_p(S0+S)={:.2f}$".format(r_disc_n[0]),
              fontsize=18)
    
    axt0.text(-0.3,np.log10(70-48), r'$\bf r_s(All)={:.2f}$'.format(rs_n[0]),fontsize=18)
    axt0.text(-0.3,np.log10(70-56), r'$r_s(E+ES)={:.2f}$'.format(rs_E_n[0]),
              fontsize=18)
    axt0.text(-0.3,np.log10(70-60), r"$r_s(S0+S)={:.2f}$".format(rs_disc_n[0]),
              fontsize=18)
    
    axt0.set_ylabel(r"$\log_{10}(R_\mathrm{e,Sph})$", fontsize=22)
    axt0.set_xlabel(r"$\log_{10}(n_\mathrm{Sph})$", fontsize=22)

    axt0.set_xlim(-0.4,1.2)
    axt0.set_ylim(-1.19,2.27)
    
    axt0.legend(loc="lower right")
    axt0.tick_params(axis='x',labelsize=15)
    axt0.tick_params(axis='y',labelsize=15)
    # The Mag vs mu_0 plot
    axt1 = plt.subplot(gs[1,0])

    # Check the dimension of the input array  
    if mu.shape[0] != len(label):
        axt1.plot(n,mu,'ko',ms=10)
    elif mu.shape[0] == len(label):
        for i in range(len(mu)):
            axt1.plot(n[i],mu[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i]) 

    # 3-sigma
    std_mu0 = np.std(mu0)
    print('median',np.median(mu0),'3sigma',3*std_mu0)
    
    # Calculate the Pearson correlation for the sample and its subgroups
    r_mu = stats.pearsonr(Sersic_n,
                          mu0)
    r_E_mu = stats.pearsonr(n[0],mu[0])
    r_disc_mu = stats.pearsonr(np.concatenate((n[1],n[2])),\
                     np.concatenate((mu[1],mu[2])))
 
    # Calculate the Spearsman correlation for the sample and its subgroups
    rs_mu = stats.spearmanr(Sersic_n,mu0)
    rs_E_mu = stats.spearmanr(n[0],mu[0])
    rs_disc_mu = stats.spearmanr(np.concatenate((n[1],n[2])),\
                     np.concatenate((mu[1],mu[2])))
      
    # show the text on screen    
    axt1.text(0.5,2.53, r"$\bf r_p(All)={:.2f}$".format(r_mu[0]),fontsize=18)
    axt1.text(0.5,2.53+1.0, r"$r_p(E+ES)={:.2f}$".format(r_E_mu[0]),fontsize=18)
    axt1.text(0.5,2.53+2.0, r"$r_p(S0+S)={:.2f}$".format(r_disc_mu[0]),fontsize=18)
        
    axt1.text(0.5,2.53+3.0, r"$\bf r_s(All)={:.2f}$".format(rs_mu[0]),fontsize=18)
    axt1.text(0.5,2.53+4.0, r"$r_s(E+ES)={:.2f}$".format(rs_E_mu[0]),fontsize=18)
    axt1.text(0.5,2.53+5.0, r"$r_s(S0+S)={:.2f}$".format(rs_disc_mu[0]),fontsize=18)

    axt1.set_xlabel(r"$\mathrm{Sersic}~n_\mathrm{Sph}$", fontsize=22)
    axt1.set_ylabel(r"$ \mu_\mathrm{0,Sph}~(\rm mag~arcsec^{-2})$", fontsize=22)
    axt1.set_xlim(-0.5,11.3)
    axt1.set_ylim(0.91,18.8)

    axt1.tick_params(axis='x',labelsize=15)
    axt1.tick_params(axis='y',labelsize=15)
    # The Mag vs R_e
    axt2 = plt.subplot(gs[1,1],sharey=axt1)

    # Check the dimension of the input array  
    if Re.shape[0] != len(label):
        axt2.plot(Re,mu,'ko',ms=10)
    elif Re.shape[0] == len(label):
        for i in range(len(Re)):
            axt2.plot(Re[i],mu[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i])   
    
    # Calculate the Pearson correlation for the sample and its subgroups
    r_R = stats.pearsonr(Re_kpc,mu0)
    r_E_R = stats.pearsonr(Re[0],mu[0])
    r_disc_R = stats.pearsonr(
        np.concatenate((Re[1],Re[2])),\
                     np.concatenate((mu[1],
                                     mu[2])))
    # Calculate the Spearsman correlation for the sample and its subgroups
    rs_R = stats.spearmanr(Re_kpc,mu0)
    rs_E_R = stats.spearmanr(Re[0],mu[0])
    rs_disc_R = stats.spearmanr(np.concatenate((Re[1],Re[2])),\
                     np.concatenate((mu[1],
                                     mu[2])))
        
    #show text on screen
    axt2.text(0.12,2.53, r"$\bf r_p(All)={:.2f}$".format(r_R[0]),fontsize=18)
    axt2.text(0.12,2.53+1, r"$r_p(E+ES)={:.2f}$".format(r_E_R[0]),
              fontsize=18)
    axt2.text(0.12,2.53+2, r"$r_p(S0+S)={:.2f}$".format(r_disc_R[0]),
              fontsize=18)
        
    axt2.text(0.12,2.53+3, r"$\bf r_s(All)={:.2f}$".format(rs_R[0]),fontsize=18)
    axt2.text(0.12,2.53+4, r"$r_s(E+ES)={:.2f}$".format(rs_E_R[0]),
              fontsize=18)
    axt2.text(0.12,2.53+5, r"$r_s(S0+S)={:.2f}$".format(rs_disc_R[0]),
              fontsize=18)
  
    axt2.set_xlabel(r"$ R_\mathrm{e,Sph}~(\rm kpc)$", fontsize=22)
    axt2.set_xscale('log')
    axt2.set_xlim(0.09,120)
    axt2.set_ylim(0.91,18.8)
    axt2.tick_params(axis='x',labelsize=15)
    axt2.tick_params(axis='y',labelsize=15)
    
    #plt.setp(axt1.get_yticklabels(), visible=False)
    plt.setp(axt2.get_yticklabels(), visible=False)
    
    plt.gca().invert_yaxis()
    plt.tight_layout()

    axt2.tick_params(axis='x',labelsize=15)

    plt.show()
    return (fig)#,*popt_n,*popt_mu)
#%%
xlim = [1e9,3e12]
ylim = [0.08,167]

#xlim = [3e7,10e12]
#ylim = [0.08,240]
# plot size-mass diagram
def plot_dexter_sample_T11(mass, size, size_err,mass_err, 
                           A, scale='log',alpha=0.65,
                           colour='#a5200b',label="This work",marker='o',s=70):
    A.scatter(mass, size,marker='o',c=colour,label=label, 
              s =s, alpha=alpha)
    A.errorbar(mass, size, yerr = size_err, ms =10,
                  xerr = mass_err*mass, ls='none',linewidth=4, 
                  color = colour,
                  ecolor= colour, capsize=0, alpha=alpha, marker=marker)

    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale(scale)
    A.set_yscale(scale)
    plt.tight_layout()
    
#%%
from scipy.optimize import curve_fit
import linmix as linmix

def plot_sizemass_z0():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    
    func_S = SAna.AnalyticFunctions.size_mass_Graham_equ(10**mass0,
                                                         -9.4, -14.3, 
                                                         2.0/3, -29.5)
    func_B = SAna.AnalyticFunctions.size_mass_Graham_equ(10**mass0,
                                                         -7.5, -20, 
                                                         0.6, -29.7,ML=8,A=4) #M_sun=5.09
    #print("func_S",func_S)
    
    ax0.plot(mass0,np.log10(func_S),'g-',lw = 4)
    ax0.plot(mass0,np.log10(func_B),'g--',lw = 4)
    
    #input the data
    mass1, mass2, mass3, mass4 = np.log10(E_T11_K), np.log10(Savorgnan_mass_36), \
                        np.log10(Davis_mass_36), np.log10(Sahu_mass_36)
        
    R1, R2, R3, R4 = np.log10(Re_kpc), np.log10(Savorgnan_size_eq_kpc), \
                        np.log10(Davis_size_eq_kpc), np.log10(Sahu_size_eq_kpc)
    
    ax0.scatter(mass1, R1, s=100,
                           alpha = 0.8,marker='o',color='#a5200b',
                          label=r"$\rm Hon~et~al.~(2022)$")
    ax0.scatter(mass2, R2, s=100,
                           alpha = 0.8, marker='o', color='#b940c8',
                           label=r"$\rm Savorgnan~&~Graham~(2016)$")
    ax0.scatter(mass3, R3, s=100,
                           alpha = 0.8, marker='o', color='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$")
    ax0.scatter(mass4, R4, s=100,
                           alpha = 0.8, marker='o', color='#e1a000',
                          label=r"$\rm Sahu~et~al.~(2019)$")
    
    #combine them all together
    E_all  = np.concatenate((mass1,mass2,mass3,mass4))
    Re_all = np.concatenate((R1,R2,R3,R4))
    
    # run fit
    lm = linmix.LinMix(E_all, Re_all, K=2)
    lm.run_mcmc(silent=True)
    
    for i in range(0, len(lm.chain), 25):
        xs = mass0
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        ax0.plot(xs, ys, color='r', alpha=0.02,lw=3)
        
    alpha = np.median(lm.chain[:]['alpha'])
    beta = np.median(lm.chain[:]['beta'])
    sigsqr = np.median(lm.chain[:]['sigsqr'])
    
    alpha_std = np.median(lm.chain[:]['alpha'].std())
    beta_std = np.median(lm.chain[:]['beta'].std())
    
    print("FIT para",10**alpha*(1e10)**beta,beta,sigsqr)
    
    #median line
    ys = alpha + xs * beta
    
    # vertical 
    ys_v = alpha + E_all * beta
    # calculate the root mean square
    delta_rms1 = np.sqrt(sum((ys_v - Re_all)**2)/len(ys_v))
    
    ax0.plot(xs, ys, linestyle='dashdot', color='k',lw=3)
    
    # run fit again, flip
    lm2 = linmix.LinMix(Re_all, E_all, K=2)
    lm2.run_mcmc(silent=True)
    
    ys2 = mass0
    for i in range(0, len(lm2.chain), 25):
        xs2 = (ys2- lm2.chain[i]['alpha']) / lm2.chain[i]['beta']
        ax0.plot(ys2, xs2, color='cyan', alpha=0.02,lw=3)
        
    alpha2 = np.median(lm2.chain[:]['alpha'])
    beta2 = np.median(lm2.chain[:]['beta'])
    sigsqr2 = np.median(lm2.chain[:]['sigsqr'])
    
    beta2_std=np.std(1/lm2.chain[:]['beta'])
    alpha2_std=np.std(-1.0*lm2.chain[:]['alpha']/lm2.chain[:]['beta'])
    print("FIT para",alpha2,beta2,sigsqr2)
    
    # median line
    xs2 = (ys2 - alpha2)/ beta2
    
    # vertical
    xs2_v = (E_all - alpha2)/ beta2
    # calculate the root mean square
    delta_rms2 = np.sqrt(sum((xs2_v - Re_all)**2)/len(xs2_v))
    
    ax0.plot(ys2, xs2,'--', color='k',lw=3)
        
    # Bisector line
    m_C, k_C = find_bisector_line(1/beta2,beta,-1.0*alpha2/beta2,alpha)
    
    xs3 = mass0
    ys3 = m_C * xs3 + k_C
    
    # vertical
    ys3_v = m_C * E_all + k_C
    # calculate the root mean square
    delta_rms3 = np.sqrt(sum((ys3_v - Re_all)**2)/len(ys3_v))


    plt.plot(xs3, ys3, '-', color='k',lw=3) # Note that x and y are reversed

    # Text
    # X as independent variable
    ax0.text(9.05,1.8+0.2,
             r"$R_{\rm e,Sph}(M_{\rm *,Sph})$",
             fontsize=12)
    ax0.text(9.05,1.8,
             r"$S={}\pm{},int.={}\pm{}, \Delta_{}={}$".format(round(beta,2),
                                                              round(beta_std,2),
                                                              round(alpha,2),
                                                              round(alpha_std,2), 
                                                              "{rms}",
                                                              round(delta_rms1,2)),
             fontsize=12)

    ax0.hlines(1.8+0.2+0.05,9.9,10.2,linestyle='dashdot',color='k',lw=4)
    
    # X as independent variable
    ax0.text(9.05,1.8-0.2,
             r"$M_{*,Sph}(R_{\rm e,Sph})$",fontsize=12)
    ax0.text(9.05,1.8-0.4,r"$S={}\pm{},int.={}\pm{}, \Delta_{}={}$".format(round(1/beta2,2),
                                                             round(beta2_std,2),
                                                             round(-1.0*alpha2/beta2,2),
                                                             round(abs(alpha2_std),2),
                                                             "{rms}",
                                                             round(delta_rms2,2)
                                                             ),
             fontsize=12)
    ax0.hlines(1.8-0.2+0.05,9.9,10.2,linestyle='dashed',color='k',lw=4)

    ax0.text(9.05,1.8-0.6,
             r"$Bisector~line$",fontsize=12)
    ax0.text(9.05,1.8-0.8,r"$S={}, int.={}, \Delta_{}={}$".format(round(m_C,2),
                                                                  round(k_C,2), 
                                                                  "{rms}", 
                                                                  round(delta_rms3,2)),
             fontsize=12)
    ax0.hlines(1.8-0.6+0.05,9.9,10.2,linestyle='solid',color='k',lw=4)

    #ax0.text(9.5,1.3-0.2,r"$S={}\pm{}$".format(),fontsize=10)
    #ax0.text(9.5,1.3-0.2,r"$b={:.2f}$".format(beta),fontsize=18)

    ax0.set_ylabel(r"$ \mathrm{log_{10}}(R_\mathrm{e,Sph}~\rm(kpc))$", fontsize=16)
    ax0.set_xlabel(r"$ \mathrm{log_{10}}(M_\mathrm{*,Sph} / \rm M_\mathrm{\odot})$", fontsize=16)
    
    plt.legend(fontsize = 10,loc="lower right")
    plt.tight_layout()
    
    ax0.set_xlim(left = np.log10(xlim[0]), right = np.log10(xlim[1]))
    ax0.set_ylim(bottom = np.log10(ylim[0]), top = np.log10(ylim[1]))

    plt.show()


def plot_sizemass_Sersic_vs_core():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    plot_dexter_sample_T11(E_T11_K_bulge, Re_kpc_bulge,Re_kpc_err_bulge,mass_err_bulge,ax0,colour='k',label=r"$\rm S\'{e}rsic$")
    plot_dexter_sample_T11(E_T11_K_corebulge, Re_kpc_corebulge,Re_kpc_err_corebulge,mass_err_corebulge,ax0,label=r"$\rm Core-S\'{e}rsic$")
    ax0.set_ylabel("$ R_\mathrm{e,sph}$ (kpc)", fontsize=16)
    ax0.set_xlabel(r"$ M_\mathrm{*,sph} / \rm M_\mathrm{\odot} $", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()

    plt.show()

Dist0 = np.linspace(0,120,2000)


def plot_sizemass_morph():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    
    #R_corebulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'],Re_kpc)
    #mass_corebulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'], E_T11_K)
    ms0 = 10

    E_all, Re_all = np.concatenate(E_combine_morph), np.concatenate(Re_kpc_combine_morph)
    E_all_SG16, Re_all_SG16 = np.concatenate(E_combine_morph_SG16), np.concatenate(Re_kpc_combine_morph_SG16)
    E_all_D19, Re_all_D19 = np.concatenate(E_combine_morph_D19), np.concatenate(Re_kpc_combine_morph_D19)
    E_all_S19, Re_all_S19 = np.concatenate(E_combine_morph_S19), np.concatenate(Re_kpc_combine_morph_S19)
    E_all_A3, Re_all_A3 = np.concatenate(E_combine_morph_A3), np.concatenate(Re_kpc_combine_morph_A3)
    
    #ALL the mass and size across 4 studies
    #E_total = np.concatenate((E_all,E_all_SG16,E_all_D19,E_all_S19))
    #Re_total = np.concatenate((Re_all,Re_all_SG16,Re_all_D19,Re_all_S19))    
    
    print('correlation_size_mass', stats.pearsonr(np.log10(E_all), np.log10(Re_all)))

    E_total = np.concatenate((E_all,E_all_A3))
    Re_total = np.concatenate((Re_all,Re_all_A3))
    
    E_total = np.log10(E_total)
    Re_total = np.log10(Re_total)
    
    print('correlation_size_mass_total', stats.pearsonr((E_total), (Re_total)))

    # plot the respective scatter
    # plot the green line
    C = SSort.group_to_bin1D( E_total, Re_total, 
                             [[9.0,9.6,10.2,10.8,11.4],
                              [9.6,10.2,10.8,11.4,12.0]])
    
    # note that in here, y is mass, x is radius    
    ax0.plot(SSort.trim_value(np.nan_to_num(C['median_x']),0),
             SSort.trim_value(np.nan_to_num(C['median_y']),0),
             "o",color='#46f12b',ms=16,zorder=50)  # turn the nan to 0 and trim it
    ax0.errorbar(SSort.trim_value(np.nan_to_num(C['median_x']),0),
                 SSort.trim_value(np.nan_to_num(C['median_y']),0), 
                 yerr = SSort.trim_value(np.nan_to_num(C['std_y']),0),
                #xerr = SSort.trim_value(np.nan_to_num(C['std_x']),0),
                 c='#46f12b',linestyle='dashed',lw=5,elinewidth=5,zorder=50)
   # ax0.plot(mass_barbulge,R_barbulge,'d',ms=14,alpha=0.7)
    
    ax0.plot(np.log10(E_combine_morph[0]), np.log10(Re_kpc_combine_morph[0]),
             'o',ms= ms0,color='k',label = "E+ES",alpha=0.8)
    ax0.plot(np.log10(E_combine_morph[1]), np.log10(Re_kpc_combine_morph[1]),
             'o',ms= ms0,color="#d20d0d",label="S0",alpha=0.8)
    ax0.plot(np.log10(E_combine_morph[2]), np.log10(Re_kpc_combine_morph[2]),
             'o',ms= ms0,color="#ebb800",label="S",alpha=0.8)

    #ax0.plot(E_combine_morph_SG16[0], Re_kpc_combine_morph_SG16[0],'o',ms= ms0,color='k')
    #ax0.plot(E_combine_morph_SG16[1], Re_kpc_combine_morph_SG16[1],'o',ms= ms0,color="#d20d0d")
    #ax0.plot(E_combine_morph_SG16[2], Re_kpc_combine_morph_SG16[2],'o',ms= ms0,color="#ebb800")
    
    #ax0.plot(E_combine_morph_D19[0], Re_kpc_combine_morph_D19[0],'o',ms= ms0,color='k')
    #ax0.plot(E_combine_morph_D19[1], Re_kpc_combine_morph_D19[1],'o',ms= ms0,color="#d20d0d")
    #ax0.plot(E_combine_morph_D19[2], Re_kpc_combine_morph_D19[2],'o',ms= ms0,color="#ebb800")
    
    #ax0.plot(E_combine_morph_S19[0], Re_kpc_combine_morph_S19[0],'o',ms= ms0,color='k')
    #ax0.plot(E_combine_morph_S19[1], Re_kpc_combine_morph_S19[1],'o',ms= ms0,color="#d20d0d")
    #ax0.plot(E_combine_morph_S19[2], Re_kpc_combine_morph_S19[2],'o',ms= ms0,color="#ebb800")    
    
    
    ax0.plot(np.log10(E_combine_morph_A3[0]), np.log10(Re_kpc_combine_morph_A3[0]),
             'o',ms= ms0,color='k',alpha=0.8)
    ax0.plot(np.log10(E_combine_morph_A3[1]), np.log10(Re_kpc_combine_morph_A3[1]),
             'o',ms= ms0,color='#d20d0d',alpha=0.8)
    ax0.plot(np.log10(E_combine_morph_A3[2]), np.log10(Re_kpc_combine_morph_A3[2]),
             'o',ms= ms0,color='#ebb800',alpha=0.8)

    ax0.set_xlim(left = np.log10(xlim[0]), right = np.log10(xlim[1]))
    ax0.set_ylim(bottom = np.log10(ylim[0]), top = np.log10(ylim[1]))
    #ax0.set_xscale('log')
    #ax0.set_yscale('log')
    ax0.set_ylabel(r"$ \log_{10}(R_\mathrm{e,Sph}~(\rm kpc))$", fontsize=16)
    ax0.set_xlabel(r"$ \log_{10}(M_\mathrm{*,sph} / \rm M_\mathrm{\odot})$", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()

    plt.show()

def plot_densmass_morph():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    ax0.scatter(E_1Sersic,E_1Sersic/(np.pi* Re_kpc_1Sersic**2),
                           alpha = 0.4,color='#45b847',
                          label=r"$\rm This~work~(1Sersic)$",s=130)
    ax0.plot(E_combine_morph[0], E_combine_morph[0]/(np.pi*Re_kpc_combine_morph[0]**2),
                           'o',color='k',label = "E+ES",ms=10)
    ax0.plot(E_combine_morph[1], E_combine_morph[1]/(np.pi*Re_kpc_combine_morph[1]**2),
                           'o',color="#d20d0d",label="S0",ms=10)
    ax0.plot(E_combine_morph[2], E_combine_morph[2]/(np.pi*Re_kpc_combine_morph[2]**2),
                           'o',color="#ebb800",label="S",ms=10)
    
    func_cS = SAna.AnalyticFunctions.size_mass_Graham_equ((mass0,np.repeat(np.average(MLR),len(mass0))),
                                                       -2.55158269, -20.94292143,0.10998085, -23.71212895)
    func_S = SAna.AnalyticFunctions.size_mass_Graham_equ((mass0,np.repeat(np.average(MLR),len(mass0))),
                                                       -9, -17.85,6.12, -113)

    #ax0.plot(mass0,func,'r-',label='Graham (2019)',lw = 4)
    ax0.plot(mass0,mass0/(np.pi*func_cS**2),'r-',label='Graham (2019)',lw = 4)
    ax0.plot(mass0,mass0/(np.pi*func_S**2),'r--',label='Graham (2019)',lw = 4)
    ax0.plot(mass0,mass0/(np.pi*func_S**2)+mass0/(np.pi*func_cS**2),'k--',label='Graham (2019)',lw = 4)
   #ax0.plot(mass_barbulge,R_barbulge,'d',ms=14,alpha=0.7)
    
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S_old, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='g',label="S_dustcorr_old",marker='s')
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='y',label="S_dustcorr",marker='s')
    #SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
    #                                                  "Barro et al. 2013", 
    #                                                  alpha0=0,AX=ax0)
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel(r"$ \Sigma_\mathrm{Sph}$", fontsize=16)
    ax0.set_xlabel(r"$ M_\mathrm{*,Sph} / \rm M_\mathrm{\odot} (RC15)$", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()

    plt.show()

def plot_Mag_Re_morph(R,Mag):
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    ax0.plot(R[0], Mag[0],'o',ms = 10, color='k',label = "E+ES")
    ax0.plot(R[1], Mag[1],'o',ms = 10, color="#d20d0d",label="S0")
    ax0.plot(R[2], Mag[2],'o',ms = 10, color="#ebb800",label="S")

    ax0.set_xlabel(r"$ R_\mathrm{e,sph}~(\rm kpc)$", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()
    ax0.set_xscale('log')
    plt.gca().invert_yaxis()

    plt.show()
    
def draw_condense_line(X, AX,errorbar=True,color='b'):
    AX.plot(SSort.trim_value(np.nan_to_num(X['median_y']),0),
             SSort.trim_value(np.nan_to_num(X['median_x']),0),
             "o--",lw=6,
            color=color,ms=12)  # turn the nan to 0 and trim it
    if errorbar == True:
        AX.errorbar(SSort.trim_value(np.nan_to_num(X['median_y']),0),
                 SSort.trim_value(np.nan_to_num(X['median_x']),0), 
                 xerr = SSort.trim_value(np.nan_to_num(X['std_y']),0),
                 c=color,linestyle='dashed',lw=4,elinewidth=4)
    
def plot_Mag_Re_morph_diverse(R1,R2,R3,Mag):
    fig = plt.figure(figsize=(6.4, 4.8))
    ax0 = plt.subplot()
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    ax0.plot(R1[0], Mag[0],'o',ms = 10, color='k',label = "z=0.05")
    ax0.plot(R1[1], Mag[1],'o',ms = 10, color="k")
    ax0.plot(R1[2], Mag[2],'o',ms = 10, color="k")

    ax0.plot(R2[0], Mag[0],'o',ms = 10, color='#d20d0d',label = "z=0.5")
    ax0.plot(R2[1], Mag[1],'o',ms = 10, color="#d20d0d")
    ax0.plot(R2[2], Mag[2],'o',ms = 10, color="#d20d0d")
    
    ax0.plot(R3[0], Mag[0],'o',ms = 10, color='#286ba9',label = "z=0.95")
    ax0.plot(R3[1], Mag[1],'o',ms = 10, color="#286ba9")
    ax0.plot(R3[2], Mag[2],'o',ms = 10, color="#286ba9")
    
    R1_all, R2_all, R3_all = np.concatenate(R1), np.concatenate(R2), np.concatenate(R3)
    Mag_all = np.concatenate(Mag)
    
    C1 = SSort.group_to_bin1D(Mag_all, R1_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    C2 = SSort.group_to_bin1D(Mag_all, R2_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    C3 = SSort.group_to_bin1D(Mag_all, R3_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    
    #print('C1',C1)
    
    # note that in here, y is mass, x is radius
    #print("C",C['median_x'],np.array(C['std_y'])/1e10)
    draw_condense_line(C1,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C2,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C3,ax0,errorbar=False,color='#46f12b')

    ax0.set_xlabel(r"$ R_\mathrm{z,sph}~\rm(kpc)$ ", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=16)

    ax0.set_xscale('log')

    plt.gca().invert_yaxis()
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.show()

def plot_Mag_mu_morph_diverse(mu1,mu2,mu3,Mag):
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    #ax0.plot(mu1[0], Mag[0],'o',ms = 10, color='k',label = "z=0.05")
    #ax0.plot(mu1[1], Mag[1],'o',ms = 10, color="k")
    #ax0.plot(mu1[2], Mag[2],'o',ms = 10, color="k")

    ax0.plot(mu2[0], Mag[0],'o',ms = 10, color='#d20d0d',label = "z=0.5")
    ax0.plot(mu2[1], Mag[1],'o',ms = 10, color="#d20d0d",)
    ax0.plot(mu2[2], Mag[2],'o',ms = 10, color="#d20d0d",)
    
    ax0.plot(mu3[0], Mag[0],'o',ms = 10, color='#286ba9',label = "z=0.95")
    ax0.plot(mu3[1], Mag[1],'o',ms = 10, color="#286ba9")
    ax0.plot(mu3[2], Mag[2],'o',ms = 10, color="#286ba9")
    
    ax0.plot(mu0_combine_morph[0], Mag[0],'o',ms = 10, color='#6833da',label = "z=0.0")
    ax0.plot(mu0_combine_morph[1], Mag[1],'o',ms = 10, color="#6833da")
    ax0.plot(mu0_combine_morph[2], Mag[2],'o',ms = 10, color="#6833da")
    
    mu1_all, mu2_all, mu3_all = np.concatenate(mu1), np.concatenate(mu2), np.concatenate(mu3)
    Mag_all = np.concatenate(Mag)
    
    # grup data points by bins
    C1 = SSort.group_to_bin1D(Mag_all, mu1_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])
    C2 = SSort.group_to_bin1D(Mag_all, mu2_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])
    C3 = SSort.group_to_bin1D(Mag_all, mu3_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                 [-24,-23,-22,-21,-20,-19,-18]])

    # Draw condense line
    # note that in here, y is mass, x is radius
    #print("C",C['median_x'],np.array(C['std_y'])/1e10)
    draw_condense_line(C1,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C2,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C3,ax0,errorbar=False,color='#46f12b')
    
    
    
    ax0.set_xlabel(r"$\mu_\mathrm{z,Sph}$ ", fontsize=18)
    ax0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=18)
    ax0.set_xlim(35,-0.5)
    plt.legend(loc="lower right")
    plt.gca().invert_yaxis()

    plt.tight_layout()

    plt.show()
    
def plot_Mag_diverse_2plot(R1,R2,R3,mu1,mu2,mu3,Mag):
    fig = plt.figure(figsize=(6.4*2, 5.8))
    gs = gridspec.GridSpec(ncols=2, nrows=1,
                               hspace=0.0, wspace=0.0)
    
    # Mag - Sersic n
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(R1[0], Mag[0],'o', ms = 10, color='k', label = "z=0.05")
    ax0.plot(R1[1], Mag[1],'o', ms = 10, color="k")
    ax0.plot(R1[2], Mag[2],'o', ms = 10, color="k")

    ax0.plot(R2[0], Mag[0],'o', ms = 10, color='#d20d0d', label = "z=0.50")
    ax0.plot(R2[1], Mag[1],'o', ms = 10, color="#d20d0d")
    ax0.plot(R2[2], Mag[2],'o', ms = 10, color="#d20d0d")
    
    ax0.plot(R3[0], Mag[0],'o', ms = 10, color='#286ba9', label = "z=0.95")
    ax0.plot(R3[1], Mag[1],'o', ms = 10, color="#286ba9")
    ax0.plot(R3[2], Mag[2],'o', ms = 10, color="#286ba9")
    
    R1_all, R2_all, R3_all = np.concatenate(R1), np.concatenate(R2), np.concatenate(R3)
    Mag_all = np.concatenate(Mag)
    C1_R = SSort.group_to_bin1D(Mag_all, R1_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    C2_R = SSort.group_to_bin1D(Mag_all, R2_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    C3_R = SSort.group_to_bin1D(Mag_all, R3_all, [[-25,-24,-23,-22,-21,-20,-19],
                                                [-24,-23,-22,-21,-20,-19,-18]])
    draw_condense_line(C1_R,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C2_R,ax0,errorbar=False,color='#46f12b')
    draw_condense_line(C3_R,ax0,errorbar=False,color='#46f12b')
   
    ax0.set_xlabel(r"$R_\mathrm{z,Sph}$ ", fontsize=18)
    ax0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=18)
    ax0.set_xscale('log')

    ax0.legend(loc="lower right")
    ax0.tick_params(axis='y',labelsize=15)
    ax0.tick_params(axis='x',labelsize=15)

    # Mag - Sigma0
    ax1 = plt.subplot(gs[1],sharey=ax0)
    
    ax1.plot(mu2[0], Mag[0],'o', ms = 10, color='#d20d0d', label = "z=0.50")
    ax1.plot(mu2[1], Mag[1],'o', ms = 10, color="#d20d0d",)
    ax1.plot(mu2[2], Mag[2],'o', ms = 10, color="#d20d0d",)
    
    ax1.plot(mu3[0], Mag[0],'o', ms = 10, color='#286ba9', label = "z=0.95")
    ax1.plot(mu3[1], Mag[1],'o', ms = 10, color="#286ba9")
    ax1.plot(mu3[2], Mag[2],'o', ms = 10, color="#286ba9")
    
    ax1.plot(mu0_combine_morph[0], Mag[0],'o', ms = 10, color='#6833da', 
             label = "z=0.0")
    ax1.plot(mu0_combine_morph[1], Mag[1],'o', ms = 10, color="#6833da")
    ax1.plot(mu0_combine_morph[2], Mag[2],'o', ms = 10, color="#6833da")
    
    mu1_all, mu2_all, mu3_all = np.concatenate(mu1), np.concatenate(mu2), np.concatenate(mu3)
    Mag_all = np.concatenate(Mag)
    C1_mu = SSort.group_to_bin1D(Mag_all, mu1_all, 
                                 [[-25,-24,-23,-22,-21,-20,-19],
                                  [-24,-23,-22,-21,-20,-19,-18]])
    C2_mu = SSort.group_to_bin1D(Mag_all, mu2_all, 
                                 [[-25,-24,-23,-22,-21,-20,-19],
                                  [-24,-23,-22,-21,-20,-19,-18]])
    C3_mu = SSort.group_to_bin1D(Mag_all, mu3_all, 
                                 [[-25,-24,-23,-22,-21,-20,-19],
                                  [-24,-23,-22,-21,-20,-19,-18]])

    draw_condense_line(C1_mu,ax1,errorbar=False,color='#46f12b')
    draw_condense_line(C2_mu,ax1,errorbar=False,color='#46f12b')
    draw_condense_line(C3_mu,ax1,errorbar=False,color='#46f12b')
    
    ax1.set_xlabel(r"$\mu_\mathrm{z,Sph}$", fontsize=18)
    ax1.set_xlim(34.7,-0.5)

    ax1.legend(loc="lower right")
    ax1.tick_params(axis='x',labelsize=15)

    plt.gca().invert_yaxis()
    plt.setp(ax1.get_yticklabels(), visible=False)

    #plt.tight_layout()

    plt.show()

    return fig
    
def plot_Dens_mass_morph(dens,mass):
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    
    ax0.plot(dens[0], mass[0],'o',ms = 10, color='k',label = "E+ES")
    ax0.plot(dens[1], mass[1],'o',ms = 10, color="#d20d0d",label="S0")
    ax0.plot(dens[2], mass[2],'o',ms = 10, color="#ebb800",label="S")

    ax0.set_xlabel(r"$\Sigma_\mathrm{e,sph}~(\rm kpc)$ ", fontsize=16)
    ax0.set_ylabel(r"$\M_*/\mathrm{M_{\odot}}$", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()
    ax0.set_xscale('log')
    plt.show()
    

def plot_size_Mag_z0():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    
    M_tot = Abs_sph_mag
    Mag00 = np.linspace(-16,-25,50)    

    func_S2 = SAna.AnalyticFunctions.size_mass_Graham_equ(Mag00,
                                                       -8.8, -17.89,2.6, -59+5.2)
    #func_S = SAna.AnalyticFunctions.size_mass_Graham_equ((mass1,np.repeat(1.0,len(mass1))),
    #                                                   -8.75, -17.85,6.12, -93)
    #ax0.plot(Mag00,func_cS,'b-',label='Graham (2019)_g',lw = 4)
    #ax0.plot(Mag00,func_S,'b--',label='Graham (2019)_g',lw = 4)

    #ax0.plot(func_cS2,Mag00,'r-',label='Graham (2019)_f',lw = 4)
    ax0.plot(func_S2,Mag00,'r--',label='Graham (2019)_f',lw = 4)
    
    ax0.plot(Re_kpc,Abs_sph_mag,'o',alpha = 0.4,color='#a5200b',
                          label=r"$\rm This~work$",ms=10)
    
    ax0.set_xlabel(r"$R_\mathrm{e,Sph}\rm ~(kpc)$ ", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_{Sph}$", fontsize=16)
    
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.ylim(0,100)
    #plt.xlim(-26,-16)
    ax0.set_xscale('log')
    plt.show()
    

def plot_size_mass_comparison():
    #fig = plt.figure(figsize=(6.4, 5.8))
    
    fig = plt.figure(figsize=(6.4, 16.9))
    gs = gridspec.GridSpec(ncols=1, nrows=3,
                               hspace=0.0, wspace=0.0)
    ax0 = plt.subplot(gs[0])

    ax0.scatter(nsa_mass_R15BC, nsa_Re, marker='x', color = 'grey', 
              label = r"", s = 16,alpha = 0.2)   
    # plot my line
    func3 = 10**(0.88*mass0-9.15)
    ax0.plot(10**mass0,func3,'k-', lw = 3,label=r'',zorder=50)

    # plot Shen 2003
    ax0.plot(10**mass_shen03,Shen2003,color='purple', linestyle= 'solid', 
             label=r"",lw = 3,zorder=50)
    
    # plot Graham curve
    func_S = SAna.AnalyticFunctions.size_mass_Graham_equ(10**mass0,
                                                         -9.4, -14.3, 
                                                         2.0/3, -29.5)
    ax0.plot(10**mass0,func_S,'g-', label=r'',zorder=50,lw = 3)
    
    # plot Nedkova 2021
    ax0.plot(10**mass_Nedkova21,Nedkova2021_lowz,color='#f2b425', 
             linestyle='solid',label=r"",lw = 3,zorder=60)

    # Baldry 2021
    ax0.plot(10**mass_Baldry21,Baldry2020,'r-',label=r'', lw=3,zorder =60)
    
    # plot Lange 2015 (M>10^10)
    ax0.plot(10**mass_Lange15,Lange2015_ETG_morph, color='blue', 
             linestyle = "solid",
             label=r"",lw = 3,zorder=30) 
    
    #ZOE Capperli
    ax0.plot(10**mass0,ZOE,"-",color='#676d84',label=r"",lw = 3,zorder=50)

    E_all, Re_all = np.concatenate((E_combine_morph)), np.concatenate((Re_kpc_combine_morph)) 
    E_all2, Re_all2 = np.concatenate((E_combine_morph_A3)), np.concatenate((Re_kpc_combine_morph_A3))

    ax0.plot(E_all,Re_all, 'o',color='k', ms = 5,alpha =0.9)
    ax0.plot(E_all2,Re_all2, 'o',color='k', ms = 5,alpha =0.9)
    
    ax0.text(1.3e9,40,r"$\rm Local~Galaxies$", weight = "bold",color="k",
             fontsize=22,zorder=60)

    ax0.text(1.0e11,0.6,r"$\rm Nedkova+21~$", weight = "bold",
             color="#f2b425",fontsize=16,zorder=60)
    ax0.text(1.0e11,0.4,r"$\rm Baldry+21~$", weight = "bold",
             color="r",fontsize=16,zorder=60)
    ax0.text(1.0e11,0.28,r"$\rm Shen+03~$", weight = "bold",
             color="purple",fontsize=16,zorder=60)
    ax0.text(1.0e11,0.2,r"$\rm Lange+15~$", weight = "bold",
             color="b",fontsize=16,zorder=60)
    ax0.text(1.0e11,0.14,r"$\rm Graham+03~$", weight = "bold",
             color="g",fontsize=16,zorder=60)
    ax0.text(1.0e11,0.1,r"$\rm Cappellari+13~(ZOE)$", weight = "bold",
             color="#676d84",fontsize=16,zorder=60)

    ax0.set_xlim(left = xlim[0], right = xlim[1])
    ax0.set_ylim(bottom = ylim[0], top = ylim[1])
       
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel(r"$ R_\mathrm{e}~(\rm kpc)$", fontsize=16)
    
    #ax0.legend(loc="upper right",prop={'size': 14})
    # plot bulge
    ax1 = plt.subplot(gs[1],sharex=ax0)
    ax1.plot(10**mass0,func3,'k-', lw = 3,label=r'',zorder=20)

    func_B = SAna.AnalyticFunctions.size_mass_Graham_equ(10**mass0,
                                                         -7.5, -20, 
                                                         0.6, -29.7,ML=8,A=4) #M_sun=5.09
    
    
    ax1.plot(E_all,Re_all, 'o',color='k', ms = 5,alpha =0.9)
    ax1.plot(E_all2,Re_all2, 'o',color='k', ms = 5,alpha =0.9)
    
    ax1.plot(10**mass_Mendez_Abreu21, Mendez_Abreu2021_bulge_highmass, '-', 
             color='#113ee4',lw = 3 ,  #dashes=(20, 10),
             label='Mendez_Abreu2021_bulge_highmass',
             zorder=70)
    #ax1.plot(10**mass0, Mendez_Abreu2021_bulge_lowmass, linestyle="dashdot", color='#113ee4',lw = 4 )
    
    ax1.plot(10**mass_Dimauro18, Dimauro2018_bulge_lowz1, linestyle="-",
             color='grey',lw = 3, #dashes=(20, 10),
             label='Dimauro2018_bulge_lowz1',
             zorder=60)

    
    ax1.plot(10**mass0,func_B,'g--',lw = 3,dashes=(20, 10))
    ax1.plot(10**mass_Lange16,Lange2016_final_bulge,'-',color='cyan',
             lw = 3)#,dashes=(20, 10))
    
    ax1.text(1.3e9,40,r"$\rm Local~bulges$", weight = "bold",color="k",
             fontsize=20,zorder=60)

    ax1.text(1.0e11,0.28,r"$\rm M\'endez-Abreu+21$", weight = "bold",color="b",
             fontsize=16,zorder=60)
    ax1.text(1.0e11,0.2,r"$\rm Dimauro+18$", weight = "bold",color="grey",
             fontsize=16,zorder=60)
    ax1.text(1.0e11,0.14,r"$\rm Lange+16$", weight = "bold",color="cyan",
             fontsize=16,zorder=60)
    ax1.text(1.0e11,0.1,r"$\rm Graham+13$", weight = "bold",color="g",
             fontsize=16,zorder=60)

    ax1.set_xlim(left = xlim[0], right = xlim[1])
    ax1.set_ylim(bottom = ylim[0], top = ylim[1])
       
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$ R_\mathrm{e}~(\rm kpc)$", fontsize=16)
    #ax1.legend()
    # high-z
    ax2 = plt.subplot(gs[2],sharex=ax0)
    
    # plot my line
    func3 = SAna.AnalyticFunctions.size_mass_powerlaw(10**mass0,0.51/(1e10**0.79),0.79)
    ax2.plot(10**mass0,func3,'k-', lw = 3,label=r'',zorder=45)

    # plot Lange 2016 Bulge
    #ax0.plot(10**mass0,Lange2016_E, '--', lw = 3,
    #         label=r"$\rm E~in~Lange~et~al.~(2016)$")
    #ax0.plot(10**mass0,Lange2016_ETG_bulge, '--', lw =3,
    #         label=r"$\rm ETG~bulge~in~Lange~et~al.~(2016)$")
    #ax0.plot(10**mass0,Lange2016_E_M1e10, '--', lw=3,
    #         label=r"$ \mathrm{E}(M_{*}/\mathrm{M_{\odot}} \gtrsim 10^{10})\rm ~in~Lange~et~al.~(2016)$")
    #ax0.plot(10**mass0,Lange2016_final_bulge, '--', lw=3,
    #         label=r"$ \rm Final~spheroids~in~Lange~et~al.~(2016)$") 
    
    
    # plot Lange 2015 (M>10^10)
    #ax0.plot(10**mass0,Lange2015_ETG_morph,'b--',label=r"$\rm ETGs~in~Lange~et~al.~(2015)$",lw = 5,zorder=50)

    # plot Saracco 2018
    ax2.plot(10**mass1,Saracco2018_high,'-',color='#156769',
             label=r"",lw = 3,zorder=50)
    ax2.plot(10**mass2,Saracco2018_low,'-',color='#156769',
             label=r"",lw = 3,zorder=50)
    
    # plot Nedkova 2021
    ax2.plot(10**mass_Nedkova21,Nedkova2021_highz2,'-',color='#f2b425',
             label=r"",lw = 3, zorder=50,alpha=0.9)
   # ax2.plot(10**mass0,Nedkova2021_highz3,linestyle='dashdot',color='#f2b425',label=r"$\rm QGs~(1.5<z<2.0)~in~Nedkova~et~al.~(2021)$",lw = 3)

    ax2.plot(vdWel_highz_mass,vdWel_highz_size*0.76,'o',ms=3,color='r',
             label=r"$\rm QGs~(z\sim1.25)~in~vdWel~et~al.~(2014)$",lw=3)    
    ax2.plot(vdWel_highz_mass2,vdWel_highz_size2*0.76,'o',ms=3,color='#f0a10f',
                          label=r"$\rm QGs~(z\sim2.25)~in~vdWel~et~al.~(2014)$",
                          zorder = 70,lw=3)

    ax2.plot(E_all,Re_all, 'o',color='k', ms = 5,alpha =0.9)
    ax2.plot(E_all2,Re_all2, 'o',color='k', ms = 5,alpha =0.9)
    
    ax2.text(1.3e9,40,r"$\rm High-z $", weight = "bold",color="k",
             fontsize=20,zorder=60)    
    ax2.text(1.3e9,40-17,r"$\rm Quiescent~Galaxies$", weight = "bold",color="k",
                          fontsize=20,zorder=60)

    ax2.text(1.0e11,0.14,r"$\rm Nedkova+16$", weight = "bold",color="#f2b425",
             fontsize=16,zorder=60)
    ax2.text(1.0e11,0.1,r"$\rm Saracco+17$", weight = "bold",color="#156769",
             fontsize=16,zorder=60)

    ax2.set_xlim(left = xlim[0], right = xlim[1])
    ax2.set_ylim(bottom = ylim[0], top = ylim[1])
       
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$ R_\mathrm{e}~(\rm kpc)$", fontsize=16)
    ax2.set_xlabel(r"$ M_\mathrm{*} / \rm M_\mathrm{\odot}$", fontsize=16)
    
    ax2.legend(loc="upper right",prop={'size': 10})
    plt.tight_layout()
    
    #twin0=ax0.twinx()
    #
    ##twin0.scatter([],[],label=r"$\rm SDSS~galaxies $", color ='grey', 
    ##           marker ="x",s=50)
    #
    #twin0.plot([],[],label=r"$\rm This~work $", color ='#2a3236',  marker="",
    #           linestyle ="-",lw=7)
    #twin0.scatter([],[],label=r"$\rm Our~Spheroids$", color ='k', 
    #           marker ="o",s=50)
    #
    #twin0.legend(loc='upper center',fontsize=13, bbox_to_anchor=(0.5, 1.1),
    #      fancybox=True, shadow=False, ncol=3)
    #plt.setp(twin0.get_yticklabels(), visible=False)

    plt.show()
    return fig



def plot_3D_para():
    fig = plt.figure(figsize=(6.4, 5.8))

    ax0 = plt.axes(projection='3d')
    
    xdata = np.concatenate(n_combine_morph)
    ydata = np.concatenate(mu0_combine_morph)
    zdata = np.concatenate(Re_kpc_combine_morph)
    
    #print(len(xdata))
    
    cdata = np.concatenate(Mag_combine_morph)
    ax0.scatter3D(np.log10(xdata), ydata, np.log10(zdata), c=cdata, cmap='Greens')
    
    ax0.set_xlabel('Sersic n ')
    ax0.set_ylabel('mu0')
    ax0.set_zlabel('Re')
    plt.show()

    return ax0

def calculate_sigma_0(Dist,mu,s,MLR,M_o):
    
    DM = 25+5*np.log10(Dist) # Dist in Mpc
    #s is in pc arcsec^-1
    A = mu-DM-M_o -2.5*np.log10(1/(s**2))-2.5*np.log10(MLR)
    
    Sigma_R = 10**(A/-2.5)
    
    return Sigma_R

def calculate_average_mu(mu_z,Sersic_n,bnz,frac_z):
    
    frac_z = np.repeat(frac_z,len(mu_z))
    print('nnn',len(Sersic_n))
    term_a = mu_z
    term_b = -2.5*np.log10(((
            2*frac_z*Sersic_n*np.exp(bnz))/((
                bnz)**(2*Sersic_n)))*gamma(2*Sersic_n))

    average_mu = term_a + term_b
    
    return average_mu, term_a, term_b

def plot_mass_n_sigma0_2plot(mass1,n1,sigma1,
                             mass2,n2,sigma2,
                             mass3=[],n3=[],sigma3=[],
                             mass4=[],n4=[],sigma4=[]):
    
    mass1, mass2, mass3, mass4 = np.log10(mass1), np.log10(mass2), np.log10(mass3), np.log10(mass4)
    sigma1, sigma2, sigma3, sigma4 = np.log10(sigma1), np.log10(sigma2), np.log10(sigma3), np.log10(sigma4)
    n1, n2, n3, n4 = np.log10(n1), np.log10(n2), np.log10(n3), np.log10(n4)

    # seperate Sigma
    sigma1_morph = seperate_morph_simple(sigma1,index, morph)
    sigma2_morph = seperate_morph_simple(sigma2,index_SG16, Savorgnan_morph)
    sigma3_morph = seperate_morph_simple(sigma3,index_D19, Davis_morph)
    sigma4_morph = seperate_morph_simple(sigma4,index_S19, Sahu_morph)

    # seperate mass 
    mass1_morph = seperate_morph_simple(mass1,index, morph)
    mass2_morph = seperate_morph_simple(mass2,index_SG16, Savorgnan_morph)
    mass3_morph = seperate_morph_simple(mass3,index_D19, Davis_morph)
    mass4_morph = seperate_morph_simple(mass4,index_S19, Sahu_morph)

    # seperate n
    n1_morph = seperate_morph_simple(n1,index, morph)
    n2_morph = seperate_morph_simple(n2,index_SG16, Savorgnan_morph)
    n3_morph = seperate_morph_simple(n3,index_D19, Davis_morph)
    n4_morph = seperate_morph_simple(n4,index_S19, Sahu_morph)
    
    # reassemble
    sigma_disc = np.concatenate((sigma1_morph[1], sigma1_morph[2],
                                            sigma2_morph[1], sigma2_morph[2],
                                            sigma3_morph[1], sigma3_morph[2], 
                                            sigma4_morph[1], sigma4_morph[2]))
    sigma_E = np.concatenate((sigma1_morph[0], sigma2_morph[0], 
                                         sigma3_morph[0], sigma4_morph[0]))
 
    mass_disc = np.concatenate((mass1_morph[1], mass1_morph[2],
                                            mass2_morph[1], mass2_morph[2],
                                            mass3_morph[1], mass3_morph[2], 
                                            mass4_morph[1], mass4_morph[2]))

    mass_E = np.concatenate((mass1_morph[0], mass2_morph[0], 
                                         mass3_morph[0], mass4_morph[0]))
    
    n_disc = np.concatenate((n1_morph[1], n1_morph[2],
                                            n2_morph[1], n2_morph[2],
                                            n3_morph[1], n3_morph[2], 
                                            n4_morph[1], n4_morph[2]))

    n_E = np.concatenate((n1_morph[0], n2_morph[0], 
                                         n3_morph[0], n4_morph[0]))
    
    # plotting
    fig = plt.figure(figsize=(8*2, 6.4))
    gs = gridspec.GridSpec(ncols=2, nrows=1,
                               hspace=0.2, wspace=0.0) 
    
    axt0 = plt.subplot(gs[0])
    axt0.plot(n1,mass1,'o', ms = 10,
              color='#a5200b',label=r"$\rm Hon~et~al.~(2022)$")
    axt0.plot(n2,mass2,'o', ms = 10,
              color='#b940c8',label=r"$\rm Savorgnan~&~Graham~(2016)$")
    axt0.plot(n3,mass3,'o', ms = 10,
              color='#2e417b',label=r"$\rm Davis~et~al.~(2019)$") 
    axt0.plot(n4,mass4,'o', ms = 10,
              color='#e1a000',label=r"$\rm Sahu~et~al.~(2019)$")
    
    # Pearson Correlation
    r_all = stats.pearsonr((np.concatenate((n1, n2, n3, n4))),
                   (np.concatenate((mass1,mass2,mass3,mass4))))[0]
    
    r_E = stats.pearsonr(n_E,mass_E)[0]

    
    r_disc = stats.pearsonr(n_disc,mass_disc)[0]
    
    # Spearman Correlation
    rs_all = stats.spearmanr((np.concatenate((n1, n2, n3, n4))),
                   (np.concatenate((mass1,mass2,mass3,mass4))))[0]
    
    rs_E = stats.spearmanr(n_E,mass_E)[0]

    
    rs_disc = stats.spearmanr(n_disc,mass_disc)[0]   
    
    #    np.log10(np.concatenate((n_combine_morph[1],n_combine_morph[2],
    #                    n_combine_morph_A3[1],n_combine_morph_A3[2]))),
    #                np.log10(np.concatenate((
    #                    E_combine_morph[1],E_combine_morph[2],
    #                                E_combine_morph_A3[1],E_combine_morph_A3[2]
    #                                ))))[0]
    axt0.text(-0.12,11.2,r"$\bf S0+S$",color='red',fontsize=18)
    axt0.text(0.42,12.1,r"$\bf E+ES$",color='black',fontsize=18)
    
    axt0.text(np.log10(5),np.log10(1e10), r"$\bf r_p(All)={:.2f}$".format(r_all),fontsize=18)
    axt0.text(np.log10(5),np.log10(1e10)-0.2, r"$ r_p(E+ES)={:.2f}$".format(r_E),fontsize=18)
    axt0.text(np.log10(5),np.log10(1e10)-0.4, r"$ r_p(S0+S)={:.2f}$".format(r_disc),fontsize=18)
    
    axt0.text(np.log10(5),np.log10(2e9), r"$\bf r_s(All)={:.2f}$".format(rs_all),fontsize=18)
    axt0.text(np.log10(5),np.log10(2e9)-0.2, r"$ r_s(E+ES)={:.2f}$".format(rs_E),fontsize=18)
    axt0.text(np.log10(5),np.log10(2e9)-0.4, r"$r_s(S0+S)={:.2f}$".format(rs_disc),fontsize=18)

    SPlot.Plot2D.confidence_ellipse(n_E, mass_E, 
                                    axt0, n_std=2.0,facecolor='None',
                                    edgecolor='black')
    
    SPlot.Plot2D.confidence_ellipse(n_disc, mass_disc, 
                                    axt0, n_std=2.0,facecolor='None',
                                    edgecolor='red')
    
    axt0.set_ylabel(r"$\mathrm{log_{10}}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    axt0.set_xlabel(r"$\mathrm{log_{10}}(\mathrm{Sersic}~n_\mathrm{Sph})$", fontsize=18)
    #axt0.set_xscale('log')
    #axt0.set_yscale('log')
    
    axt0.set_xlim(left=np.log10(0.4),right=np.log10(16))
    axt0.set_ylim(bottom=np.log10(3e8),top=np.log10(4e12))
    axt0.legend(loc="upper left",fontsize=10)

    # Mass vs Sigma0 
    axt1 = plt.subplot(gs[1])

    axt1.plot(sigma1,mass1,'o', ms = 10,
              color='#a5200b',label=r"$\rm Hon~et~al.~(2022)$")
    axt1.plot(sigma2,mass2,'o', ms = 10,
              color='#b940c8',label=r"$\rm Savorgnan~&~Graham~(2016)$")
    axt1.plot(sigma3,mass3,'o', ms = 10,
              color='#2e417b',label=r"$\rm Davis~et~al.~(2019)$")
    axt1.plot(sigma4,mass4,'o', ms = 10,
              color='#e1a000',label=r"$\rm Sahu~et~al.~(2019)$")    
    
    # Pearson Correlation
    r_all = stats.pearsonr(np.concatenate((sigma1, sigma2, 
                                                    sigma3, sigma4)),
                   np.concatenate((mass1,mass2,mass3,mass4)))[0]
    
    r_E = stats.pearsonr(sigma_E,mass_E)[0]
    r_disc = stats.pearsonr(sigma_disc,mass_disc)[0]
    
    # Spearman Correlation
    rs_all = stats.spearmanr(np.concatenate((sigma1, sigma2, 
                                                    sigma3, sigma4)),
                   np.concatenate((mass1,mass2,mass3,mass4)))[0]
    
    rs_E = stats.spearmanr(sigma_E,mass_E)[0]
    rs_disc = stats.spearmanr(sigma_disc,mass_disc)[0]

    axt1.text(2.37,11.0,r"$\bf S0+S$",color='red',fontsize=18)
    axt1.text(3.0,12,r"$\bf E+ES$",color='black',fontsize=18)
    
    axt1.text(np.log10(3.2e7),np.log10(1e10), r"$\bf r_p(All)={:.2f}$".format(r_all),fontsize=18)
    axt1.text(np.log10(3.2e7),np.log10(1e10)-0.2, r"$ r_p(E+ES)={:.2f}$".format(r_E),fontsize=18)
    axt1.text(np.log10(3.2e7),np.log10(1e10)-0.4, r"$r_p(S0+S)={:.2f}$".format(r_disc),fontsize=18)

    axt1.text(np.log10(3.2e7),np.log10(2e9), r"$\bf r_s(All)={:.2f}$".format(rs_all),fontsize=18)
    axt1.text(np.log10(3.2e7),np.log10(2e9)-0.2, r"$ r_s(E+ES)={:.2f}$".format(rs_E),fontsize=18)
    axt1.text(np.log10(3.2e7),np.log10(2e9)-0.4, r"$r_s(S0+S)={:.2f}$".format(rs_disc),fontsize=18)
    
    #axt1.set_ylabel(r"$\mathrm{log_{10}}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    axt1.set_xlabel(r"$\mathrm{log_{10}}(\Sigma_\mathrm{0,Sph})$", fontsize=18)
    #axt1.set_xscale('log')
    #axt1.set_yscale('log')
  
    axt1.set_xlim(left=np.log10(200),right=np.log10(3e10))
    axt1.set_ylim(bottom=np.log10(3e8),top=np.log10(4e12))
    
    #axt2 = plt.subplot(gs[2])

    #axt2.set_xlabel(r"$\mathrm{log_{10}}(R_\mathrm{e,Sph})$", fontsize=18)

       
    SPlot.Plot2D.confidence_ellipse(sigma_E, mass_E, 
                                    axt1, n_std=2.0,facecolor='None',
                                    edgecolor='black')
    
    SPlot.Plot2D.confidence_ellipse(sigma_disc, mass_disc, 
                                    axt1, n_std=2.0,facecolor='None',
                                    edgecolor='red')

    plt.setp(axt1.get_yticklabels(), visible=False)
    #plt.setp(axt2.get_yticklabels(), visible=False)
    #plt.tight_layout()
    #axt1.contour((sigma1), (mass1), h, 20, cmap='RdGy');

    plt.show()
 
    return fig

from scipy.special import gamma
from scipy.special import gammainc


def plot_average_mu_2plot(mu1,n1,mass1,bnz,z=0.5):
    
    term_a = mu1
    
    term_b = -2.5*np.log10(((2*z*n1*np.exp(bnz))/((bnz)**(2*n1)))*gamma(2*n1))
    
    #print(term_b)
    fig = plt.figure(figsize=(6.4, 6.4*2.5))
    gs = gridspec.GridSpec(ncols=1, nrows=2,
                               hspace=0.2, wspace=0.0) 
    
    r_a = stats.pearsonr(term_a,np.log10(mass1))
    r_b = stats.pearsonr(term_b,np.log10(mass1))
    r_combine = stats.pearsonr(term_a+term_b,np.log10(mass1))

    print('ra,rb,rcombine',r_a,r_b,r_combine)
    
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(term_b,np.log10(mass1),'bo')
    ax0.plot(term_a,np.log10(mass1),'ro')
    
    ax0.set_xlabel(
        r"$(b)-2.5\mathrm{log_{10}}(\frac{2ze^{b_{n,z}}}{(b_{n,z})^{2n}}\Gamma(2n))$", 
        fontsize=18)
    
    ax0.set_ylabel(
        r"$\mathrm{log_{10}}(M_\mathrm{*,Sph}/\rm M_{\odot})$", 
        fontsize=18)
    
    #ax0.set_yscale('log')
    
    ax1 = plt.subplot(gs[1])

    
    ax1.plot(term_a+term_b,np.log10(mass1),'ko')
    ax1.set_xlabel(
        r"$\langle\Sigma\rangle_{z}$", 
        fontsize=18)
    
    ax1.set_ylabel(
        r"$\mathrm{log_{10}}(M_\mathrm{*,Sph}/\rm M_{\odot})$", 
        fontsize=18)
    plt.show()

    #Plot mass vs 
    return None

def plot_corr_trend(Mass,R_z,Sigma_z):
    #calculate r_p
    frac_z = [5,10,20,30,40,50,60,70,80,90]
    frac_z0 = [0,5,10,20,30,40,50,60,70,80,90,95]
    frac_z1 = [5,10,20,30,40,50,60,70,80,90,95]
    
    r_p_list, r_p_err_list = [], []
    r_s_list, r_s_err_list = [], []
    for i in range(len(Sigma_z_all)):
        r_p = stats.pearsonr(np.log10(Sigma_z_all[i]),np.log10(Mass))
        r_s = stats.spearmanr(np.log10(Sigma_z_all[i]),np.log10(Mass))
        
        r_p_list.append(r_p[0])
        r_s_list.append(r_s[0])
        
        r_p_err_list.append(r_p[1])
    
    r_p_list2, r_p_err_list2 = [], []
    r_s_list2, r_s_err_list2 = [], []
    for i in range(len(R_z)):
        print("R_z[i]",i,np.median(R_z[i]))
        r_p = stats.pearsonr(np.log10(R_z[i]),np.log10(Mass))
        r_s = stats.spearmanr(np.log10(R_z[i]),np.log10(Mass))
        
        r_p_list2.append(r_p[0])
        r_s_list2.append(r_s[0])
        
        r_p_err_list2.append(r_p[1])
        
    r_p_list3, r_p_err_list3 = [], []
    r_p_list4, r_p_err_list4 = [], []
    r_p_list5, r_p_err_list5 = [], []
    
    r_s_list3, r_s_err_list3 = [], []
    r_s_list4, r_s_err_list4 = [], []
    r_s_list5, r_s_err_list5 = [], []
    
    for i in range(len(mu_z_all)):
        term_a = mu_z_all[i]
        term_b = -2.5*np.log10(((
            2*frac_z1[i]*Sersic_n_master*np.exp(bnz_all[i]))/((
                bnz_all[i])**(2*Sersic_n_master)))*gamma(2*Sersic_n_master))

        average_mu = term_a + term_b
        average_sigma = calculate_sigma_0(Dist_master,average_mu,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
        print(len(term_a),len(term_b),len(average_sigma))
        
        #print('average_sigma',average_sigma)
        r_p = stats.pearsonr(np.log10(average_sigma),np.log10(Mass))
        r_p_b = stats.pearsonr(term_b,np.log10(Mass))
        r_p_a = stats.pearsonr(term_a,np.log10(Mass))
        
        r_s = stats.spearmanr(np.log10(average_sigma),np.log10(Mass))
        r_s_b = stats.spearmanr(term_b,np.log10(Mass))
        r_s_a = stats.pearsonr(term_a,np.log10(Mass))

        r_p_list3.append(r_p[0])
        r_s_list3.append(r_s[0])

        r_p_err_list3.append(r_p[1])
        
        r_p_list4.append(r_p_b[0])
        r_s_list4.append(r_s_b[0])
        
        r_p_list5.append(r_p_a[0])
        r_s_list5.append(r_s_a[0])

    #calculate r_s
    
    #plot panels
    fig = plt.figure(figsize=(6.4,4*2))
    gs = gridspec.GridSpec(ncols=1, nrows=2,
                               hspace=0.0, wspace=0.0) 
    ax0 = plt.subplot(gs[0])
    #frac_z = [5,10,20,30,40,50,60,70,80,90,95,100]
    
    ax0.plot(frac_z0,r_p_list,'o-',label=r'$\log_{10}(M_\mathrm{*,Sph})$-$\log_{10}(\Sigma_z)$')
    ax0.plot(frac_z,r_p_list2,'o-',label=r'$\log_{10}(M_\mathrm{*,Sph})$-$\log_{10}(R_z)$')
    ax0.plot(frac_z1,r_p_list3,'o-',label=r'$\log_{10}(M_\mathrm{*,Sph})$-$\log_{10}(\langle\Sigma\rangle_z)$')
    ax0.plot(frac_z1,r_p_list5,'o-',label=r'$\log_{10}(M_\mathrm{*,Sph})$-$\mu_z$')
    ax0.plot(frac_z1,r_p_list4,'o-',label=r'$\log_{10}(M_\mathrm{*,Sph})$-$-2.5\log_{10}(B_z(n))$')

    
    ax0.hlines(0,0,100,linestyle='dashed',color='k',lw=3)
    
    ax0.legend(loc="center right",fontsize=12)
    
    #ax0.set_xlabel(
    #    r"$\mathrm{fraction}~z$", 
    #    fontsize=18)
    
    ax0.set_ylabel(
        r"$r_\mathrm{p}$", 
        fontsize=18)
    
    #ax0.set_yscale('log')
    
    ax1 = plt.subplot(gs[1],sharex=ax0)

    ax1.plot(frac_z0,r_s_list,'o-',label=r'$M_*$-$\Sigma_z$')
    ax1.plot(frac_z,r_s_list2,'o-',label=r'$M_*$-$R_z$')
    ax1.plot(frac_z1,r_s_list3,'o-',label=r'$M_*$-$\langle\Sigma\rangle_z$')
    ax1.plot(frac_z1,r_s_list5,'o-',label=r'$M_*$-$term~a$')
    ax1.plot(frac_z1,r_s_list4,'o-',label=r'$M_*$-$term~b$')

    ax1.hlines(0,0,100,linestyle='dashed',color='k',lw=3)

    ax1.set_xlabel(
        r"$\mathrm{fraction}~z$", 
        fontsize=18)
    
    ax1.set_ylabel(
        r"$r_s$", 
        fontsize=18)
    plt.setp(ax0.get_xticklabels(), visible=False)

    plt.show()
    return fig


def plot_Sigma_z():
    fig = plt.figure(figsize=(6.4,4.8))
    
    A1 = calculate_average_mu(mu_05,Sersic_n_master,
                            bn_05,0.05)
    A2 = calculate_average_mu(mu_50,Sersic_n_master,
                            bn_50,0.5)
    A3 = calculate_average_mu(mu_95,Sersic_n_master,
                            bn_95,0.95) 
    
    
    A1_1 = calculate_sigma_0(Dist_master,A1[0],scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
    A2_1 = calculate_sigma_0(Dist_master,A2[0],scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
    A3_1 = calculate_sigma_0(Dist_master,A3[0],scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
    #A3[0][147] = np.median(A3[0])
    print("Before",A3[0][147])
    A3_2= SSort.replace_outliers(A3[0],1.0,condition="<",replace='median')
    print("After",A3_2[147])

    #plt.plot(np.log10(A1[0]),np.log10(mass_master),'ro',label="z=0.05")
    #plt.plot(np.log10(A2[0]),np.log10(mass_master),'bo',label='z=0.5')
    #plt.plot(np.log10(A3[0]),np.log10(mass_master),'ko',label='z=0.80')
    
    C1_A1 = SSort.group_to_bin1D(np.log10(mass_master), np.log10(A1_1), 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    C2_A2 = SSort.group_to_bin1D(np.log10(mass_master), np.log10(A2_1), 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    C3_A3 = SSort.group_to_bin1D(np.log10(mass_master), np.log10(A3_1), 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    
    #print('C1',C1)
    plt.plot(np.log10(A1_1),np.log10(mass_master),'ko',label="z=0.05")
    plt.plot(np.log10(A2_1),np.log10(mass_master),'ro',label='z=0.5')
    plt.plot(np.log10(A3_1),np.log10(mass_master),'bo',label='z=0.95')
    
    # note that in here, y is mass, x is radius
    #print("C",C['median_x'],np.array(C['std_y'])/1e10)
    draw_condense_line(C1_A1,plt,errorbar=False,color='#46f12b')
    draw_condense_line(C2_A2,plt,errorbar=False,color='#46f12b')
    draw_condense_line(C3_A3,plt,errorbar=False,color='#46f12b')


    print('check')
    print(stats.pearsonr(np.nan_to_num(np.log10(A1[0])),np.log10(mass_master)))
    print(stats.pearsonr(np.nan_to_num(np.log10(A2[0])),np.log10(mass_master)))
    print(stats.pearsonr(np.nan_to_num(np.log10(A3[0])),np.log10(mass_master)))
    
    print(stats.pearsonr(np.nan_to_num(np.log10(A1_1)),np.log10(mass_master)))
    print(stats.pearsonr(np.nan_to_num(np.log10(A2_1)),np.log10(mass_master)))
    print(stats.pearsonr(np.nan_to_num(np.log10(A3_1)),np.log10(mass_master)))
    #print(A3[0])
   # 
    #boolArr = (A3[0] <1.0)
    #print(
    #    np.where(boolArr))
    plt.xlim(-0.8, 5.5)
    plt.ylim(8.8, 12.5)
        
    plt.xlabel(r"$\log_{10}(\langle\Sigma\rangle_z/\rm M_{\odot}pc^{-2})$", fontsize=18)
    plt.ylabel(r"$\log_{10}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    plt.legend()
    plt.show()


def plot_term_b(n1):
    fig = plt.figure(figsize=(6.4,4.8))
    
    z1, z2, z3 = 0.05, 0.5, 0.95
    bnz1, bnz2, bnz3 = bn_05, bn_50, bn_95
    
    A1_1_B = -2.5*np.log10(((2*z1*n1*np.exp(bnz1))/((bnz1)**(2*n1)))*gamma(2*n1))
    A2_1_B = -2.5*np.log10(((2*z2*n1*np.exp(bnz2))/((bnz2)**(2*n1)))*gamma(2*n1))
    A3_1_B = -2.5*np.log10(((2*z3*n1*np.exp(bnz3))/((bnz3)**(2*n1)))*gamma(2*n1))
    #A3[0][147] = np.median(A3[0])


    #plt.plot(np.log10(A1[0]),np.log10(mass_master),'ro',label="z=0.05")
    #plt.plot(np.log10(A2[0]),np.log10(mass_master),'bo',label='z=0.5')
    #plt.plot(np.log10(A3[0]),np.log10(mass_master),'ko',label='z=0.80')
    
    C1_B1 = SSort.group_to_bin1D(np.log10(mass_master), A1_1_B, 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    C2_B2 = SSort.group_to_bin1D(np.log10(mass_master), A2_1_B, 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    C3_B3 = SSort.group_to_bin1D(np.log10(mass_master), A3_1_B, 
                                 [[9.5,10,10.5,11,11.5],
                                                [10,10.5,11,11.5,12]])
    
    #print('C1',C1)
    plt.plot(A1_1_B,np.log10(mass_master),'ko',label="z=0.05")
    plt.plot(A2_1_B,np.log10(mass_master),'ro',label='z=0.50')
    plt.plot(A3_1_B,np.log10(mass_master),'bo',label='z=0.95')
    
    # note that in here, y is mass, x is radius
    #print("C",C['median_x'],np.array(C['std_y'])/1e10)
    draw_condense_line(C1_B1,plt,errorbar=False,color='#46f12b')
    draw_condense_line(C2_B2,plt,errorbar=False,color='#46f12b')
    draw_condense_line(C3_B3,plt,errorbar=False,color='#46f12b')
    
    print("np.log10(A2_1_B)",A2_1_B/-2.5)
    print('rp_Bn',stats.pearsonr(np.nan_to_num(np.log10(A2_1_B/-2.5)),np.log10(mass_master)))
    print('rs_Bn',stats.spearmanr(np.nan_to_num(np.log10(A2_1_B/-2.5)),np.log10(mass_master)))
    
    
    plt.xlim(-5.4,0.5)
    plt.ylim(8.6,12.2)
    
    plt.xlabel(r"$-2.5\log_{10}(B_z(n))$", fontsize=18)
    plt.ylabel(r"$\log_{10}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    plt.legend(loc='lower left')
    plt.show()
    return fig

def plot_Bn_n(n1):
    fig = plt.figure(figsize=(6.4,4.8))
    
    z1, z2, z3 = 0.05, 0.5, 0.95
    bnz1, bnz2, bnz3 = bn_05, bn_50, bn_95
    
    A1_1_B = np.log10(((2*z1*n1*np.exp(bnz1))/((bnz1)**(2*n1)))*gamma(2*n1))
    A2_1_B = np.log10(((2*z2*n1*np.exp(bnz2))/((bnz2)**(2*n1)))*gamma(2*n1))
    A3_1_B = np.log10(((2*z3*n1*np.exp(bnz3))/((bnz3)**(2*n1)))*gamma(2*n1))
    
    plt.plot(10**np.log10(n1), 10**A1_1_B,'ko',label="z=0.05")
    plt.plot(10**np.log10(n1), 10**A2_1_B,'ro',label="z=0.5")
    plt.plot(10**np.log10(n1), 10**A3_1_B,'bo',label="z=0.95")

   # draw_condense_line(C2_B2,plt,errorbar=False,color='#46f12b')
   
    plt.ylim(0,40)
    plt.xlim(-0.2,12)
    
    #plt.yscale('log')
    #plt.xscale('log')
    plt.ylabel(r"$B_z(n)$", fontsize=18)
    plt.xlabel(r"$n_\mathrm{Sph}$", fontsize=18)
    plt.legend(loc='upper left')
    return fig

def plot_Ltot_Rz(n1):
    fig = plt.figure(figsize=(6.4,4.8))

    L_tot1
#Re_kpc_05

    plt.plot(Re_kpc_05, 10**A1_1_B,'ko',label="z=0.05")
    plt.plot(Re_kpc_05, 10**A2_1_B,'ro',label="z=0.5")
    plt.plot(Re_kpc_95, 10**A3_1_B,'bo',label="z=0.95")


    plt.ylabel(r"$L_{tot}(n)$", fontsize=18)
    plt.xlabel(r"$n_\mathrm{Sph}$", fontsize=18)
    plt.legend(loc='upper left')
    return fig


def plot_n_mass_z0():
    fig = plt.figure(figsize=(6.4,5.8))
    
    z0=0.5
    n1,n2,n3,n4 = np.log10(Sersic_n), np.log10(Savorgnan_Sersic_n), np.log10(Sahu_Sersic_n), np.log10(Davis_Sersic_n)

    mass1, mass2, mass3, mass4 = np.log10(E_T11_K), np.log10(Savorgnan_mass_36), \
                         np.log10(Sahu_mass_36), np.log10(Davis_mass_36)
    
    plt.scatter(mass1, n1, s=100,
                           alpha = 0.8,marker='o',color='#a5200b',
                          label=r"$\rm Hon~et~al.~(2022)$")
    plt.scatter(mass2, n2, s=100,
                           alpha = 0.8, marker='o', color='#b940c8',
                           label=r"$\rm Savorgnan~&~Graham~(2016)$")
    plt.scatter(mass4, n4, s=100,
                          alpha = 0.8, marker='o', color='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$")
    plt.scatter(mass3, n3, s=100,
                           alpha = 0.8, marker='o', color='#e1a000',
                         label=r"$\rm Sahu~et~al.~(2019)$")
    
    
    #combine them all together
    E_all  = np.concatenate((mass1,mass2,mass3,mass4))
    n_all = np.concatenate((n1,n2,n3,n4))
    
    # print correlation coefficient
    print('rp_n',stats.pearsonr(E_all,n_all))
    print('rs_n',stats.spearmanr(E_all,n_all))

    print('rp_n',stats.pearsonr(n_all,E_all))
    print('rs_n',stats.spearmanr(n_all,E_all))
    # run fit
    lm = linmix.LinMix(E_all, n_all, K=2)
    lm.run_mcmc(silent=True)
    
    for i in range(0, len(lm.chain), 25):
        xs = mass0
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        plt.plot(xs, ys, color='r', alpha=0.02,lw=3)
        
    alpha = np.median(lm.chain[:]['alpha'])
    beta = np.median(lm.chain[:]['beta'])
    sigsqr = np.median(lm.chain[:]['sigsqr'])
    
    alpha_std = np.median(lm.chain[:]['alpha'].std())
    beta_std = np.median(lm.chain[:]['beta'].std())
    
    # vertical 
    ys_v = alpha + E_all * beta
    # horizontal
    xs_h = (n_all - alpha)/beta
    # calculate the root mean square
    delta_rms1 = np.sqrt(sum((ys_v - n_all)**2)/len(ys_v))
    delta_rms1_x = np.sqrt(sum((xs_h - E_all)**2)/len(xs_h))
    
    print("chain")
    print(alpha,beta)
    print("FIT_para",alpha,beta,sigsqr)
    
    ys = alpha + xs * beta
    plt.plot(xs, ys,linestyle='dashdot', color='k',lw=3)
    
    # run fit again
    lm2 = linmix.LinMix(n_all, E_all, K=2) # x and y flipped
    lm2.run_mcmc(silent=True)
    
    for i in range(0, len(lm2.chain), 25):
        ys2 = mass0
        xs2 = (ys2 - lm2.chain[i]['alpha'])/ lm2.chain[i]['beta']
        plt.plot(ys2,xs2, color='cyan', alpha=0.02,lw=3)
        
    alpha2 = np.median(lm2.chain[:]['alpha'])
    beta2 = np.median(lm2.chain[:]['beta'])
    sigsqr2 = np.median(lm2.chain[:]['sigsqr'])
    
    alpha2_std = np.std(-1.0*lm2.chain[:]['alpha']/lm2.chain[:]['beta'])
    beta2_std = np.std(1/lm2.chain[:]['beta'])

    print("chain")
    print(alpha2,beta2)
    print("chain",lm2.chain[0],lm2.chain[-1])
    print("FIT_para",alpha2,beta2,sigsqr2)
    print("FIT_para_std",np.median(lm2.chain[:]['alpha'].std()),
          np.median(lm2.chain[:]['alpha'].std()))
        
    # vertical 
    xs_v2 = (E_all - alpha2)/ beta2 # this is y (n)
    # horizontal 
    ys_h2 = alpha2 + beta2* n_all # this is x (mass)
    # calculate the root mean square
    delta_rms2 = np.sqrt(sum((ys_h2- E_all)**2)/len(ys_h2))
    delta_rms2_y = np.sqrt(sum((xs_v2- n_all)**2)/len(xs_v2))
    
    xs2 = (ys2 - alpha2)/ beta2
    plt.plot(ys2, xs2, '--', color='k',lw=3) # Note that x and y are reversed
    
    
    m_C, k_C = find_bisector_line(1/beta2,beta,-1.0*alpha2/beta2,alpha)
    
    xs3 = mass0
    ys3 = m_C * xs3 + k_C
    
    # vertical 
    ys_v3 = k_C + E_all * m_C
    # hortizontal
    xs_h3 = (n_all - k_C) / m_C
    
    # calculate the root mean square
    delta_rms3 = np.sqrt(sum((ys_v3- n_all)**2)/len(ys_v3))
    delta_rms3_x = np.sqrt(sum((xs_h3- E_all)**2)/len(xs_h3))
    
    plt.plot(xs3, ys3, '-', color='k',lw=3) # Note that x and y are reversed

    print("k_C, m_C", k_C, m_C)
    print("delta_rms1_x",delta_rms1_x, 'delta_rms2_y', delta_rms2_y,
          'delta_rms3_x', delta_rms3_x)
    
    
    plt.hlines(1.36+0.01,9.78,10.15,linestyle='dashdot',color='k',lw=4)
    plt.hlines(1.36-0.22+0.01,9.78,10.15,linestyle='dashed',color='k',lw=4)
    plt.hlines(1.36-0.46+0.01,9.78,10.15,linestyle='solid',color='k',lw=4)
    
    plt.text(9.0,1.36,
             r"$n_{\rm Sph}(M_{\rm *,Sph}))$",fontsize=12)
    plt.text(9.0,1.36-0.08,r"$S={}\pm{},int.={}\pm{},$".format(round(beta,2),
                                                        round(beta_std,2),
                                                        round(alpha,2),
                                                        round(alpha_std,2)),
             fontsize=12)
    #plt.text(9.0,1.36-0.16,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format(
    #    "{rms}",round(delta_rms1,2),"{rms}",round(delta_rms1_x,2)),
    #         fontsize=12)
    plt.text(9.0,1.36-0.16,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms1,2)),
             fontsize=12)
    
    # X as independent variable
    plt.text(9.0,1.36-0.22,
             r"$M_{*,Sph}(n_{\rm Sph})$",fontsize=12)
    plt.text(9.0,1.36-0.30,r"$S={}\pm{},int.={}\pm{}, $".format(round(1/beta2,2),
                                                             round(beta2_std,2),
                                                             round(alpha2,2),
                                                             round(alpha2_std,2)),
             fontsize=12)
    #plt.text(9.0,1.36-0.38,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format(
    #    "{rms}",round(delta_rms2_y,2),"{rms}",round(delta_rms2,2)),
    #         fontsize=12)
    plt.text(9.0,1.36-0.38,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms2_y,2)),
             fontsize=12)
    
    plt.text(9.0,1.36-0.46,
             r"$Bisector~line$",fontsize=12)
    plt.text(9.0,1.36-0.54,r"$S={},int.={}, $".format(round(m_C,2),
                                                   round(k_C,2)),
             fontsize=12)
    #plt.text(9.0,1.36-0.62,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format("{rms}",
    #                                               round(delta_rms3,2),
    #                                               "{rms}",
    #                                               round(delta_rms3_x,2)),
    #         fontsize=12)
    plt.text(9.0,1.36-0.62,r"$\Delta_{}={}$".format("{rms}",
                                                   round(delta_rms3,2)),
             fontsize=12)
    
    plt.ylabel(r"$\log_{10}(n_\mathrm{Sph})$", fontsize=18)
    plt.xlabel(r"$\log_{10}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    plt.xlim(8.9,12.2)
    plt.ylim(-0.35,1.45)
    #plt.legend(loc='lower right', fontsize=10)
    plt.show()
    return fig

def find_bisector_line(m_A,m_B,k_A,k_B):
    x0 = (k_A-k_B) / (m_B-m_A)
    y0 = ((m_B*k_A)-(m_A*k_B))/(m_B-m_A)
    
    print('x0,y0',x0,y0)
    #Get x-intercept for A("
    q_A = (-1.0*k_A)/m_A
    
    #x_c= (x0-q_A) + (x0-0)
    #y_c = (y0-0) + (y0-k_B)
    
    x_c=  q_A
    y_c = k_B
    
    m_C = np.tan((np.arctan(m_B)+np.arctan(m_A))/2)
    k_C = y0-m_C*x0
    #calculate the bisector line slope and intercept
    
    return m_C,k_C

def plot_Bn_mass_z1(bn_50):
    fig = plt.figure(figsize=(6.4,5.8))
    
    z0=0.5
    
    n1,n2,n3,n4 = Sersic_n, Savorgnan_Sersic_n, Sahu_Sersic_n, Davis_Sersic_n
    bnz1, bnz2, bnz3, bnz4 = bn_50[0:101],bn_50[101:144], bn_50[170:203],bn_50[144:170]
    
    z1 = np.repeat(np.pi, len(bnz1))/np.exp(bnz1)
    z2 = np.repeat(np.pi, len(bnz2))/np.exp(bnz2)
    z3 = np.repeat(np.pi, len(bnz3))/np.exp(bnz3)
    z4 = np.repeat(np.pi, len(bnz4))/np.exp(bnz4)
    
    B1=np.log10(((2*z0*n1*np.repeat(np.pi, len(bnz1)))/((bnz1)**(2*n1)))*gamma(2*n1))
    B2=np.log10(((2*z0*n2*np.repeat(np.pi, len(bnz2)))/((bnz2)**(2*n2)))*gamma(2*n2))
    B3=np.log10(((2*z0*n3*np.repeat(np.pi, len(bnz3)))/((bnz3)**(2*n3)))*gamma(2*n3))
    B4=np.log10(((2*z0*n4* np.repeat(np.pi, len(bnz4)))/((bnz4)**(2*n4)))*gamma(2*n4))

    mass1, mass2, mass3, mass4 = np.log10(E_T11_K), np.log10(Savorgnan_mass_36), \
                         np.log10(Sahu_mass_36), np.log10(Davis_mass_36)
    
    plt.scatter(mass1, B1, s=100,
                           alpha = 0.8,marker='o',color='#a5200b',
                          label=r"$\rm Hon~et~al.~(2022)$")
    plt.scatter(mass2, B2, s=100,
                           alpha = 0.8, marker='o', color='#b940c8',
                           label=r"$\rm Savorgnan~&~Graham~(2016)$")
    plt.scatter(mass4, B4, s=100,
                          alpha = 0.8, marker='o', color='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$")
    plt.scatter(mass3, B3, s=100,
                           alpha = 0.8, marker='o', color='#e1a000',
                         label=r"$\rm Sahu~et~al.~(2019)$")
    
    
    #combine them all together
    E_all  = np.concatenate((mass1,mass2,mass3,mass4))
    B_all = np.concatenate((B1,B2,B3,B4))
    
    # print correlation coefficient
    print('rp_Bn',stats.pearsonr(E_all,B_all))
    print('rs_Bn',stats.spearmanr(E_all,B_all))
  
    plt.ylabel(r"$\log_{10}(B_\mathrm{e,Sph}(n_\mathrm{Sph}))$", fontsize=18)
    plt.xlabel(r"$\log_{10}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    #plt.xlim(8.5,12.5)
    #plt.ylim(0.1,0.9)
    #plt.legend(loc='lower right', fontsize=10)
    plt.show()
    
def plot_Bn_mass_z0(bn_50):
    fig = plt.figure(figsize=(6.4,5.8))
    
    z0=0.5
    
    n1,n2,n3,n4 = Sersic_n, Savorgnan_Sersic_n, Sahu_Sersic_n, Davis_Sersic_n
    bnz1, bnz2, bnz3, bnz4 = bn_50[0:101],bn_50[101:144], bn_50[170:203],bn_50[144:170]
    
    z1 = np.repeat(np.pi, len(bnz1))/np.exp(bnz1)
    z2 = np.repeat(np.pi, len(bnz2))/np.exp(bnz2)
    z3 = np.repeat(np.pi, len(bnz3))/np.exp(bnz3)
    z4 = np.repeat(np.pi, len(bnz4))/np.exp(bnz4)
    
    B1=np.log10(((2*z0*n1*np.exp(bnz1))/((bnz1)**(2*n1)))*gamma(2*n1))
    B2=np.log10(((2*z0*n2*np.exp(bnz2))/((bnz2)**(2*n2)))*gamma(2*n2))
    B3=np.log10(((2*z0*n3*np.exp(bnz3))/((bnz3)**(2*n3)))*gamma(2*n3))
    B4=np.log10(((2*z0*n4*np.exp(bnz4))/((bnz4)**(2*n4)))*gamma(2*n4))
   # A3_1_B = -2.5*np.log10(((2*z3*n1*np.exp(bnz3))/((bnz3)**(2*n1)))*gamma(2*n1))

    mass1, mass2, mass3, mass4 = np.log10(E_T11_K), np.log10(Savorgnan_mass_36), \
                         np.log10(Sahu_mass_36), np.log10(Davis_mass_36)
    
    plt.scatter(mass1, B1, s=100,
                           alpha = 0.8,marker='o',color='#a5200b',
                          label=r"$\rm Hon~et~al.~(2022)$")
    plt.scatter(mass2, B2, s=100,
                           alpha = 0.8, marker='o', color='#b940c8',
                           label=r"$\rm Savorgnan~&~Graham~(2016)$")
    plt.scatter(mass4, B4, s=100,
                          alpha = 0.8, marker='o', color='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$")
    plt.scatter(mass3, B3, s=100,
                           alpha = 0.8, marker='o', color='#e1a000',
                         label=r"$\rm Sahu~et~al.~(2019)$")
    
    
    #combine them all together
    E_all  = np.concatenate((mass1,mass2,mass3,mass4))
    B_all = np.concatenate((B1,B2,B3,B4))
    
    # print correlation coefficient
    print('rp_Bn',stats.pearsonr(E_all,B_all))
    print('rs_Bn',stats.spearmanr(E_all,B_all))

    print('rp_Bn',stats.pearsonr(B_all,E_all))
    print('rs_Bn',stats.spearmanr(B_all,E_all))

    # run fit
    lm = linmix.LinMix(E_all, B_all, K=2)
    lm.run_mcmc(silent=True)
    
    for i in range(0, len(lm.chain), 25):
        xs = mass0
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        plt.plot(xs, ys, color='r', alpha=0.02,lw=3)
        
    alpha = np.median(lm.chain[:]['alpha'])
    beta = np.median(lm.chain[:]['beta'])
    sigsqr = np.median(lm.chain[:]['sigsqr'])
    
    alpha_std = np.median(lm.chain[:]['alpha'].std())
    beta_std = np.median(lm.chain[:]['beta'].std())
    
    # vertical 
    ys_v = alpha + E_all * beta
    #horizontal
    xs_h = (B_all - alpha)/beta
    
    # calculate the root mean square
    delta_rms1 = np.sqrt(sum((ys_v - B_all)**2)/len(ys_v))
    delta_rms1_x = np.sqrt(sum((xs_h - E_all)**2)/len(xs_h))
    
    print("chain")
    print(alpha,beta)
    print("FIT_para",alpha,beta,sigsqr)
    
    ys = alpha + xs * beta
    plt.plot(xs, ys,linestyle='dashdot', color='k',lw=3)
    
    # run fit again
    lm2 = linmix.LinMix(B_all, E_all, K=2)
    lm2.run_mcmc(silent=True)
    
    for i in range(0, len(lm2.chain), 25):
        ys2 = mass0
        xs2 = (ys2 - lm2.chain[i]['alpha'])/ lm2.chain[i]['beta']
        plt.plot(ys2,xs2, color='cyan', alpha=0.02,lw=3)
        
    alpha2 = np.median(lm2.chain[:]['alpha'])
    beta2 = np.median(lm2.chain[:]['beta'])
    sigsqr2 = np.median(lm2.chain[:]['sigsqr'])
    
    ##flip it
    #alpha2 = -1.0*alpha2/beta2
    #beta2 = 1 / beta2
    alpha2_std = np.std(-1.0*lm2.chain[:]['alpha']/lm2.chain[:]['beta'])
    beta2_std = np.std(1/lm2.chain[:]['beta'])

    print("chain")
    print(alpha2,beta2)
    print("chain",lm2.chain[0],lm2.chain[-1])
    print("FIT_para",alpha2,beta2,sigsqr2)
    print("FIT_para_std",np.median(lm2.chain[:]['alpha'].std()),
          np.median(lm2.chain[:]['alpha'].std()))
        
    # vertical 
    xs_v2 = (E_all - alpha2)/ beta2 # This is y (B)
    # horizontal
    ys_h2 = alpha2 + B_all*beta2 # This is x (mass)
    
    # calculate the root mean square
    delta_rms2 = np.sqrt(sum((ys_h2- E_all)**2)/len(ys_h2))
    delta_rms2_y = np.sqrt(sum((xs_v2- B_all)**2)/len(xs_v2))
    
    xs2 = (ys2 - alpha2)/ beta2
    #ys2 =  alpha2+ beta2*xs2

    plt.plot(ys2,xs2, '--', color='k',lw=3) # Note that x and y are reversed
    
    
    m_C, k_C = find_bisector_line(1/beta2,beta,-1.0*alpha2/beta2,alpha)
    
    xs3 = mass0
    ys3 = m_C * xs3 + k_C
    
    # vertical 
    ys_v3 = k_C + E_all * m_C
    # horizontal
    xs_h3 = (B_all - k_C) / m_C
    
    # calculate the root mean square
    delta_rms3 = np.sqrt(sum((ys_v3- B_all)**2)/len(ys_v3))
    delta_rms3_x = np.sqrt(sum((xs_h3- E_all)**2)/len(xs_h3))
    
    plt.plot(xs3, ys3, '-', color='k',lw=3) # Note that x and y are reversed

    print("k_C, m_C", k_C, m_C)
    plt.hlines(0.8+0.04+0.01,9.65,10.0,linestyle='dashdot',color='k',lw=4)
    plt.hlines(0.8-0.08+0.01,9.65,10.0,linestyle='dashed',color='k',lw=4)
    plt.hlines(0.8-0.20+0.01,9.65,10.0,linestyle='solid',color='k',lw=4)
    
    plt.text(8.6,0.8+0.04,
             r"$B_{\rm e,Sph}(M_{\rm *,Sph}))$",fontsize=12)
    plt.text(8.6,0.8,r"$S={}\pm{},int.={}\pm{},$".format(round(beta,2),
                                                        round(beta_std,2),
                                                        round(alpha,2),
                                                        round(alpha_std,2)),
             fontsize=12)
    #plt.text(8.6,0.8-0.04,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format(
    #    "{rms}",round(delta_rms1,2),"{rms}",round(delta_rms1_x,2)),
    #         fontsize=12)
    plt.text(8.6,0.8-0.04,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms1,2)),
             fontsize=12)

    # X as independent variable
    plt.text(8.6,0.8-0.08,
             r"$M_{*,Sph}(B_{\rm e,Sph})$",fontsize=12)
    plt.text(8.6,0.8-0.12,r"$S={}\pm{},int.={}\pm{}, $".format(round(1/beta2,2),
                                                             round(beta2_std,2),
                                                             round(alpha2,2),
                                                             round(alpha2_std,2)),
             fontsize=12)
   # plt.text(8.6,0.8-0.16,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format(
   #     "{rms}",round(delta_rms2_y,2),"{rms}",round(delta_rms2,2)),
   #          fontsize=12)
    plt.text(8.6,0.8-0.16,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms2_y,2)),
             fontsize=12)
    
    plt.text(8.6,0.8-0.20,
             r"$Bisector~line$",fontsize=12)
    plt.text(8.6,0.8-0.24,r"$S={},int.={}, $".format(round(m_C,2),
                                                   round(k_C,2)),
             fontsize=12)
    #plt.text(8.6,0.8-0.28,r"$\Delta_{}={}(y),~\Delta_{}={}(x)$".format(
    #    "{rms}",round(delta_rms3,2),"{rms}",round(delta_rms3_x,2)),
    #         fontsize=12)
    plt.text(8.6,0.8-0.28,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms3,2)),
             fontsize=12)
    
    plt.ylabel(r"$\log_{10}(B_\mathrm{e,Sph}(n_\mathrm{Sph}))$", fontsize=18)
    plt.xlabel(r"$\log_{10}(M_\mathrm{*,Sph}/\rm M_{\odot})$", fontsize=18)
    plt.xlim(8.5,12.5)
    #plt.ylim(0.1,0.9)
    plt.legend(loc='lower right', fontsize=10)
    plt.show()
    return fig

def find_bisector_line(m_A,m_B,k_A,k_B):
    x0 = (k_A-k_B) / (m_B-m_A)
    y0 = ((m_B*k_A)-(m_A*k_B))/(m_B-m_A)
    
    print('x0,y0',x0,y0)
    #Get x-intercept for A("
    q_A = (-1.0*k_A)/m_A
    
    #x_c= (x0-q_A) + (x0-0)
    #y_c = (y0-0) + (y0-k_B)
    
    x_c=  q_A
    y_c = k_B
    
    m_C = np.tan((np.arctan(m_B)+np.arctan(m_A))/2)
    k_C = y0-m_C*x0
    #calculate the bisector line slope and intercept
    
    return m_C,k_C

def plot_n_mu0_z0():
    fig = plt.figure(figsize=(6.4,5.8))
    
    z0=0.5
    n1,n2,n3,n4 = Sersic_n, Savorgnan_Sersic_n, Sahu_Sersic_n, Davis_Sersic_n
    mu1,mu2,mu3,mu4 =  mu0,Savorgnan_mu_0,Sahu_mu_0, Davis_mu_0
    Sersic_n0 = np.linspace(0, 12,num=20)
    
    plt.scatter(n1, mu1, s=100,
                           alpha = 0.8,marker='o',color='#a5200b',
                          label=r"$\rm Hon~et~al.~(2022)$")
    plt.scatter(n2, mu2, s=100,
                           alpha = 0.8, marker='o', color='#b940c8',
                           label=r"$\rm Savorgnan~&~Graham~(2016)$")
    plt.scatter(n4, mu4, s=100,
                          alpha = 0.8, marker='o', color='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$")
    plt.scatter(n3, mu3, s=100,
                           alpha = 0.8, marker='o', color='#e1a000',
                         label=r"$\rm Sahu~et~al.~(2019)$")
    
    
    #combine them all together
    mu_all  = np.concatenate((mu1,mu2,mu3,mu4))
    n_all = np.concatenate((n1,n2,n3,n4))
    
    # print correlation coefficient
    print('rp_n',stats.pearsonr(n_all,mu_all))
    print('rs_n',stats.spearmanr(n_all,mu_all))

    # run fit
    lm = linmix.LinMix(n_all, mu_all, K=2)
    lm.run_mcmc(silent=True)
    
    for i in range(0, len(lm.chain), 25):
        xs = Sersic_n0
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        plt.plot(xs, ys, color='r', alpha=0.02,lw=3)
        
    alpha = np.median(lm.chain[:]['alpha'])
    beta = np.median(lm.chain[:]['beta'])
    sigsqr = np.median(lm.chain[:]['sigsqr'])
    
    alpha_std = np.median(lm.chain[:]['alpha'].std())
    beta_std = np.median(lm.chain[:]['beta'].std())
    
    # vertical 
    ys_v = alpha + n_all * beta
    # horizontal
    xs_h = (mu_all - alpha)/beta
    # calculate the root mean square
    delta_rms1 = np.sqrt(sum((ys_v - mu_all)**2)/len(ys_v))
    delta_rms1_x = np.sqrt(sum((xs_h - n_all)**2)/len(xs_h))
    
    print("chain")
    print(alpha,beta)
    print("FIT_para",alpha,beta,sigsqr)
    
    ys = alpha + xs * beta
    plt.plot(xs, ys,linestyle='dashdot', color='k',lw=3)
    
    # run fit again
    lm2 = linmix.LinMix(mu_all, n_all, K=2) # x and y flipped
    lm2.run_mcmc(silent=True)
    
    for i in range(0, len(lm2.chain), 25):
        ys2 = Sersic_n0
        xs2 = (ys2 - lm2.chain[i]['alpha'])/ lm2.chain[i]['beta']
        plt.plot(ys2,xs2, color='cyan', alpha=0.02,lw=3)
        
    alpha2 = np.median(lm2.chain[:]['alpha'])
    beta2 = np.median(lm2.chain[:]['beta'])
    sigsqr2 = np.median(lm2.chain[:]['sigsqr'])
    
    alpha2_std = np.std(-1.0*lm2.chain[:]['alpha']/lm2.chain[:]['beta'])
    beta2_std = np.std(1/lm2.chain[:]['beta'])

    print("chain")
    print(alpha2,beta2)
    print("chain",lm2.chain[0],lm2.chain[-1])
    print("FIT_para",alpha2,beta2,sigsqr2)
    print("FIT_para_std",np.median(lm2.chain[:]['alpha'].std()),
          np.median(lm2.chain[:]['alpha'].std()))
        
    # vertical 
    xs_v2 = (n_all - alpha2)/ beta2 # this is y (n)
    # horizontal 
    ys_h2 = alpha2 + beta2* mu_all # this is x (mass)
    # calculate the root mean square
    delta_rms2 = np.sqrt(sum((ys_h2- n_all)**2)/len(ys_h2))
    delta_rms2_y = np.sqrt(sum((xs_v2- mu_all)**2)/len(xs_v2))
    
    xs2 = (ys2 - alpha2)/ beta2
    plt.plot(ys2, xs2, '--', color='k',lw=3) # Note that x and y are reversed
    
    m_C, k_C = find_bisector_line(1/beta2,beta,-1.0*alpha2/beta2,alpha)
    
    xs3 = Sersic_n0
    ys3 = m_C * xs3 + k_C
    
    # vertical 
    ys_v3 = k_C + n_all * m_C
    # hortizontal
    xs_h3 = (n_all - k_C) / m_C
    
    # calculate the root mean square
    delta_rms3 = np.sqrt(sum((ys_v3- mu_all)**2)/len(ys_v3))
    delta_rms3_x = np.sqrt(sum((xs_h3- n_all)**2)/len(xs_h3))
    
    plt.plot(xs3, ys3, '-', color='k',lw=3) # Note that x and y are reversed

    print("k_C, m_C", k_C, m_C)
    print("delta_rms1_x",delta_rms1_x, 'delta_rms2_y', delta_rms2_y,
          'delta_rms3_x', delta_rms3_x)
    
    
    plt.hlines(7.0+0.08,2.61,3.66,linestyle='dashdot',color='k',lw=4)
    plt.hlines(7.0-2.1+0.08,2.61,3.66,linestyle='dashed',color='k',lw=4)
    plt.hlines(7.0-4.2+0.08,2.61,3.66,linestyle='solid',color='k',lw=4)
    
    plt.text(0.33,7.0,
             r"$\mu_{\rm 0, Sph}(n_{\rm Sph}))$",fontsize=12)
    plt.text(0.33,7.0-0.7,r"$S={}\pm{},int.={}\pm{},$".format(round(beta,2),
                                                        round(beta_std,2),
                                                        round(alpha,2),
                                                        round(alpha_std,2)),
             fontsize=12)

    plt.text(0.33,7.0-1.4,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms1,2)),
             fontsize=12)
    
    # X as independent variable
    plt.text(0.33,7.0-2.1,
             r"$n_{Sph}(\mu_{\rm Sph})$",fontsize=12)
    plt.text(0.33,7.0-2.8,r"$S={}\pm{},int.={}\pm{}, $".format(round(1/beta2,2),
                                                             round(beta2_std,2),
                                                             round(alpha2,2),
                                                             round(alpha2_std,2)),
             fontsize=12)
    plt.text(0.33,7.0-3.5,r"$\Delta_{}={}$".format(
        "{rms}",round(delta_rms2_y,2)),
             fontsize=12)
    
    plt.text(0.33,7.0-4.2,
             r"$Bisector~line$",fontsize=12)
    plt.text(0.33,7.0-4.9,r"$S={},int.={}, $".format(round(m_C,2),
                                                   round(k_C,2)),
             fontsize=12)

    plt.text(0.33,7.0-5.6,r"$\Delta_{}={}$".format("{rms}",
                                                   round(delta_rms3,2)),
             fontsize=12)
    
    plt.xlabel(r"$n_\mathrm{Sph}$", fontsize=18)
    plt.ylabel(r"$\mu_\mathrm{0,Sph}/\rm magarcsec^{-2}$", fontsize=18)
    plt.xlim(-0.05,11.0)
    plt.ylim(0.61,18.71)
    #plt.legend(loc='lower right', fontsize=10)
    plt.show()
    return fig


def plot_para_comparison():
    #Pioneer comparison
    fig = plt.figure(figsize=(6.4*3, 5.2))
    gs = gridspec.GridSpec(ncols=3, nrows=1,
                               hspace=0, wspace=0.0) 
    axt0 = plt.subplot(gs[0])   
    # plot the scatter point     
    axt0.scatter(Sersic_n_master,Abs_mag_master,
                 facecolors='k', edgecolors='none', s = 10,
                 label=r"$\rm This~work$")
    axt0.scatter(Baldry_Sersic_n,Baldry_sph_Mag,
                 facecolors='r', edgecolors='none', s = 20,
                 label=r"$\rm Baldry~et~al.~(2009)$")
    axt0.scatter(Khosroshahi_Sersic_n,Khosroshahi_sph_mag,
                 facecolors='g', edgecolors='none', s = 40,
                 label=r"$\rm Khosroshahi~et~al.~(2000)$") 
    axt0.scatter(Laurikainen_Sersic_n, Laurikainen_sph_mag,
                 facecolors='blue', edgecolors='none', s = 20,
                 label=r"$\rm Laurikainen~et~al.~(2010)$")
    
    axt0.set_ylabel(r"$\mathfrak{M}_\mathrm{Sph}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n_\mathrm{Sph} $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
    
    axt0.tick_params(axis='x',labelsize=15)
    axt0.tick_params(axis='y',labelsize=15)
    
    axt1 = plt.subplot(gs[1],sharey=axt0)   
    # plot the scatter point     
    axt1.scatter(mu0_master,Abs_mag_master,
                 facecolors='k', edgecolors='none', s = 10,
                 label=r"$\rm $")
    axt1.scatter(Baldry_mu0,Baldry_sph_Mag,
                 facecolors='r', edgecolors='none', s = 20,
                 label=r"$\rm $")
    axt1.scatter(Khosroshahi_mu0,Khosroshahi_sph_mag,
                 facecolors='g', edgecolors='none', s = 20,
                 label=r"$\rm $")
    axt1.scatter(Laurikainen_mu0, Laurikainen_sph_mag,
                 facecolors='blue', edgecolors='none', s = 20,
                 label=r"$\rm Laurikainen~et~al.~(2010)$")
    
    axt1.set_xlabel(r"$ \mu_\mathrm{0,Sph}~(\rm mag~arcsec^{-2})$", fontsize=22)
    axt1.set_xlim(19.2,-0.5)
    axt1.tick_params(axis='x',labelsize=15)
    
    axt2 = plt.subplot(gs[2],sharey=axt0)   
    # plot the scatter point     
    axt2.scatter(Re_kpc_master,Abs_mag_master,
                 facecolors='k', edgecolors='none', s = 30,
                 label=r"$\rm This~work$",zorder=40)
    axt2.scatter(Baldry_R_kpc,Baldry_sph_Mag,
                 facecolors='r', edgecolors='none', s = 20,
                 label=r"$\rm Baldry~et~al.~(2009)$")
    axt2.scatter(Khosroshahi_R_kpc,Khosroshahi_sph_mag,
                 facecolors='g', edgecolors='none', s = 20,
                 label=r"$\rm Khosroshahi~et~al.~(2000)$")
    axt2.scatter(Laurikainen_R_kpc, Laurikainen_sph_mag,
                 facecolors='blue', edgecolors='none', s = 20,
                 label=r"$\rm Laurikainen~et~al.~(2010)$")
    
    axt2.set_xlabel(r"$ R_\mathrm{e,Sph}~(\rm kpc)$", fontsize=22)
    axt2.set_xscale('log')
    axt2.set_xlim(0.09,120)
    axt2.legend(loc="lower right")

    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.setp(axt2.get_yticklabels(), visible=False)
    
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
    plt.tight_layout()

    axt2.tick_params(axis='x',labelsize=15)

    plt.show()
    return fig

def plot_mu0_n(n,mu,xlabel,ylabel):
    fig = plt.figure(figsize=(6.4, 5.2))
    
    plt.plot(n,mu,'o')
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    
#%%
def run_Linmix_bisector_fit(x,y,x0):
    # run fit
    lm = linmix.LinMix(x, y, K=2)
    lm.run_mcmc(silent=True)
    
    for i in range(0, len(lm.chain), 25):
        xs = x0
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        plt.plot(xs, ys, color='r', alpha=0.02,lw=3)
        
    alpha = np.median(lm.chain[:]['alpha'])
    beta = np.median(lm.chain[:]['beta'])
    sigsqr = np.median(lm.chain[:]['sigsqr'])
    
    alpha_std = np.median(lm.chain[:]['alpha'].std())
    beta_std = np.median(lm.chain[:]['beta'].std())
    
    # vertical 
    ys_v = alpha + x * beta
    # horizontal
    xs_h = (y - alpha)/beta
    # calculate the root mean square
    delta_rms1 = np.sqrt(sum((ys_v - y)**2)/len(ys_v))
    delta_rms1_x = np.sqrt(sum((xs_h - x)**2)/len(xs_h))
    
    print("chain")
    print(alpha,beta)
    print("FIT_para",alpha,beta,sigsqr)
    
    ys = alpha + xs * beta
    plt.plot(xs, ys,linestyle='dashdot', color='k',lw=3)
    
    # run fit again
    lm2 = linmix.LinMix(y, x, K=2) # x and y flipped
    lm2.run_mcmc(silent=True)
    
    for i in range(0, len(lm2.chain), 25):
        ys2 = x0
        xs2 = (ys2 - lm2.chain[i]['alpha'])/ lm2.chain[i]['beta']
        plt.plot(ys2,xs2, color='cyan', alpha=0.02,lw=3)
        
    alpha2 = np.median(lm2.chain[:]['alpha'])
    beta2 = np.median(lm2.chain[:]['beta'])
    sigsqr2 = np.median(lm2.chain[:]['sigsqr'])
    
    alpha2_std = np.std(-1.0*lm2.chain[:]['alpha']/lm2.chain[:]['beta'])
    beta2_std = np.std(1/lm2.chain[:]['beta'])

    print("chain")
    print(alpha2,beta2)
    print("chain",lm2.chain[0],lm2.chain[-1])
    print("FIT_para",alpha2,beta2,sigsqr2)
    print("FIT_para_std",np.median(lm2.chain[:]['alpha'].std()),
          np.median(lm2.chain[:]['alpha'].std()))
        
    # vertical 
    xs_v2 = (x - alpha2)/ beta2 # this is y (n)
    # horizontal 
    ys_h2 = alpha2 + beta2* y # this is x (mass)
    # calculate the root mean square
    delta_rms2 = np.sqrt(sum((ys_h2- x)**2)/len(ys_h2))
    delta_rms2_y = np.sqrt(sum((xs_v2- y)**2)/len(xs_v2))
    
    xs2 = (ys2 - alpha2)/ beta2
    plt.plot(ys2, xs2, '--', color='k',lw=3) # Note that x and y are reversed
    
    m_C, k_C = find_bisector_line(1/beta2,beta,-1.0*alpha2/beta2,alpha)
    
    xs3 = x0
    ys3 = m_C * xs3 + k_C
    
    # vertical 
    ys_v3 = k_C + x * m_C
    # hortizontal
    xs_h3 = (x - k_C) / m_C
    
    # calculate the root mean square
    delta_rms3 = np.sqrt(sum((ys_v3- y)**2)/len(ys_v3))
    delta_rms3_x = np.sqrt(sum((xs_h3- x)**2)/len(xs_h3))
    
    plt.plot(xs3, ys3, '-', color='k',lw=3) # Note that x and y are reversed

    print("k_C, m_C", k_C, m_C)
    print("delta_rms1_x",delta_rms1_x, 'delta_rms2_y', delta_rms2_y,
          'delta_rms3_x', delta_rms3_x)
    
    return (xs,ys,xs2,ys2,xs3,ys3), (alpha,beta,alpha_std,beta_std,delta_rms1),\
        (alpha2,beta2,alpha2_std,beta2_std,delta_rms2), (m_C, k_C,delta_rms3)


#%% Execution Area

def ScS_type_generate(file_name):
    table = SRead.read_list(file_name)
    ScStype= []
    for i in range(len(table)):
        if "Bulge" in table[i]:
            ScStype.append("S")
        elif "CoreBulge" in table[i]:
            ScStype.append("cS")
    return ScStype

ScStype = ScS_type_generate("/home/dexter/SphProject/F_Gal_bundle_equvi_cpt")

# extrapolate the central surface brightness
R_gen = np.linspace(0,300,300*2)

# produce the stacked radial profile figure, as well as the mu0 
#stack = plot_stack_surface_brightness_profile(R_gen)
#list_mu0_extrapolate(stack[1])
#new_mu_list = stack[1]

# create morphlogy dependent parameters arrays
# For multi-cpt decomposition in H21
n_combine_morph, mu0_combine_morph, Mag_combine_morph, E_combine_morph, \
    Re_kpc_combine_morph, Re_kpc_err_combine, mass_err_combine, \
        Abs_total_mag_combine_morph = seperation_morph(
    index, morph, sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag)
        
# For multi-cpt decomposition in SG16
ZSG16 = np.zeros(len(Savorgnan_mass_36))

n_combine_morph_SG16, mu0_combine_morph_SG16, Mag_combine_morph_SG16, \
    E_combine_morph_SG16, Re_kpc_combine_morph_SG16, Re_kpc_err_combine_SG16, \
        mass_err_combine_SG16, Abs_total_mag_combine_morph_SG16 = seperation_morph(
    index_SG16, Savorgnan_morph, ZSG16, Savorgnan_size_eq_kpc, ZSG16, ZSG16,ZSG16, 
                     ZSG16, ZSG16, ZSG16, ZSG16, np.full(len(ZSG16),1/ars),
                     ZSG16, ZSG16, ZSG16, ZSG16, Savorgnan_mass_36,
                     ZSG16, ZSG16,ZSG16)

#Savorgnan_size_eq_kpc = Savorgnan_data[:,8]
#Savorgnan_mass_36 = 10**Savorgnan_data[:,10]

# For multi-cpt decomposition in D19
ZD19 = np.zeros(len(Davis_mass_36))

n_combine_morph_D19, mu0_combine_morph_D19, Mag_combine_morph_D19, \
    E_combine_morph_D19, Re_kpc_combine_morph_D19, Re_kpc_err_combine_D19, \
        mass_err_combine_D19, Abs_total_mag_combine_morph_D19 =seperation_morph(
    index_D19, Davis_morph, ZD19, Davis_size_eq_kpc, ZD19, ZD19,ZD19, 
                     ZD19, ZD19, ZD19, ZD19, np.full(len(ZD19),1/ars),
                     ZD19, ZD19, ZD19, ZD19, Davis_mass_36,
                     ZD19, ZD19,ZD19)
# For multi-cpt decomposition in S19
ZS19 = np.zeros(len(Sahu_mass_36))

n_combine_morph_S19, mu0_combine_morph_S19, Mag_combine_morph_S19, \
    E_combine_morph_S19, Re_kpc_combine_morph_S19, Re_kpc_err_combine_S19,\
        mass_err_combine_S19, Abs_total_mag_combine_morph_S19 = seperation_morph(
    index_S19, Sahu_morph, ZS19, Sahu_size_eq_kpc, ZS19, ZS19,ZS19, 
                     ZS19, ZS19, ZS19, ZS19, np.full(len(ZS19),1/ars),
                     ZS19, ZS19, ZS19, ZS19, Sahu_mass_36,
                     ZS19, ZS19,ZS19)

# for multi-cpt decompiosition in all 3
ZALL3 = np.zeros(len(ALL3_mass))
print("len(ALL3_mass))",len(ALL3_mass))

n_combine_morph_A3, mu0_combine_morph_A3, Mag_combine_morph_A3, \
    E_combine_morph_A3, Re_kpc_combine_morph_A3, Re_kpc_err_combine_A3, \
        mass_err_combine_A3, Abs_total_mag_combine_morph_A3 = seperation_morph(
    index_A3, ALL3_morph, ALL3_Mag, ALL3_size_eq_kpc, ALL3_Sersic_n, ALL3_mu_e,ALL3_gal_Mag, 
                     ALL3_name, ALL3_mu_0, ZALL3, ZALL3, np.full(len(ZALL3),1/ars),
                     ZALL3, ZALL3, ALL3_Mag, ZALL3, ALL3_mass,
                     ZALL3, ZALL3, ALL3_gal_Mag)


# For B+D decomposition in H21
n_combine_morph_BD, mu0_combine_morph_BD, Mag_combine_morph_BD, E_combine_morph_BD, Re_kpc_combine_morph_BD, Re_kpc_err_combine_BD, mass_err_combine_BD, Abs_total_mag_combine_morph_BD = seperation_morph(
    index, morph, sph_mag_BD, Re_BD, Sersic_n_BD, mu_e_BD,total_mag_BD, 
                     name, mu0_BD, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag_BD, elle, E_BD,
                     MLR, mass_uerr,Abs_total_mag_BD)

# For 1 Sersic decomposition in H21
n_combine_morph_1Sersic, mu0_combine_morph_1Sersic, Mag_combine_morph_1Sersic, E_combine_morph_1Sersic, Re_kpc_combine_morph_1Sersic, Re_kpc_err_combine_1Sersic, mass_err_combine_1Sersic, Abs_total_mag_combine_morph_1Sersic = seperation_morph(
    index, morph, sph_mag_1Sersic, Re_1Sersic, Sersic_n_1Sersic, mu_e_1Sersic,total_mag_1Sersic, 
                     name, mu0_1Sersic, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag_1Sersic, elle, E_1Sersic,
                     MLR, mass_uerr,Abs_total_mag_1Sersic)

n_combine_morph, mu0_combine_morph, Mag_combine_morph, E_combine_morph, Re_kpc_combine_morph, Re_kpc_err_combine, mass_err_combine, Abs_total_mag_combine_morph =seperation_morph(
    index, morph, sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag)


# plot 3 plots
#CC = plot_n_mu0_Mag_3plot(n_combine_morph,mu0_combine_morph, Re_kpc_combine_morph, Mag_combine_morph,
#                         Sersic_n_corebulge,
#                         mu0_corebulge, Re_kpc_corebulge, 
#                         Abs_sph_mag_corebulge,
#                         label=[r"$\rm E+ES$",r"$\rm S0$", r"$\rm S$"])

#plot_n_mu0_Mag_3Lplot(n_combine_morph,mu0_combine_morph, Re_kpc_combine_morph, 
#                         label=[r"$\rm E+ES$",r"$\rm S0$", r"$\rm S$"],fit_instruc=0)

print('Corrlation',
      stats.pearsonr(np.log10(Sersic_n),Abs_sph_mag),
      stats.pearsonr(mu0,Abs_sph_mag),
      stats.pearsonr(np.log10(Re_kpc),Abs_sph_mag))

print('Corrlation_E+ES',
      stats.pearsonr(np.log10(n_combine_morph[0]),Abs_total_mag_combine_morph[0]),
      stats.pearsonr(mu0_combine_morph[0],Abs_total_mag_combine_morph[0]),
      stats.pearsonr(np.log10(Re_kpc_combine_morph[0]),Abs_total_mag_combine_morph[0]))

print('Corrlation_S0+S',
      stats.pearsonr(np.log10(np.concatenate((n_combine_morph[1],n_combine_morph[2]))),
                     np.concatenate((Abs_total_mag_combine_morph[1],Abs_total_mag_combine_morph[2]))),
      stats.pearsonr(np.concatenate((mu0_combine_morph[1],mu0_combine_morph[2])),
                     np.concatenate((Abs_total_mag_combine_morph[1],Abs_total_mag_combine_morph[2]))),
      stats.pearsonr(np.log10(np.concatenate((Re_kpc_combine_morph[1],Re_kpc_combine_morph[2]))),
                     np.concatenate((Abs_total_mag_combine_morph[1],Abs_total_mag_combine_morph[2]))))



#CCC =  plot_n_mu0_Mag_3plot(n_combine_morph_A3,mu0_combine_morph_A3, 
#                            Re_kpc_combine_morph_A3, Mag_combine_morph_A3,
#                         0,
#                         0, 0, 
#                         0,                     
#                     label=[r"$\rm E+ES$",r"$\rm S0$", r"$\rm S$"])


print('Corrlation_ALL3',
      stats.pearsonr(np.log10(ALL3_Sersic_n),ALL3_Mag),
      stats.pearsonr(ALL3_mu_0,ALL3_Mag),
      stats.pearsonr(np.log10(ALL3_size_eq_kpc),ALL3_Mag))

#print('morph, B+D')
#DD = plot_n_mu0_Mag_2plot(n_combine_morph_BD,mu0_combine_morph_BD,
#                     Mag_combine_morph_BD,label=[r"$\rm E+ES$",
#                                        r"$\rm S0$", r"$\rm S$"])
#print('morph, 1-Sersic')
#E = plot_n_mu0_Mag_2plot(n_combine_morph_1Sersic,mu0_combine_morph_1Sersic,
#                     Abs_total_mag_combine_morph_1Sersic,label=[r"$\rm E+ES$",
#                                        r"$\rm S0$", r"$\rm S$"])

#%%
## plot size-mass relation
# plot_sizemass_morph()
# plot_sizemass_z0()


## Plot comparison with other studies
#plot_size_mass_comparison()

#%%
#calculate bn Re and mu 10 and 90 percent
#Re_kpc_10_combine_morph, Re_kpc_90_combine_morph, mu0_10_combine_morph, mu0_50_combine_morph, mu0_90_combine_morph = get_bn_Re_mu_1090()

#Read the pickle file
bn_dict = SRead.read_list("bn_dict.pkl")
Re_kpc_10_combine_morph =  bn_dict["Re_kpc_10_combine_morph"]
Re_kpc_90_combine_morph = bn_dict["Re_kpc_90_combine_morph"]
mu0_10_combine_morph = bn_dict["mu0_10_combine_morph"]
mu0_50_combine_morph = bn_dict["mu0_50_combine_morph"]
mu0_90_combine_morph = bn_dict["mu0_90_combine_morph"]

# plot the characteristic radius
#plot_Mag_Re_morph_diverse(Re_kpc_10_combine_morph,Re_kpc_combine_morph,Re_kpc_90_combine_morph,Mag_combine_morph)  
#plot_Mag_mu_morph_diverse(mu0_combine_morph,mu0_50_combine_morph,mu0_90_combine_morph,Mag_combine_morph)

#plot_Mag_diverse_2plot(Re_kpc_10_combine_morph,Re_kpc_combine_morph,
 #                         Re_kpc_90_combine_morph, mu0_combine_morph,
 #                         mu0_50_combine_morph,mu0_90_combine_morph, 
 #                         Mag_combine_morph)


print('BH',
      stats.pearsonr(np.log10(ALL3_mass),np.log10(ALL3_BH_mass)),
      stats.pearsonr(np.log10(ALL3_gal_mass),np.log10(ALL3_BH_mass)))

#%%
#plot the 3D sitribution of the structral parameters
#plot_3D_para()
#%% 
# calculate the sigma value
Sigma0 = calculate_sigma_0(D,mu0,scale_other,MLR,np.repeat(4.53,len(MLR)))
Sigma0_SG16 = calculate_sigma_0(Savorgnan_dist,Savorgnan_mu_0,Savorgnan_scale,Savorgnan_MLR,np.repeat(4.53,len(Savorgnan_MLR)))
Sigma0_D19 = calculate_sigma_0(Davis_dist,Davis_mu_0,Davis_scale,Davis_MLR,np.repeat(4.53,len(Davis_MLR)))
Sigma0_S19 = calculate_sigma_0(Sahu_dist,Sahu_mu_0,Sahu_scale,Sahu_MLR,np.repeat(4.53,len(Sahu_MLR)))


Sigma0_A3 = calculate_sigma_0(ALL3_dist,ALL3_mu_0,ALL3_scale,ALL3_MLR,ALL3_Mo)

# Plot the sigma0 relations
#plot_mass_n_sigma0_2plot(E_T11_K,Sersic_n,Sigma0,
#                         Savorgnan_mass_36,Savorgnan_Sersic_n,Sigma0_SG16,
#                         Davis_mass_36,Davis_Sersic_n,Sigma0_D19,
#                         Sahu_mass_36,Sahu_Sersic_n,Sigma0_S19)
                      
bn_50 =  bn_dict['bn_50']

bn_50 = np.delete(bn_50, [43,73])

mu_50 = mu0 + (2.5*bn_50)/np.log(10)

#%% average sigma
# plot term a and term b + average Sigma0
#plot_average_mu_2plot(mu_50, Sersic_n, E_T11_K, bn_50)   

#make all data into one array
Sersic_n_master = np.concatenate((Sersic_n,Savorgnan_Sersic_n,Davis_Sersic_n,Sahu_Sersic_n))
mass_master = np.concatenate((E_T11_K,Savorgnan_mass_36,Davis_mass_36,Sahu_mass_36))
Sigma0_master = np.concatenate((Sigma0,Sigma0_SG16,Sigma0_D19,Sigma0_S19))

Dist_master = np.concatenate((D,Savorgnan_dist,Davis_dist,Sahu_dist))
mu0_master = np.concatenate((mu0,Savorgnan_mu_0,Davis_mu_0,Sahu_mu_0))
scale_master = np.concatenate((scale_other,Savorgnan_scale,Davis_scale,Sahu_scale))
MLR_master = np.concatenate((MLR,Savorgnan_MLR,Davis_MLR,Sahu_MLR))

Re_kpc_master = np.concatenate((Re_kpc, Savorgnan_size_eq_kpc, 
                               Davis_size_eq_kpc, Sahu_size_eq_kpc))
morph_master = np.concatenate((morph,Savorgnan_morph,Davis_morph,Sahu_morph))

Abs_mag_master = np.concatenate((Abs_sph_mag,Savorgnan_Mag,Davis_Mag,Sahu_Mag))
# Get the bn and store it. One time use only
index_all = np.arange(len(morph_master))
#bn_dic_BIG = get_bn_Re_mu_1090_BIG(Sersic_n_master,Re_kpc_master,mu0_master,index_all,morph_master)

## Get the Sigma  with different fraction z
bn_dict_BIG = SRead.read_list("bn_dict_BIG.pkl")
#
#bn_05 = np.delete(bn_dict_BIG['bn_05'], [43,73])
#bn_10 = np.delete(bn_dict_BIG['bn_10'], [43,73])
#bn_20 = np.delete(bn_dict_BIG['bn_20'], [43,73])
#bn_30 = np.delete(bn_dict_BIG['bn_30'], [43,73])
#bn_40 = np.delete(bn_dict_BIG['bn_40'], [43,73])
#bn_50 = np.delete(bn_dict_BIG['bn_50'], [43,73])
#bn_60 = np.delete(bn_dict_BIG['bn_60'], [43,73])
#bn_70 = np.delete(bn_dict_BIG['bn_70'], [43,73])
#bn_80 = np.delete(bn_dict_BIG['bn_80'], [43,73])
#bn_90 = np.delete(bn_dict_BIG['bn_90'], [43,73])
#bn_95 = np.delete(bn_dict_BIG['bn_95'], [43,73])
#bn_100 = np.delete(bn_dict_BIG['bn_100'], [43,73])
#
#Re_kpc_05 = np.delete(bn_dict_BIG['Re_kpc_05'],[43,73])
#Re_kpc_10 = np.delete(bn_dict_BIG['Re_kpc_10'],[43,73])
#Re_kpc_20 = np.delete(bn_dict_BIG['Re_kpc_20'],[43,73])
#Re_kpc_30 = np.delete(bn_dict_BIG['Re_kpc_30'],[43,73])
#Re_kpc_40 = np.delete(bn_dict_BIG['Re_kpc_40'],[43,73])
#Re_kpc_50 = np.delete(bn_dict_BIG['Re_kpc_50'],[43,73])
#Re_kpc_60 = np.delete(bn_dict_BIG['Re_kpc_60'],[43,73])
#Re_kpc_70 = np.delete(bn_dict_BIG['Re_kpc_70'],[43,73])
#Re_kpc_80 = np.delete(bn_dict_BIG['Re_kpc_80'],[43,73])
#Re_kpc_90 = np.delete(bn_dict_BIG['Re_kpc_90'],[43,73])
##Re_kpc_95 = np.delete(bn_dict_BIG['Re_kpc_95'],[43,73])
##Re_kpc_100 = np.delete(bn_dict_BIG['Re_kpc_100'],[43,73])
#
#mu_05 = np.delete(bn_dict_BIG['mu_05'],[43,73])
#mu_10 = np.delete(bn_dict_BIG['mu_10'],[43,73])
#mu_20 = np.delete(bn_dict_BIG['mu_20'],[43,73])
#mu_30 = np.delete(bn_dict_BIG['mu_30'],[43,73])
#mu_40 = np.delete(bn_dict_BIG['mu_40'],[43,73])
#mu_50 = np.delete(bn_dict_BIG['mu_50'],[43,73])
#mu_60 = np.delete(bn_dict_BIG['mu_60'],[43,73])
#mu_70 = np.delete(bn_dict_BIG['mu_70'],[43,73])
#mu_80 = np.delete(bn_dict_BIG['mu_80'],[43,73])
#mu_90 = np.delete(bn_dict_BIG['mu_90'],[43,73])
#mu_95 = np.delete(bn_dict_BIG['mu_95'],[43,73])
#mu_100 = np.delete(bn_dict_BIG['mu_100 (0, (1, 1))'],[43,73])

bn_05 = bn_dict_BIG['bn_05']
bn_10 = bn_dict_BIG['bn_10']
bn_20 = bn_dict_BIG['bn_20']
bn_30 = bn_dict_BIG['bn_30']
bn_40 = bn_dict_BIG['bn_40']
bn_50 = bn_dict_BIG['bn_50']
bn_60 = bn_dict_BIG['bn_60']
bn_70 = bn_dict_BIG['bn_70']
bn_80 = bn_dict_BIG['bn_80']
bn_90 = bn_dict_BIG['bn_90']
bn_95 = bn_dict_BIG['bn_95']
bn_100 = bn_dict_BIG['bn_100']

Re_kpc_05 = bn_dict_BIG['Re_kpc_05']
Re_kpc_10 = bn_dict_BIG['Re_kpc_10']
Re_kpc_20 = bn_dict_BIG['Re_kpc_20']
Re_kpc_30 = bn_dict_BIG['Re_kpc_30']
Re_kpc_40 = bn_dict_BIG['Re_kpc_40']
Re_kpc_50 = bn_dict_BIG['Re_kpc_50']
Re_kpc_60 = bn_dict_BIG['Re_kpc_60']
Re_kpc_70 = bn_dict_BIG['Re_kpc_70']
Re_kpc_80 = bn_dict_BIG['Re_kpc_80']
Re_kpc_90 = bn_dict_BIG['Re_kpc_90']
#Re_kpc_95 = np.delete(bn_dict_BIG['Re_kpc_95'],[43,73])
#Re_kpc_100 = np.delete(bn_dict_BIG['Re_kpc_100'],[43,73])

mu_05 = bn_dict_BIG['mu_05']
mu_10 = bn_dict_BIG['mu_10']
mu_20 = bn_dict_BIG['mu_20']
mu_30 = bn_dict_BIG['mu_30']
mu_40 = bn_dict_BIG['mu_40']
mu_50 = bn_dict_BIG['mu_50']
mu_60 = bn_dict_BIG['mu_60']
mu_70 = bn_dict_BIG['mu_70']
mu_80 = bn_dict_BIG['mu_80']
mu_90 = bn_dict_BIG['mu_90']
mu_95 = bn_dict_BIG['mu_95']
mu_100 = bn_dict_BIG['mu_100']

Sigma05 = calculate_sigma_0(Dist_master,mu_05,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma10 = calculate_sigma_0(Dist_master,mu_10,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma20 = calculate_sigma_0(Dist_master,mu_20,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma30 = calculate_sigma_0(Dist_master,mu_30,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma40 = calculate_sigma_0(Dist_master,mu_40,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma50 = calculate_sigma_0(Dist_master,mu_50,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma60 = calculate_sigma_0(Dist_master,mu_60,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma70 = calculate_sigma_0(Dist_master,mu_70,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma80 = calculate_sigma_0(Dist_master,mu_80,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma90 = calculate_sigma_0(Dist_master,mu_90,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma95 = calculate_sigma_0(Dist_master,mu_95,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
Sigma100 = calculate_sigma_0(Dist_master,mu_100,scale_master, 
                            MLR_master,np.repeat(4.53,len(MLR_master)))
# List them all up
Sigma_z_all = [Sigma0_master,Sigma05,Sigma10,Sigma20,Sigma30,Sigma40,Sigma50,
             Sigma60,Sigma70,Sigma80,Sigma90,Sigma95]
R_z_all = [Re_kpc_05,Re_kpc_10,Re_kpc_20,
           Re_kpc_30,Re_kpc_40,Re_kpc_50,
             Re_kpc_60,Re_kpc_70,Re_kpc_80,
             Re_kpc_90]
mu_z_all = [mu_05,mu_10,mu_20,
            mu_30,mu_40,mu_50,
             mu_60,mu_70,mu_80,
             mu_90,mu_95]

bnz_all = [bn_05,bn_10,bn_20,
            bn_30,bn_40,bn_50,
             bn_60,bn_70,bn_80,
             bn_90,bn_95]

# plot 3 plots
#CC = plot_n_mu0_Mag_3plot(n_combine_morph,mu0_combine_morph, Re_kpc_combine_morph, Mag_combine_morph,
#                         Sersic_n_corebulge,
#                         mu0_corebulge, Re_kpc_corebulge, 
#                         Abs_sph_mag_corebulge,
#                         label=[r"$\rm E+ES$",r"$\rm S0$", r"$\rm S$"])

#plot_n_mu0_Mag_3Lplot(n_combine_morph,mu0_combine_morph, Re_kpc_combine_morph, 
#                         label=[r"$\rm E+ES$",r"$\rm S0$", r"$\rm S$"],fit_instruc=0)

# plot the correlation by fraction z diagram

#plot_Sigma_z()

#plot_term_b(Sersic_n_master)
#plot_corr_trend(mass_master,R_z_all,Sigma_z_all)

#plot_n_mass_z0()
#plot_Bn_mass_z0(bn_50)
#plot_Bn_mass_z1(bn_50)
#plot_Bn_n(Sersic_n_master)

# Pioneer comparison
#plot_para_comparison()

#plot_mu0_n(Sersic_n_master,mu0_master,'n','mu0')
#plot_mu0_n(np.log10(Sersic_n_master),np.log10(Re_kpc_master),'n','Re')

print('rp_Bn',stats.pearsonr(Sersic_n_master,mu0_master))    
print('rs_Bn',stats.spearmanr(Sersic_n_master,mu0_master))

print('rp_Bn',stats.pearsonr(np.log10(Sersic_n_master),np.log10(Re_kpc_master)))
print('rs_Bn',stats.spearmanr(np.log10(Sersic_n_master),np.log10(Re_kpc_master)))

# other scaling relations
plot_n_mu0_z0()
