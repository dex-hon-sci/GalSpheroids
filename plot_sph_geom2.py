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

from scipy.optimize import curve_fit

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0
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
print('len(sph_mag_1Sersic)',len(sph_mag_1Sersic))
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

#print(Abs_sph_mag_dustCorr,Abs_disc_mag_dustCorr)

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
MLR_dust = ML_select_T11_dustcorr
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
ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale
scale = D* ars

scale_lerr, scale_uerr = (D-D_lerr)*ars, (D+D_uerr)*ars

# Calculate Re with multi-cpt, B+D, and 1-Sersic
Re_kpc = Re* scale
Re_kpc_BD = Re_BD* scale
Re_kpc_1Sersic = Re_1Sersic* scale

Re_kpc_lerr, Re_kpc_uerr = abs(Re_kpc - Re* scale_lerr) , abs(Re* scale_uerr - Re_kpc)
Re_kpc_err =[Re_kpc_lerr, Re_kpc_uerr]

print(Sersic_n[26])

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

    return Re_kpc_10_combine_morph,Re_kpc_90_combine_morph, mu0_10_combine_morph, mu0_50_combine_morph, mu0_90_combine_morph
#K correction and extinction for sph mag
sph_mag = sph_mag - i_EXT - i_kcorr
#%%
# Read Sahu, Davis, and Savorgnan

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
#Savorgnan_mass_T11 = 10**Savorgnan_mass_36

Davis_name = Davis_data_n[:,0]
Davis_size_eq_kpc = Davis_data[:,8]
Davis_mass_36 = Davis_data[:,10]

Davis_mass_T11 = 10**(0.88 * Davis_mass_36+1.02)
#Davis_mass_T11 = 10**Davis_mass_36

Sahu_name = Sahu_data_n[:,0]
Sahu_size_eq_kpc = Sahu_data[:,8]
Sahu_mass_36 = Sahu_data[:,10]

Sahu_mass_T11 = 10**(0.88 * Sahu_mass_36+1.02)
#Sahu_mass_T11 = 10**Sahu_mass_36


## End in reading base data#################################################

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

##############End Seperation 1#################################################

####Start Seperation 2 ########################################################
index = np.arange(103)

def seperate_morph_simple(A,index, morph):
    morph_dict = SSort.morph_str_selection(index, morph)
    
    A_E = SSort.cherry_pick(morph_dict['E'], A)
    A_S0 = SSort.cherry_pick(morph_dict['S0'], A)
    A_S = SSort.cherry_pick(morph_dict['S'], A)
    
    A_combine_morph = np.array([A_E,A_S0,A_S])
    return A_combine_morph


def seperation_morph(sph_mag, Re, Sersic_n, mu_e,total_mag, 
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
    Abs_sph_mag_dustCorr_S = SSort.cherry_pick(morph_dict['S'],Abs_sph_mag_dustCorr)
    
    E_T11_K_E = SSort.cherry_pick(morph_dict['E'], E_T11_K)
    E_T11_K_S0 = SSort.cherry_pick(morph_dict['S0'], E_T11_K)
    E_T11_K_S = SSort.cherry_pick(morph_dict['S'], E_T11_K)

    E_T11_K_dustcorr_E = SSort.cherry_pick(morph_dict['E'], E_T11_K_dustCorr)
    E_T11_K_dustcorr_S0 = SSort.cherry_pick(morph_dict['S0'], E_T11_K_dustCorr)
    E_T11_K_dustcorr_S = SSort.cherry_pick(morph_dict['S'], E_T11_K_dustCorr)

    E_T11_K_dustcorr_E_old = SSort.cherry_pick(morph_dict['E'], E_T11_K_dustCorr_old)   
    E_T11_K_dustcorr_S0_old = SSort.cherry_pick(morph_dict['S0'], E_T11_K_dustCorr_old)
    E_T11_K_dustcorr_S_old = SSort.cherry_pick(morph_dict['S'], E_T11_K_dustCorr_old)

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
############End reading 

#%%
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
                
                b = SAna.get_bn(n) # get_bn is the exact function used in profiler
                # obtain mu_0 by equ7 from Graham2005
                new_mu = mu_e-2.5*b/np.log(10)
                print(name[i],new_mu,Re[i],n,sph_mag[i]) # list mu0
                new_mu_list.append(new_mu) # save the new mu_0
                
                #list name, first element of from the radial SB, and the mu_0 
                #from equ 7
                #print(name[i],line[0])#, new_mu
                pass
            elif identifier == "CoreBulge": #if the spheroid is a Core Sersic bulge
                para = bundle[i][j+1]
                
                # calculate the mu_e of the core Sersic funtion
                Re, n = para[1], para[2]
                mu_e = SAna.AnalyticFunctions.simple_mu_core_sersic(Re,*list(para))

                # put in the core Sersic fit parameter mu_e, r_e , n_ser
                line = SAna.AnalyticFunctions.mu_sersic_func(r,mu_e,para[1],para[2])
                plt.plot(r,line,"k-",lw= 3)
                new_mu_list.append(line[0])

                print(name[i],line[0], Re[i],n,sph_mag[i]) #list mu0
                pass
    plt.gca().invert_yaxis()
    #plt.ylim(np.log10(24),np.log(10))
    plt.ylim(24,12)
    plt.xlim(0,120)
    
    plt.tight_layout()
    plt.show()
    
    return fig,new_mu_list


#%%
def list_mu0_extrapolate(new_mu_list):
    """
    List the name and mu0
    mu0 is an extrapolation to the centre

    """
    for i in range(len(name)):
        print(name[i],new_mu_list[i])
        
#%%
def plot_n_mu0_Mag_2plot(n,mu,Mag,label=[],fit_instruc=0):
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
    n_line = np.linspace(-5, 20, 100)
    mu_line = np.linspace(-5, 30, 100)

    #Fitting the relevant data
    popt_mu,pcov_mu = [], []
    popt_n,pcov_n = [], []

    for i in range(len(mu)):
        popt_mu_ind,pcov_mu_ind = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                mu[i], Mag[i],p0=[6.12, -117])
        popt_n_ind,pcov_n_ind = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                              np.log10(n[i]), Mag[i],p0=[-9, -17.5])  
        
        popt_mu.append(popt_mu_ind)
        pcov_mu.append(pcov_mu_ind)
        popt_n.append(popt_n_ind)
        pcov_n.append(pcov_n_ind)
        
    # Fitting the whole data set
    # restich them together
    n_inone = n[1]#np.concatenate(n)
    mu_inone= mu[1]#np.concatenate(mu)
    Mag_inone= Mag[1]#np.concatenate(Mag)
        
    print('size', len(n_inone),len(mu_inone),len(Mag_inone))
    
    # Fit all data in one
    popt_mu_inone,pcov_mu_inone = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                mu_inone, Mag_inone,p0=[2.62, -60],#p0=[.63,-28.4])#,
                                bounds = ((2.6,-62),(2.7,-59)), method='trf',
                                max_nfev = 100000)                           
    popt_n_inone,pcov_n_inone = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                              np.log10(n_inone), Mag_inone,p0=[-9, -17.5],
                              bounds = ((-9.6,-17.9),(-8.8,-17.0)), method='trf',
                              max_nfev = 100000)  
    

#-8.75, -17.85,6.12, -93
#-9, -17.85,6.12, -113
    print('===========================================')
    print("popt_mu_inone",*popt_mu_inone)
    print("popt_n_inone", *popt_n_inone)
    print('-------------------------------------------')
    #print("1",mu0[fit_instruc],Mag[fit_instruc])
    #print("2",np.log10(n[fit_instruc]),Mag[fit_instruc])

    #Plotting
    fig = plt.figure(figsize=(12, 4.8))
    gs = gridspec.GridSpec(ncols=2, nrows=1,
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
            axt0.plot(n_line, 
              SAna.AnalyticFunctions.linear_func_1D(np.log10(n_line),*popt_n[i]),
              color=colour[i],linestyle='dashed',lw=3,label='in_one')
    
    axt0.plot(n_line, 
              SAna.AnalyticFunctions.linear_func_1D(np.log10(n_line),*popt_n_inone),
              color='b',linestyle='solid',lw=3,label='in_one')
    
    
    #print(S_bar[0]['PrimBar_index'])
    
    #plot bar bulge
    n_barbulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'],Sersic_n)
    Mag_barbulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'], Abs_sph_mag)   
    mu_barbulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'], mu0)


   # axt0.plot(n_barbulge,Mag_barbulge,'d',ms=14,alpha=0.6)
        
    axt0.plot(n_line,SAna.AnalyticFunctions.linear_func_1D(np.log10(n_line),-8.5, -17.85),
              'b--',lw=4,label='best')

    #axt0.legend(loc="lower right")
    axt0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
   
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
            axt1.plot(mu_line, 
              SAna.AnalyticFunctions.linear_func_1D(mu_line,*popt_mu[i]),
              color=colour[i],linestyle='dashed',lw=3)
            
    print("Linear fit mu",*popt_mu)
    print("Linear fit n",*popt_n)
    print('===========================================')

    axt1.plot(mu_line, 
              SAna.AnalyticFunctions.linear_func_1D(mu_line,popt_mu_inone[0],popt_mu_inone[1]+5.2),
              color='b',linestyle='solid',lw=3,label='in_one')   
    
 
    axt1.plot(mu_line,SAna.AnalyticFunctions.linear_func_1D(mu_line,6.12, -93),
              'b--',lw=4, label='best')

    #axt1.plot(mu_barbulge,Mag_barbulge,'d',ms=14,alpha=0.6,label='bulge in barred gal')

    axt1.legend(loc="lower right")
    axt1.set_xlabel(r"$ \mu_\mathrm{0,i}$", fontsize=22)
    axt1.set_xlim(19.2,-0.5)

    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
    plt.show()
    return (fig)#,*popt_n,*popt_mu)
#%%
def plot_n_mu0_Mag_3plot(n,mu,Re,Mag,label=[],fit_instruc=0):
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
    fig = plt.figure(figsize=(6.4*3, 4.8))
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
            
    axt0.scatter(Sersic_n_corebulge,Abs_sph_mag_corebulge,
                 facecolors='none', edgecolors='r', s = 500)

    #axt0.legend(loc="lower right")
    axt0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
   
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
    axt1.scatter(mu0_corebulge,Abs_sph_mag_corebulge,
                 facecolors='none', edgecolors='r', s = 500)
    
    axt1.set_xlabel(r"$ \mu_\mathrm{0,i}~(\rm mag~arcsec^{-2})$", fontsize=22)
    axt1.set_xlim(19.2,-0.5)
    
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
    
    axt2.legend(loc="lower right")
    axt2.set_xlabel(r"$ R_\mathrm{e,Sph}~(\rm kpc)$", fontsize=22)
    axt2.set_xscale('log')
    axt2.set_xlim(0.09,120)
 
    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.setp(axt2.get_yticklabels(), visible=False)
    
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
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

#size_mass_Graham_equ(Mass,ML,n,A,B)
#print("fitting", curve_fit(SAna.AnalyticFunctions.size_mass_Graham_equ,(E_T11_K,MLR),Re_kpc))

#mass_gen = np.linspace()
#MLR_gen = np.linspace()
#Sersic_n_gen = np.

mass0 = np.linspace(2e7,0.5e13,50)
# Shen 2003 size-mass relatio
Shen2003 = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,0.347e-5,0.56)

# Lange 2016 size-mass relation from table 1
#Lange2016_E = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,7.4e-5,0.44)
#Lange2016_E_M1e10 = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,1.382,0.643)
#Lange2016_ETG_bulge = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,1.836,0.267)
Lange2016_E = SAna.AnalyticFunctions.size_mass_powerlaw(mass0/1e10,2.114,0.329)
Lange2016_E_M1e10 = SAna.AnalyticFunctions.size_mass_powerlaw(mass0/1e10,1.382,0.643)
Lange2016_ETG_bulge = SAna.AnalyticFunctions.size_mass_powerlaw(mass0/1e10,1.836,0.267)
Lange2016_final_bulge = SAna.AnalyticFunctions.size_mass_powerlaw(mass0/1e10,2.063,0.263)

def plot_sizemass_z0():
    fig = plt.figure(figsize=(6.4, 5.8))
    ax0 = plt.subplot()
    
    M_tot0 = 4.53 - 2.5*np.log10(mass0/np.repeat(np.median(MLR),len(mass0)))  

    ax0.plot(mass0,Shen2003,'-',label=r"$\rm Shen~et~al.~(2003)$",lw = 4)
    
    ax0.plot(mass0,Lange2016_E, '--', label=r"$\rm E~in~Lange~et~al.~(2016)$")
    ax0.plot(mass0,Lange2016_ETG_bulge, '--', label=r"$\rm ETG~bulge~in~Lange~et~al.~(2016)$")
    ax0.plot(mass0,Lange2016_E_M1e10, '--', label=r"$ \mathrm{E}(M_{*}/\mathrm{M_{\odot}} \gtrsim 10^{10})\rm ~in~Lange~et~al.~(2016)$")
    ax0.plot(mass0,Lange2016_final_bulge, '--', label=r"$ \rm Final~spheroids~in~Lange~et~al.~(2016)$")


    plot_dexter_sample_T11(E_BD, Re_kpc_BD,0,0,ax0,
                           alpha = 0.6,colour='#2bcbc9',
                          label=r"$\rm This~work~(BD)$",s=130)
    plot_dexter_sample_T11(E_1Sersic, Re_kpc_1Sersic,0,0,ax0,
                           alpha = 0.4,colour='#45b847',
                          label=r"$\rm This~work~(1Sersic)$",s=130)
    
    plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0,
                           alpha = 0.4,colour='#a5200b',
                          label=r"$\rm This~work$")
    
    #pp, ppv = curve_fit(SAna.AnalyticFunctions.size_mass_Graham_equ,(E_T11_K,MLR),Re_kpc)
    #func = SAna.AnalyticFunctions.size_mass_Graham_equ((mass0,np.repeat(np.average(MLR),len(mass0))),
    #                                                   -4.35029071, -19.13050166,0.59445763, -29.299029)
    func_cS = SAna.AnalyticFunctions.size_mass_Graham_equ(M_tot0,
                                                       -2.55158269, -20.94292143,0.10998085, -23.71212895)
    func_S = SAna.AnalyticFunctions.size_mass_Graham_equ(M_tot0,
                                                       -8.75, -17.85,6.12, -93)
#-9, -17.85,6.12, -113

    func_cS2 = SAna.AnalyticFunctions.size_mass_Graham_equ(M_tot0,
                                                       -2.55158269, -20.94292143,0.10998085, -23.71212895)
    func_S2 = SAna.AnalyticFunctions.size_mass_Graham_equ(M_tot0,
                                                       -8.8, -17.89,2.6, -59+5.2)
    #ax0.plot(mass0,func,'r-',label='Graham (2019)',lw = 4)
    ax0.plot(mass0,func_cS,'b-',label='Graham (2019)_g',lw = 4)
    ax0.plot(mass0,func_S,'b--',label='Graham (2019)_g',lw = 4)

    ax0.plot(mass0,func_cS2,'r-',label='Graham (2019)_f',lw = 4)
    ax0.plot(mass0,func_S2,'r--',label='Graham (2019)_f',lw = 4)
    
  
    # Single power law fit
    pp_3power, ppv_3power= curve_fit(SAna.AnalyticFunctions.linear_func_1D,np.log10(E_T11_K),np.log10(Re_kpc))    
    func3 = SAna.AnalyticFunctions.size_mass_powerlaw(mass0,10**pp_3power[1],pp_3power[0])
    ax0.plot(mass0,func3,'-',label='1-power')
    
    # two power law fit
    pp_2power, ppv_2power= curve_fit(SAna.AnalyticFunctions.size_mass_2powerlaw,E_T11_K,
                                     Re_kpc,p0 = [0.001,1.17,0.169,6.2e9],
                                     bounds = ((0.0009,1.1,0.14,6.19e9),(0.0011,1.18,0.17,6.21e9)),
                                     max_nfev = 100000)

                                     #np.log10(Re_kpc),p0 = [0.01,0.74,0.3,0.3e10],
                                     #bounds = ((0.009,0.70,0.29,0.29e10),(0.011,0.75,0.31,0.31e10)),
                                     #max_nfev = 100000)
                                     #bounds = ((0.009,0.7,0.029,0.3e10),(0.12,0.8,0.031,0.4e10)),max_nfev = 10000)
    print(*pp_2power)
    func2 = SAna.AnalyticFunctions.size_mass_2powerlaw(mass0,*pp_2power)
    ax0.plot(mass0,func2,'-',label='2-power',lw = 4)
      
    func22 =SAna.AnalyticFunctions.size_mass_2powerlaw(mass0,0.09,0.77,0.18,2.17e10)
    ax0.plot(mass0,func22,'-',label='2-power Lange 2015')
   
    
   #print('func',func)E_T11_K
    #plot_dexter_sample_T11(E_T11_K_E, Re_kpc_E,Re_kpc_err_E,mass_err_E,ax0,colour='k',alpha = 0.1, label=r"$\rm This~work$")
    #plot_dexter_sample_T11(E_T11_K_S0, Re_kpc_S0,Re_kpc_err_S0,mass_err_S0,ax0,colour='k',alpha = 0.1, label=r"$\rm This~work$")
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='g',alpha = 0.1,label=r"$\rm S~(Dust corrected)$")


    plot_dexter_sample_T11(Savorgnan_mass_36, Savorgnan_size_eq_kpc,0,0,ax0, 
                           alpha = 0.8,colour='#b940c8',
                           label=r"$\rm Savorgnan~et~al.~(2016)$", s=100)
    plot_dexter_sample_T11(Davis_mass_36, Davis_size_eq_kpc, 0, 0,ax0, 
                           alpha = 0.8,colour='#2e417b',
                           label=r"$\rm Davis~et~al.~(2019)$", s=100)
    plot_dexter_sample_T11(Sahu_mass_36, Sahu_size_eq_kpc, 0, 0,ax0,
                           alpha = 0.8,colour='#e1a000',
                          label=r"$\rm Sahu~et~al.~(2019)$", s =100)

    ax0.set_ylabel("$ R_\mathrm{e,Sph}$ (kpc)", fontsize=16)
    ax0.set_xlabel(r"$ M_\mathrm{*,Sph} / \rm M_\mathrm{\odot} (T11)$", fontsize=16)
    plt.legend(fontsize = 12.5,loc="lower right")
    plt.tight_layout()

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
    #plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
    
    #R_barbulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'],Re_kpc)
    #mass_barbulge = SSort.cherry_pick(S_bar[0]['PrimBar_index'], E_T11_K)
    
    plot_dexter_sample_T11(E_1Sersic, Re_kpc_1Sersic,0,0,ax0,
                           alpha = 0.4,colour='#45b847',
                          label=r"$\rm This~work~(1Sersic)$",s=130)
    plot_dexter_sample_T11(E_combine_morph[0], Re_kpc_combine_morph[0],Re_kpc_err_combine[0],mass_err_combine[0],ax0,colour='k',label = "E+ES")
    plot_dexter_sample_T11(E_combine_morph[1], Re_kpc_combine_morph[1],Re_kpc_err_combine[1],mass_err_combine[1],ax0,colour="#d20d0d",label="S0")
    plot_dexter_sample_T11(E_combine_morph[2], Re_kpc_combine_morph[2],Re_kpc_err_combine[2],mass_err_combine[2],ax0,colour="#ebb800",label="S")
    
    
   # ax0.plot(mass_barbulge,R_barbulge,'d',ms=14,alpha=0.7)
    
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S_old, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='g',label="S_dustcorr_old",marker='s')
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='y',label="S_dustcorr",marker='s')
    SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
                                                      "Barro et al. 2013", 
                                                      alpha0=0,AX=ax0)
    ax0.set_ylabel("$ R_\mathrm{e,sph}$ (kpc)", fontsize=16)
    ax0.set_xlabel(r"$ M_\mathrm{*,sph} / \rm M_\mathrm{\odot} (T11)$", fontsize=16)
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
   # ax0.plot(mass_barbulge,R_barbulge,'d',ms=14,alpha=0.7)
    
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S_old, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='g',label="S_dustcorr_old",marker='s')
    #plot_dexter_sample_T11(E_T11_K_dustcorr_S, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='y',label="S_dustcorr",marker='s')
    #SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", 
    #                                                  "Barro et al. 2013", 
    #                                                  alpha0=0,AX=ax0)
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel("$ \Sigma_\mathrm{Sph}$ (kpc)", fontsize=16)
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

    ax0.set_xlabel("$ R_\mathrm{e,sph}$ (kpc)", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=16)
    plt.legend(loc="lower right")
    plt.tight_layout()
    ax0.set_xscale('log')
    plt.gca().invert_yaxis()

    plt.show()
    
def plot_Mag_Re_morph_diverse(R1,R2,R3,Mag):
    fig = plt.figure(figsize=(6.4, 5.8))
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

    ax0.set_xlabel(r"$ R_\mathrm{z,sph}~\rm(kpc)$ ", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=16)

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
    
    ax0.plot(mu0_combine_morph[0], Mag[0],'o',ms = 10, color='#6833da',label = "z=0")
    ax0.plot(mu0_combine_morph[1], Mag[1],'o',ms = 10, color="#6833da")
    ax0.plot(mu0_combine_morph[2], Mag[2],'o',ms = 10, color="#6833da")
    
    ax0.set_xlabel(r"$\mu_\mathrm{z,Sph}$ ", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=16)
    ax0.set_xlim(35,-0.5)
    plt.legend(loc="lower right")
    plt.gca().invert_yaxis()

    plt.tight_layout()

    plt.show()
    
    
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

    #func_cS = SAna.AnalyticFunctions.size_mass_Graham_equ((mass1,np.repeat(1.0,len(mass1))),
    #                                                   -2.55158269, -20.94292143,0.10998085, -23.71212895)
    #func_S = SAna.AnalyticFunctions.size_mass_Graham_equ((mass1,np.repeat(1.0,len(mass1))),
    #                                                   -8.75, -17.85,6.12, -93)

    #func_cS2 = SAna.AnalyticFunctions.size_mass_Graham_equ(Mag00,
    #                                                  -2.55158269, -20.94292143,0.10998085, -23.71212895)
    func_S2 = SAna.AnalyticFunctions.size_mass_Graham_equ(Mag00,
                                                       -8.8, -17.89,2.6, -59+5.2)

    #ax0.plot(Mag00,func_cS,'b-',label='Graham (2019)_g',lw = 4)
    #ax0.plot(Mag00,func_S,'b--',label='Graham (2019)_g',lw = 4)

    #ax0.plot(func_cS2,Mag00,'r-',label='Graham (2019)_f',lw = 4)
    ax0.plot(func_S2,Mag00,'r--',label='Graham (2019)_f',lw = 4)
    
    ax0.plot(Re_kpc,Abs_sph_mag,'o',alpha = 0.4,color='#a5200b',
                          label=r"$\rm This~work$",ms=10)
    
    ax0.set_xlabel(r"$R_\mathrm{e,Sph}\rm ~(kpc)$ ", fontsize=16)
    ax0.set_ylabel(r"$\mathfrak{M}_{i,Sph}$", fontsize=16)
    
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.ylim(0,100)
    #plt.xlim(-26,-16)
    ax0.set_xscale('log')
    plt.show()
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

def list_para():
    for i in range(len(name)):
        print(name[i],total_mag[i])#new_mu_list[i],Re[i],Sersic_n[i],sph_mag[i],ScStype[i])

#list_para()
# plot the Mag vs n and mu0 plot
#plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[r"$type~1$",r"$type~2$"])

n_combine_morph, mu0_combine_morph, Mag_combine_morph, E_combine_morph, Re_kpc_combine_morph, Re_kpc_err_combine, mass_err_combine, Abs_total_mag_combine_morph =seperation_morph(
    sph_mag, Re, Sersic_n, mu_e,total_mag, 
                     name, mu0, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag, elle, E_T11_K,
                     MLR, mass_uerr,Abs_total_mag)

n_combine_morph_BD, mu0_combine_morph_BD, Mag_combine_morph_BD, E_combine_morph_BD, Re_kpc_combine_morph_BD, Re_kpc_err_combine_BD, mass_err_combine_BD, Abs_total_mag_combine_morph_BD = seperation_morph(
    sph_mag_BD, Re_BD, Sersic_n_BD, mu_e_BD,total_mag_BD, 
                     name, mu0_BD, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag_BD, elle, E_BD,
                     MLR, mass_uerr,Abs_total_mag_BD)

n_combine_morph_1Sersic, mu0_combine_morph_1Sersic, Mag_combine_morph_1Sersic, E_combine_morph_1Sersic, Re_kpc_combine_morph_1Sersic, Re_kpc_err_combine_1Sersic, mass_err_combine_1Sersic, Abs_total_mag_combine_morph_1Sersic = seperation_morph(
    sph_mag_1Sersic, Re_1Sersic, Sersic_n_1Sersic, mu_e_1Sersic,total_mag_1Sersic, 
                     name, mu0_1Sersic, mag_g_kcorr, mag_i_kcorr, D,
                     D_lerr, D_uerr, Abs_sph_mag_1Sersic, elle, E_1Sersic,
                     MLR, mass_uerr,Abs_total_mag_1Sersic)

#A = plot_n_mu0_Mag_2plot(n_combine,mu0_combine,
#                         total_mag_combine,label=[r"$\rm Core-S\'{e}rsic$",
#                                        r"$\rm S\'{e}rsic$"])
print('cS vs S, multi')
#B = plot_n_mu0_Mag_2plot(n_combine,mu0_combine,
#                     Mag_combine,label=[r"$\rm Core-S\'{e}rsic$",
#                                        r"$\rm S\'{e}rsic$"])

print('morph, multi')
#C = plot_n_mu0_Mag_2plot(n_combine_morph,mu0_combine_morph,
#                     Mag_combine_morph,label=[r"$\rm E+ES$",
#                                        r"$\rm S0$", r"$\rm S$"])
CC = plot_n_mu0_Mag_3plot(n_combine_morph,mu0_combine_morph, Re_kpc_combine_morph,
                     Mag_combine_morph,label=[r"$\rm E+ES$",
                                        r"$\rm S0$", r"$\rm S$"])
print('morph, B+D')

#DD = plot_n_mu0_Mag_2plot(n_combine_morph_BD,mu0_combine_morph_BD,
#                     Mag_combine_morph_BD,label=[r"$\rm E+ES$",
#                                        r"$\rm S0$", r"$\rm S$"])
print('morph, 1-Sersic')

#E = plot_n_mu0_Mag_2plot(n_combine_morph_1Sersic,mu0_combine_morph_1Sersic,
#                     Abs_total_mag_combine_morph_1Sersic,label=[r"$\rm E+ES$",
#                                        r"$\rm S0$", r"$\rm S$"])
print('====================')


# D = plot_n_mu0_Mag_2plot(n_combine_ELtype,mu0_combine_ELtype,
#                     Mag_combine_ELtype,label=[r"$\rm ETG$", r"$\rm LTG$"])
plot_Mag_Re_morph(Re_kpc_combine_morph,Mag_combine_morph)
#plot_Mag_Re_morph(Re_kpc_combine_morph_BD,Abs_total_mag_combine_morph)



plot_sizemass_morph()
#plot_densmass_morph()

#
#plot_sizemass_z0()

#calculate bn Re and mu 10 and 90 percent
Re_kpc_10_combine_morph, Re_kpc_90_combine_morph, mu0_10_combine_morph, mu0_50_combine_morph, mu0_90_combine_morph = get_bn_Re_mu_1090()
plot_Mag_Re_morph_diverse(Re_kpc_10_combine_morph,Re_kpc_combine_morph,Re_kpc_90_combine_morph,Mag_combine_morph)  
plot_Mag_mu_morph_diverse(mu0_10_combine_morph,mu0_50_combine_morph,mu0_90_combine_morph,Mag_combine_morph)

#plot_size_Mag_z0()

#plot_dustvsnodust()
#plot_dustvsnodust_g()
#plot_dustvsnodust_i()
#plot_dustvsnodust_BD_g()
#plot_dustvsnodust_BD_i()

#plot_dustvsnodust_BD_gi()


#plot_sphmag_g()
#plot_discmag_g()
#plot_sphmag_i()
#plot_discmag_i()