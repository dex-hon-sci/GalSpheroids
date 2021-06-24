#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:14:07 2021

@author: dexter

This is a script dedicated to
visualise the spheroid parameters.

This script produce the following plots:
0) stacked radial profile, as well as listing the Sersic mu0 (Done) 
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

bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
bundle = SRead.read_list(bundle_name)

name = geom_file_n[:,0]
mu0 = geom_file[:,6]

# Read the data from the galaxy bundle
sph_mag = SRead.grab_mag("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"])
Re = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 1) 
Sersic_n = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 2) 
mu_e = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 0) 

core_sersic_mu_p = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["CoreBulge"], 0)

total_mag = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_cpt")

# Get the distance from the parent sample
D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")
D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt",
    dtype = "str")

K_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin_all_Kcorr.dat")

mag_g, mag_i = D0_all_table[:,11], D0_all_table[:,10] #*Galaxy colours
mag_g_kcorr, mag_i_kcorr = K_table[:,19], K_table[:,18] #Galaxy colours with K correction

D, D_lerr, D_uerr = D0_all_table[:,29], D0_all_table[:,30], D0_all_table[:,31]

morph = D0_all_table_n[:,-1]

# Calculate the absoulte magnitude and the stellar mass
ML_select_T11 = SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Taylor11_MassRatio
M = SPlot.MassCalculation(sph_mag, D, 4.53,mag_g_kcorr,mag_i_kcorr)
Abs_sph_mag = M.cal_abs_mag()
E_T11_K = M.cal_Mass(ML_select_T11)

# Calculate the error of the Re and stellar mass
MLR = ML_select_T11
MLR_e = 10**0.1    
mag_e = 0.3 #magnitude error

mass_uerr = np.sqrt(((mag_e/2.5)**2)+((2*D_uerr/(D*np.log(10)))**2)+((MLR_e/(MLR*np.log(10)))**2))
mass_err = mass_uerr

# calculate the scale
ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale
scale = D* ars

scale_lerr, scale_uerr = (D-D_lerr)*ars, (D+D_uerr)*ars
Re_kpc = Re* scale

Re_kpc_lerr, Re_kpc_uerr = abs(Re_kpc - Re* scale_lerr) , abs(Re* scale_uerr - Re_kpc)
Re_kpc_err =[Re_kpc_lerr, Re_kpc_uerr]

## End in reading base data#################################################

# Read in the data from the others
core = SRead.grab_parameter_whole(bundle_name, ["CoreBulge"])
#######Seperation 1###########################################################
#core vs not core seperation here
S = SSort.seperator_label_generic(bundle, ["Bulge","CoreBulge"])

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

# Extract the info from Core Bulge only
sph_corebulge_mag = list(SRead.grab_mag(S[1]["CoreBulge"], ["CoreBulge"]))
Re_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 1))
Sersic_n_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 2)) 
mu_e_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 0))
total_mag_corebulge = list(SRead.grab_total_mag(S[1]["CoreBulge"]))

# cherry pick info based on the index list (output list)
name_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], name)
mu0_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mu0)

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

Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge = abs(Re_kpc_corebulge - Re_corebulge* scale_lerr_corebulge) , abs(Re_corebulge* scale_uerr_corebulge - Re_kpc_corebulge)
Re_kpc_err_corebulge =[Re_kpc_lerr_corebulge, Re_kpc_uerr_corebulge]

# put the two array together  
n_combine = np.array([Sersic_n_bulge,Sersic_n_corebulge])
mu0_combine = np.array([mu0_bulge,mu0_corebulge])
Mag_combine = np.array([Abs_sph_mag_bulge,Abs_sph_mag_corebulge])


##############End Seperation 1#################################################

####Start Seperation 2 ########################################################
index = np.arange(103)

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

mag_g_E = SSort.cherry_pick(morph_dict['E'], mag_g)
mag_g_S0 = SSort.cherry_pick(morph_dict['S0'], mag_g)
mag_g_S = SSort.cherry_pick(morph_dict['S'], mag_g)
mag_i_E = SSort.cherry_pick(morph_dict['E'], mag_i)
mag_i_S0 = SSort.cherry_pick(morph_dict['S0'], mag_i)
mag_i_S = SSort.cherry_pick(morph_dict['S'], mag_i)

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

E_T11_K_E = SSort.cherry_pick(morph_dict['E'], E_T11_K)
E_T11_K_S0 = SSort.cherry_pick(morph_dict['S0'], E_T11_K)
E_T11_K_S = SSort.cherry_pick(morph_dict['S'], E_T11_K)

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

# put the Three array together  
n_combine_morph = np.array([Sersic_n_E,Sersic_n_S0,Sersic_n_S])
mu0_combine_morph = np.array([mu0_E,mu0_S0,mu0_S])
Mag_combine_morph = np.array([Abs_sph_mag_E,Abs_sph_mag_S0, Abs_sph_mag_S])

##############End Seperation 2#################################################

####Start Seperation 3 ########################################################
# Seperate the sample into LTG and ETG
sph_mag_ETG = sph_mag_E + sph_mag_S0
sph_mag_LTG = sph_mag_S

Re_ETG = Re_E + Re_S0 
Re_LTG = Re_S 

Sersic_n_ETG = Sersic_n_E + Sersic_n_S0 
Sersic_n_LTG = Sersic_n_S 

mu_e_ETG = mu_e_E + mu_e_S0 
mu_e_LTG = mu_e_S 

total_mag_ETG = total_mag_E + total_mag_S0 
total_mag_LTG = total_mag_S 

name_ETG = name_E+name_S0 
name_LTG = name_S 

mu0_ETG = mu0_E + mu0_S0 
mu0_LTG = mu0_S 

mag_g_ETG = mag_g_E + mag_g_S0 
mag_g_LTG = mag_g_S 
mag_i_ETG = mag_i_E + mag_i_S0
mag_i_LTG = mag_i_S 

mag_g_kcorr_ETG = mag_g_kcorr_E + mag_g_kcorr_S0
mag_g_kcorr_LTG = mag_g_kcorr_S 
mag_i_kcorr_ETG = mag_i_kcorr_E + mag_i_kcorr_S0 
mag_i_kcorr_LTG = mag_i_kcorr_S 

D_ETG = D_E + D_S0 
D_LTG = D_S 

D_lerr_ETG = D_lerr_E + D_lerr_S0
D_lerr_LTG  = D_lerr_S 
D_uerr_ETG = D_uerr_E + D_uerr_S0 
D_uerr_LTG = D_uerr_S 

Abs_sph_mag_ETG = Abs_sph_mag_E + Abs_sph_mag_S0
Abs_sph_mag_LTG = Abs_sph_mag_S

E_T11_K_ETG = E_T11_K_E + E_T11_K_S0
E_T11_K_LTG = E_T11_K_S 

MLR_ETG = MLR_E + MLR_S0
MLR_LTG = MLR_S 

mass_uerr_ETG = mass_uerr_E +mass_uerr_S0 
mass_uerr_LTG = mass_uerr_S 
mass_err_ETG = np.array(mass_uerr_ETG)
mass_err_LTG = np.array(mass_uerr_LTG)

scale_ETG = np.array(D_ETG)*ars
scale_LTG = np.array(D_LTG)*ars

scale_lerr_ETG, scale_uerr_ETG = (np.array(D_ETG)-np.array(D_lerr_ETG))*ars, (np.array(D_ETG)+np.array(D_uerr_ETG))*ars
scale_lerr_LTG, scale_uerr_LTG = (np.array(D_LTG)-np.array(D_lerr_LTG))*ars, (np.array(D_LTG)+np.array(D_uerr_LTG))*ars

Re_kpc_ETG = Re_ETG* scale_ETG
Re_kpc_LTG = Re_LTG* scale_LTG

Re_kpc_lerr_ETG, Re_kpc_uerr_ETG = abs(Re_kpc_ETG - Re_ETG* scale_lerr_ETG) , abs(Re_ETG* scale_uerr_ETG - Re_kpc_ETG)
Re_kpc_lerr_LTG, Re_kpc_uerr_LTG = abs(Re_kpc_LTG - Re_LTG* scale_lerr_LTG) , abs(Re_LTG* scale_uerr_LTG - Re_kpc_LTG)

Re_kpc_err_ETG =[Re_kpc_lerr_ETG, Re_kpc_uerr_ETG]
Re_kpc_err_LTG =[Re_kpc_lerr_LTG, Re_kpc_uerr_LTG]

# put the Three array together  
n_combine_ELtype = np.array([Sersic_n_ETG,Sersic_n_LTG])
mu0_combine_ELtype = np.array([mu0_ETG,mu0_LTG])
Mag_combine_ELtype = np.array([Abs_sph_mag_ETG,Abs_sph_mag_LTG])
##############End Seperation 3#################################################

# investigation
# M_i > -20
# weird core-Sersic gal, NGC4382, NGC4636
# COre+ disk nothing special??
Q = SSort.selection_generic(name_corebulge, Abs_sph_mag_corebulge, np.repeat(-20,len(name_corebulge)),
                            direction="high")
Q2 = SSort.selection_generic(sph_corebulge_mag, Abs_sph_mag_corebulge, np.repeat(-20,len(name_corebulge)),
                            direction="high")
print(Q,Q2)

print(D[12],D[14])

# Find 3.2,4.6, 4e10 5.2e10 S0

# weird S gal, high bugle size and mass NGC3270, paticular high size. nuclear cpt
Q3 = SSort.zone_in_2D(np.array([list(E_T11_K_S),list(Re_kpc_S)]),[4e10,3.7],[5.2e10,4.6])

print(Q3)
print(name[Q3["index"]])

#weird dim S0, NGC3665 nuclear **Dusty disk
Q4 = SSort.zone_in_2D([Sersic_n_S0,Abs_sph_mag_S0],[0,-18],[1,-16])

print(Q4,name[Q4["index"]])
############End reading 

#%%
def plot_stack_surface_brightness_profile(r):
    bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
    bundle=SRead.read_list(bundle_name)
    name = SRead.grab_name(bundle_name)
    
    fig = plt.figure()
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
                mu_e, n = para[0], para[2]
                b = SAna.get_bn(n) # get_bn is the exact function used in profiler
                # obtain mu_0 by equ7 from Graham2005
                new_mu = mu_e-2.5*b/np.log(10)
                #print(name[i],new_mu) # list mu0
                new_mu_list.append(new_mu) # save the new mu_0
                
                #list name, first element of from the radial SB, and the mu_0 
                #from equ 7
                print(name[i],line[0],new_mu) 
                pass
            elif identifier == "CoreBulge": #if the spheroid is a Core Sersic bulge
                para = bundle[i][j+1]
                line = SAna.AnalyticFunctions.mu_core_sersic_func(r,*list(para))
                plt.plot(r,line,"k-",lw= 3)
                new_mu_list.append(line[0])

                #print(name[i],line[0]) #list mu0
                pass
    plt.gca().invert_yaxis()
    #plt.ylim(np.log10(24),np.log(10))
    plt.ylim(24,12)
    plt.xlim(0,120)
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
def plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[],fit_instruc=0):
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
    mu_line = np.linspace(5, 30, 100)

    #Fitting the relevant data
    popt_mu,pcov_mu = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                mu0[fit_instruc], Mag[fit_instruc])
    popt_n,pcov_n = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                              np.log10(n[fit_instruc]), Mag[fit_instruc])
  
    #Plotting
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=2, nrows=1,
                               hspace=0, wspace=0.0) 
    # Define colour and markers
    markers = ["o","X", "s","+","*"]
    colour = ["k","#d20d0d","#ebb800","#0f920f","#0f4592"]
    
    # The Mag vs Sersic n plot
    axt0 = plt.subplot(gs[0])
    axt0.plot(10**n_line, 
              SAna.AnalyticFunctions.linear_func_1D(n_line,*popt_n),'b--',lw=4)
    # Check the dimension of the input array 
    if n.shape[0] != len(label):
        axt0.plot(n,Mag,'ko',ms=10)
    elif n.shape[0] == len(label):
        for i in range(len(n)):
            axt0.plot(n[i],Mag[i], linestyle="None",
                      marker=markers[i], color = colour[i],
                      ms=10,label=label[i])
    
    #axt0.legend(loc="lower right")
    axt0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
   
    # The Mag vs mu_0 plot
    axt1 = plt.subplot(gs[1],sharey=axt0)
    axt1.plot(mu_line, 
              SAna.AnalyticFunctions.linear_func_1D(mu_line,*popt_mu),'b--',lw=4)
   
    # Check the dimension of the input array  
    if mu0.shape[0] != len(label):
        axt1.plot(mu0,Mag,'ko',ms=10)
    elif mu0.shape[0] == len(label):
        for i in range(len(mu0)):
            axt1.plot(mu0[i],Mag[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i])    
    print("Linear fit mu",*popt_mu)
    print("Linear fit n",*popt_n)
    
    axt1.legend(loc="lower right")
    axt1.set_xlabel(r"$ \mu_\mathrm{0,i}$", fontsize=22)
    axt1.set_xlim(19.2,12.2)

    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
    plt.show()
    return (fig,*popt_n,*popt_mu)

#%%
xlim = [3e8,1.3e12]
ylim = [0.08,167]
# plot size-mass diagram
def plot_dexter_sample_T11(mass, size, size_err,mass_err, 
                           A, scale='log',alpha=0.65,colour='#a5200b',label="This work"):
    A.scatter(mass, size,marker='o',c=colour,label=label, 
              s =70, alpha=0.7)
    A.errorbar(mass, size, yerr = size_err, 
                  xerr = mass_err*mass, ls='none',linewidth=4, 
                  color = colour,
                  ecolor= colour, capsize=0, alpha=alpha, marker='o')

    
    A.set_xlim(left = xlim[0], right = xlim[1])
    A.set_ylim(bottom = ylim[0], top = ylim[1])
       
    A.set_xscale(scale)
    A.set_yscale(scale)
    
#%% Execution Area

# extrapolate the central surface brightness
R_gen = np.linspace(0,300,300*2)

# produce the stacked radial profile figure, as well as the mu0 
stack = plot_stack_surface_brightness_profile(R_gen)
list_mu0_extrapolate(stack[1])
# plot the Mag vs n and mu0 plot
#plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[r"$type~1$",r"$type~2$"])

#plot_n_mu0_Mag_2plot(Sersic_n,mu0,Abs_sph_mag)
B = plot_n_mu0_Mag_2plot(n_combine,mu0_combine,
                     Mag_combine,label=[r"$\rm S\'{e}rsic$",
                                        r"$\rm Core-S\'{e}rsic$"])
#,label=[r"$type~1$",r"$type~2$"])


C = plot_n_mu0_Mag_2plot(n_combine_morph,mu0_combine_morph,
                     Mag_combine_morph,label=[r"$\rm E$",
                                        r"$\rm S0$", r"$\rm S$"])

D = plot_n_mu0_Mag_2plot(n_combine_ELtype,mu0_combine_ELtype,
                     Mag_combine_ELtype,label=[r"$\rm ETG$", r"$\rm LTG$"])

fig = plt.figure()
ax0 = plt.subplot()
#plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
plot_dexter_sample_T11(E_T11_K_bulge, Re_kpc_bulge,Re_kpc_err_bulge,mass_err_bulge,ax0,colour='k',label=r"$\rm S\'{e}rsic$")
plot_dexter_sample_T11(E_T11_K_corebulge, Re_kpc_corebulge,Re_kpc_err_corebulge,mass_err_corebulge,ax0,label=r"$\rm Core-S\'{e}rsic$")
ax0.set_ylabel("$ R_{e,sph}$ (kpc)", fontsize=16)
ax0.set_xlabel(r"$ M_{*,sph} / \rm M_{\odot} $", fontsize=16)
plt.legend(loc="lower right")
plt.show()


fig = plt.figure()
ax0 = plt.subplot()
#plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
plot_dexter_sample_T11(E_T11_K_E, Re_kpc_E,Re_kpc_err_E,mass_err_E,ax0,colour='k',label = "E")
plot_dexter_sample_T11(E_T11_K_S0, Re_kpc_S0,Re_kpc_err_S0,mass_err_S0,ax0,colour='r',label="S0")
plot_dexter_sample_T11(E_T11_K_S, Re_kpc_S,Re_kpc_err_S,mass_err_S,ax0,colour='b',label="S")
ax0.set_ylabel("$ R_{e,sph}$ (kpc)", fontsize=16)
ax0.set_xlabel(r"$ M_{*,sph} / \rm M_{\odot} $", fontsize=16)
plt.legend(loc="lower right")
plt.show()

fig = plt.figure()
ax0 = plt.subplot()
#plot_dexter_sample_T11(E_T11_K, Re_kpc,Re_kpc_err,mass_err,ax0)
plot_dexter_sample_T11(E_T11_K_ETG, Re_kpc_ETG,Re_kpc_err_ETG,mass_err_ETG,ax0,colour='k',label = "ETG")
plot_dexter_sample_T11(E_T11_K_LTG, Re_kpc_LTG,Re_kpc_err_LTG,mass_err_LTG,ax0,colour='r',label="LTG")
ax0.set_ylabel("$ R_{e,sph}$ (kpc)", fontsize=16)
ax0.set_xlabel(r"$ M_{*,sph} / \rm M_{\odot} $", fontsize=16)
plt.legend(loc="lower right")
plt.show()