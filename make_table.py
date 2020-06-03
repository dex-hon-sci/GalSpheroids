#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 03:02:36 2020

@author: dexter
"""


import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import numpy as np
import pickle

M_sun = 4.53

#/home/dexter/result/stat/Into_selection_all
#/home/dexter/result/distance/
#vel_disp_list_all_mag.txt

# Name & RA (2000) & Dec (2000) & $Dist$  & $g$ & $i$ & scale  & Seeing & $\sigma_{velocity}$ & $morph$ & $Mag_{i-band}$ & \log $(M_*/M_\odot)_{total}$\\

###############



    
#Bin4cpt = SRead.read_list("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")
name = SRead.grab_name("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")
distance_file = SRead.grab_dist(
    "/home/dexter/result/distance/dir_dist_list.txt",
    "/home/dexter/result/distance/dist_list_bin4.txt")

vdis_file = SRead.read_table(
    "/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt")
vdis_file2 = SRead.read_table(
    "/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt",dtype='str')

seeing_file1 = SRead.read_table(
    "/home/dexter/result/seeing_correction/seeing_list.txt",dtype='str')
seeing_file2 = SRead.read_table(
    "/home/dexter/result/seeing_correction/seeing_list.txt")

RA = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,6])
DEC = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,7])

Dist = distance_file["Dist"]
Dist_err = distance_file["Dist_err"]
mag_g = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,10]))
mag_i = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,9]))
scale = distance_file["Scale"]

print('---------')
seeing = SRead.extract_match(name, seeing_file1[:,0] ,seeing_file2[:,2])
print('---------')

v_dis = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,3])
v_dis_err = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,4])


morph = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file2[:,8])


Mag_i = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,12]))

ML_select = SPlot.MLRelationIband(mag_g,mag_i).Taylor11_MassRatio
#M_4 = SPlot.MassCalculation(mag_i, Dist, M_sun,mag_g,mag_i)



Dist_final = SSort.str_zipping(np.around(Dist,decimals=2), 
                               np.around(Dist_err,decimals=2), 
                               zip_symbol="$\pm$")
v_dis_final = SSort.str_zipping(np.around(v_dis, decimals=2),
                          np.around(v_dis_err, decimals=2),
                          zip_symbol="$\pm$")


M_sun_array = np.full(len(Mag_i), M_sun)

mass = ML_select*(10**((M_sun_array-Mag_i)/2.5)) 

table = {"Name": name,
         "RA": np.around(RA, decimals=2),
         "DEC": np.around(DEC, decimals=2),
         "Dist": Dist_final,
         "mag_g": np.around(mag_g, decimals=2) ,
         "mag_i": np.around(mag_i, decimals=2) ,
         "Scale": np.around(scale, decimals=2) ,
         "seeing": np.around(seeing, decimals=2),
         "v_dis": v_dis_final,
         "morph": morph,
         "Mag_i": np.around(Mag_i, decimals=2),
         "mass": np.around(mass/1e10, decimals=2)}

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table1_bin4_Tmass", 'wb') as f:
    pickle.dump(table, f)

#SRead.convert_dict_ascii("Gal_table1_bin3","Gal_table1_bin3.txt")

##################################
Bundle = "/home/dexter/result/Gal_bundle_equvi_bin4_cpt"
distance_file = SRead.grab_dist(
    "/home/dexter/result/distance/dir_dist_list.txt",
    "/home/dexter/result/distance/dist_list_bin4.txt")
name = SRead.grab_name(Bundle)

mu_e= SRead.grab_parameter(Bundle, 
                     ["Bulge", "CoreBulge"], 0)

scale = distance_file["Scale"]


R_e= SRead.grab_parameter(Bundle, 
                     ["Bulge", "CoreBulge"], 1)*scale

n= SRead.grab_parameter(Bundle, 
                     ["Bulge", "CoreBulge"], 2)

mag_sph = SRead.grab_mag(Bundle, ["Bulge", "CoreBulge"])

mag_g = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,10]))
mag_i = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,9]))


Dist = distance_file["Dist"]
Dist_err = distance_file["Dist_err"]

Mag_sph=mag_sph-25-5*np.log10(Dist)  

Mag_sph_plus = (mag_sph-25-5*np.log10(Dist+Dist_err))-Mag_sph
Mag_sph_minus = (mag_sph-25-5*np.log10(Dist-Dist_err))-Mag_sph

print(Dist_err)

print(Dist+Dist_err, Dist-Dist_err)

print(Mag_sph_plus-Mag_sph_minus)

ML_select_T = SPlot.MLRelationIband(mag_g,mag_i).Taylor11_MassRatio
ML_select_Z = SPlot.MLRelationIband(mag_g,mag_i).Zibetti09_MassRatio
ML_select_R = SPlot.MLRelationIband(mag_g,mag_i).Roediger15BC03_MassRatio
ML_select_I = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio


Taylor_mass = ML_select_T*(10**((4.53-Mag_sph)/2.5))
Zibetti_mass = ML_select_Z*(10**((4.53-Mag_sph)/2.5))
Roediger_mass = ML_select_R*(10**((4.53-Mag_sph)/2.5))
Into_mass = ML_select_I*(10**((4.53-Mag_sph)/2.5))


table = {"Name": name,
         "mu_e": np.around(mu_e, decimals=2),
         "R_e": np.around(R_e, decimals=2),
         "n": np.around(n, decimals=2),
         "mag_sph": np.around(mag_sph, decimals=2),
         "Mag_sph": np.around(Mag_sph, decimals=2),
         "Taylor_mass": np.around(np.log10(Taylor_mass), decimals=2),
         "Zibetti_mass": np.around(np.log10(Zibetti_mass), decimals=2), 
         "Roediger_mass": np.around(np.log10(Roediger_mass), decimals=2),
         "Into_mass": np.around(np.log10(Into_mass), decimals=2)} 

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table2_bin4", 'wb') as f:
    pickle.dump(table, f)

#SRead.convert_dict_ascii("Gal_table2_bin2","Gal_table2_bin4.txt")
#Mag_i = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")


#ML_select_4 = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
#M_4 = SPlot.MassCalculation(Mag_i, Dist, M_sun,mag_g,mag_i)
#mass_4 = M_4.cal_Mass(ML_select_4)


#   Name & $\mu_e$ & $R_e$ & $n$  & $mag_{Sph}$ & $Abs Mag_{Sph}$ & \log$(M_*/M_\odot)_{Taylor}$  & log$(M_*/M_\odot)_{Into}$ & log$(M_*/M_\odot)$ & log$(M_*/M_\odot)$  \\
