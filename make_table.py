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
vdis_file = SRead.read_table("/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt")
vdis_file2 = SRead.read_table("/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt",dtype='str')
distance_file = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list_bin4.txt")
seeing_file1 = SRead.read_table("/home/dexter/result/seeing_correction/seeing_list.txt",dtype='str')
seeing_file2 = SRead.read_table("/home/dexter/result/seeing_correction/seeing_list.txt")

RA = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,6])
DEC = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,7])

Dist = distance_file["Dist"]
Dist_err = distance_file["Dist_err"]
mag_g = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,9]))
mag_i = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,10]))
scale = distance_file["Scale"]
seeing = SRead.extract_match(name, seeing_file1[:,0] ,seeing_file2[:,2])

v_dis = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,3])
v_dis_err = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,4])


morph = SRead.extract_match(name, vdis_file2[:,0] ,vdis_file2[:,8])


Mag_i = np.array(SRead.extract_match(name, vdis_file2[:,0] ,vdis_file[:,12]))


ML_select_4 = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
#M_4 = SPlot.MassCalculation(mag_i, Dist, M_sun,mag_g,mag_i)

print(ML_select_4)
M_sun_array = np.full(len(Mag_i), M_sun)
print(M_sun_array)
print(Mag_i)
mass_4 = ML_select_4*(10**((M_sun_array-Mag_i)/2.5)) 

print(mass_4)
table4 = {"Name": name,
         "RA": np.around(RA, decimals=2),
         "DEC": np.around(DEC, decimals=2),
         "Dist":np.around(Dist,decimals=2),
         "Dist_err":np.around(Dist_err,decimals=2),
         "mag_g": np.around(mag_g, decimals=2) ,
         "mag_i": np.around(mag_i, decimals=2) ,
         "Scale": np.around(scale, decimals=2) ,
         "seeing": np.around(seeing, decimals=2),
         "v_dis": np.around(v_dis, decimals=2) ,
         "v_dis_err": np.around(v_dis_err, decimals=2),
         "morph": morph,
         "Mag_i": np.around(Mag_i, decimals=2),
         "mass": np.around(mass_4/1e10, decimals=2)}

value = list(table4.values())
key = list(table4.keys())

for i in range(len(list(table4.keys()))):
    print(key[i],len(value[i]))
    
with open("Gal_table1_bin4", 'wb') as f:
    pickle.dump(table4, f)

print(table4)

SRead.convert_dict_ascii("Gal_table1_bin4","Gal_table1_bin4_2.txt")
#A = SRead.read_list("Gal_vdis_dict_equvi")

#mag_g,mag_i = A["mag_g"], A["mag_i"]

#####second table#####

#Mag_i = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")


#ML_select_4 = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
#M_4 = SPlot.MassCalculation(Mag_i, Dist, M_sun,mag_g,mag_i)
#mass_4 = M_4.cal_Mass(ML_select_4)


#   Name & $\mu_e$ & $R_e$ & $n$  & $mag_{Sph}$ & $Abs Mag_{Sph}$ & \log$(M_*/M_\odot)_{Taylor}$  & log$(M_*/M_\odot)_{Into}$ & log$(M_*/M_\odot)$ & log$(M_*/M_\odot)$  \\
