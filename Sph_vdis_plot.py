#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:16:32 2020

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import matplotlib.pyplot as plt
import numpy as np

M_sun = 4.53


# creating the velocity dispersion dictionary and calculating the stellar mass of the spheroid 
## If you need further investigation on the distance list and the velocity dispersion sources
# Just let me know, I can send the tables to you
#%% 
A = SSort.vdis_match("Gal_bundle_equvi_cpt", 
                 "./velocity_disp/vel_disp_list_all_mag.txt", 
                 ["./distance/dir_dist_list.txt","./distance/dist_list.txt"],
                 "Gal_vdis_dict_equvi")

B = SSort.LTG_ETG_seperator("Gal_bundle_equvi_cpt",
                            "Gal_bundle_equvi_cpt_LTG",
                            "Gal_bundle_equvi_cpt_ETG") 

# removing and adding the exception, NGC4334 is a LTG
BBE =SSort.remove_item("Gal_bundle_equvi_cpt_ETG",["NGC4334"])
BBL = SSort.add_item("Gal_bundle_equvi_cpt", "Gal_bundle_equvi_cpt_LTG", ["NGC4334"])

BL = SSort.vdis_match("Gal_bundle_equvi_cpt_LTG", 
                 "./velocity_disp/vel_disp_list_all_mag.txt", 
                 ["./distance/dir_dist_list.txt","./distance/dist_list.txt"],
                 "Gal_vdis_dict_equvi_LTG")

BE = SSort.vdis_match("Gal_bundle_equvi_cpt_ETG", 
                 "./velocity_disp/vel_disp_list_all_mag.txt", 
                 ["./distance/dir_dist_list.txt","./distance/dist_list.txt"],
                 "Gal_vdis_dict_equvi_ETG")


#%%
#calculate the stellar mass, base on two different M/L realtion
ML_select_1 = SPlot.MLRelationIband(A["mag_g"],A["mag_i"]).Into13_MassRatio
M_1 = SPlot.MassCalculation(A["sph_mag"], A["Dist"], M_sun,A["mag_g"],A["mag_i"])
mass_1 = M_1.cal_Mass(ML_select_1)

T_M_1 = SPlot.MassCalculation(A["total_mag"], A["Dist"], M_sun,A["mag_g"],A["mag_i"])
T_mass_1 = T_M_1.cal_Mass(ML_select_1)

ML_select_2 = SPlot.MLRelationIband(A["mag_g"],A["mag_i"]).Taylor11_MassRatio
M_2 = SPlot.MassCalculation(A["sph_mag"], A["Dist"], M_sun,A["mag_g"],A["mag_i"])
mass_2 = M_2.cal_Mass(ML_select_2)

T_M_2 = SPlot.MassCalculation(A["total_mag"], A["Dist"], M_sun,A["mag_g"],A["mag_i"])
T_mass_2 = T_M_2.cal_Mass(ML_select_2)


##

ML_select_BL_1 = SPlot.MLRelationIband(BL["mag_g"],BL["mag_i"]).Into13_MassRatio
M_BL_1 = SPlot.MassCalculation(BL["sph_mag"], BL["Dist"], M_sun,BL["mag_g"],BL["mag_i"])
mass_BL_1 = M_BL_1.cal_Mass(ML_select_BL_1)

T_M_BL_1 = SPlot.MassCalculation(BL["total_mag"], BL["Dist"], M_sun,BL["mag_g"],BL["mag_i"])
T_mass_BL_1 = T_M_BL_1.cal_Mass(ML_select_BL_1)

ML_select_BE_1 = SPlot.MLRelationIband(BE["mag_g"],BE["mag_i"]).Into13_MassRatio
M_BE_1 = SPlot.MassCalculation(BE["sph_mag"], BE["Dist"], M_sun,BE["mag_g"],BE["mag_i"])
mass_BE_1 = M_BE_1.cal_Mass(ML_select_BE_1)

T_M_BE_1 = SPlot.MassCalculation(BE["total_mag"], BE["Dist"], M_sun,BE["mag_g"],BE["mag_i"])
T_mass_BE_1 = T_M_BE_1.cal_Mass(ML_select_BE_1)
###

ML_select_BL_2 = SPlot.MLRelationIband(BL["mag_g"],BL["mag_i"]).Taylor11_MassRatio
M_BL_2 = SPlot.MassCalculation(BL["sph_mag"], BL["Dist"], M_sun,BL["mag_g"],BL["mag_i"])
mass_BL_2 = M_BL_2.cal_Mass(ML_select_BL_2)

T_M_BL_2 = SPlot.MassCalculation(BL["total_mag"], BL["Dist"], M_sun,BL["mag_g"],BL["mag_i"])
T_mass_BL_2 = T_M_BL_2.cal_Mass(ML_select_BL_2)

ML_select_BE_2 = SPlot.MLRelationIband(BE["mag_g"],BE["mag_i"]).Taylor11_MassRatio
M_BE_2 = SPlot.MassCalculation(BE["sph_mag"], BE["Dist"], M_sun,BE["mag_g"],BE["mag_i"])
mass_BE_2 = M_BE_2.cal_Mass(ML_select_BE_2)

T_M_BE_2 = SPlot.MassCalculation(BE["total_mag"], BE["Dist"], M_sun,BE["mag_g"],BE["mag_i"])
T_mass_BE_2 = T_M_BE_2.cal_Mass(ML_select_BE_2)
#%%
#Plot the sigma-mass diagram
mass_B_1 = [mass_BL_1,mass_BE_1]
mass_B_2 = [mass_BL_2,mass_BE_2]

T_mass_B_1 = [T_mass_BL_1,T_mass_BE_1]
T_mass_B_2 = [T_mass_BL_2,T_mass_BE_2]

vdis_B = [BL["vdis"],BE["vdis"]]
colour_B = [ "ro", "bo"]
legend_B = ["LTGs", "ETGs"]



SPlot.ShowcaseIndi.vdis_mass_plot(mass_1, A["vdis"], "ro", "All")
SPlot.ShowcaseIndi.vdis_mass_plot(T_mass_1, A["vdis"], "bo", "All")

plt.show()


plt.plot(mass_1, A["vdis"], "ro", label = "Sph mass")
plt.plot(T_mass_1, A["vdis"], "bo",label = "Total mass")
plt.xlabel("$M_*/M_{\odot}$",fontsize=16)
plt.ylabel("$\sigma$ (km/s)",fontsize=16)
            
plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
plt.show()


SPlot.ShowcaseIndi.vdis_mass_plot(mass_B_1, vdis_B, colour_B, legend_B)
plt.title("Velocity dispersion vs spheroid mass (Into2013)")

plt.show()

SPlot.ShowcaseIndi.vdis_mass_plot(mass_B_2, vdis_B, colour_B, legend_B)
plt.title("Velocity dispersion vs spheroid mass (Taylor2011)")

plt.show()


SPlot.ShowcaseIndi.vdis_mass_plot(T_mass_B_1, vdis_B, colour_B, legend_B)
plt.title("Velocity dispersion vs total stellar mass (Into2013)")
plt.show()
SPlot.ShowcaseIndi.vdis_mass_plot(T_mass_B_2, vdis_B, colour_B, legend_B)
plt.title("Velocity dispersion vs total stellar mass (Taylor2011)")

plt.show()
#%%

print(BE["total_mag"])
# saving in the new dictionary files
A["Into_mass"] = mass_1
A["Taylor_mass"] = mass_2

BL["Into_mass"] = mass_BL_1
BL["Taylor_mass"] = mass_BL_2

BE["Into_mass"] = mass_BE_1
BE["Into_mass"] = mass_BE_2

print(list(A.keys()))

print(BL.keys())

print(BE.keys())

with open("Gal_vdis_dict_equvi_mass", 'wb') as f:
    pickle.dump(A, f)
    
with open("Gal_vdis_dict_equvi_LTG_mass", 'wb') as f:
    pickle.dump(BL, f)
    
with open("Gal_vdis_dict_equvi_ETG_mass", 'wb') as f:
    pickle.dump(BE, f)