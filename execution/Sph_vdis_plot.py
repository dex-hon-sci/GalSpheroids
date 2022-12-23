#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 18:16:32 2020

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import pickle

import matplotlib.pyplot as plt
import numpy as np

M_sun = 4.53


# creating the velocity dispersion dictionary and calculating the stellar mass of the spheroid 
## If you need further investigation on the distance list and the velocity dispersion sources
# Just let me know, I can send the tables to you
#%% 
A = SSort.vdis_match("./galaxy_bundle/Gal_bundle_equvi_cpt", 
                 "../result/velocity_disp/vel_disp_list_all_mag.txt", 
                 ["../result/distance/dir_dist_list.txt","../result/distance/dist_list.txt"],
                 "./galaxy_bundle/Gal_vdis_dict_equvi")

B = SSort.LTG_ETG_seperator("./galaxy_bundle/Gal_bundle_equvi_cpt",
                            "./galaxy_bundle/Gal_bundle_equvi_cpt_LTG",
                            "./galaxy_bundle/Gal_bundle_equvi_cpt_ETG") 

# removing and adding the exception, NGC4334 is a LTG
#NGC3126 (Sd) ETG->LTG
#NGC4015 E
#NGC4224 SA 
#NGC4235 SA
#NGC5318 S0
#NGC5363 S0-a 
#UGC8736 Sc

###

#IC947 L->E
#NGC2824 L->E
#NGC2962 SAB0 L->E
#NGC2968 Sa
#NGC3106 L->E
#NGC3189 SAB(rs)0/a
#NGC3658 SA0
#NGC4429 SA0^+(r) 
#NGC4461 SB0^+(s)
#NGC4477 SB0(s):?
#NGC4570 S0edge-on
#NGC4578 SA0^0(r)
#NGC4596 SB0^+(r
#NGC4608 SB0^0(r)
#NGC4643 SB(rs)0/a 
#NGC3658 SA0
#NGC4429 SA0^+(r) 
#NGC4461 SB0^+(s)
#NGC4477 SB0(s):?
#NGC4570 S0edge-on
#NGC4578 SA0^0(r)
#NGC4596 SB0^+(r
#NGC4608 SB0^0(r)
#NGC4643 SB(rs)0/a 

#NGC5174 Sd
#NGC5473 SAB0^-(s)
#UGC8564 S0

SSort.remove_item("./galaxy_bundle/Gal_bundle_equvi_cpt_ETG",["NGC3126"
                                                              ,"NGC4224",
                                                              "NGC4235",
                                                              "UGC8736"])
SSort.add_item("./galaxy_bundle/Gal_bundle_equvi_cpt",
               "./galaxy_bundle/Gal_bundle_equvi_cpt_LTG", ["NGC3126",
                                                            "NGC4224",
                                                            "NGC4235",
                                                            "UGC8736"])

SSort.remove_item("./galaxy_bundle/Gal_bundle_equvi_cpt_LTG",["IC00947",
                                                              "NGC2824",
                                                              "NGC2962",
                                                              "NGC3106", 
                                                              "NGC3189",  
                                                              "NGC3658", 
                                                              "NGC4429", 
                                                              "NGC4461","NGC4477", 
                                                              "NGC4570", "NGC4578", 
                                                              "NGC4596","NGC4608",
                                                              "NGC4643", "NGC5473", 
                                                              "UGC8564"])
SSort.add_item("./galaxy_bundle/Gal_bundle_equvi_cpt", 
               "./galaxy_bundle/Gal_bundle_equvi_cpt_ETG", ["IC00947","NGC2824",
                                                            "NGC2962","NGC3106", 
                                                            "NGC3189",  "NGC3658",
                                                            "NGC4429", "NGC4461",
                                                            "NGC4477", "NGC4570", 
                                                            "NGC4578", "NGC4596",
                                                            "NGC4608","NGC4643", 
                                                            "NGC5473", "UGC8564"])


BL = SSort.vdis_match("./galaxy_bundle/Gal_bundle_equvi_cpt_LTG", 
                 "../result/velocity_disp/vel_disp_list_all_mag.txt", 
                 ["../result/distance/dir_dist_list.txt","../result/distance/dist_list.txt"],
                 "./galaxy_bundle/Gal_vdis_dict_equvi_LTG")

BE = SSort.vdis_match("./galaxy_bundle/Gal_bundle_equvi_cpt_ETG", 
                 "../result/velocity_disp/vel_disp_list_all_mag.txt", 
                 ["../result/distance/dir_dist_list.txt","../result/distance/dist_list.txt"],
                 "./galaxy_bundle/Gal_vdis_dict_equvi_ETG")


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
SPlot.ShowcaseIndi.show_name(mass_BL_1, BL["vdis"],BL["name"],size=9)
SPlot.ShowcaseIndi.show_name(mass_BE_1, BE["vdis"],BE["name"],size=9)

plt.title("Velocity dispersion vs spheroid mass (Into2013)")

plt.show()

SPlot.ShowcaseIndi.vdis_mass_plot(mass_B_2, vdis_B, colour_B, legend_B)
SPlot.ShowcaseIndi.show_name(mass_BL_2, BL["vdis"],BL["name"],size=9)
SPlot.ShowcaseIndi.show_name(mass_BE_2, BE["vdis"],BE["name"],size=9)

plt.title("Velocity dispersion vs spheroid mass (Taylor2011)")

plt.show()


SPlot.ShowcaseIndi.vdis_mass_plot(T_mass_B_1, vdis_B, colour_B, legend_B)
SPlot.ShowcaseIndi.show_name(T_mass_BL_1, BL["vdis"],BL["name"],size=9)
SPlot.ShowcaseIndi.show_name(T_mass_BE_1, BE["vdis"],BE["name"],size=9)
plt.title("Velocity dispersion vs total stellar mass (Into2013)")
plt.show()

SPlot.ShowcaseIndi.vdis_mass_plot(T_mass_B_2, vdis_B, colour_B, legend_B)
SPlot.ShowcaseIndi.show_name(T_mass_BL_2, BL["vdis"],BL["name"],size=9)
SPlot.ShowcaseIndi.show_name(T_mass_BE_2, BE["vdis"],BE["name"],size=9)
plt.title("Velocity dispersion vs total stellar mass (Taylor2011)")

plt.show()
#%%

# saving in the new dictionary files
A["Into_mass"] = mass_1
A["Taylor_mass"] = mass_2

BL["Into_mass"] = mass_BL_1
BL["Taylor_mass"] = mass_BL_2

BE["Into_mass"] = mass_BE_1
BE["Taylor_mass"] = mass_BE_2

print(list(A.keys()))

print(BL.keys())

print(BE.keys())

print((BL["name"]))
print((BE["name"]))

with open("Gal_vdis_dict_equvi_mass", 'wb') as f:
    pickle.dump(A, f)
    
with open("Gal_vdis_dict_equvi_LTG_mass", 'wb') as f:
    pickle.dump(BL, f)
    
with open("Gal_vdis_dict_equvi_ETG_mass", 'wb') as f:
    pickle.dump(BE, f)
    
    
SRead.convert_dict_ascii("Gal_vdis_dict_equvi_mass","Gal_vdis_dict_equvi_mass.dat")
SRead.convert_dict_ascii("Gal_vdis_dict_equvi_LTG_mass","Gal_vdis_dict_equvi_LTG_mass.dat")
SRead.convert_dict_ascii("Gal_vdis_dict_equvi_ETG_mass","Gal_vdis_dict_equvi_ETG_mass.dat")
