#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:02:25 2020

@author: dexter
"""

import SphRead as SRead
import SphSort as SSort
import SphPlot as SPlot


import matplotlib.pyplot as plt

name = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list.txt")["Gal_name"]
DD = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list.txt")["Dist"]
scale = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list.txt")["Scale"]

#total mag difference
mcpt_total_mag = SRead.grab_total_mag("Gal_bundle_equvi_cpt")

twocpt_total_mag = SRead.grab_total_mag("Gal_bundle_BD_equvi_cpt")


SPlot.ShowcaseCompare2.plot_compare_total_mag("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",sub =True, div = False, label = ["mulit \, cpt", "2cpt"])
plt.show()
SPlot.ShowcaseCompare2.plot_compare_total_mag("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",sub =False, div = True, label = ["mulit \, cpt", "2cpt"])
plt.show()

#Sersic index difference
mcpt_n = SRead.grab_parameter("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"],2)
twocpt_n = SRead.grab_parameter("Gal_bundle_BD_equvi_cpt", ["Bulge","CoreBulge"],2)

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],2, para_name="n", label = ["mulit \, cpt", "2cpt"])
plt.show()

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],2, sub =False, div = True, para_name="n", label = ["mulit \, cpt", "2cpt"])
plt.show()

#vRe plots
mcpt_Re = SRead.grab_parameter("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"],1)*scale
twocpt_Re = SRead.grab_parameter("Gal_bundle_BD_equvi_cpt", ["Bulge","CoreBulge"],1)*scale

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],1, para_name="R_e(kpc)", label = ["mulit \, cpt", "2cpt"])
plt.show()

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],1, sub =False, div = True, para_name="R_e(kpc)", label = ["mulit \, cpt", "2cpt"])
plt.show()


##mu plots

mcpt_mu = SRead.grab_parameter("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"],0)
twocpt_Rmu = SRead.grab_parameter("Gal_bundle_BD_equvi_cpt", ["Bulge","CoreBulge"],0)

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],0, para_name="mu", label = ["mulit \, cpt", "2cpt"])
plt.show()

SPlot.ShowcaseCompare2.plot_compare_feature_para("Gal_bundle_equvi_cpt","Gal_bundle_BD_equvi_cpt",["Bulge","CoreBulge"],0, sub =False, div = True, para_name="mu", label = ["mulit \, cpt", "2cpt"])
plt.show()




R_e = SRead.grab_parameter("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_kpc = R_e * scale

R_e_BD = SRead.grab_parameter("Gal_bundle_BD_equvi_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_kpc_BD = R_e_BD * scale


mag = SRead.grab_mag("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"])
mag_BD = SRead.grab_mag("Gal_bundle_BD_equvi_cpt", ["Bulge","CoreBulge"])

#mag2_total = SRead.grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Total_mag"],0)

A = SRead.read_list("Gal_vdis_dict_equvi")

mag_g,mag_i = A["mag_g"], A["mag_i"]

ML_select = SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
M = SPlot.MassCalculation(mag, DD, 4.53,mag_g,mag_i)
E = M.cal_Mass(ML_select)

M_BD = SPlot.MassCalculation(mag_BD, DD, 4.53,mag_g,mag_i)
E_BD = M_BD.cal_Mass(ML_select)


#SPlot.ShowcaseCompare2.plot_scat_arrow(E_BD,Re_kpc_BD,[],'bo',"B+D", E,Re_kpc,[],'ro',"Multi")
#plt.show()


SPlot.ShowcaseCompare2.plot_seperation_generic(Re_kpc_BD, Re_kpc, 1, para_name="", 
                                colour1="blue", colour2="orange", name=[], 
                                label=[1,2])
plt.show()


SPlot.ShowcaseCompare2.plot_seperation_generic(mcpt_n, twocpt_n, 1, para_name="", 
                                colour1="blue", colour2="orange", name=[], 
                                label=[1,2])
plt.show()
