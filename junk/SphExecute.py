#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:38:57 2020

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import matplotlib.pyplot as plt
import numpy as np
mass = np.linspace(2e8,0.5e13,2000)



#import astropy.cosmology.core as A
###########initial parameters###############
c, H0=    299792.458, 68.0        #speed of light in km/s, Hubble constant in(km/s)/Mpc
M_sun = 4.53


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=68, Om0=0.3) 

dc = cosmo.comoving_distance(0.025)
sc = cosmo.kpc_proper_per_arcmin(0.025) / 60




############################################

#fig, ax = plt.subplots()

#plt.show()
#SRead.lookup_bundle("Gal_bundle_equvi_bin2_cpt","NGC2918")

#print("lookup", SRead.lookup_bundle("Gal_bundle_equvi_bin2_cpt","NGC2918"))
#print("Para",SRead.grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1))
#print("mag",SRead.grab_mag("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"]))

DD = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list.txt")["Dist"]


#input_list, dist = "Gal_bundle_equvi_cpt",DD


SPlot.ShowcaseIndi.plot_hist_percentage("Gal_bundle_equvi_cpt",DD)
plt.show()

name2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Gal_name"]
DD2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Dist"]
scale2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Scale"]

CCC2 = np.genfromtxt("Ctable_Into_bin2_final.txt" , dtype='float')

R_e_2 = SRead.grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag2 = SRead.grab_mag("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"])


AbsMag_g_bin2 = CCC2[:,6]
AbsMag_i_bin2 = CCC2[:,8]


Re_kpc2 = R_e_2 * scale2

Mag2 = mag2 - 25 -5*np.log10(DD2)
Lum2 = 10**((Mag2-(4.53))/(-2.5))
##########################

name3 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Gal_name"]
DD3 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Dist"]
scale3 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Scale"]

CCC3 = np.genfromtxt("Ctable_Into_bin3_final.txt" , dtype='float')

R_e_3 = SRead.grab_parameter("Gal_bundle_equvi_bin3_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag3 = SRead.grab_mag("Gal_bundle_equvi_bin3_cpt", ["Bulge","CoreBulge"])


AbsMag_g_bin3 = CCC3[:,6]
AbsMag_i_bin3 = CCC3[:,8]

#Mass_bulge_Into2 = cal_Mass(mag2,DD2,ML_relation_Iband(AbsMag_g_bin2,AbsMag_i_bin2).Into13_MassRatio,M_sun)
Re_kpc3 = R_e_3 * scale3

Mag3 = mag3 - 25 -5*np.log10(DD3)
Lum3 = 10**((Mag2-(4.53))/(-2.5))

##########################

name4 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Gal_name"]
DD4 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Dist"]
scale4 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Scale"]

CCC4 = np.genfromtxt("Ctable_Into_bin4_final.txt" , dtype='float')

R_e_4 = SRead.grab_parameter("Gal_bundle_equvi_bin4_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag4 = SRead.grab_mag("Gal_bundle_equvi_bin4_cpt", ["Bulge","CoreBulge"])


AbsMag_g_bin4 = CCC4[:,6]
AbsMag_i_bin4 = CCC4[:,8]

#Mass_bulge_Into2 = cal_Mass(mag2,DD2,ML_relation_Iband(AbsMag_g_bin2,AbsMag_i_bin2).Into13_MassRatio,M_sun)
Re_kpc4 = R_e_4 * scale4

Mag4 = mag4 - 25 -5*np.log10(DD4)
Lum4 = 10**((Mag4-(4.53))/(-2.5))
################################
ML_select2 = SPlot.MLRelationIband(AbsMag_g_bin2,AbsMag_i_bin2).Into13_MassRatio
M2 = SPlot.MassCalculation(mag2, DD2, M_sun,AbsMag_g_bin2,AbsMag_i_bin2)
E2 = M2.cal_Mass(ML_select2)

ML_select3 = SPlot.MLRelationIband(AbsMag_g_bin3,AbsMag_i_bin3).Into13_MassRatio
M3 = SPlot.MassCalculation(mag3, DD3, M_sun,AbsMag_g_bin3,AbsMag_i_bin3)
E3 = M3.cal_Mass(ML_select3)

ML_select4 = SPlot.MLRelationIband(AbsMag_g_bin4,AbsMag_i_bin4).Into13_MassRatio
M4 = SPlot.MassCalculation(mag4, DD4, M_sun,AbsMag_g_bin4,AbsMag_i_bin4)
E4 = M4.cal_Mass(ML_select4)



SPlot.ShowcaseIndi.Mass_Re_plot(E2,Re_kpc2,name2,"ro","test",1.0)
SPlot.ShowcaseIndi.Mass_Re_plot(E3,Re_kpc3,name3,"bo","test",1.0)
SPlot.ShowcaseIndi.Mass_Re_plot(E4,Re_kpc4,name4,"ko","test",1.0)
plt.show()



SPlot.ShowcaseCompare2.plot_compare_rms('Gal_bundle_major_cpt','Gal_bundle_BD_major_cpt')
plt.show()


index2 = np.linspace(0,len(Mag2)-1,len(Mag2))
index3 = np.linspace(0+len(Mag2),len(Mag3)-1+len(Mag2),len(Mag3))
index4 = np.linspace(0+len(Mag2)+len(Mag3),len(Mag4)-1+len(Mag2)+len(Mag3),len(Mag4))

delta2 = Mag2-AbsMag_i_bin2
delta3 = Mag3-AbsMag_i_bin3
delta4 = Mag4-AbsMag_i_bin4

index = np.concatenate((np.concatenate((index2, index3)),index4))
delta = np.concatenate((np.concatenate((delta2, delta3)),delta4))

plt.plot(index2,delta2,'ro')
plt.plot(index3,delta3,'bo')
plt.plot(index4,delta4,'ko')

plt.text(50, np.average(delta)+0.005, r"$\langle \Delta Abs Mag_i \rangle$ = %s $\pm$ %s" %(round(np.average(delta),3),round(np.std(delta),2)),fontsize=14, color="#453816")

plt.hlines(np.average(delta),0,103, color='k', linestyle="dashed", linewidth=3)

plt.hlines(0,0,103, color='k', linestyle="dashed", linewidth=3)

plt.show()