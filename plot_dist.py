#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 01:30:40 2020

@author: dexter
"""

from astropy.cosmology import FlatLambdaCDM

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)


CCC2 = np.genfromtxt("/home/dexter/result/Ctable_Into_bin2_final.txt" , dtype='float')
CCC3 = np.genfromtxt("/home/dexter/result/Ctable_Into_bin3_final.txt" , dtype='float')
CCC4 = np.genfromtxt("/home/dexter/result/Ctable_Into_bin4_final.txt" , dtype='float')

zdist2 = CCC2[:,1]
zdist3 = CCC3[:,1]
zdist4 = CCC4[:,1]

dc2 = cosmo.comoving_distance(zdist2).value
sc2 = cosmo.kpc_proper_per_arcmin(zdist2).value / 60

dc3 = cosmo.comoving_distance(zdist3).value
sc3 = cosmo.kpc_proper_per_arcmin(zdist3).value / 60

dc4 = cosmo.comoving_distance(zdist4).value
sc4 = cosmo.kpc_proper_per_arcmin(zdist4).value / 60

D2= SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt",
                    "/home/dexter/result/distance/dist_list_bin2.txt")
D3= SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt",
                    "/home/dexter/result/distance/dist_list_bin3.txt")
D4= SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt",
                    "/home/dexter/result/distance/dist_list_bin4.txt")

name2, name3, name4 = D2["Gal_name"], D3["Gal_name"], D4["Gal_name"]
DD2, DD3, DD4 = D2["Dist"], D3["Dist"], D4["Dist"]
DD2_err, DD3_err, DD4_err = D2["Dist_err"], D3["Dist_err"], D4["Dist_err"]
scale2,scale3,scale4 = D2["Scale"], D3["Scale"],D4["Scale"]

#################

SPlot.ShowcaseCompare2.plot_distdist(DD2, scale2, dc2, sc2, name2, 100.0, DD_err = DD2_err)
SPlot.ShowcaseCompare2.plot_distdist(DD3, scale3, dc3, sc3, name3, 65.0, DD_err = DD3_err)
SPlot.ShowcaseCompare2.plot_distdist(DD4, scale4, dc4, sc4, name4, 43.6, DD_err = DD4_err)
###############
total_mag2 = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin2_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin3_cpt")
total_mag4 = SRead.grab_total_mag("/home/dexter/result/Gal_bundle_equvi_bin4_cpt")


vdis_file = SRead.read_table("/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt")
vdis_file2 = SRead.read_table("/home/dexter/result/velocity_disp/vel_disp_list_all_mag.txt",dtype='str')


mag_g2 = np.array(SRead.extract_match(name2, vdis_file2[:,0] ,vdis_file[:,10]))
mag_i2 = np.array(SRead.extract_match(name2, vdis_file2[:,0] ,vdis_file[:,9]))

mag_g3 = np.array(SRead.extract_match(name3, vdis_file2[:,0] ,vdis_file[:,10]))
mag_i3 = np.array(SRead.extract_match(name3, vdis_file2[:,0] ,vdis_file[:,9]))

mag_g4 = np.array(SRead.extract_match(name4, vdis_file2[:,0] ,vdis_file[:,10]))
mag_i4 = np.array(SRead.extract_match(name4, vdis_file2[:,0] ,vdis_file[:,9]))


ML_select2 = SPlot.MLRelationIband(mag_g2,mag_i2).Into13_MassRatio
M2 = SPlot.MassCalculation(total_mag2, DD2, 4.53,mag_g2,mag_i2)
E2 = M2.cal_Mass(ML_select2)

ML_select3 = SPlot.MLRelationIband(mag_g3,mag_i3).Into13_MassRatio
M3 = SPlot.MassCalculation(total_mag3, DD3, 4.53,mag_g3,mag_i3)
E3 = M3.cal_Mass(ML_select3)

ML_select4 = SPlot.MLRelationIband(mag_g4,mag_i4).Into13_MassRatio
M4 = SPlot.MassCalculation(total_mag4, DD4, 4.53,mag_g4,mag_i4)
E4 = M4.cal_Mass(ML_select4)


#############


M2_o = SPlot.MassCalculation(mag_i2, DD2, 4.53,mag_g2,mag_i2)
E2_o = M2_o.cal_Mass(ML_select2)

M3_o = SPlot.MassCalculation(mag_i3, DD3, 4.53,mag_g3,mag_i3)
E3_o = M3_o.cal_Mass(ML_select3)

M4_o = SPlot.MassCalculation(mag_i4, DD4, 4.53,mag_g4,mag_i4)
E4_o = M4_o.cal_Mass(ML_select4)


##############
SPlot.ShowcaseCompare2.plot_seperation_generic(E2, E2_o, 4e11, para_name="$M_*/M_\odot$", 
                                colour1="blue", colour2="orange", name=name2, 
                                label=['new','old'])  

SPlot.ShowcaseCompare2.plot_seperation_generic(E3, E3_o, 2e11, para_name="$M_*/M_\odot$", 
                                colour1="blue", colour2="orange", name=name3, 
                                label=['new','old'])  

SPlot.ShowcaseCompare2.plot_seperation_generic(E4, E4_o, 1e11, para_name="$M_*/M_\odot$", 
                                colour1="blue", colour2="orange", name=name4, 
                                label=['new','old'])  


#############
fig, ax = plt.subplots()

plt.plot(DD2,E2_o,'r^')
plt.plot(DD3,E3_o,'b^')
plt.plot(DD4,E4_o,'k^')
#plt.plot(DD2,E2,'ro')
#plt.plot(DD3,E3,'bo')
#plt.plot(DD4,E4,'ko')
plt.hlines(4e11,0,100,linestyle="solid")
plt.hlines(2e11,0,65,linestyle="solid")
plt.hlines(1e11,0,43.6,linestyle="solid")
plt.hlines(0.5e11,0,43.6,linestyle="solid")

plt.vlines(120,0,1e12,linestyle="solid")
plt.vlines(100,0,1e12,linestyle="solid")
plt.vlines(65,0,4e11,linestyle="solid")
plt.vlines(43.6,0,2e11,linestyle="solid")


plt.show()


