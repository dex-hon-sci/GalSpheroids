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
def plot_dist_difference():
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


############# new distance, new i-band magnitude

ML_select2 = SPlot.MLRelationIband(mag_g2,mag_i2).Into13_MassRatio
M2 = SPlot.MassCalculation(total_mag2, DD2, 4.53,mag_g2,mag_i2)
E2 = M2.cal_Mass(ML_select2)

ML_select3 = SPlot.MLRelationIband(mag_g3,mag_i3).Into13_MassRatio
M3 = SPlot.MassCalculation(total_mag3, DD3, 4.53,mag_g3,mag_i3)
E3 = M3.cal_Mass(ML_select3)

ML_select4 = SPlot.MLRelationIband(mag_g4,mag_i4).Into13_MassRatio
M4 = SPlot.MassCalculation(total_mag4, DD4, 4.53,mag_g4,mag_i4)
E4 = M4.cal_Mass(ML_select4)

############# old distance, old i-band magnitude

M2_o1 = SPlot.MassCalculation(mag_i2, DD2, 4.53,mag_g2,mag_i2)
E2_o1 = M2_o1.cal_Mass(ML_select2)

M3_o1 = SPlot.MassCalculation(mag_i3, DD3, 4.53,mag_g3,mag_i3)
E3_o1 = M3_o1.cal_Mass(ML_select3)

M4_o1 = SPlot.MassCalculation(mag_i4, DD4, 4.53,mag_g4,mag_i4)
E4_o1 = M4_o1.cal_Mass(ML_select4)


############# old distance, old i-band magnitude

M2_o = SPlot.MassCalculation(mag_i2, dc2, 4.53,mag_g2,mag_i2)
E2_o = M2_o.cal_Mass(ML_select2)

M3_o = SPlot.MassCalculation(mag_i3, dc3, 4.53,mag_g3,mag_i3)
E3_o = M3_o.cal_Mass(ML_select3)

M4_o = SPlot.MassCalculation(mag_i4, dc4, 4.53,mag_g4,mag_i4)
E4_o = M4_o.cal_Mass(ML_select4)


##############

def plot_mass_difference():

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
from astropy.io import fits
#############
master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass.txt"

#hdul = fits.open(master_file)
#data = np.array(hdul[1].data)

master = SRead.read_table(master_file)

dist_all = master[:,0]
Into_mass_all = master[:,1]

fig, ax = plt.subplots()

#plt.plot(dc2,E2_o,'r^',label='pre-process-bin1', ms=12,alpha=0.5)
##plt.plot(dc3,E3_o,'b^',label='pre-process-bin2', ms=12,alpha=0.5)
#plt.plot(dc4,E4_o,'k^',label='pre-process-bin3',ms=12,alpha=0.5)

plt.plot(DD2,E2_o1,'ro',label='new Dist, SDSS mag-bin1', ms=12,alpha=0.2)
plt.plot(DD3,E3_o1,'bo',label='new Dist, SDSS mag-bin2', ms=12,alpha=0.2)
plt.plot(DD4,E4_o1,'ko',label='new Dist, SDSS mag-bin3', ms=12,alpha=0.2)


plt.plot(DD2,E2,'r^',label='new Dist, profiler mag-bin1', ms=12,alpha=0.5)
plt.plot(DD3,E3,'b^',label='new Dist, profiler mag-bin2', ms=12,alpha=0.5)
plt.plot(DD4,E4,'k^',label='new Dist, profiler mag-bin3', ms=12,alpha=0.5)

from scipy import stats
import statistics
print(np.median(mag_i2))
print(np.median(total_mag2))
print(np.median(dc2))
print(np.median(DD2))
print(np.median(E2_o))
print(np.median(E2_o1))
print(np.median(E2))
print(np.median(ML_select2))

print(np.median(E2-E2_o1))

print(np.median(total_mag2-mag_i2))
print(np.median(DD2-dc2))

print(DD2-dc2)

plt.hlines(np.average(E2_o1),0,120,linestyle="dashed",color="red")
plt.hlines(np.average(E3_o1),0,120,linestyle="dashed",color="blue")
plt.hlines(np.average(E4_o1),0,120,linestyle="dashed",color="black")

plt.plot(dist_all,Into_mass_all,'x',color='grey' , alpha=0.9, label='All')


plt.hlines(4e11,0,100,linestyle="solid")
plt.hlines(2e11,0,65,linestyle="solid")
plt.hlines(1e11,0,43.6,linestyle="solid")
plt.hlines(0.5e11,0,43.6,linestyle="solid")

plt.vlines(120,0,1e15,linestyle="solid")
plt.vlines(100,0,1e15,linestyle="solid")
plt.vlines(65,0,4e11,linestyle="solid")
plt.vlines(43.6,0,2e11,linestyle="solid")
plt.ylim(top =1e12 , bottom = 1e8)
plt.xlim(left=0, right = 130)



x_edge1,y_edge1= [0,0,100,100], [4e11,1e13,1e13,4e11]
x_edge2,y_edge2= [0,0,65,65], [2e11,4e11,4e11,2e11]
x_edge3,y_edge3= [0,0,43.6,43.6], [1e11,2e11,2e11,1e11]

plt.fill(x_edge1,y_edge1, alpha=0.3, color='green',
                  label='selection')
plt.fill(x_edge2,y_edge2, alpha=0.2, color='green',
                  label='selection')
plt.fill(x_edge3,y_edge3, alpha=0.1, color='green',
                  label='selection')

plt.ylabel("$M_*/M_\odot(Into Mass)$",fontsize=16)
plt.xlabel("$Distance(H_0 = 68, \Omega_m =0.3)/Mpc$",fontsize=16)
plt.yscale( 'log' )

plt.legend()
plt.show()

