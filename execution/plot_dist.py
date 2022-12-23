#!/usr/bin/envhon3
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





############### reading old data###############

#plot_dist_difference()

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

#######################################################
################reading the new set of data############
#######################################################
#d##the cosmicflow distance
d1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/cosmicflow_dist_bag3_Bin1.dat")
d2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/cosmicflow_dist_bag3_Bin2.dat")
d3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/cosmicflow_dist_bag3_Bin3.dat")

d1 = d1_table[:,7] * (75.0/75)
d2 = d2_table[:,7] * (75.0/75)
d3 = d3_table[:,7] * (75.0/75)

d1_uerr = (d1_table[:,9] - d1_table[:,7])*(75.0/75)
d2_uerr = (d2_table[:,9] - d2_table[:,7])*(75.0/75)
d3_uerr = (d3_table[:,9] - d3_table[:,7])*(75.0/75)

d1_lerr = (d1_table[:,7]- d1_table[:,5]) * (75.0/75)
d2_lerr = (d2_table[:,7]- d2_table[:,5]) * (75.0/75)
d3_lerr = (d3_table[:,7]- d3_table[:,5]) * (75.0/75)


print(np.size(d1),np.size(d1_uerr),np.size(d1_lerr))
print(np.size(d2),np.size(d2_uerr),np.size(d2_lerr))
print(np.size(d3),np.size(d3_uerr),np.size(d3_lerr))

#DD## mould and redshift independent measurement
DDD1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Mould_dist_bag3_Bin1.dat")
DDD2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Mould_dist_bag3_Bin2.dat")
DDD3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/Mould_dist_bag3_Bin3.dat")

DDD1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1.dat",
    dtype='str')
DDD2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2.dat",
    dtype='str')
DDD3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3.dat",
    dtype='str')

ddc1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1.dat")
ddc2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2.dat")
ddc3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3.dat")

d_name1 = DDD1_table_n[:,4]
d_name2 = DDD2_table_n[:,4]
d_name3 = DDD3_table_n[:,4]

DDD1 = DDD1_table[:,1] / 68.0
DDD2 = DDD2_table[:,1] / 68.0
DDD3 = DDD3_table[:,1] / 68.0

DDD1_err = DDD1_table[:,2] / 68.0
DDD2_err = DDD2_table[:,2] / 68.0
DDD3_err = DDD3_table[:,2] / 68.0

#dc## Willick2007 redshift , h68

ddc1 = ddc1_table[:,2]
ddc2 = ddc2_table[:,2]
ddc3 = ddc3_table[:,2]

corr_dist1, corr_E1 = ddc1_table[:,12],  ddc1_table[:,13]
corr_dist2, corr_E2 = ddc2_table[:,12],  ddc2_table[:,13]
corr_dist3, corr_E3 = ddc3_table[:,12],  ddc3_table[:,13]


d_spec_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",
    dtype="float")
d_spec_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",
    dtype="str")

d_spec = d_spec_table[:,3] 
d_spec_name = d_spec_table_n[:,0] 
d_spec_err = d_spec_table[:,4]


d_spec_index, d_spec_d, d_spec_d_err = [],[],[]

decision_table1 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin1.dat", dtype='str')
decision_table2 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin2.dat", dtype='str')
decision_table3 = SRead.read_table(
    "/home/dexter/result/stat/completeness/dist_judgement_Bin3.dat", dtype='str')

decision1 = decision_table1[:,1]
decision2 = decision_table2[:,1]
decision3 = decision_table3[:,1]



def plot_dist_difference2():
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD1,ddc1,d1,d_name1,110, 
                                                 DD_err=DDD1_err,dc_err=None, 
                                                 d_err= [d1_lerr,d1_uerr],                                                 
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = decision1)
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD2,ddc2,d2,d_name2,75, 
                                                 DD_err=DDD2_err,dc_err=None, 
                                                 d_err=[d2_lerr,d2_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name, 
                                                 d_spec_err = d_spec_err,
                                                 decision = decision2)
    
    SPlot.ShowcaseCompare2.plot_distdist_3points(DDD3,ddc3,d3,d_name3,45, 
                                                 DD_err=DDD3_err,dc_err=None, 
                                                 d_err=[d3_lerr,d3_uerr],
                                                 d_spec = d_spec, 
                                                 d_spec_name = d_spec_name , 
                                                 d_spec_err = d_spec_err,
                                                 decision = decision3)
    plt.show()
    
plot_dist_difference2()
########################################################
########################################################
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

###############################################

master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"
broad_cut_file="/home/dexter/result/stat/completeness/diagonal_selection_bag3_2.dat"

master = SRead.read_table(master_file)
broad_cut = SRead.read_table(broad_cut_file)
broad_cut_str = SRead.read_table(broad_cut_file,dtype="str")


dc = master[:,3]
initial_mass = master[:,4]

name_b = broad_cut_str[:,4]
dc_b = broad_cut[:,2]
initial_mass_b = broad_cut[:,3]


D = SRead.grab_dist(
    "/home/dexter/result/stat/completeness/dir_dist_list_broad3_2.txt",                    
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_2.dat",
              name_index = 4, vel_index =6, vel_err_index = 7, scale_index = 8,
              dir_name_index = 0, dir_dist_index = 3, dir_dist_err_index =4, 
              dir_method_index = 5, dir_method_flag_index = 7)

name1 = D["Gal_name"]
DD1 = D["Dist"]
DD1_err = D["Dist_err"]
scale1 = D["Scale"]


mag_g, mag_i = broad_cut[:,10],broad_cut[:,9]

ML_select= SPlot.MLRelationIband(mag_g,mag_i).Into13_MassRatio
M = SPlot.MassCalculation(mag_i, DD1, 4.53,mag_g,mag_i)
E = M.cal_Mass(ML_select)

#for i in range(np.size(name1)):
#    print(name1[i], DD1[i],E[i])

#############

#Bad_list = SRead.read_table("/home/dexter/result/stat/completeness/bad_list2.txt")

#Bad_list = Bad_list[]
#Bad_list_dist = Bad_list[]

###############################################

from scipy import stats
import statistics
import matplotlib.gridspec as gridspec

fig = plt.figure()

gs = gridspec.GridSpec(2,1) 

axs1 = plt.subplot(gs[1])
axs0 = plt.subplot(gs[0])

#plt.plot(dc2,E2_o,'r^',label='pre-process-bin1', ms=12,alpha=0.5)
#plt.plot(dc3,E3_o,'b^',label='pre-process-bin2', ms=12,alpha=0.5)
#plt.plot(dc4,E4_o,'k^',label='pre-process-bin3',ms=12,alpha=0.5)


SPlot.ShowcaseIndi.show_name(DD1,E,name_b)

####################################################################

axs1.hlines(2e12,0,110,linestyle="solid")
axs1.hlines(5e11,0,110,linestyle="solid")
axs1.hlines(2e11,0,75,linestyle="solid")
axs1.hlines(1e11,0,45,linestyle="solid")

axs1.vlines(110,5e11,2e12,linestyle="solid")
axs1.vlines(75,2e11,5e11,linestyle="solid")
axs1.vlines(45,1e11,2e11,linestyle="solid")

x_edge1,y_edge1= [0,0,110,110], [5e11,2e12,2e12,5e11]
x_edge2,y_edge2= [0,0,75,75], [2e11,5e11,5e11,2e11]
x_edge3,y_edge3= [0,0,45,45], [1e11,2e11,2e11,1e11]

##################################################################


#plt.plot(DD2,E2,'r^',label='new Dist, profiler mag-bin1', ms=12,alpha=0.5)
#plt.plot(DD3,E3,'b^',label='new Dist, profiler mag-bin2', ms=12,alpha=0.5)
#plt.plot(DD4,E4,'k^',label='new Dist, profiler mag-bin3', ms=12,alpha=0.5)

axs0.plot(dc,initial_mass,'x',color='grey', alpha=0.7, label='All')
axs0.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)
axs0.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)
#axs0.plot(dc,E,'x',color='green', alpha=0.4, label='New Distance')
axs0.plot(dc_b,initial_mass_b,'x',color='green', alpha=0.4, label='Broad Selection')


axs1.plot(DD1,E,'x',color='#a51a74', alpha=0.4, label='New Distance')

axs1.plot(corr_dist1,corr_E1,'ro',label='Bin1 samples', ms=12,alpha=0.5)
axs1.plot(corr_dist2,corr_E2,'bo',label='Bin2 samples', ms=12,alpha=0.5)
axs1.plot(corr_dist3,corr_E3,'ko',label='Bin3 samples', ms=12,alpha=0.5)

axs1.plot(dc, SPlot.SelectionCut(initial_mass, dc).
         parent_sample_cut(),linestyle="solid", color="blue",
         linewidth=2)
axs1.vlines(115,4.17e11,1e15,linestyle="solid", color="blue",
           linewidth=2)

axs1.fill(x_edge1,y_edge1, alpha=0.3, color='green')
axs1.fill(x_edge2,y_edge2, alpha=0.2, color='green')
axs1.fill(x_edge3,y_edge3, alpha=0.1, color='green')
axs0.set_ylim(top =5e12 , bottom = 1e9)
axs0.set_xlim(left=0, right = 120)


axs0.set_ylabel("$M_*/M_\odot(Into Mass)$",fontsize=24)
axs0.set_xlabel("$Distance(H_0 = 68, \Omega_m =0.3)/Mpc$",fontsize=24)
axs0.set_yscale( 'log' )

axs1.set_ylabel("$M_*/M_\odot(Into Mass)$",fontsize=24)
axs1.set_xlabel("$Distance(H_0 = 68, \Omega_m =0.3)/Mpc$",fontsize=24)
axs1.set_yscale( 'log' )

axs1.set_ylim(top =5e13 , bottom = 1e10)
axs1.set_xlim(left=0, right = 120)

axs0.legend()
axs1.legend()
plt.show()

#4.6e11  111.5

#2.1039e1	75.47

#1.35e11	46.2


#plt.hlines(np.average(E2_o1),0,120,linestyle="dashed",color="red")
#plt.hlines(np.average(E3_o1),0,120,linestyle="dashed",color="blue")
#plt.hlines(np.average(E4_o1),0,120,linestyle="dashed",color="black")

#plt.plot(dc,initial_mass,'x',color='grey', alpha=0.7, label='All')

####################################################################
#
#plt.hlines(1.12e12,0,108,linestyle="solid")
#plt.hlines(4.6e11,0,108,linestyle="solid")
#plt.hlines(2.104e11,0,74.5,linestyle="solid")
#plt.hlines(1.18e11,0,44.7,linestyle="solid")
#
#plt.vlines(108,4.6e11,1.12e12,linestyle="solid")
#plt.vlines(74.5,2.104e11,4.6e11,linestyle="solid")
#plt.vlines(44.7,1.18e11,2.104e11,linestyle="solid")
#
#x_edge1,y_edge1= [0,0,108,108], [4.6e11,1.12e12,1.12e12,4.6e11]
#x_edge2,y_edge2= [0,0,74.5,74.5], [2.104e11,4.6e11,4.6e11,2.104e11]
#x_edge3,y_edge3= [0,0,44.7,44.7], [1.18e11,2.104e11,2.104e11,1.18e11]

##################################################################



#Bag = SPlot.SelectionCut(Into_mass_all, dist_all).selection_subsample(master)
#print(Bag["bag"])
#print("len(Bag)", len(Bag["index"]))


#SRead.pickle_save(Bag, "diagonal_selection_bag")

#SRead.convert_list_ascii(Bag["bag"],"diagonal_selection_bag.dat")




fig, ax = plt.subplots()

plt.plot(E2,mag_g2-mag_i2,'ro')
plt.plot(E3,mag_g3-mag_i3,'bo')
plt.plot(E4,mag_g4-mag_i4,'ko')
plt.xscale( 'log' )

plt.show()
