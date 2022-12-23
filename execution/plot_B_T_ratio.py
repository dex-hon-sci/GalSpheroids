#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 21:10:56 2021

@author: dexter
"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

import matplotlib.gridspec as gridspec

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0

#input cpt bundle
bundle_name0= "/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt"
bundle = SRead.read_list(bundle_name0)

#input morph info
morph_file0 = "/home/dexter/result/stat/completeness/morph_list2.txt"

#morph_source = SRead.read_table(morph_file)

#input mag and dist info
D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_4.txt")


K_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_BinV_Kcorr_EXT.dat")

K_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_BinV_Kcorr_EXT.dat", 
    dtype='str')

morph_file_n = SRead.read_table(morph_file0,dtype='str')
# the corrected mag g and i, Kcorrection+EXTINCTIOn
mag_g, mag_i = K_table[:,10], K_table[:,9]
g_EXT, i_EXT = K_table[:,23], K_table[:,24]
g_kcorr, i_kcorr = K_table[:,25], K_table[:,26]

mag_g_corr, mag_i_corr = mag_g-g_kcorr, mag_i-i_kcorr

Dist, D_lerr, D_uerr = D0_all_table[:,29], D0_all_table[:,30], D0_all_table[:,31]

ES_index = [33,93]

morph_dict = SSort.morph_str_selection(np.arange(0,103,1),  morph_file_n[:,1])

#get core sersic
S = SSort.seperator_label_generic(bundle, ["Bulge","CoreBulge"])

sph_corebulge_mag = list(SRead.grab_mag(S[1]["CoreBulge"], ["CoreBulge"]))
Re_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 1))
Sersic_n_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 2)) 
mu_e_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 0))
total_mag_corebulge = list(SRead.grab_total_mag(S[1]["CoreBulge"]))

mag_g_EXT_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], g_EXT))
mag_i_EXT_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], i_EXT))
mag_g_kcorr_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], g_kcorr))
mag_i_kcorr_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], i_kcorr))

mag_g_corr_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], mag_g_corr))
mag_i_corr_corebulge = np.array(SSort.cherry_pick(S[1]["CoreBulge_index"], mag_i_corr))

Dist_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], Dist)

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale_corebulge = np.array(Dist_corebulge)* ars
Re_kpc_corebulge = Re_corebulge* scale_corebulge

sph_corebulge_mag = np.array([sph_corebulge_mag])
total_mag_corebulge = np.array([total_mag_corebulge])

sph_corebulge_mag = sph_corebulge_mag -mag_i_kcorr_corebulge -mag_i_EXT_corebulge
total_mag_corebulge = total_mag_corebulge -mag_i_kcorr_corebulge -mag_i_EXT_corebulge

#Read Méndez-Abreu 
Range_data = SRead.read_table("/home/dexter/result/stat/completeness/MendezAbreu_2017_morph_range.dat")

mean_BT = Range_data[:,4]
mean_mass = Range_data[:,1]

high_BT = Range_data[:,5]
high_mass = Range_data[:,2]

low_BT = Range_data[:,6]
low_mass = Range_data[:,3]

def plot_B_T_ratio(morph_file,bundle_name,num_E):
    
    morph_source_t = SRead.read_table(morph_file,dtype="str")
    morph_name, morph = morph_source_t[:,0], morph_source_t[:,1]

    fig = plt.figure(figsize=(6.4,4.8))

    gs = gridspec.GridSpec(1,1) 
    ax = plt.subplot(gs[0])
    S = SPlot.ShowcaseIndi.cpt_to_total_by_type_plot(bundle_name, morph_name, 
                                                 morph, num_E=num_E, AX =ax)
    plt.tight_layout()

    return S

def plot_B_T_mass(bundle_name): # A function per request by the referee
    
    # input: Bundle, what component, morphology type
    name = SRead.grab_name(bundle_name)
    cpt_mag = SRead.grab_mag(bundle_name, ["Bulge","CoreBulge"])
    total_mag = SRead.grab_total_mag(bundle_name)
    #K correction and extinction for sph mag
    cpt_mag = cpt_mag - i_EXT - i_kcorr
    total_mag = total_mag - i_EXT - i_kcorr
    
    # calculate flux ratio
    mag_ratio = 10**((cpt_mag-total_mag) / -2.5)
    #mag_ratio = 10**(cpt_mag)/10**total_mag
    mag_ratio = np.log10(mag_ratio)
    
    # B/T ratio for corebulge
    mag_ratio_corebulge = 10**((sph_corebulge_mag-total_mag_corebulge) / -2.5)
    mag_ratio_corebulge = np.log10(mag_ratio_corebulge)

    #calculate masses
    MLR_RC15 = SPlot.MLRelationIband(mag_g_corr,mag_i_corr).Roediger15BC03_MassRatio
    M = SPlot.MassCalculation(cpt_mag, Dist, 4.53,mag_g_corr,mag_i_corr)
    
    MLR_RC15_corebulge = SPlot.MLRelationIband(mag_g_corr_corebulge,
                                               mag_i_corr_corebulge).\
        Roediger15BC03_MassRatio
    M_corebulge = SPlot.MassCalculation(sph_corebulge_mag, Dist_corebulge, 4.53,
                                        mag_g_corr_corebulge,mag_i_corr_corebulge)

    mass = M.cal_Mass(MLR_RC15)
    mass_corebulge = M_corebulge.cal_Mass(MLR_RC15_corebulge)

    # separate masses
    name_E = np.array(SSort.cherry_pick(morph_dict['E'], name))
    mass_E = np.array(SSort.cherry_pick(morph_dict['E'], mass))
    mag_ratio_E = np.array(SSort.cherry_pick(morph_dict['E'], mag_ratio))
    
    name_S0 = np.array(SSort.cherry_pick(morph_dict['S0'], name))
    mass_S0 = np.array(SSort.cherry_pick(morph_dict['S0'], mass))
    mag_ratio_S0 = np.array(SSort.cherry_pick(morph_dict['S0'], mag_ratio))
    
    name_S = np.array(SSort.cherry_pick(morph_dict['S'], name))
    mass_S = np.array(SSort.cherry_pick(morph_dict['S'], mass))
    mag_ratio_S = np.array(SSort.cherry_pick(morph_dict['S'], mag_ratio))

    mass_ES= np.array([mass[ES_index[0]],mass[ES_index[1]]])
    mag_ratio_ES= np.array([mag_ratio[ES_index[0]],mag_ratio[ES_index[1]]])

    # plotting
    fig = plt.figure(figsize=(6.4,4.8))
    gs = gridspec.GridSpec(1,1) 
    ax = plt.subplot(gs[0])
        
    ax.scatter(mass_corebulge, 10**mag_ratio_corebulge,
                 facecolors='none', edgecolors='r', s = 500, 
                 label=r'$\rm core-S\'ersic$')
    # make the scatter plot for each morphology
    ax.plot(mass_E,10**mag_ratio_E,'s',color='#1b872a',ms=10,label=r"$\rm E$")
    ax.plot(mass_ES,10**mag_ratio_ES,'s',color="#ccab05",ms=10,label=r"$\rm ES$")
    
    ax.plot(mass_S0,10**mag_ratio_S0,'o',color="#d20d0d",ms=10,label=r"$\rm S0$")
    ax.plot(mass_S,10**mag_ratio_S,'o',color="k",ms=10,label=r"$\rm S$")
    
    #plot the range of Méndez-Abreu 2017
    #ax.errorbar(10**mean_mass, mean_BT, xerr = [10**mean_mass-10**low_mass,10**high_mass-10**mean_mass], 
    #              yerr = [mean_BT-low_BT,high_BT-mean_BT],ls='none',linewidth=4,
    #              color = '#a5200b',
    #              ecolor='#a5200b',capsize=0,
    #              alpha=1.0, marker='o') 

    # Find outliers lower than our mass limit    
    low = SSort.selection_generic(mag_ratio, mass, np.repeat(6.7e10, len(mass)))
    print(len(low['index']))
    
    # find outliers in S 2e11
    out1 = SSort.selection_generic(mag_ratio_S, mass_S, np.repeat(6e11, len(mass_S)),direction='high')
    print(out1)
    print(out1['index'])
    
    out2 = SSort.selection_generic(mag_ratio_E, mass_E, np.repeat(np.log10(0.8), len(mass_E)),direction='low',axis='x')
    print(out2)
    print(name_E[out2['index']])
    
    out3 = SSort.selection_generic(mag_ratio_S0, mass_S0, np.repeat(6e11, len(mass_S0)),direction='high')
    
    #show all name
    #SPlot.ShowcaseIndi.show_name(mass, 10**mag_ratio, name)
    
    def show_name_BT():
        SPlot.ShowcaseIndi.show_name(mass_E[out2['index']], 10**mag_ratio_E[out2['index']], name_E[out2['index']])
    
        # print each outlier's name
        for i in range(len(out1['index'])):
            print(name_S[out1['index'][i]])
            SPlot.ShowcaseIndi.show_name(mass_S[out1['index']], 10**mag_ratio_S[out1['index']], name_S[out1['index']])
        for j in range(len(out3['index'])):
            SPlot.ShowcaseIndi.show_name(mass_S0[out3['index']], 10**mag_ratio_S0[out3['index']], name_S0[out3['index']])
    
    show_name_BT()

    ax.set_ylabel(r"$B/T$",fontsize=18)
    ax.set_xlabel(r"$M_\mathrm{*,Sph} / \rm M_{\odot} (\rm RC15)$",fontsize=18)
    
    ax.set_xscale("log")
    
    plt.ylim(-0.1,1.1)
    #plt.xlim(2e10,2e12)
    plt.legend(loc="upper left", fontsize=12)
    plt.tight_layout()

    return fig

    #A.plot(mass_E, Re_kpc_E ,marker='s',color='#1b872a',label=r'$\rm E$', 
    #          ms =10, alpha=1.0,linestyle="None")   
    ## Plot ES galaxies
    #A.plot(mass_ES, Re_kpc_ES,marker='s',color='#ccab05',label=r'$\rm ES$', 
    #          ms =10, alpha=1.0,linestyle="None")

plot_B_T_ratio(morph_file0,bundle_name0,num_E=11.0)    #

bundle_2 = "/home/dexter/SphProject/F_Gal_bundle_equvi_V"

B = SRead.grab_parameter(bundle_2, ["BrokenExp"], 2)
C = SRead.grab_parameter(bundle_2, ["BrokenExp"], 3)

D = SRead.grab_parameter(bundle_2, ["CoreSersic"], 0)
D2 = SRead.grab_parameter(bundle_name0, ["CoreBulge"], 0)

#print(SRead.read_list(bundle_2))
#print(len(B))
#print(D)
#print(D2)

print('===============')

bundle_name3 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt"
morph_file3 = "/home/dexter/result/stat/completeness/morph_list_Bin3.txt"


plot_B_T_ratio(morph_file3,bundle_name3, num_E=3.0)
plt.text(0,-1.8,r"$\rm Bin~3$",fontsize = 26)

print('===============')

bundle_name2 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt"
morph_file2 = "/home/dexter/result/stat/completeness/morph_list_Bin2.txt"

plot_B_T_ratio(morph_file2,bundle_name2, num_E=3.0)
plt.text(0,-1.8,r"$\rm Bin~2$",fontsize = 26)

print('===============')


bundle_name1 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt"
morph_file1 = "/home/dexter/result/stat/completeness/morph_list_Bin1.txt"

plot_B_T_ratio(morph_file1,bundle_name1,num_E=5.0)
plt.text(0,-1.8,r"$\rm Bin~1$",fontsize = 26)

print('===============')

plot_B_T_mass(bundle_name0)

#print(B-C)

#+ 14 (typeII)
#- 17 (typeIII)

#fig = plt.figure()
#
#gs2 = gridspec.GridSpec(1,1) 
#
#ax2 = plt.subplot(gs2[0])
##
#S2 = SPlot.ShowcaseIndi.cpt_to_total_by_type_plot(bundle_name, morph_name, 
      #                                           morph, cpt=["nucDisk","Disk"],
      #                                           AX =ax2)
