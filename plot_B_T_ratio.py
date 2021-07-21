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



bundle_name0= "/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt"

morph_file0 = "/home/dexter/result/stat/completeness/morph_list.txt"

#morph_source = SRead.read_table(morph_file)


def plot_B_T_ratio(morph_file,bundle_name):
    
    morph_source_t = SRead.read_table(morph_file,dtype="str")
    morph_name, morph = morph_source_t[:,0], morph_source_t[:,1]

    fig = plt.figure(figsize=(6.4,4.8))

    gs = gridspec.GridSpec(1,1) 
    ax = plt.subplot(gs[0])
    S = SPlot.ShowcaseIndi.cpt_to_total_by_type_plot(bundle_name, morph_name, 
                                                 morph, AX =ax)
    plt.tight_layout()

plot_B_T_ratio(morph_file0,bundle_name0)    #

bundle_2 = "/home/dexter/SphProject/F_Gal_bundle_equvi_V"

B = SRead.grab_parameter(bundle_2, ["BrokenExp"], 2)
C = SRead.grab_parameter(bundle_2, ["BrokenExp"], 3)

D = SRead.grab_parameter(bundle_2, ["CoreSersic"], 0)
D2 = SRead.grab_parameter(bundle_name0, ["CoreBulge"], 0)

#print(SRead.read_list(bundle_2))
#print(len(B))
print(D)
print(D2)


bundle_name3 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt"
morph_file3 = "/home/dexter/result/stat/completeness/morph_list_Bin3.txt"


#plot_B_T_ratio(morph_file3,bundle_name3)
#plt.text(0,-1.8,r"$\rm Bin~3$",fontsize = 26)



bundle_name2 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt"
morph_file2 = "/home/dexter/result/stat/completeness/morph_list_Bin2.txt"

#plot_B_T_ratio(morph_file2,bundle_name2)
#plt.text(0,-1.8,r"$\rm Bin~2$",fontsize = 26)


bundle_name1 = "/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt"
morph_file1 = "/home/dexter/result/stat/completeness/morph_list_Bin1.txt"

#plot_B_T_ratio(morph_file1,bundle_name1)
#plt.text(0,-1.8,r"$\rm Bin~1$",fontsize = 26)


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
