#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 03:02:36 2020

@author: dexter
"""
#/home/dexter/result/stat/Into_selection_all
#/home/dexter/result/distance/
#vel_disp_list_all_mag.txt

# Name & RA (2000) & Dec (2000) & $Dist$  & $g$ & $i$ & scale  & Seeing & $\sigma_{velocity}$ & $morph$ & $Mag_{i-band}$ & \log $(M_*/M_\odot)_{total}$\\
    
DD = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list.txt")["Dist"]
scale = SRead.grab_dist("/home/dexter/result/distance/dir_dist_list.txt","/home/dexter/result/distance/dist_list.txt")["Scale"]


A = SRead.read_list("Gal_vdis_dict_equvi")

mag_g,mag_i = A["mag_g"], A["mag_i"]

#####second table#####

#   Name & $\mu_e$ & $R_e$ & $n$  & $mag_{Sph}$ & $Abs Mag_{Sph}$ & \log$(M_*/M_\odot)_{Taylor}$  & log$(M_*/M_\odot)_{Into}$ & log$(M_*/M_\odot)$ & log$(M_*/M_\odot)$  \\
