#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 23:59:09 2020

@author: dexter
"""

##Select from SphExecute.py
import SphRead as SRead
import SphSort as SSort
import SphPlot as SPlot

import numpy as np

c, H0=    299792.458, 100.0       #speed of light in km/s, Hubble constant in(km/s)/Mpc
M_sun = 4.53

A = SSort.read_list("Gal_vdis_dict_equvi")

mag_g, mag_i = A["mag_g"],A["mag_i"]

name2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Gal_name"]
DD2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Dist"]
scale2 = SRead.grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Scale"]

CCC2 = np.genfromtxt("Ctable_Into_bin2_final.txt" , dtype='float')


R_e_2 = SRead.grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag2 = SRead.grab_mag("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"])
mag2_total = SRead.grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Total_mag"],0)

Re_kpc2 = R_e_2 * scale2

Mag2 = mag2 - 25 -5*np.log10(DD2) # with new distance system, NED, hyperleda
Lum2 = 10**((Mag2-(4.53))/(-2.5))

#BD for bin2

name2BD = SRead.grab_parameter("Gal_bundle_BD_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1) 
R_e_2_BD = SRead.grab_parameter("Gal_bundle_BD_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag2_BD = SRead.grab_mag("Gal_bundle_BD_equvi_bin2_cpt", ["Bulge","CoreBulge"])

Re_kpc2_BD = R_e_2_BD * scale2

Mag2_BD = mag2_BD - 25 -5*np.log10(DD2)
Lum2_BD = 10**((Mag2_BD-(4.53))/(-2.5))



ML_select2 = SPlot.MLRelationIband(AbsMag_g_bin2,AbsMag_i_bin2).Into13_MassRatio
M2 = SPlot.MassCalculation(mag2, DD2, M_sun,AbsMag_g_bin2,AbsMag_i_bin2)
M2_total_old = SPlot.MassCalculation(mag2_total, dc2, M_sun,AbsMag_g_bin2,AbsMag_i_bin2)

E2 = M2.cal_Mass(ML_select2)
E2_total_old = M2_total_old.cal_Mass(ML_select2)

M2_BD = SPlot.MassCalculation(mag2_BD, DD2, M_sun,AbsMag_g_bin2,AbsMag_i_bin2)
E2_BD = M2_BD.cal_Mass(ML_select2)
