#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 18:35:04 2021

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

import numpy as np
import pickle

from astropy.cosmology import FlatLambdaCDM

name = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_BD_equvi_cpt")


total_mag = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_cpt")

total_mag_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_BD_equvi_cpt")


total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt")

total_mag1_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt")
total_mag2_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt")
total_mag3_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3_cpt", ["Bulge","CoreBulge"])

table = {"namw": name, 
         "total_mag": total_mag,
         "total_mag_BD": total_mag_BD} 

value = list(table.values())
key = list(table.keys())

for i in range(len(list(table.keys()))):
    print(key[i],len(value[i]))
#    
with open("Gal_mag_BD_and_multi", 'wb') as f:
    pickle.dump(table, f)

SRead.convert_dict_ascii("Gal_mag_BD_and_multi","Gal_mag_BD_and_multi.txt")