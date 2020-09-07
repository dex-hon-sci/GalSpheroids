#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 19:37:20 2020

@author: dexter
"""

import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot

#mass = np.linspace(2e8,0.5e13,2000)

#fig, ax = plt.subplots()

#plt.show()

override_list_maj = \
["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source", "NGC4570", 2, "Disk",\
 "UGC8736",2,"Disk", "NGC4583", 5,"IntDisk","NGC5382", 5,"IntDisk","NGC4772", 14 ,"Point Source", "NGC4845", 14, "Point Source",
 "NGC5375",11 , "Point Source"]

override_list_equ = \
["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source",\
 "UGC8736",2,"Disk","NGC4583", 5,"IntDisk", "NGC5382", 5,"IntDisk","NGC4772", 14 ,"Point Source", "NGC4845", 14, "Point Source",
 "NGC5375",11 , "Point Source"]

SRead.run_list("F_equvi_list_Bin1.txt","F_Gal_bundle_equvi_Bin1",True)
SRead.run_list("F_equvi_list_Bin2.txt","F_Gal_bundle_equvi_Bin2",True)
SRead.run_list("F_equvi_list_Bin3.txt","F_Gal_bundle_equvi_Bin3",True)


SRead.run_list("F_BD_equvi_list_Bin1.txt","F_Gal_bundle_BD_equvi_Bin1",True)
SRead.run_list("F_BD_equvi_list_Bin2.txt","F_Gal_bundle_BD_equvi_Bin2",True)
SRead.run_list("F_BD_equvi_list_Bin3.txt","F_Gal_bundle_BD_equvi_Bin3",True)


C2,C3,C4 = SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin1'), SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin2'),SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin3')
A2,A3,A4 = SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin1'), SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin2'),SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin3')

#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)

SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin1',C2,'F_Gal_bundle_equvi_Bin1_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin2',C3,'F_Gal_bundle_equvi_Bin2_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin3',C4,'F_Gal_bundle_equvi_Bin3_cpt',override_list_equ)


SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin1',A2,'F_Gal_bundle_BD_equvi_Bin1_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin2',A3,'F_Gal_bundle_BD_equvi_Bin2_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin3',A4,'F_Gal_bundle_BD_equvi_Bin3_cpt',[])

###############
SRead.run_list("F_BD_major_list.txt","F_Gal_bundle_BD_major", False)
SRead.run_list("F_BD_equvi_list.txt","F_Gal_bundle_BD_equvi", True)

SRead.run_list("F_major_list.txt","F_Gal_bundle_major",False)
SRead.run_list("F_equvi_list.txt","F_Gal_bundle_equvi",True)

Am, Ae =SSort.cpt_seperator_demo('F_Gal_bundle_BD_major'), SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi')
B, C = SSort.cpt_seperator_demo('F_Gal_bundle_major'), SSort.cpt_seperator_demo('F_Gal_bundle_equvi')

SSort.cpt_classifier_demo('F_Gal_bundle_BD_major',Am,'F_Gal_bundle_BD_major_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi',Ae,'F_Gal_bundle_BD_equvi_cpt',[])

SSort.cpt_classifier_demo('F_Gal_bundle_major',B,'F_Gal_bundle_major_cpt',override_list_maj)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi',C,'F_Gal_bundle_equvi_cpt',override_list_equ)
