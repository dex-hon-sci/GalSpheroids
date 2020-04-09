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
 "UGC8736",2,"Disk" ]

override_list_equ = \
["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source",\
 "UGC8736",2,"Disk" ]

SRead.run_list("equvi_list_bin2.txt","Gal_bundle_equvi_bin2",True)
SRead.run_list("equvi_list_bin3.txt","Gal_bundle_equvi_bin3",True)
SRead.run_list("equvi_list_bin4.txt","Gal_bundle_equvi_bin4",True)


SRead.run_list("BD_equvi_list_bin2.txt","Gal_bundle_BD_equvi_bin2",True)
SRead.run_list("BD_equvi_list_bin3.txt","Gal_bundle_BD_equvi_bin3",True)
SRead.run_list("BD_equvi_list_bin4.txt","Gal_bundle_BD_equvi_bin4",True)


C2,C3,C4 = SSort.cpt_seperator_demo('Gal_bundle_equvi_bin2'), SSort.cpt_seperator_demo('Gal_bundle_equvi_bin3'),SSort.cpt_seperator_demo('Gal_bundle_equvi_bin4')
A2,A3,A4 = SSort.cpt_seperator_demo('Gal_bundle_BD_equvi_bin2'), SSort.cpt_seperator_demo('Gal_bundle_BD_equvi_bin3'),SSort.cpt_seperator_demo('Gal_bundle_BD_equvi_bin4')

#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)

SSort.cpt_classifier_demo('Gal_bundle_equvi_bin2',C2,'Gal_bundle_equvi_bin2_cpt',override_list_equ)
SSort.cpt_classifier_demo('Gal_bundle_equvi_bin3',C3,'Gal_bundle_equvi_bin3_cpt',override_list_equ)
SSort.cpt_classifier_demo('Gal_bundle_equvi_bin4',C4,'Gal_bundle_equvi_bin4_cpt',override_list_equ)


SSort.cpt_classifier_demo('Gal_bundle_BD_equvi_bin2',A2,'Gal_bundle_BD_equvi_bin2_cpt',[])
SSort.cpt_classifier_demo('Gal_bundle_BD_equvi_bin3',A3,'Gal_bundle_BD_equvi_bin3_cpt',[])
SSort.cpt_classifier_demo('Gal_bundle_BD_equvi_bin4',A4,'Gal_bundle_BD_equvi_bin4_cpt',[])

###############
SRead.run_list("BD_major_list.txt","Gal_bundle_BD_major", False)
SRead.run_list("BD_equvi_list.txt","Gal_bundle_BD_equvi", True)

SRead.run_list("major_list_2.txt","Gal_bundle_major",False)
SRead.run_list("equvi_list_2.txt","Gal_bundle_equvi",True)

Am, Ae =SSort.cpt_seperator_demo('Gal_bundle_BD_major'), SSort.cpt_seperator_demo('Gal_bundle_BD_equvi')
B, C = SSort.cpt_seperator_demo('Gal_bundle_major'), SSort.cpt_seperator_demo('Gal_bundle_equvi')

SSort.cpt_classifier_demo('Gal_bundle_BD_major',Am,'Gal_bundle_BD_major_cpt',[])
SSort.cpt_classifier_demo('Gal_bundle_BD_equvi',Ae,'Gal_bundle_BD_equvi_cpt',[])

SSort.cpt_classifier_demo('Gal_bundle_major',B,'Gal_bundle_major_cpt',override_list_maj)
SSort.cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)
