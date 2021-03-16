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

    
###############
#Run based on Bins equvi
SRead.run_list("F_equvi_list_Bin1V.txt","F_Gal_bundle_equvi_Bin1V",True)
SRead.run_list("F_equvi_list_Bin2V.txt","F_Gal_bundle_equvi_Bin2V",True)
SRead.run_list("F_equvi_list_Bin3V.txt","F_Gal_bundle_equvi_Bin3V",True)


SRead.run_list("F_BD_equvi_list_Bin1V.txt","F_Gal_bundle_BD_equvi_Bin1V",True)
SRead.run_list("F_BD_equvi_list_Bin2V.txt","F_Gal_bundle_BD_equvi_Bin2V",True)
SRead.run_list("F_BD_equvi_list_Bin3V.txt","F_Gal_bundle_BD_equvi_Bin3V",True)


C2,C3,C4 = SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin1V'), SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin2V'),SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin3V')
A2,A3,A4 = SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin1V'), SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin2V'),SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_Bin3V')

#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)

SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin1V',C2,'F_Gal_bundle_equvi_Bin1V_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin2V',C3,'F_Gal_bundle_equvi_Bin2V_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin3V',C4,'F_Gal_bundle_equvi_Bin3V_cpt',override_list_equ)


SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin1V',A2,'F_Gal_bundle_BD_equvi_Bin1V_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin2V',A3,'F_Gal_bundle_BD_equvi_Bin2V_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_Bin3V',A4,'F_Gal_bundle_BD_equvi_Bin3V_cpt',[])

###############
#Run based on Bins major
SRead.run_list("F_major_list_Bin1V.txt","F_Gal_bundle_major_Bin1V",False)
SRead.run_list("F_major_list_Bin2V.txt","F_Gal_bundle_major_Bin2V",False)
SRead.run_list("F_major_list_Bin3V.txt","F_Gal_bundle_major_Bin3V",False)


SRead.run_list("F_BD_major_list_Bin1V.txt","F_Gal_bundle_BD_major_Bin1V",False)
SRead.run_list("F_BD_major_list_Bin2V.txt","F_Gal_bundle_BD_major_Bin2V",False)
SRead.run_list("F_BD_major_list_Bin3V.txt","F_Gal_bundle_BD_major_Bin3V",False)


C2_m,C3_m,C4_m = SSort.cpt_seperator_demo('F_Gal_bundle_major_Bin1V'), SSort.cpt_seperator_demo('F_Gal_bundle_major_Bin2V'),SSort.cpt_seperator_demo('F_Gal_bundle_major_Bin3V')
A2_m,A3_m,A4_m = SSort.cpt_seperator_demo('F_Gal_bundle_BD_major_Bin1V'), SSort.cpt_seperator_demo('F_Gal_bundle_BD_major_Bin2V'),SSort.cpt_seperator_demo('F_Gal_bundle_BD_major_Bin3V')

#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)

SSort.cpt_classifier_demo('F_Gal_bundle_major_Bin1V',C2_m,'F_Gal_bundle_major_Bin1V_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_major_Bin2V',C3_m,'F_Gal_bundle_major_Bin2V_cpt',override_list_equ)
SSort.cpt_classifier_demo('F_Gal_bundle_major_Bin3V',C4_m,'F_Gal_bundle_major_Bin3V_cpt',override_list_equ)


SSort.cpt_classifier_demo('F_Gal_bundle_BD_major_Bin1V',A2_m,'F_Gal_bundle_BD_major_Bin1V_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_major_Bin2V',A3_m,'F_Gal_bundle_BD_major_Bin2V_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_major_Bin3V',A4_m,'F_Gal_bundle_BD_major_Bin3V_cpt',[])

###############
#Run all
SRead.run_list("F_BD_major_list.txt","F_Gal_bundle_BD_major_V", False)
SRead.run_list("F_BD_equvi_list.txt","F_Gal_bundle_BD_equvi_V", True)

SRead.run_list("F_major_list.txt","F_Gal_bundle_major_V",False)
SRead.run_list("F_equvi_list.txt","F_Gal_bundle_equvi_V",True)

Am, Ae =SSort.cpt_seperator_demo('F_Gal_bundle_BD_major_V'), SSort.cpt_seperator_demo('F_Gal_bundle_BD_equvi_V')
B, C = SSort.cpt_seperator_demo('F_Gal_bundle_major_V'), SSort.cpt_seperator_demo('F_Gal_bundle_equvi_V')

SSort.cpt_classifier_demo('F_Gal_bundle_BD_major_V',Am,'F_Gal_bundle_BD_major_V_cpt',[])
SSort.cpt_classifier_demo('F_Gal_bundle_BD_equvi_V',Ae,'F_Gal_bundle_BD_equvi_V_cpt',[])

SSort.cpt_classifier_demo('F_Gal_bundle_major_V',B,'F_Gal_bundle_major_V_cpt',override_list_maj)
SSort.cpt_classifier_demo('F_Gal_bundle_equvi_V',C,'F_Gal_bundle_equvi_V_cpt',override_list_equ)

print(len(SRead.grab_name("F_Gal_bundle_equvi_V_cpt")),
      len(SRead.grab_name("F_Gal_bundle_major_V_cpt")),
      len(SRead.grab_name("F_Gal_bundle_equvi_Bin1V_cpt")),
      len(SRead.grab_name("F_Gal_bundle_major_Bin1V_cpt"))
      )