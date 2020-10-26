#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 16:49:02 2020

@author: dexter
"""

markers = ["E0","E1","E2","E3","EAS","EABS","EBS","SA0","SAB0", "SB0", "SAa", 
           "SABa", "SBa", "SAb","SABb", "SBb", "SAc","SABc", "SBb","SAc","SABc",
           "SBc","SAd","SABd","SBd","SAm","SABm","SBm"]

markers_group = ["E","ES","S0","Sa","Sb","Sc","Sd","Sm"]



#read the data


marker_morph_dict = {"E0":"./morph/E0.jpg",
               "E1":"./morph/E1-2.jpg",
               "E2":"./morph/E1-2.jpg",
               "E3":"./morph/E3.jpg",
               "EAS":"./morph/EAS.jpg",
               "EABS":"./morph/EABS.jpg",
               "EBS":"./morph/EBS.jpg",
               "SA0":"./morph/SA0.jpg",
               "SAB0":"./morph/SAB0.jpg",
               "SB0":"./morph/SB0.jpg",
               "SAa":"./morph/SAa.jpg",
               "SABa":"./morph/SABa.jpg",
               "SBa":"./morph/SBa.jpg",
               "SAb":"./morph/SAb.jpg",
               "SABb":"./morph/SABb.jpg",
               "SBb":"./morph/SBb.jpg",
               "SAc":"./morph/SAc.jpg",
               "SABc":"./morph/SABc.jpg",
               "SBc":"./morph/SBc.jpg",
               "SAd":"./morph/SAd.jpg",    
               "SABd":"./morph/SABd.jpg",
               "SBd":"./morph/SBd.jpg",
               "SAm":"./morph/SAm.jpg",
               "SABm":"./morph/SABm.jpg",
               "SBm":"./morph/SBm.jpg"
               }



#plot the size mass relation transision


