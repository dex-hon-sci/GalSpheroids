#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 00:04:02 2020

@author: dexter
"""


import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import astro_func_list as func

import numpy as np

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

#box=[[8.1,8.2],
#     [8.2,8.3],
#     [8.3,8.4],
#     [8.4,8.5],
#     [8.5,8.6],
#     [8.6,8.7],
#     [8.7,8.8],
#     [8.8,8.9],
#     [8.9,9.0],
#     [9.0,9.1],
#     [9.1,9.2],  
#     [9.2,9.3],
#     [9.3,9.4],
#     [9.4,9.5],
#     [9.5,9.6],
#     [9.6,9.7],
#     [9.7,9.8],
#     [9.8,9.9],
#     [9.9,10.0],
#     [10.0,10.1],
#     [10.1,10.2],
#     [10.2,10.3],
#     [10.3,10.4],
#     [10.4,10.5],
#     [10.5,10.6],
#     [10.6,10.7],
#     [10.7,10.8],
#     [10.8,10.9],
#     [10.9,11.0],
#     [11.0,11.3],
#     [11.3,11.4],
#     [11.4,11.5],
#     [11.5,11.6]]

# calculate the volume
D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180)-np.cos(np.pi/2)))

print('voll',voll)

box=[[8.0,8.3],[8.3,8.6],[8.6,8.9],
     [8.9,9.2],[9.2,9.5],[9.5,9.8],
     [9.8,10.1],[10.1,10.4],[10.4,10.7],
     [10.7,11.0],[11.0,11.3],[11.3,11.6],
     [11.6,11.9],[11.9,12.2],[12.2,12.5],[12.5,12.8]]

#box=[[7,8],[8,9],[9,10],[10,11],[11,12]]
#box =[[9.9,10.5],[10.5,11.0],[11.0,11.6]]

#old volume
#volume = [29526.108,97704.819,357422.506]

volume = voll

#volume = [357422.506,357422.506,357422.506]

##########################################
# input, horizontal bins
# "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt"
# input, vertical bins

# Calculate the local Sph mass
D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V.txt")

mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])

Re_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"], 1) #get Re

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale

scale1 = D1* ars
scale2 = D2* ars
scale3 = D3* ars

Re_1_kpc = Re_1* scale1
Re_2_kpc = Re_2* scale2
Re_3_kpc = Re_3* scale3

# Calculate mass
ML_select1_IP13 = SPlot.MLRelationIband(mag_g1,mag_i1).Into13_MassRatio
ML_select1_R15BC = SPlot.MLRelationIband(mag_g1,mag_i1).Roediger15BC03_MassRatio
ML_select1_Z09 = SPlot.MLRelationIband(mag_g1,mag_i1).Zibetti09_MassRatio
ML_select1_T11 = SPlot.MLRelationIband(mag_g1,mag_i1).Taylor11_MassRatio

M1 = SPlot.MassCalculation(sph_mag1, D1, 4.53,mag_g1,mag_i1)

E1_IP13 = M1.cal_Mass(ML_select1_IP13)
E1_R15BC = M1.cal_Mass(ML_select1_R15BC)
E1_Z09 = M1.cal_Mass(ML_select1_Z09)
E1_T11 = M1.cal_Mass(ML_select1_T11)


ML_select2_IP13 = SPlot.MLRelationIband(mag_g2,mag_i2).Into13_MassRatio
ML_select2_R15BC = SPlot.MLRelationIband(mag_g2,mag_i2).Roediger15BC03_MassRatio
ML_select2_Z09 = SPlot.MLRelationIband(mag_g2,mag_i2).Zibetti09_MassRatio
ML_select2_T11 = SPlot.MLRelationIband(mag_g2,mag_i2).Taylor11_MassRatio

M2 = SPlot.MassCalculation(sph_mag2, D2, 4.53,mag_g2,mag_i2)

E2_IP13 = M2.cal_Mass(ML_select2_IP13)
E2_R15BC = M2.cal_Mass(ML_select2_R15BC)
E2_Z09 = M2.cal_Mass(ML_select2_Z09)
E2_T11 = M2.cal_Mass(ML_select2_T11)


ML_select3_IP13 = SPlot.MLRelationIband(mag_g3,mag_i3).Into13_MassRatio
ML_select3_R15BC = SPlot.MLRelationIband(mag_g3,mag_i3).Roediger15BC03_MassRatio
ML_select3_Z09 = SPlot.MLRelationIband(mag_g3,mag_i3).Zibetti09_MassRatio
ML_select3_T11 = SPlot.MLRelationIband(mag_g3,mag_i3).Taylor11_MassRatio

M3 = SPlot.MassCalculation(sph_mag3, D3, 4.53,mag_g3,mag_i3)

E3_IP13 = M3.cal_Mass(ML_select3_IP13)
E3_R15BC = M3.cal_Mass(ML_select3_R15BC)
E3_Z09 = M3.cal_Mass(ML_select3_Z09)
E3_T11 = M3.cal_Mass(ML_select3_T11)


array1_2kpc = np.repeat(2.0,len(E1_R15BC))
array2_2kpc = np.repeat(2.0,len(E2_R15BC))
array3_2kpc = np.repeat(2.0,len(E3_R15BC))

print(array1_2kpc)
print(np.size(array1_2kpc),np.size(array2_2kpc),np.size(array3_2kpc))
print(np.size(Re_1_kpc),np.size(Re_2_kpc),np.size(Re_3_kpc))

mass1 = np.log10(E1_R15BC)
mass2 = np.log10(E2_R15BC)
mass3 = np.log10(E3_R15BC)

E1_2kpc = SSort.selection_generic(10**mass1, Re_1_kpc, array1_2kpc, direction="down")['bag_x']
E2_2kpc = SSort.selection_generic(10**mass2, Re_2_kpc, array2_2kpc, direction="down")['bag_x']
E3_2kpc = SSort.selection_generic(10**mass3, Re_3_kpc, array3_2kpc, direction="down")['bag_x']

print(SSort.selection_generic(E1_R15BC, Re_1_kpc, array1_2kpc, direction="down")['bag_y']
      , E1_2kpc)


mass1_2kpc = np.log10(E1_2kpc)
mass2_2kpc = np.log10(E2_2kpc)
mass3_2kpc = np.log10(E3_2kpc)

##### calculate the mass function (GAMA)
Kalvin14_all1 ={'M_star':[],
                 'M_star_lerr':[],
                 'M_star_uerr':[],
                 'alpha':[],
                 'alpha_lerr':[],
                 'alpha_uerr':[],
                 'phi_0':[],
                 'phi_0_lerr':[],
                 'phi_0_uerr':[],
                 'line_style':[],
                 'line_colour':[],
                 'label':[],
                 'satu':[]}

Kalvin14_morph ={'M_star':[10**(10.94), 10**(10.25),10**(10.09)],
                 'M_star_lerr':[0.10,0.07,0.15],
                 'M_star_uerr':[0.18,0.03,0.09],
                 'alpha':[-0.79,0.87,-0.01],
                 'alpha_lerr':[0.13,0.23,0.31],
                 'alpha_uerr':[0.23,0.15,0.26],
                 'phi_0':[ 0.85*1e-3,2.38*1e-3,3.57 * 1e-3],
                 'phi_0_lerr':[0.27,0.83,0.81],
                 'phi_0_uerr':[0.49,0.27,0.63],
                 'line_style':['dashdot','dashdot','dashdot'],
                 'line_colour':['#ee9119','#a2579b','#58906d'],
                 'label':['E','S0-Sa','Sab-Scd'],
                 'satu':[0.2,0.2,0.2]}

Kalvin14_dichotomy ={'M_star':[10.6,10.7],
                 'M_star_lerr':[0.05,0.23],
                 'M_star_uerr':[0.08,0.07],
                 'alpha':[-0.27,-1.37],
                 'alpha_lerr':[0.16,0.11],
                 'alpha_uerr':[0.20,0.04],
                 'phi_0':[3.96,0.42],
                 'phi_0_lerr':[1.05,0.42],
                 'phi_0_uerr':[0.37,0.14],
                 'line_style':['solid','solid'],
                 'line_colour':['#b65050','#616a9a'],
                 'label':['Sph','Disk'],
                 'satu':[0.2,0.2]}


M = np.linspace(10**8, 10**11.9, 10**3)

M_star = 10**(10.64)
alpha1, alpha2 = -0.43, -1.53
phi1_0, phi2_0 = 4.18 *1e-3, 0.74*1e-3

Phi1 = func.Schechter_func(M, alpha1, M_star, phi1_0)
Phi2 = func.Schechter_func(M, alpha2, M_star, phi2_0)
Phi = Phi1 + Phi2


import matplotlib.pyplot as plt

M_star = Kalvin14_morph["M_star"]
alpha = Kalvin14_morph["alpha"]
phi_0 = Kalvin14_morph["phi_0"]
colour = Kalvin14_morph['line_colour']
label = Kalvin14_morph['label']
line_style = Kalvin14_morph['line_style']        


##Ploting##################################

V1_V = volume[2] - volume[1]
V2_V = volume[1] - volume[0]
V3_V = volume[0]
fig, ax = plt.subplots()

print('V1','V2','V3',volume[2],volume[1],volume[0])

print('V1_V','V2_V','V3_V',V1_V,V2_V,V3_V)

#plt.plot(M,Phi,color="black", linestyle="solid", lw = 3,alpha=0.6)

#Plot GAMA mass function
#for i in range(len(M_star)):
#    Phi = func.Schechter_func(M, alpha[i], M_star[i], phi_0[i])
#    plt.plot(M,Phi, color=colour[i]
#             , label =label[i],
#             linestyle=line_style[i], lw=3,alpha=0.6)#
#
#    plt.xscale( 'log' )
#    plt.yscale( 'log' )
#
#plt.xlim(10**7.9,10**11.9)
#plt.ylim(2*10**-6,10**-2 )

#mass1 = np.log10(SRead.read_list("Gal_table1_bin2_Tmass")["mass"]*1e10)
#mass2 = np.log10(SRead.read_list("Gal_table1_bin3_Tmass")["mass"]*1e10)
#mass3 = np.log10(SRead.read_list("Gal_table1_bin4_Tmass")["mass"]*1e10)

#print(np.log10(SRead.read_list("Gal_table1_bin2_Tmass")["mass"]*1e10))
#print(np.log10(SRead.read_list("Gal_table1_bin3_Tmass")["mass"]*1e10))
#print(np.log10(SRead.read_list("Gal_table1_bin4_Tmass")["mass"]*1e10))

#nu_dens1, mid_pt1 = SPlot.ShowcaseIndi.mass_function_plot(mass3, box, V3_V, 
#                                                          colour='#2a3236',
#                                                          label="Bin3")
#nu_dens2, mid_pt2 = SPlot.ShowcaseIndi.mass_function_plot(mass2, box, V2_V, 
#                                                          colour='#0b5786',
#                                                          label="Bin2")
#nu_dens3, mid_pt3 = SPlot.ShowcaseIndi.mass_function_plot(mass1, box, V1_V, 
#                                                          colour='#a5200b',
#                                                          label="Bin1")

nu_dens1_t, mid_pt1_t = SPlot.ShowcaseIndi.mass_function_plot(mass3, box, V3_V, 
                                                          colour='#2a3236',
                                                          label="Bin 3",
                                                          trim=False)
nu_dens2_t, mid_pt2_t = SPlot.ShowcaseIndi.mass_function_plot(mass2, box, V2_V, 
                                                          colour='#0b5786',
                                                          label="Bin 2",
                                                          trim=False)
nu_dens3_t, mid_pt3_t = SPlot.ShowcaseIndi.mass_function_plot(mass1, box, V1_V, 
                                                          colour='#a5200b',
                                                          label="Bin 1",
                                                          trim=False)


nu_dens_t_sum = nu_dens1_t +nu_dens2_t +nu_dens3_t

nu_dens1_2kpc, mid_pt1_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass3_2kpc, box, V3_V, 
                                                          colour='#2a3236',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)
nu_dens2_2kpc, mid_pt2_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass2_2kpc, box, V2_V, 
                                                          colour='#0b5786',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)
nu_dens3_2kpc, mid_pt3_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass1_2kpc, box, V1_V, 
                                                          colour='#a5200b',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)

# create mass function for sph Re<2kpc
nu_dens_2kpc = nu_dens1_2kpc +nu_dens2_2kpc +nu_dens3_2kpc

#set limit for each bin
Bin1_limit, Bin2_limit, Bin3_limit = 2.11e11,6.70e10,2.41e10
Bin1_sigma, Bin2_sigma, Bin3_sigma = Bin1_limit*0.6, Bin2_limit*0.6, Bin3_limit*0.6

Bin1_l,Bin1_u = Bin1_limit - Bin1_sigma, Bin1_limit + Bin1_sigma
Bin2_l,Bin2_u = Bin2_limit - Bin2_sigma, Bin2_limit + Bin2_sigma
Bin3_l,Bin3_u = Bin3_limit - Bin3_sigma, Bin3_limit + Bin3_sigma

#Bin1_shade_x, Bin1_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin1_l,Bin1_u, Bin1_u,Bin1_l]
#Bin2_shade_x, Bin2_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin2_l,Bin2_u, Bin2_u,Bin2_l]
#Bin3_shade_x, Bin3_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin3_l,Bin3_u, Bin3_u,Bin3_l]

high_end =5e12

Bin1_shade_y, Bin1_shade_x = [3.1e-6,3.6e-6,3.6e-6,3.1e-6], [Bin1_limit, Bin1_limit, high_end, high_end]
Bin2_shade_y, Bin2_shade_x = [4.0e-6,4.6e-6,4.6e-6,4.0e-6], [Bin2_limit, Bin2_limit, high_end, high_end]
Bin3_shade_y, Bin3_shade_x = [5.0e-6,5.9e-6,5.9e-6,5.0e-6], [Bin3_limit, Bin3_limit, high_end, high_end]


#Bin1_shade_y, Bin1_shade_x = [3.1e-3,3.6e-6,3.6e-6,3.1e-3], [Bin1_limit, Bin1_limit, high_end, high_end]
#Bin2_shade_y, Bin2_shade_x = [4.0e-3,4.6e-6,4.6e-6,4.0e-3], [Bin2_limit, Bin2_limit, high_end, high_end]
#Bin3_shade_y, Bin3_shade_x = [5.0e-3,5.9e-6,5.9e-6,5.0e-3], [Bin3_limit, Bin3_limit, high_end, high_end]

ax.fill(Bin1_shade_x, Bin1_shade_y, '#a5200b',alpha=0.7)
ax.fill(Bin2_shade_x, Bin2_shade_y, '#0b5786',alpha=0.7)
ax.fill(Bin3_shade_x, Bin3_shade_y, '#2a3236',alpha=0.7)


ax.plot( mid_pt1_t , nu_dens_t_sum,"o--",label="Sum",linewidth=7)
ax.plot(mid_pt1_t,nu_dens_2kpc, 'o--',label="<2kpc",linewidth=7)

#ax.vlines(x=Bin1_limit,ymin=1e-2,ymax=1e-6,color='#a5200b',lw=5, alpha=0.7)
#ax.vlines(x=Bin2_limit,ymin=1e-2,ymax=1e-6,color='#0b5786',lw=5, alpha=0.7)
#ax.vlines(x=Bin3_limit,ymin=1e-2,ymax=1e-6,color='#2a3236',lw=5, alpha=0.7)

ax.set_xlabel(r"$M_*/ \rm M_{\odot}(IP13)$",fontsize=16)

#plt.grid(True)
ax.set_ylim(1e-6,1e-2)
ax.set_xlim(1e8,6e12)
plt.legend()
plt.show()