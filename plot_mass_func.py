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
import astropy as astropy
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck13, z_at_value
from astropy import units as u
from astropy.coordinates import Distance

import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0

# Box collection
box_01=[[8.1,8.2], [8.2,8.3], [8.3,8.4], [8.4,8.5], [8.5,8.6],
     [8.6,8.7], [8.7,8.8], [8.8,8.9], [8.9,9.0], [9.0,9.1],
     [9.1,9.2], [9.2,9.3], [9.3,9.4], [9.4,9.5], [9.5,9.6],
     [9.6,9.7], [9.7,9.8], [9.8,9.9], [9.9,10.0], [10.0,10.1],
     [10.1,10.2], [10.2,10.3], [10.3,10.4], [10.4,10.5], [10.5,10.6],
     [10.6,10.7], [10.7,10.8], [10.8,10.9], [10.9,11.0], [11.0,11.3],
     [11.3,11.4], [11.4,11.5], [11.5,11.6]]

box_03=[[8.0,8.3],[8.3,8.6],[8.6,8.9],
     [8.9,9.2],[9.2,9.5],[9.5,9.8],
     [9.8,10.1],[10.1,10.4],[10.4,10.7],
     [10.7,11.0],[11.0,11.3],[11.3,11.6],
     [11.6,11.9],[11.9,12.2],[12.2,12.5],[12.5,12.8]]

box_02 = [[8.0,8.2],[8.2,8.4],[8.4,8.6],[8.6,8.8],[8.8,9.0],[9.0,9.2],[9.2,9.4],
          [9.4,9.6],[9.6,9.8],[9.8,10.0],[10.0,10.2],[10.2,10.4],[10.4,10.6],
          [10.6,10.8],[10.8,11.0],[11.0,11.2],[11.2,11.4],[11.4,11.6],[11.6,11.8],
          [11.8,12.0]]

box_10=[[7,8],[8,9],[9,10],[10,11],[11,12]]
#box =[[9.9,10.5],[10.5,11.0],[11.0,11.6]]

box_05=[[8.0,8.5],[8.5,9.0],[9.0,9.5],[9.5,10.0],
     [10.0,10.5],[10.5,11.0],[11.0,11.5],[11.5,12.0]]

box = box_03


# calculate the volume
D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180)-np.cos(np.pi/2)))

volume = voll

V1_V = volume[2] - volume[1]
V2_V = volume[1] - volume[0]
V3_V = volume[0]

#V1_V = volume[2] 
#V2_V = volume[1]
#V3_V = volume[0]


print('voll',voll)



#old volume
#volume = [29526.108,97704.819,357422.506]


#volume = [357422.506,357422.506,357422.506]


def plot_solid_circles():
    Solid_Bin1_mass = 10**np.array([np.average(box[12]),np.average(box[13])])
    Solid_Bin1_nudens = np.array([3,2])/V1_V/0.3

    Solid_Bin2_mass = 10**np.array([np.average(box[10])])
    Solid_Bin2_nudens = np.array([6])/V2_V/0.3

    Solid_Bin3_mass = 10**np.array([np.average(box[9]), np.average(box[11])])
    Solid_Bin3_nudens = np.array([9,1])/V3_V/0.3
    
    ax.plot(Solid_Bin3_mass, Solid_Bin3_nudens, "o", color='#2a3236',label=r"$\rm Bin~3$",ms=14)
    ax.plot(Solid_Bin2_mass, Solid_Bin2_nudens, "o", color='#0b5786',label=r"$\rm Bin~2$",ms=14)
    ax.plot(Solid_Bin1_mass, Solid_Bin1_nudens, "o",color='#a5200b',label=r"$\rm Bin~1$",ms=14)

def plot_half_circles():
    HalfSolid_Bin1_mass = 10**np.array([np.average(box[11])])
    HalfSolid_Bin1_nudens = np.array([3])/V1_V/0.3

    HalfSolid_Bin2_mass = 10**np.array([np.average(box[9])])
    HalfSolid_Bin2_nudens = np.array([4])/V2_V/0.3

    HalfSolid_Bin3_mass = 10**np.array([np.average(box[7]), np.average(box[8])])
    HalfSolid_Bin3_nudens = np.array([9,13])/V3_V/0.3
    
    ax.plot(HalfSolid_Bin3_mass, HalfSolid_Bin3_nudens, "o", fillstyle = "bottom",color='#2a3236',ms=14, alpha=0.3)
    ax.plot(HalfSolid_Bin2_mass, HalfSolid_Bin2_nudens, "o", fillstyle = "bottom",color='#0b5786',ms=14, alpha=0.3)
    ax.plot(HalfSolid_Bin1_mass, HalfSolid_Bin1_nudens, "o", fillstyle = "bottom",color='#a5200b',ms=14, alpha=0.3)

##########################################
# input, horizontal bins
# "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt"
# input, vertical bins

# Calculate the local Sph mass
D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_2.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_2.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_2.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1V_4_2.txt",dtype="str")
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2V_4_2.txt",dtype="str")
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3V_4_2.txt",dtype="str")


mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

# calculate the redshift with these distance
cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)

#z_1 = z_at_value(Planck13.distmod,100*u.Mpc)

# convert to Luminosity distance

z_1,z_2,z_3 = [], [], []
for i in range(len(D1)):
    z_1.append(Distance(D1[i],unit=u.Mpc).compute_z(cosmology=cosmo))
for i in range(len(D2)):
    z_2.append(Distance(D2[i],unit=u.Mpc).compute_z(cosmology=cosmo))
for i in range(len(D3)):
    z_3.append(Distance(D3[i],unit=u.Mpc).compute_z(cosmology=cosmo))
z_1,z_2,z_3 = np.array(z_1), np.array(z_2), np.array(z_3)

print("D1",D1)
print("z1",z_1)
print("DL1",D1/(1+z_1))

#calculate the Cosmological dimming correction 

I1_corr_factor= (1+z_1)**(-4)
I2_corr_factor= (1+z_2)**(-4)
I3_corr_factor= (1+z_3)**(-4)

# extract the total magnitude
total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1V_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2V_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3V_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1V_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2V_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3V_cpt", ["Bulge","CoreBulge"])

# magntiude correction (cosmological dimming)
Extinction_g1 = list(D0_Bin1_table[:,-4])
Extinction_g2 = list(D0_Bin2_table[:,-4])
Extinction_g3 = list(D0_Bin3_table[:,-4])

Extinction_i1 = list(D0_Bin1_table[:,-3])
Extinction_i2 = list(D0_Bin2_table[:,-3])
Extinction_i3 = list(D0_Bin3_table[:,-3])

sph_mag1 = sph_mag1 - Extinction_i1
sph_mag2 = sph_mag2 - Extinction_i2
sph_mag3 = sph_mag3 - Extinction_i3

mag_g1, mag_i1 = mag_g1 - Extinction_g1, mag_i1 - Extinction_i1
mag_g2, mag_i2 = mag_g2 - Extinction_g2, mag_i2 - Extinction_i2
mag_g3, mag_i3 = mag_g3 - Extinction_g3, mag_i3 - Extinction_i3

# extinction correction
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

# make array of 2 kpc
array1_2kpc = np.repeat(2.0,len(E1_R15BC))
array2_2kpc = np.repeat(2.0,len(E2_R15BC))
array3_2kpc = np.repeat(2.0,len(E3_R15BC))

print(np.size(array1_2kpc),np.size(array2_2kpc),np.size(array3_2kpc))
print(np.size(Re_1_kpc),np.size(Re_2_kpc),np.size(Re_3_kpc))

# Choose mass
mass1 = np.log10(E1_R15BC)
mass2 = np.log10(E2_R15BC)
mass3 = np.log10(E3_R15BC)

# select for the data pount that has R_e < 2 kpc
E1_2kpc = SSort.selection_generic(10**mass1, Re_1_kpc, array1_2kpc, direction="low")['bag_x']
E2_2kpc = SSort.selection_generic(10**mass2, Re_2_kpc, array2_2kpc, direction="low")['bag_x']
E3_2kpc = SSort.selection_generic(10**mass3, Re_3_kpc, array3_2kpc, direction="low")['bag_x']


print("Re_2kpc",SSort.selection_generic(E1_R15BC, Re_1_kpc, array1_2kpc, direction="low",axis="y")['bag_y']
      , E1_2kpc)

# define the mass of data points with R_e < 2 kpc
mass1_2kpc = np.log10(E1_2kpc)
mass2_2kpc = np.log10(E2_2kpc)
mass3_2kpc = np.log10(E3_2kpc)

print("mass1_2kpc",mass1_2kpc)

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
fig, ax = plt.subplots(figsize=(12.2,6.4))

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

nu_dens3_t, mid_pt3_t = SPlot.ShowcaseIndi.mass_function_plot(mass3, box, V3_V, 
                                                          colour='#2a3236',
                                                          label="",
                                                          trim=False)
nu_dens2_t, mid_pt2_t = SPlot.ShowcaseIndi.mass_function_plot(mass2, box, V2_V, 
                                                          colour='#0b5786',
                                                          label="",
                                                          trim=False)
nu_dens1_t, mid_pt1_t = SPlot.ShowcaseIndi.mass_function_plot(mass1, box, V1_V, 
                                                          colour='#a5200b',
                                                          label="",
                                                          trim=False)


nu_dens_t_sum = (nu_dens1_t*V1_V*0.3 +nu_dens2_t*V2_V*0.3 +nu_dens3_t*V3_V*0.3)/volume[2]/0.3


nu_dens3_2kpc, mid_pt3_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass3_2kpc, box, V3_V, 
                                                          colour='#2a3236',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)
nu_dens2_2kpc, mid_pt2_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass2_2kpc, box, V2_V, 
                                                          colour='#0b5786',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)
nu_dens1_2kpc, mid_pt1_2kpc = SPlot.ShowcaseIndi.mass_function_plot(mass1_2kpc, box, V1_V, 
                                                          colour='#a5200b',
                                                          label="",
                                                          trim=False,
                                                          plot_yes=False)



# create mass function for sph Re<2kpc
nu_dens_2kpc = (nu_dens1_2kpc*V1_V*0.3 +nu_dens2_2kpc*V2_V*0.3 +nu_dens3_2kpc*V3_V*0.3)/volume[2]/0.3

#set limit for each bin
Bin1_limit, Bin2_limit, Bin3_limit = 3.4e11, 1.3e11, 6.7e10
Bin1_limit_az, Bin2_limit_az, Bin3_limit_az = 3.4e11*0.42, 1.3e11*0.36, 6.7e10*0.24
#Bin1_sigma, Bin2_sigma, Bin3_sigma = Bin1_limit*0.3, Bin2_limit*0.3, Bin3_limit*0.3

#Bin1_l,Bin1_u = Bin1_limit - Bin1_sigma, Bin1_limit + Bin1_sigma
#Bin2_l,Bin2_u = Bin2_limit - Bin2_sigma, Bin2_limit + Bin2_sigma
#Bin3_l,Bin3_u = Bin3_limit - Bin3_sigma, Bin3_limit + Bin3_sigma

#Bin1_shade_x, Bin1_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin1_l,Bin1_u, Bin1_u,Bin1_l]
#Bin2_shade_x, Bin2_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin2_l,Bin2_u, Bin2_u,Bin2_l]
#Bin3_shade_x, Bin3_shade_y = [1e-6,1e-2,1e-2,1e-6],[Bin3_l,Bin3_u, Bin3_u,Bin3_l]

high_end =5e12

Bin1_shade_y, Bin1_shade_x = [3.1e-6,3.6e-6,3.6e-6,3.1e-6], [Bin1_limit, Bin1_limit, high_end, high_end]
Bin2_shade_y, Bin2_shade_x = [4.0e-6,4.6e-6,4.6e-6,4.0e-6], [Bin2_limit, Bin2_limit, high_end, high_end]
Bin3_shade_y, Bin3_shade_x = [5.0e-6,5.9e-6,5.9e-6,5.0e-6], [Bin3_limit, Bin3_limit, high_end, high_end]

Bin1_shade_az_y, Bin1_shade_az_x = [3.1e-6,3.6e-6,3.6e-6,3.1e-6], [Bin1_limit_az, Bin1_limit_az, high_end, high_end]
Bin2_shade_az_y, Bin2_shade_az_x = [4.0e-6,4.6e-6,4.6e-6,4.0e-6], [Bin2_limit_az, Bin2_limit_az, high_end, high_end]
Bin3_shade_az_y, Bin3_shade_az_x = [5.0e-6,5.9e-6,5.9e-6,5.0e-6], [Bin3_limit_az, Bin3_limit_az, high_end, high_end]

#Bin1_shade_y, Bin1_shade_x = [3.1e-3,3.6e-6,3.6e-6,3.1e-3], [Bin1_limit, Bin1_limit, high_end, high_end]
#Bin2_shade_y, Bin2_shade_x = [4.0e-3,4.6e-6,4.6e-6,4.0e-3], [Bin2_limit, Bin2_limit, high_end, high_end]
#Bin3_shade_y, Bin3_shade_x = [5.0e-3,5.9e-6,5.9e-6,5.0e-3], [Bin3_limit, Bin3_limit, high_end, high_end]

ax.fill(Bin1_shade_x, Bin1_shade_y, '#a5200b',alpha=0.7)
ax.fill(Bin2_shade_x, Bin2_shade_y, '#0b5786',alpha=0.7)
ax.fill(Bin3_shade_x, Bin3_shade_y, '#2a3236',alpha=0.7)

ax.fill(Bin1_shade_az_x, Bin1_shade_az_y, '#a5200b',alpha=0.3)
ax.fill(Bin2_shade_az_x, Bin2_shade_az_y, '#0b5786',alpha=0.3)
ax.fill(Bin3_shade_az_x, Bin3_shade_az_y, '#2a3236',alpha=0.3)


print("nu_dens_t_sum",nu_dens_t_sum)
print("mid_pt1_t",mid_pt1_t)

nu_dens_t_sum[nu_dens_t_sum==0] = np.nan
nu_dens_2kpc[nu_dens_2kpc==0] = np.nan

ax.plot(mid_pt1_t, nu_dens_t_sum,"o--",label=r"$n_\mathrm{Sph}$",linewidth=5)
ax.plot(mid_pt1_t,nu_dens_2kpc, 'o--',label=r"$<~2~\rm kpc$",linewidth=5)

plot_solid_circles()
plot_half_circles()


#ax.vlines(x=Bin1_limit,ymin=1e-2,ymax=1e-6,color='#a5200b',lw=5, alpha=0.7)
#ax.vlines(x=Bin2_limit,ymin=1e-2,ymax=1e-6,color='#0b5786',lw=5, alpha=0.7)
#ax.vlines(x=Bin3_limit,ymin=1e-2,ymax=1e-6,color='#2a3236',lw=5, alpha=0.7)

ax.set_xlabel(r"$M_\mathrm{*,Sph}/ \rm M_{\odot}(RC15)$",fontsize=16)

#plt.grid(True)
ax.set_ylim(2e-6,6e-3)
ax.set_xlim(4e8,2e12)
plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
plt.show()