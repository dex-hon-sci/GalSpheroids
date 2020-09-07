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



box=[[8.0,8.3],[8.3,8.6],[8.6,8.9],
     [8.9,9.2],[9.2,9.5],[9.5,9.8],
     [9.8,10.1],[10.1,10.4],[10.4,10.7],
     [10.7,11.0],[11.0,11.3],[11.3,11.6],
     [11.6,11.9],[11.9,12.0]]

#box=[[7,8],[8,9],[9,10],[10,11],[11,12]]
#box =[[9.9,10.5],[10.5,11.0],[11.0,11.6]]
volume = [29526.108,97704.819,357422.506]

#volume = [357422.506,357422.506,357422.506]

#####
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
                 'line_style':['solid','solid','solid'],
                 'line_colour':['#b13322','#a2579b','#58906d'],
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
                 'line_style':['solid','solid','solid'],
                 'line_colour':['#b65050','#616a9a'],
                 'label':['Sph','Disk'],
                 'satu':[0.2,0.2]}
#############
M = np.linspace(10**8, 10**11.6, 10**3)


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

fig, ax = plt.subplots()


plt.plot(M,Phi,color="black", linestyle="solid")

for i in range(len(M_star)):
    Phi = func.Schechter_func(M, alpha[i], M_star[i], phi_0[i])
    plt.plot(M,Phi, color=colour[i]
             , label =label[i],
             linestyle=line_style[i])

    plt.xscale( 'log' )
    plt.yscale( 'log' )

plt.xlim(10**7.9,10**11.6)
plt.ylim(3*10**-5,8*10**-2 )


##########################################

D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt")


mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]


total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3_cpt", ["Bulge","CoreBulge"])


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


mass1 = np.log10(E1_T11)
mass2 = np.log10(E2_T11)
mass3 = np.log10(E3_T11)

###########################################

#mass1 = np.log10(SRead.read_list("Gal_table1_bin2_Tmass")["mass"]*1e10)
#mass2 = np.log10(SRead.read_list("Gal_table1_bin3_Tmass")["mass"]*1e10)
#mass3 = np.log10(SRead.read_list("Gal_table1_bin4_Tmass")["mass"]*1e10)

#print(np.log10(SRead.read_list("Gal_table1_bin2_Tmass")["mass"]*1e10))
#print(np.log10(SRead.read_list("Gal_table1_bin3_Tmass")["mass"]*1e10))
#print(np.log10(SRead.read_list("Gal_table1_bin4_Tmass")["mass"]*1e10))
print('----------------------------------------------------------------')
#print(np.log10(SRead.read_list("Gal_table2_bin2")["Taylor_mass"]*1e10))
print('----------------------------------------------------------------')
SPlot.ShowcaseIndi.mass_function_plot(mass3, box, volume[0], colour="black",label="mass_bin3")
SPlot.ShowcaseIndi.mass_function_plot(mass2, box, volume[1], colour="blue",label="mass_bin2")
SPlot.ShowcaseIndi.mass_function_plot(mass1, box, volume[2], colour="red",label="mass_bin1")

A = np.array([316227766.01683793,3162277660.1683793, 31622776601.683792, 316227766016.83795])
B = np.array([2,17+4+4,23+20+22,4+7])/((357422.506))
#plt.plot(A,B*1,'-o')

AA = np.array([177827941.00389227, 354813389.23357606, 707945784.3841373, 1412537544.6227555, 2818382931.264449, 5623413251.903491, 11220184543.019653, 22387211385.683376, 44668359215.09635, 79432823472.42822, 141253754462.27554, 281838293126.44495])

C = np.array([0.00011289443679245951+3.411636567622456e-05+0.0,
              0.0019192054254718116+0.00013646546270489824+3.730412357786259e-05,
              0.002596572046226569+0.0006823273135244912+0.00020517267967824424,    
              0.0+0.00013646546270489824+6.528221626125954e-05])

C1 = np.array([0.00011289443679245951, 0.0, 0.0, 0.00011289443679245951, 0.0003386833103773786, 0.0007902610575472167, 0.0014676276783019738, 0.0005644721839622977, 0.0009031554943396761, 0.00011289443679245951, 0.0, 0.0])
C2 = np.array([0.0, 0.0, 3.411636567622456e-05, 0.0, 3.411636567622456e-05, 0.00010234909702867368, 0.00020469819405734735, 0.00020469819405734735, 0.0001705818283811228, 0.00010234909702867368, 0.00010234909702867368, 3.411636567622456e-05])
C3 = np.array([0.0, 0.0, 0.0, 0.0, 1.8652061788931295e-05, 9.326030894465648e-06, 2.7978092683396946e-05, 6.528221626125954e-05, 5.595618536679389e-05, 5.595618536679389e-05, 2.7978092683396946e-05, 2.7978092683396946e-05])

CC = C1+C2+C3


#plt.plot(A,C,'-o')
#plt.plot(AA,CC,'-o')

AAA = np.array([141253754.46227497, 177827941.00389227, 223872113.85683468, 281838293.12644494, 354813389.23357606, 446683592.15096164, 562341325.1903491, 707945784.3841403, 891250938.1337441, 1122018454.3019652, 1412537544.6227496, 1778279410.0389228, 2238721138.568347, 2818382931.264449, 3548133892.3357606, 4466835921.509617, 5623413251.903491, 7079457843.841402, 8912509381.33744, 11220184543.019653, 14125375446.227497, 17782794100.38923, 22387211385.683468, 28183829312.644493, 35481338923.357605, 44668359215.09617, 56234132519.034904, 70794578438.41402, 89125093813.3744, 141253754462.27554, 223872113856.8347, 281838293126.44495, 354813389233.57605])

CC1 = np.array([0.0, 0.0003386833103773785, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003386833103773785, 0.0, 0.0, 0.0, 0.0003386833103773785, 0.0003386833103773785, 0.000677366620754757, 0.0010160499311321356, 0.000677366620754757, 0.0016934165518868928, 0.0010160499311321356, 0.0016934165518868928, 0.0010160499311321356, 0.0003386833103773785, 0.0003386833103773785, 0.0003386833103773785, 0.001354733241509514, 0.0010160499311321356, 0.0003386833103773785, 0.0, 0.0, 0.0, 0.0, 0.0])
CC2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00010234909702867368, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00010234909702867368, 0.0, 0.0, 0.000307047291086021, 0.0, 0.0004093963881146947, 0.00020469819405734735, 0.00020469819405734735, 0.000307047291086021, 0.0, 0.00020469819405734735, 0.0, 0.000307047291086021, 0.000307047291086021, 0.0, 0.000307047291086021, 0.00010234909702867368, 0.0, 0.0] )
CC3 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.7978092683396943e-05, 2.7978092683396943e-05, 0.0, 0.0, 0.0, 0.0, 2.7978092683396943e-05, 0.0, 5.5956185366793886e-05, 2.7978092683396943e-05, 8.393427805019084e-05, 5.5956185366793886e-05, 8.393427805019084e-05, 2.7978092683396943e-05, 2.7978092683396943e-05, 0.0001398904634169847, 2.7978092683396943e-05, 8.393427805019084e-05, 0.0, 5.5956185366793886e-05, 2.7978092683396943e-05])

CCC = CC1+CC2+CC3

#plt.plot(AAA,CCC,'-o')

plt.legend()
plt.show()