#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 12:14:11 2020

@author: dexter
"""

# solve problem ->error, circle->mark number density->

from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

##
D = np.array([45,75,110])



voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180)-np.cos(np.pi/2)))

V1,V2,V3=voll[2],voll[1],voll[0]

#voll[0],voll[1],voll[2]

print("volume", V1, V2, V3)
print("1/V", 1/V1, 1/V2, 1/V3)
print("2/V", 2/V1, 2/V2, 2/V3)
print("5/V", 5/V1, 5/V2, 5/V3)

###

D0_Bin1_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt")
D0_Bin2_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt")
D0_Bin3_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt")

D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

D0_Bin1_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin1.txt", 
    dtype = 'str')
D0_Bin2_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin2.txt",
    dtype = 'str')
D0_Bin3_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW_Bin3.txt",
    dtype = 'str')

D0_all_table_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

mag_g1, mag_i1 = D0_Bin1_table[:,11], D0_Bin1_table[:,10]
mag_g2, mag_i2 = D0_Bin2_table[:,11], D0_Bin2_table[:,10]
mag_g3, mag_i3 = D0_Bin3_table[:,11], D0_Bin3_table[:,10]

D1, D1_lerr, D1_uerr = D0_Bin1_table[:,29], D0_Bin1_table[:,30], D0_Bin1_table[:,31]
D2, D2_lerr, D2_uerr = D0_Bin2_table[:,29], D0_Bin2_table[:,30], D0_Bin2_table[:,31]
D3, D3_lerr, D3_uerr = D0_Bin3_table[:,29], D0_Bin3_table[:,30], D0_Bin3_table[:,31]

############### reading result files###############
master_file="/home/dexter/result/stat/completeness/master_file_h68dist_Intomass_RADEC_2.txt"

name1 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt")
name2 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt")
name3 = SRead.grab_name("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt")

total_mag1 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin1_cpt")
total_mag2 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin2_cpt")
total_mag3 = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_Bin3_cpt")

sph_mag1 = SRead.grab_mag("F_Gal_bundle_equvi_Bin1_cpt", ["Bulge","CoreBulge"])
sph_mag2 = SRead.grab_mag("F_Gal_bundle_equvi_Bin2_cpt", ["Bulge","CoreBulge"])
sph_mag3 = SRead.grab_mag("F_Gal_bundle_equvi_Bin3_cpt", ["Bulge","CoreBulge"])


Re_1 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin1_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_2 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin2_cpt", ["Bulge","CoreBulge"], 1) #get Re
Re_3 = SRead.grab_parameter("F_Gal_bundle_equvi_Bin3_cpt", ["Bulge","CoreBulge"], 1) #get Re

ars = (4.84814e-6)*1e3 # 1arcsec = (4.84814e-6) rad ars:arcsec to rad scale


scale1 = D1* ars
scale2 = D2* ars
scale3 = D3* ars

scale1_lerr, scale1_uerr = (D1-D1_lerr)*ars, (D1+D1_uerr)*ars
scale2_lerr, scale2_uerr =  (D2-D2_lerr)*ars, (D2+D2_uerr)*ars
scale3_lerr, scale3_uerr =  (D3-D3_lerr)*ars, (D3+D3_uerr)*ars

Re_1_kpc = Re_1* scale1
Re_2_kpc = Re_2* scale2
Re_3_kpc = Re_3* scale3

Re_1_kpc_lerr, Re_1_kpc_uerr = abs(Re_1_kpc - Re_1* scale1_lerr) , abs(Re_1* scale1_uerr - Re_1_kpc)
Re_2_kpc_lerr, Re_2_kpc_uerr = abs(Re_2* scale2_lerr - Re_2_kpc) , abs(Re_2* scale2_uerr - Re_2_kpc)
Re_3_kpc_lerr, Re_3_kpc_uerr = abs(Re_3* scale3_lerr - Re_3_kpc), abs(Re_3* scale3_uerr - Re_3_kpc)


Re_1_kpc_err =[Re_1_kpc_lerr, Re_1_kpc_uerr]
Re_2_kpc_err =[Re_2_kpc_lerr, Re_2_kpc_uerr]
Re_3_kpc_err =[Re_3_kpc_lerr, Re_3_kpc_uerr]

############# calculating mass

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

################################
#calculate the mass error


MLR1 = ML_select1_IP13
MLR_e1 = 10**0.1

MLR2 = ML_select2_IP13
MLR_e2 = 10**0.1

MLR3 = ML_select3_IP13
MLR_e3 = 10**0.1

mag_e = 0.6 #magnitude error

mass_uerr1 = np.sqrt(((mag_e/2.5)**2)+((2*D1_uerr/(D1*np.log(10)))**2)+((MLR_e1/(MLR1*np.log(10)))**2))
mass_uerr2 = np.sqrt(((mag_e/2.5)**2)+((2*D2_uerr/(D2*np.log(10)))**2)+((MLR_e2/(MLR2*np.log(10)))**2))
mass_uerr3 = np.sqrt(((mag_e/2.5)**2)+((2*D3_uerr/(D3*np.log(10)))**2)+((MLR_e3/(MLR3*np.log(10)))**2))

mass_err1 = mass_uerr1
mass_err2 = mass_uerr2
mass_err3 = mass_uerr3
##########################################
# read NSA-Sloan catalog
nsa = SRead.read_table('/home/dexter/result/stat/completeness/nsa_sizemass.dat')

nsa_z = nsa[:,0]

cosmo = FlatLambdaCDM(H0=68.0, Om0=0.3)
nsa_d = cosmo.comoving_distance(nsa_z).value

nsa_scale = nsa_d * ars

nsa_Re = nsa[:,1]*nsa_scale
nsa_mass_T11 = nsa[:,2]
nsa_mass_Z09 = nsa[:,3]
nsa_mass_R15BC = nsa[:,4]
nsa_mass_IP13 = nsa[:,5]


print(np.size(nsa_mass_IP13),np.size(nsa_Re))


# read Benzanson catalog

Benzanson2015 =SRead.read_table("/home/dexter/result/Bezanson2015table.txt")
Zahid2015 = SRead.read_table("/home/dexter/result/Zahid2015table.txt")

Benzanson2015_Re = np.array(Benzanson2015[:,4])
Benzanson2015_mass = np.array(10**Benzanson2015[:,7])

Zahid2015_Re  = np.array(Zahid2015[:,2])
Zahid2015_mass  = np.array(10**Zahid2015[:,4])



#################################
# Define ploting in RDJ15
def plot_dexter_sample_Bin():
    SPlot.ShowcaseIndi.Mass_Re_plot(E1_R15BC, Re_1_kpc, yerr = Re_1_kpc_err, 
                                xerr = mass_err1*E1_IP13,
                                colour = '#a5200b',
                                name=name1,legend='Bin1',ms=8,alpha0 = 0.2, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E2_R15BC, Re_2_kpc, yerr = Re_2_kpc_err, 
                                xerr = mass_err2*E2_IP13,
                                colour = '#0b5786',
                                name=name2,legend='Bin2',ms=8,alpha0 = 0.2, lw=3)
    SPlot.ShowcaseIndi.Mass_Re_plot(E3_R15BC, Re_3_kpc, yerr = Re_3_kpc_err, 
                                xerr = mass_err3*E3_IP13,
                                colour='#2a3236',
                                name=name1,legend='Bin3',ms=8,alpha0 = 0.2, lw=3)
    
def plot_dexter_sample_all():
    
    SPlot.ShowcaseIndi.Mass_Re_plot(E1_R15BC, Re_1_kpc, yerr = Re_1_kpc_err,
                                xerr = mass_err1*E1_R15BC,
                                colour = '#a5200b',
                                name=name1,legend='This work',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E2_R15BC, Re_2_kpc, yerr = Re_2_kpc_err,
                                xerr = mass_err2*E2_R15BC,
                                colour = '#a5200b',
                                name=name2,legend='',ms=8,alpha0 = 0.4, lw=3)

    SPlot.ShowcaseIndi.Mass_Re_plot(E3_R15BC, Re_3_kpc, yerr = Re_3_kpc_err,
                                xerr = mass_err3*E3_R15BC,
                                colour='#a5200b',
                                name=name1,legend='',ms=8,alpha0 = 0.4,lw=3 )
    
    
txsep1 = 40#10*5.5
txsep2 = 65#10*6.0
txsep3 = 80#10*6.5

x,y = 8.6e8,100



circle_size = 310
def plot_Barro_cut_all():
    Bcut1 = SPlot.SelectionCut(E1_R15BC, D1).Barro13_cut()
    Bcut2 = SPlot.SelectionCut(E2_R15BC, D2).Barro13_cut()
    Bcut3 = SPlot.SelectionCut(E3_R15BC, D3).Barro13_cut()

    S1 = SSort.selection_generic(E1_R15BC, Re_1_kpc, Bcut1)
    S2 = SSort.selection_generic(E2_R15BC, Re_2_kpc, Bcut2)
    S3 = SSort.selection_generic(E3_R15BC, Re_3_kpc, Bcut3)
    
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    plt.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    plt.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    plt.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    plt.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    plt.scatter(S1["bag_x"],S1["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S2["bag_x"],S2["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
def plot_vdWel_cut_all():
    vdWcut1 = SPlot.SelectionCut(E1_R15BC, D1).vdWel14_cut()
    vdWcut2 = SPlot.SelectionCut(E2_R15BC, D2).vdWel14_cut()
    vdWcut3 = SPlot.SelectionCut(E3_R15BC, D3).vdWel14_cut()

    S1 = SSort.selection_generic(E1_R15BC, Re_1_kpc, vdWcut1)
    S2 = SSort.selection_generic(E2_R15BC, Re_2_kpc, vdWcut2)
    S3 = SSort.selection_generic(E3_R15BC, Re_3_kpc, vdWcut3)
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    plt.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    plt.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    plt.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    plt.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)
    
    plt.scatter(S1["bag_x"],S1["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S2["bag_x"],S2["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
    
def plot_vDokkum_cut_all():
    vDcut1 = SPlot.SelectionCut(E1_R15BC, D1).vDokkum15_cut()
    vDcut2 = SPlot.SelectionCut(E2_R15BC, D2).vDokkum15_cut()
    vDcut3 = SPlot.SelectionCut(E3_R15BC, D3).vDokkum15_cut()

    S1 = SSort.selection_generic(E1_R15BC, Re_1_kpc, vDcut1)
    S2 = SSort.selection_generic(E2_R15BC, Re_2_kpc, vDcut2)
    S3 = SSort.selection_generic(E3_R15BC, Re_3_kpc, vDcut3)
    
    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    plt.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    plt.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    plt.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    plt.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)
    
    plt.scatter(S1["bag_x"],S1["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S2["bag_x"],S2["bag_y"],  facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)


def plot_Damjanov_cut_all():
    Damcut1 = SPlot.SelectionCut(E1_R15BC, D1).Damjanov14_cut()
    Damcut2 = SPlot.SelectionCut(E2_R15BC, D2).Damjanov14_cut()
    Damcut3 = SPlot.SelectionCut(E3_R15BC, D3).Damjanov14_cut()

    S1 = SSort.selection_generic(E1_R15BC, Re_1_kpc, Damcut1)
    S2 = SSort.selection_generic(E2_R15BC, Re_2_kpc, Damcut2)
    S3 = SSort.selection_generic(E3_R15BC, Re_3_kpc, Damcut3)

    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    plt.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    plt.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    plt.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    plt.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    plt.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    
def plot_Graham_cut_all():
    Gcut1 = SPlot.SelectionCut(E1_R15BC, D1).Graham15_broad_cut()
    Gcut2 = SPlot.SelectionCut(E2_R15BC, D2).Graham15_broad_cut()
    Gcut3 = SPlot.SelectionCut(E3_R15BC, D3).Graham15_broad_cut()

    S1 = SSort.selection_generic(E1_R15BC, Re_1_kpc, Gcut1)
    S2 = SSort.selection_generic(E2_R15BC, Re_2_kpc, Gcut2)
    S3 = SSort.selection_generic(E3_R15BC, Re_3_kpc, Gcut3)

    n1 = np.size(S1["bag_y"]) / V1
    n2 = np.size(S2["bag_y"]) / V2
    n3 = np.size(S3["bag_y"]) / V3
    
    
    plt.text(x, y, "$n(Mpc^{-3})$",fontsize=16)
    plt.text(x, y-txsep1, "Bin1: {:.2e} (N={})".format(n1,np.size(S1["bag_y"])),fontsize=16)
    plt.text(x, y-txsep2, "Bin2: {:.2e} (N={})".format(n2,np.size(S2["bag_y"])),fontsize=16)
    plt.text(x, y-txsep3, "Bin3: {:.2e} (N={})".format(n3,np.size(S3["bag_y"])),fontsize=16)

    plt.scatter(S1["bag_x"],S1["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S2["bag_x"],S2["bag_y"], facecolors='none', edgecolors='b', s = circle_size)
    plt.scatter(S3["bag_x"],S3["bag_y"], facecolors='none', edgecolors='b', s = circle_size)

###### Ploting ##################
xlim = [3e8,1.3e12]
ylim = [0.05,1e3]

mass0 = np.linspace(2e8,0.5e13,2000)
Dist0 = np.linspace(0,120,2000)

fig, ax = plt.subplots()        




counts,xbins,ybins,image = plt.hist2d(np.log10(nsa_mass_IP13),np.log10(nsa_Re),bins=(1000,1000)
                                      , norm=LogNorm()
                                      , cmap = plt.cm.gray, alpha=0.3)

plt.plot(np.log10(E1_IP13), np.log10(Re_1_kpc),'ro',ms=14)
plt.plot(np.log10(E2_IP13), np.log10(Re_2_kpc),'bo',ms=14)
plt.plot(np.log10(E3_IP13), np.log10(Re_3_kpc),'ko',ms=14)

plt.show()

fig, ax = plt.subplots()        


SPlot.ShowcaseIndi.Mass_Re_plot(nsa_mass_IP13, nsa_Re, marker='x', colour = 'grey',
                                ms = 5, legend='SDSS galaxies',alpha0 = 0.05)


#SPlot.SelectionCut(mass0,Dist0).plot_cut()


plot_dexter_sample_Bin()
plt.show()






fig, ax = plt.subplots()        
#SPlot.SelectionCut(mass0,Dist0).plot_cut()
SPlot.ShowcaseIndi.Mass_Re_plot(Benzanson2015_mass,Benzanson2015_Re, colour = '#4277cf',
                                legend='Benzanson et al. 2015',ms=8,alpha0 = 0.2, lw=3, marker = "^")

SPlot.ShowcaseIndi.Mass_Re_plot(Zahid2015_mass,Zahid2015_Re, colour = '#c17dd5',
                                legend='Zahid et al. 2015',ms=8,alpha0 = 0.2, lw=3, marker = "^")

plot_dexter_sample_all()

plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Damjanov","Damjanov et al. 2014",alpha0=0)
plot_dexter_sample_all()
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Damjanov","")

plot_Damjanov_cut_all()
plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro", "Barro et al. 2013", alpha0=0)
plot_dexter_sample_all()
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Barro","")

plot_Barro_cut_all()
plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum","van Dokkum et al. 2015",alpha0=0)
plot_dexter_sample_all()
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vDokkum", "")

plot_vDokkum_cut_all()
plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        

SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", "van der Wel et al. 2014",alpha0=0)
plot_dexter_sample_all()
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("vdWel", "")

plot_vdWel_cut_all()
plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        

SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Graham", "Graham et al. 2015",alpha0=0)
plot_dexter_sample_all()
SPlot.SelectionCut(mass0,Dist0).plot_cut_specific("Graham", "")

plot_Graham_cut_all()
plt.ylim(ylim[0],ylim[1])
plt.xlim(xlim[0],xlim[1])
plt.show()


fig, ax = plt.subplots()        
plt.plot(E1_R15BC,mag_g1-mag_i1,'o')
plt.plot(E2_R15BC,mag_g2-mag_i2,'o')
plt.plot(E3_R15BC,mag_g3-mag_i3,'o')

plt.show()


def six_plot():
    
    return None
