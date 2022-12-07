#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:55:51 2021

@author: dexter
"""

import numpy as np
import pandas as pd

loc = "/home/dexter/Downloads/S4G_Mass_Comp/Gal_vdis_dict_equvi_mass2.csv"

df=pd.read_csv(loc)

median = 11
#df=df[df.name != 'NGC3976'] # Hon  Similar distance, magnitude seems overestimated
#df=df[df.name != 'NGC2968'] # Hon  15 Mpc higher distance

df=df[df.name != 'NGC5354'] # S4G Dist not avail
df=df[df.name != 'NGC4795'] # S4G Dist not avail
df=df[df.name != 'NGC5311'] # S4G Dist not avail

#df=df[df.name != 'NGC4321'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4461'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4596'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4725'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4754'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4795'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4826'] # 2component decomposition based Mgal not provided
#df=df[df.name != 'NGC4866'] # 2component decomposition based Mgal not provided

df=df[df.S4G_in_out == 'i']

labels=df['name']  #Hon

Hon_gal_mass=df['log_Gal_mass_Bell']  #Dexter
Hon_gal_mass_2cpt=df['BD_log_Gal_mass_Bell'] #Mgal from 2 component (B+D) decomp
Hon_gal_mass_taylorML=df['Gal_mass_Taylor_multi_cpt'] #using Taylor et al. 2011 M/L ratio on multi-cpt mag

S4G_gal_mass=df['Confirm_S4G_Mgal'] #Dexter

#Hon_ML=10**df['Log_ML_bell']  #avg 2.7

Hon_dist = df['Dist'] 
S4G_dist =df['S4G_Dist']

S4G_abs_mag=df['S4G_3.6mu_abs_mag_gal']
S4G_app_mag= S4G_abs_mag+ 5*np.log10(S4G_dist)+25  #app mag
Log_ML_S4G=S4G_gal_mass - 0.4*(6.02-S4G_abs_mag)

S4G_abs_mag_mod= S4G_app_mag- (5*np.log10(Hon_dist)+25 )

S4G_gal_mass_mod= Log_ML_S4G+0.4*(6.02-S4G_abs_mag_mod)

#print(np.mean(10**Log_ML_S4G))

DexC_Taylor_mass_mulit_cpt = df['DexC_Taylor_mass_mulit_cpt']
DexC_Roediger_mass_multi_cpt = df['DexC_Roediger_mass_multi_cpt']
DexC_Zibetti_mass_multi_cpt = df['DexC_Zibetti_mass_multi_cpt']
DexC_Into_mass_multi_cpt =  df['DexC_Into_mass_multi_cpt']

xdata=Hon_gal_mass_taylorML
ydata=S4G_gal_mass_mod

xdata1 = DexC_Taylor_mass_mulit_cpt
xdata2 = DexC_Roediger_mass_multi_cpt
xdata3 = DexC_Zibetti_mass_multi_cpt
xdata4 = DexC_Into_mass_multi_cpt
# fitting

xlim_l,xlim_u = 9.4,12.4
ylim_l,ylim_u = 9.4, 12.4

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "Times New Roman"


def compare_mag_plot(x,y):
    fig = plt.figure()

    gs = gridspec.GridSpec(1,1) 
    ax = plt.subplot(gs[0])
    
    ax.plot(np.log10(xdata1), ydata,'o',label="Taylor et al. 2011")
    ax.plot(np.log10(xdata2), ydata,'o',label="Roediger et al. 2015")
    ax.plot(np.log10(xdata3), ydata,'o',label="Zibetti et al. 2009")
    ax.plot(np.log10(xdata4), ydata,'o',label="Into et al.2013")
    
    #ax.plot(x,y,'o',ms=12)
    ax.plot([0,10,20],[0,10,20],'k')
    ax.set_xlim(xlim_l,xlim_u)
    ax.set_ylim(ylim_l,ylim_u)
    
    ax.set_ylabel(r"$M_*/\rm M_\odot~(S4G)$",fontsize=14)
    ax.set_xlabel(r"$M_*/\rm M_\odot~(Hon~et~al.)$",fontsize=14)
    return None

compare_mag_plot(xdata, ydata)
plt.legend()
plt.show()
