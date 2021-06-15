#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 16:54:05 2021

@author: dexter

This script is dedicated to 
analyzing the surface brightness (SB) profile of galaxies

This script will produce the following data:
1) The various radius for the spheroids

This script also plot:
1) The illustration of percentage covered by the calculation


"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAna
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

import matplotlib.gridspec as gridspec

from scipy.interpolate import interp1d

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

#scan, trim, info

NGC4772SB ='/home/dexter/result/NGC4772/out1_NGC4772.dat'

NGC3646SB = '/home/dexter/result/NGC3646/out1_NGC3646.dat'

NGC5382SB = '/home/dexter/result/NGC5382/out1_NGC5382.dat'


pix = SRead.read_table(NGC4772SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R = SRead.read_table(NGC4772SB)[:,0]
R_s = SRead.read_table(NGC4772SB)[:,0]*0.4
e = SRead.read_table(NGC4772SB)[:,5]
Rmax= max(R)

R_e = SAna.Isophote.circularized(R_s, e)

Rmax_e = SAna.Isophote.circularized(Rmax,e[-1])

mu = SAna.Isophote.pix_val_to_mu(pix)

f2 = interp1d(R, mu, kind='cubic')
xnew = np.linspace(0, max(R), num=800, endpoint=True)
#print("totalmag",cal_SB_mag(R,pix,Rmax,e)) #10.3

#print("NGC4772 totalmag",SAna.Isophote.cal_SB_mag(R,pix,112/0.4,e),
#      'right answer',10.3) #10.3 #no conv 10.299
#print("NGC4772 totalmag",SAna.Isophote.cal_SB_mag(xnew,f2(xnew),max(R),e),
#      'right answer',10.3) #10.3 #no conv 10.299

pix_1 = SRead.read_table(NGC3646SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R_1 = SRead.read_table(NGC3646SB)[:,0]
R_s_1 = SRead.read_table(NGC3646SB)[:,0]*0.4
e_1 = SRead.read_table(NGC3646SB)[:,5]
Rmax_1= max(R_1)

R_e_1 = SAna.Isophote.circularized(R_s_1, e_1)
Rmax_e_1 = SAna.Isophote.circularized(Rmax_1,e_1[-1])
mu_1 = SAna.Isophote.pix_val_to_mu(pix_1)

ef_1 = interp1d(R_1, e_1, kind='linear')
f2_1 = interp1d(R_1, mu_1, kind='cubic')
xnew_1 = np.linspace(0, max(R_1), num=5*len(R_1), endpoint=True)
    
fig = plt.figure()
plt.plot(R_1, e_1, 'o', xnew_1, np.nan_to_num(ef_1(xnew_1)), '-')
plt.legend(['data','cubic'], loc='best')
plt.show()

#print(np.nan_to_num(ef_1(xnew_1)))
#print(len(xnew_1),len(f2_1(xnew_1)),len(ef_1(xnew_1)))
#29
#print("NGC3646 totalmag",SAna.Isophote.cal_SB_mag(R_1,pix_1,29/0.4,e_1),
#      'right answer',11.07) #11.07 # no cov
#print("NGC3646 totalmag",SAna.Isophote.cal_SB_mag(xnew_1,f2_1(xnew_1),38.0/0.4,np.nan_to_num(ef_1(xnew_1)),step = 0.1),
#      'right answer', 11.07) #11.54 # no cov 11.4695

pix_2 = SRead.read_table(NGC5382SB)[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
R_2 = SRead.read_table(NGC5382SB)[:,0]
R_s_2 = SRead.read_table(NGC5382SB)[:,0]*0.4
e_2 = SRead.read_table(NGC5382SB)[:,5]
Rmax_2 = 46

R_e_2 = SAna.Isophote.circularized(R_s_2, e_2)
Rmax_e_2 = SAna.Isophote.circularized(Rmax_2,e_2[-1])
mu_2 = SAna.Isophote.pix_val_to_mu(pix_2)

ef_2 = interp1d(R_2, e_2, kind='linear')
f2_2 = interp1d(R_2, mu_2, kind='cubic')
xnew_2 = np.linspace(0, max(R_2), num=5*len(R_2), endpoint=True)

#42
#print("NGC5382 totalmag",SAna.Isophote.cal_SB_mag(R_2,pix_2,42.0/0.4,e_2),
#      'right answer', 11.54) #11.54 # no cov 11.4695

#print("NGC5382 totalmag",SAna.Isophote.cal_SB_mag(xnew_2,f2_2(xnew_2),42.0/0.4,np.nan_to_num(ef_2(xnew_2)),step=0.1),
#      'right answer', 11.54) #11.54 # no cov 11.4695
##
x = R_2 #np.linspace(0, 10, num=11, endpoint=True)
y = mu_2 #np.cos(-x**2/9.0)
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')

#print(len(x))

xnew = np.linspace(0, max(x), num=800, endpoint=True)

fig = plt.figure()

#print(len(xnew))
#omg num=800 is right
plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.gca().invert_yaxis()

plt.show()


nxt = 2
extsize = nxt*len(x)
xx = np.array(np.linspace(x[0],1.*nxt*x[len(x)-1],extsize))
#print(len(xx))


##

def list_centre():
    outlist = SRead.read_table("/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')
    geom_file = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat",
                               dtype='str')
    
    name = geom_file[:,0]
    
    for i in range(len(name)):
        outfile_name  = outlist[i]
        outfile = SRead.read_table(outfile_name)
    
        X0, Y0 = outfile[:,9], outfile[:,11]
        
        print(name[i],X0[0],Y0[0])
        
def list_prof_mag():
    geom_file = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat",
                               dtype='str')    
    
    name = geom_file[:,0]
    total_mag = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi")
    
    for i in range(len(name)):
        print(name[i],total_mag[i])
        
def list_mu0():
    outlist = SRead.read_table("/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')
    geom_file = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat")
    geom_file_n = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')

    name = geom_file_n[:,0]

    for i in range(len(name)):
        outfile_name  = outlist[i]
        outfile = SRead.read_table(outfile_name) 
   
        pix = outfile[:,1][0]
        mu0 = SAna.Isophote.pix_val_to_mu(pix)
        print(name[i],mu0)

#list_centre()
#list_prof_mag()
list_mu0()

def run_cal_mag():
    outlist = SRead.read_table("/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')
    geom_file = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat")
    geom_file_n = SRead.read_table("/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')

    name = geom_file_n[:,0]
    total_mag_mine_bag = []
    
    for i in range(len(name)):
        outfile_name  = outlist[i]
        outfile = SRead.read_table(outfile_name) 
   
        pix = outfile[:,1]#/SRead.read_table(NGC4772SB)[:,0]**2
        R = outfile[:,0]
        R_s = outfile[:,0]*0.4
        e = outfile[:,5]

        mu = SAna.Isophote.pix_val_to_mu(pix) 
        
        Rmax = geom_file[:,2][i]/ 0.4
        
        
        total_mag_prof = geom_file[:,5]
        R_e = SAna.Isophote.circularized(R_s, e)
        #print(len(R),len(mu),len(Rmax),len(e))
        total_mag_mine = SAna.Isophote.cal_SB_mag(R,pix,Rmax,e)
        
        print(name[i], total_mag_prof[i],total_mag_mine)
        
        total_mag_mine_bag.append(total_mag_mine)
    total_mag_mine_bag = np.array(total_mag_mine_bag)
    return name, total_mag_prof, total_mag_mine_bag
    

#A = run_cal_mag()


#SPlot.ShowcaseCompare2.plot_compare_generic(A[1], A[2], para_name="Mag", name=A[0], label=["prof","mine"])

#total_mag_BD = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_BD_equvi")

#SPlot.ShowcaseCompare2.plot_compare_generic(A[1], total_mag_BD, para_name="Mag", name=A[0], label=["multi","BD"])


#		n_ex = 7#
#		lx_ex = n_ex*len(x)#
#		x_ex = array(linspace(x[0],1.*n_ex*x[len(x)-1],4*lx_ex))
#		fit_nice_extended=construct_model.residual(fit_parameters, x_ex, n_all, sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,id_entry,td_entry,PSF_ext,fit_min,fit_max,ellipticity)
#		total_galaxy_magnitude=integrator.magnitude(x_ex,fit_nice_extended,mzp)
