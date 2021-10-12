#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 01:11:25 2020

@author: dexter

The module for Sph Project analysis.


"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import scipy as scipy
from scipy import stats

from astropy.table import Table, Column, MaskedColumn
import matplotlib.patches as mpatches

import matplotlib.pyplot as plt

from scipy.special import gammainc, k1
from math import cosh

import SphRead as SRead
import SphPlot as SPlot

__all__ = ['AnalyticFunctions','Isophote']

__author__="Dexter S.H. Hon"

#%% utility
def closest(lst, K):
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

#%% The two function from Profiler. For the sake of consistency I do not attempt to write my own functions
def gammadif(a, x):
    out = gammainc(a, x) - 0.5
    return out

def get_bn(n):
    precision = 1e-06
    a = 2.0 * n
    bn = 2.0 * n - 0.333
    if bn < 0:
        bn = 0.0
    d = gammadif(a, bn)
    inc = 1.0
    while abs(d) > precision:
        d = gammadif(a, bn)
        bnplus, bnminus = bn + inc, bn - inc
        if bnminus < 0:
            bnminus = 0.0
        dplus, dminus = gammadif(a, bnplus), gammadif(a, bnminus)
        if abs(dplus) < precision:
            bn = bnplus
            break
        elif abs(dminus) < precision:
            bn = bnminus
            break
        else:
            if abs(dplus) < abs(d):
                bn = bnplus
            elif abs(dminus) < abs(d):
                bn = bnminus
            inc = inc / 2.0
    return bn

#%%
class AnalyticFunctions(object):
    """
    """
    def __init__(self,r=None,theta = None):
        self._r = r
        self.theta = theta
        
    @property
    def r(self):
        return(self._r)
    
    def linear_func_1D(x,A,B):
    
        return A*x+B
        
    def mu_sersic_func(r, mu_e,r_e,n_Ser):
        """
        The Sersic function in surface brightness form.

        Parameters
        ----------
        r : 1D numpy array
            Radius.
        mu_e : float
            Surface brightness at r_e.
        r_e : float
            Effective half-light radius.
        n_Ser : float
            Sersic index.

        """
        b_n = get_bn(n_Ser)
        return (mu_e + 1.0857362*b_n*((r / r_e)**(1.0/n_Ser)-1))
    def simple_mu_core_sersic(r, mu_p, r_e, n_ser, r_b, al, ga):
        
        b_n = get_bn(n_ser)

        mu =  -1.0 * (mu_p - 2.5 * ga / al * np.log10(1.0 + (r_b / r) ** al) + 
                      1.0857362 * b_n * ((r ** al + r_b ** al) / r_e ** al) ** 
                      (1 / (n_ser * al)))
        
        return -1.0 * mu
    
    def mu_core_sersic_func(r, mu_p, r_e, n_ser, r_b, al, ga):
        """
        The core Sersic function in surface brightness form.

        Parameters
        ----------
        r : 1D numpy array
            Radius.
        mu_pr : float
            The surface brightness constant.
        al : float
            DESCRIPTION.
        ga : float
            DESCRIPTION.
        r_b : float
            The breaking radius.
        r_e : float
            The effective radius.
        n_Ser : TYPE
            The Sersic index.
        """
        b_n = get_bn(n_ser)
        if r[0] == 0.0:
            r[0] = 1e-10
        y = -1.0 * (mu_p - 2.5 * ga / al * np.log10(1.0 + (r_b / r) ** al) 
                    + 1.0857362 * b_n * ((r ** al + r_b ** al) / r_e ** al) ** 
                    (1 / (n_ser * al)))
        y[0] = -1.0 * (mu_p - 2.5 * ga / al * 
                       np.log10(1.0 + (r_b / (r[1] / 10.0)) ** al) + 
                       1.0857362 * b_n * (((r[1] / 10.0) ** al + r_b ** al) / 
                                          r_e ** al) ** (1 / (n_ser * al)))
        return y*-1.0
    
    def mu_exp_func(r,mu_0,h):
        """
        The exponential function in surface brightness form.

        Parameters
        ----------
        r : TYPE
            Radius.
        mu_0 : float
            surface brightness at r=0.
        h : float
            Scale length.


        """
        return  (mu_0 + 1.0857362 * (r / h))   
    
    def mu_broken_exp_func(r, mu_0, r_b, h1, h2):
        """
        The broken exponential function in surface brightness form

        Parameters
        ----------
        r : TYPE
            DESCRIPTION.
        mu_0 : TYPE
            DESCRIPTION.
        r_b : TYPE
            DESCRIPTION.
        h1 : TYPE
            DESCRIPTION.
        h2 : TYPE
            DESCRIPTION.

        Returns
        -------
        y : TYPE
            DESCRIPTION.

        """
        mu_b = mu_0 + 1.0857362 * (r_b / h1)
        i = len(r) - 1
        y = mu_0 + 1.0857362 * (r / h1)
        while r[i] > r_b:
            y[i] = mu_b + 1.0857362 * ((r[i] - r_b) / h2)
            i = i - 1
            
        return  y
    
    def mu_incl_exp_func(r, mu_0z, z0, case):
        s = np.array([0.0] * len(r))
        for i in range(0, len(r)):
            if case == 'r0':
                s[i] = -1.0 * (mu_0z - 5.0 * np.log10(1.0 / cosh(r[i] / z0)))
            if case == 'z0':
                s[i] = -1.0 * (mu_0z - 2.5 * np.log10(r[i] / z0 * k1(r[i] / z0)))

        if r[0] == 0.0:
            s[0] = -1.0 * mu_0z
        return s*-1.0
    
    def mu_ferrer_func(r, mu_0f, r_out, alpha_F, beta_F):
        """
        The Ferrer function in surface brightness form

        Parameters
        ----------
        r : 1D numpy array
            Radius.
        mu_0f : float
            Surface brightness at r=0.
        r_out : float
            The maximum radius.
        alpha_F : float
            Coefficient determining the outer slope.
        beta_F : float
            Coefficient determining the inner slope.
        """
        r = np.array(r)
        fprof = [0.0] * len(r)
        for i in range(0, len(r)):
            if r[i] >= r_out:
                fprof[i] = -99.0
            else:
                fprof = -1.0 * (mu_0f - alpha_F * np.log10(1.0 - (r / r_out) ** (2.0 - beta_F)))

        return fprof*-1.0


    def size_mass_powerlaw(M,a,b):
        """
        R_e = a*(M_*/(1e10*M_{\odot}))**b
        
        Parameters
        ----------
        M : 1D numpy array
            Stellar mass.
            
        Returns
        -------
        Effective radius.

        """
        return a*(M/1e10)**b
    
    def size_mass_twopowerlaw(M,a,b,c,Mo,M_sun):
        """
        R_e = c*(M_*/(M_{\odot})*(1+(M/Mo))**(b-a)
        
        Parameters
        ----------
        M : 1D numpy array
            Stellar mass.
            
        Returns
        -------
        Effective radius.

        """
        return c*(M/M_sun)*(1+(M/Mo))**(b-a)
    
    
#%%
from scipy.special import gamma

def b_value(n,x):
   bn=[]
   diff=[]
   b= np.linspace(0.1, 50.0, 50000) # works for most z
#   x=1/z=10 # 10, 2, 10/9    
#   b= np.linspace(0.1, 25.0, 2500)  #Re  2
#   b= np.linspace(0.1, 50.0, 25000)  #R10  10
#   b= np.linspace(0.1, 100.0, 50000)  #R90  10/9
   for j in b:        
        g2= x*gammainc(2*n,j)
        g22=np.round(g2,2)
        if g22==1.00: ## gammainc(2*n,j) is normalized by dividing with gamma(2*n)
            k=np.round(j,4)
            bn.append(k)
            dif=abs(1-g2)
            diff.append(dif)
   
   diff_min= min(diff)
   for s, d in zip(bn, diff):
       if d==diff_min:
           b_final=s
       else:
           continue     
   
   return b_final    
#%%    
class PlotProfile(AnalyticFunctions):
    
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
#%%
class Fitting():
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
#class 2Dplot(object)

#%%

from scipy.interpolate import interp1d

class Isophote(object):
    
    def pix_val_to_mu(pix_val,zp=22.5,scale=0.4):
        """
        Convert pixel value to surface brightness density

        Parameters
        ----------
        pix_val : 1D numpy array
            The pixel array.
        zp: float
            zero-point magnitude, default: 22.5 (SDSS optical standard)
        scale: float
            arcsec/pix scale, default: 0.4 (SDSS standard) 
            
        Returns
        -------
        1D array in mag arcsec-2 scale.

        """
        return zp-2.5*np.log10(pix_val/scale**2)

    def circularized(R_m,e):
        """
        Convert major axis radius into circularised radius 

        Parameters
        ----------
        R_m : 1D numpy array
            An array of radius in major axis.
        e : 1D numpy array
            An arrray of ellipicity to each radius.
            Elements need to be < 1.0

        Returns
        -------
        1D array
            circularised radius.

        """
        return R_m*np.sqrt(1-e)


    def cal_SB_mag(input_SB_x,input_SB_y,Rmax,e,step=0.1,scale=0.4):
        """
        Cacluate the SB magnitude in total numerically.
        (the result is slightly different than profiler. I suspect it is because
         Bogdan use "comp" to calculate the total mag, with different step size, 
         0.1 or so)    

        This function only operate on equ axis
        
        Parameters
        ----------
        input_SB_x : 1D list, np array
            Radius, pixel or arcsec.
        input_SB_y : 1D list, np array
            The brightness, flux, intensity.
        Rmax : float
            The maximum radius.
        e : 1D numpy array
            The ellipicity at each radius
        step: float
            The interval between each radial data point, default: 0.5 (pix)
        scale: float
            The arcsec/pixel scale, default: 0.4
        
        Returns
        -------
        float
            The total magnitude.

        """
        # Turn the input from major to equvi axis
        input_SB_x = np.nan_to_num(Isophote.circularized(input_SB_x,e))
        
        #Find the closest value to Rmax
        A = closest(input_SB_x,Rmax) 
        
        # profiler weird operation
        #nxt = 2.0
        #extsize = nxt*len(input_SB_x)
        #xx = np.array(np.linspace(input_SB_x[0],nxt*input_SB_x[len(input_SB_x)-1],int(extsize)))
        #print(xx,input_SB_x)
        
     		#total_galaxy_magnitude=integrator.magnitude(x_ex,fit_nice_extended,mzp)
        
        #yy = interp1d(input_SB_x, input_SB_y, kind='linear')
        #ee = interp1d(input_SB_x, e, kind='linear')

        #yy_o = yy(xx)
        #ee_o = ee(xx)
        
        # find the index of maximum radius
        index= np.where(input_SB_x==A)[0][:]
        #print(input_SB_y[index],int(index), Rmax*scale)
    
        dx = step #0.5 is the sma distance, 0.4 is the pixel to arcsec ratio

        #print(len(input_SB_x),len(input_SB_y))
        #for i in range(1,int(index)): #int temp
        #
        #    if input_SB_y[i] > -50.0:
        #        
        #        #dl = 2*np.pi*(1-e[i])*input_SB_x[i]*10**(input_SB_y[i]/-2.5)*dx
        ##        dl = 2*np.pi*(1-ee_o[i])*xx[i]*10**(yy_o[i]/-2.5)*dx
        #        #dl = 2*np.pi*input_SB_x[i]*10**(input_SB_y[i]/-2.5)*dx
        #        l = l +dl
        #        
        #        total_mag = -2.5*np.log10(l*scale*scale)
        #        #total_mag = input_SB_y[0:index] +2.5*np.log10(2*input_SB_x[0:index]*np.pi*0.5)
        
        
        r = input_SB_x*0.4 # x-input in arcsec
        p = input_SB_y     # y-input in SB
        # operate integration on r
        rr = r
        

        #n_ex = 7
        #lx_ex = n_ex *len(r)
        #x_ex = np.array(np.linspace(r[0],1.*n_ex*r[len(r)-1],4*lx_ex))  
        #
        #rr = x_ex
        
        # Interpolation
        pp = interp1d(r, p, kind='linear')
        ee = interp1d(r, e, kind='linear')
        
        pp_o, ee_o = pp(rr), ee(rr)
        
        # profiler code
        mzp = 22.5
        # Turn pixel value into SB
        mu = Isophote.pix_val_to_mu(pp_o)
        mu = mu - mzp
        lum = [0.0] * len(mu)
        lum = 10.0 ** (-mu / 2.5)
        dl, t = 0.0, 0.0
        for i in range(1, int(index)-1):#len(r) - 1):
            s = np.pi * (rr[i] ** 2 - rr[i-1] ** 2)
            l = 0.5 * (lum[i] + lum[i-1])
            dl = l * s
            t = t + dl

        total_mag = -2.5 * np.log10(t) + mzp
        #ax = plt.gca()
        #ax.plot(input_SB_x*scale,input_SB_y,'o')
        #ax.vlines(Rmax*scale,  max(input_SB_y), min(input_SB_y))
        #ax.fill(list(R_e)+[0,R_e[0],R_e[0]],
        #        list(mu)+[min(mu),mu[index0],mu[0]])

        #ax.invert_yaxis()
        #plt.show()

        return total_mag
    
    #WIP
    def scan_basic1D(input_array_x, input_array_y, e, Rmax,
                 percentage=1.0, start_point=0.5):
        """
        A function to scan the SB profile, and find the radius that matches the 
        percentage value.
        i.e. percentage =0.5, this function return the half-light radius

        Parameters
        ----------
        input_SB_x : 1D list, np array
            Radius, pixel or arcsec.
        input_SB_y : 1D list, np array
            The brightness, flux, intensity, or surface brightness.
        Rmax : float
            The maximum radius.
        e : 1D numpy array
            The ellipicity at each radius
        percentage : TYPE, optional
            DESCRIPTION. The default is 1.0.
        start_point : TYPE, optional
            DESCRIPTION. The default is 0.5.

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """
        # Find max
        # start with start_point, find area, 
        return None #x,y coordinate for the percentage

    # WIP
    def radius_by_percentage(x,y,e,total_mag,fraction=1.0,
                             ):
        """
        A function to calculate the radius of the galaxy at some percentage of 
        light. 
        If one 
        
        Parameters
        ----------
        input_file: str
            The name of the output file of ISOFIT.
        max_radius: float
            The maximum radius of the galaxy, given in arcsec.
        fraction: float
            The fraction of light, from the centre extending outward.
            The value should be between 0 to 1.
            
        Returns
        -------
        perc_radius: float
            The radius of the galaxy in which the amount of light matches the 
            assigned percentage parameter.

        """
        #Read inputfile, b/a e and so on
        #table, table_n = SRead.read_list(input_file), SRead.read_list(input_file,dtype ="str")
        
        #x = table[:,0]*0.4 # sma*pixscale
        #y = table[:,1] # intens
        #e = table[:,5] # ellipticity   
        
        # Turn intens pixel value to SB density mu
        yy = Isophote.pix_val_to_mu(y)
        yy = y/0.4**2 
        # Turn intens pixel value into Luminosity
        
        #calculate the total light covered within max_radius
        #if total_mag == "Auto":
        #    total_mag = Isophote.cal_SB_mag(x, y, max_radius, e)
        #elif total_mag > 1:
        #    total_mag = total_mag
            
        # Set the target magnitude
        target_lum = 10**(total_mag/-2.5) * fraction
        #loop from the centre, and numerically determine the radius in which 
        Lum = []
        #calculate the cumulative luminosity at each radius
        for i in range(2,len(x)):
            Lum.append(0.5*(yy[i]-yy[i-1])*(x[i]**2-x[i-1]**2)*np.pi)
            #print(i,Lum[i],target_lum)
            print(i, yy[i],x[i])
            print(0.5*(yy[i]-yy[i-1])*(x[i]**2-x[i-1]**2)*np.pi)
            #x percentage of light are included (stopping condition)
            if sum(Lum[0:i]) > target_lum and sum(Lum[0:i-1]) < target_lum:
                R_traget = x[i]
                break
        print("target_lum",target_lum)

        return R_traget
