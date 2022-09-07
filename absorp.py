#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 09:44:41 2022

@author: rishat
"""

import numpy as np
import os
import subprocess
from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
import lmfit as lm
import re

"""
read the absorption files. 
"""

ls = ["ls","*.dat"] 
ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
list_ = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
match = "absorptionCS_w_point_charges"


list_ = [i for i in list_ if bool(re.match(match, i)) is True]
f = open(list_[1])
lines = f.readlines()
line_ = lines[-1].split()

def read_absorp(dimention):
    
    t_E_C = np.zeros(dimention)
    for i in range(len(list_[:])):
        f = open(list_[i])
        lines = f.readlines()
        lines = lines[4:]
        for j in range(len(lines)):
            t_E_C[i, j, :] = lines[j].split()
    
    ene = t_E_C[:, :, 1::2]
    coe = t_E_C[:, :, 0::2]
    time = coe[:, :, 0]
    coe = coe[:, :, 1:]

    return {"time": time, "E": ene, "coe": coe}
    

"""
Gaussian, Lorentzian and convolution
"""

def gauss_lorentz(sigma_E, gamma_E, tEcoe, tick_size): 
    
    tick_E = np.linspace(min(tEcoe["E"].flatten()), max(tEcoe["E"].flatten()), tick_size)
    
    gauss = np.zeros((np.shape(tEcoe["E"][0, :])[0], tick_size))
    lorentz = np.zeros((np.shape(tEcoe["E"][0, :])[0], tick_size))
    gauss_lorentz = np.zeros((np.shape(tEcoe["E"][0, :])[0], tick_size))
    
    for i in range(np.shape(tEcoe["E"][0, :])[0]):
        for j in range(len(tick_E)):
            gauss[i, j] = sum(lm.models.gaussian(tEcoe["E"][:, i].flatten(), 
                                                tEcoe["coe"][:, i].flatten(), tick_E[j], sigma_E))
            
            lorentz[i, j] = sum(lm.models.lorentzian(tEcoe["E"][:, i].flatten(), 
                                                    tEcoe["coe"][:, i].flatten(), tick_E[j], gamma_E))
            
            gauss_lorentz[i , j] = sum(np.convolve(gauss[i, j], lorentz[i, j] , mode='same'))
    
    return {"tick_E" : tick_E, "gaussian" : gauss, "lorentzian" : lorentz, "gauss_lorentz_convol" : gauss_lorentz}
     
"""
Plot with Gaussian, Lorentz and Gaussian+Lorentz 
"""
    
    
def plot_gauss_lorentz_convol(tEcoe ,t_E_coe):

    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    fig, axs = plt.subplots(1,1)
    # axs[0].plot(t_E_coe["tick_E"], t_E_coe["gaussian"][0, :], "--", label = 'Gaussian')
    # axs[0].plot(t_E_coe["tick_E"], t_E_coe["lorentzian"][0, :], label = 'Lorentzian')
    # axs[0].plot(t_E_coe["tick_E"], t_E_coe["gauss_lorentz_convol"][0, :], label = 'Convolution')
    # axs[0].set_title('$\u03C3=0.1ev$')
    # axs[0].legend()

    print(np.trapz(t_E_coe["gaussian"][0, :], t_E_coe["tick_E"]))
    print(np.trapz(t_E_coe["lorentzian"][0, :], t_E_coe["tick_E"]))
    print(np.trapz(t_E_coe["gauss_lorentz_convol"][0, :], t_E_coe["tick_E"]))
    
    """
    plot contourf and save data
    """
    
    tick_E = np.array([t_E_coe["tick_E"] for i in range(len(t_E_coe["gaussian"][:, 0]))])
    tick_time = []
    for i in tEcoe["time"][0]:
        tick_ = []
        for j in range(len(t_E_coe["gaussian"][0, :])):
            tick_.append(i)
        tick_time.append(tick_)
        
    gau = axs.contourf(tick_E, tick_time, t_E_coe["gaussian"], 70, cmap='Blues')
    axs.tick_params(axis='both', which='major', labelsize=20)
    axs.set_xlabel('Energy (Ev)', fontsize = 20)
    axs.set_ylabel('Time (fs)', fontsize = 20)
    axs.set_title('$\u03C3=0.1ev$')
    cbar = plt.colorbar(gau, ax = axs)
    cbar.set_label("Relative Photon Intensity", labelpad=1, size=14)
    
    # gau = axs[0].contourf(tick_E, tick_time, t_E_coe["gaussian"], 70, cmap='Blues')
    # axs[0].tick_params(axis='both', which='major', labelsize=20)
    # axs[0].set_xlabel('Energy (Ev)', fontsize = 20)
    # axs[0].set_ylabel('Time (fs)', fontsize = 20)
    # axs[0].set_title('$\u03C3=0.1ev$')
    # cbar = plt.colorbar(gau, ax = axs[0])
    # cbar.set_label("Relative Photon Intensity", labelpad=1, size=14)
    
    
    # lo = ax2.contourf(tick_E, tick_time, t_E_coe["lorentzian"], 70, cmap='Blues')
    # ax2.tick_params(axis='both', which='major', labelsize=20)
    # ax2.set_xlabel('Energy (Ev)', fontsize=20)
    # ax2.set_ylabel('Time (fs)', fontsize=20)
    # cbar = plt.colorbar(lo,ax=ax2)
    # cbar.set_label("Relative Photon Intensity", labelpad=-1, size=14)

    # conv=ax3.contourf(tick_E, tick_time, t_E_coe["gauss_lorentz_convol"], 70, cmap='Blues')
    # ax3.tick_params(axis='both', which='major', labelsize=20)
    # ax3.set_xlabel('Energy (Ev)', fontsize=20)
    # ax3.set_ylabel('Time (fs)', fontsize=20)
    # plt.colorbar(conv,ax=ax3)

    # ax3.set_title('$\u03C3=0.1ev$')
    # ax3.legend()
    plt.show()
    
    
    return 
    

def main():
    
    tEcoe = read_absorp((len(list_), len(lines[4:]), len(line_)))
    sigma_E = gamma_E = 0.1          #FWHM ev
    tick_size = 400        #Eenergy tick size
    t_E_coe = gauss_lorentz(sigma_E, gamma_E, tEcoe, tick_size)
    plot_gauss_lorentz_convol(tEcoe, t_E_coe)
         
    return  

if __name__ == "__main__":
    main()

    
