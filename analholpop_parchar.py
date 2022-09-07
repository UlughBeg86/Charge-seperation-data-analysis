#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 08:43:50 2022

@author: rishat
"""

import numpy as np
import os
import subprocess
from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
import re
import warnings


"""
Please choose which function you want to use and update match variable accordingly
"""
functions = {"occup_state_index" : "S.log", "occupied_orbital_E" :"occupied_orbital_energies_w_point_charges", 
             "t_KE_PE_KEPE" : "E.log", "Ad_V.log" : "adiabatic_potential_energy_surface", 
             "hole_density" : "hole_populations_w_point_charges", "charge_center" : "charge_center_w_point_charges", 
             "hole_center" : "hole_center_w_point_charges",
             "partial_charges" : "partial_charges_w_point_charges", "probability_find_hole_given_state" : "P.log",
             "velocity" : "V_"}

ls = ["ls","*.log"] 
ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
match = functions["velocity"]
list_ = [i for i in list_ if bool(re.match(match[0], i)) is True]

ls = ["ls","*.dat"] 
ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
list1_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
match1 = functions["charge_center"]
list1_ = [i for i in list1_ if bool(re.match(match1[0], i)) is True]

"""
check what is given time fs 
"""
t_f = 200

"""
Read S.log file
"""

def readlog(list_, t_f):
    
    length = [];t = np.linspace(0, t_f, (t_f-0)*2)
    for file in list_:
        f = open(file)
        lines = f.readlines()
        length.append(len(lines))
    if len(t) > min(length):
        warnings.warn("hi, please check, there are short time running files!!")
    if match == "charge_center_w_point_charges" or match == "hole_center_w_point_charges":
        dtype = object
    else:
        dtype = float
        
    time_state = np.zeros((len(list_), min(length), len(lines[0].split())), dtype = dtype)
    for i in range(len(list_)):
        f = open(list_[i])
        lines = f.readlines()
        for j in range(min(length)):
            time_state[i][j] = lines[j].split() 
            # print(list_[i], lines)

    return {match : time_state}


def readlog1(list1_, t_f):
    
    length = [];t = np.linspace(0, t_f, (t_f-0)*2)
    for file in list1_:
        f = open(file)
        lines = f.readlines()
        length.append(len(lines))
    if len(t) > min(length):
        warnings.warn("hi, please check, there are short time running files!!")
    if match == "charge_center_w_point_charges" or match == "hole_center_w_point_charges":
        dtype = object
    else:
        dtype = float
        
    time_state = np.zeros((len(list1_), min(length), len(lines[0].split())), dtype = dtype)
    for i in range(len(list1_)):
        f = open(list1_[i])
        lines = f.readlines()
        for j in range(min(length)):
            time_state[i][j] = lines[j].split() 
            # print(list1_[i], lines)

    return {match1 : time_state}


def plotSlog(time_state, t_f):

    #atomic unit to fs
    time_state[match][:, :, 0] = 0.02418884254 * time_state[match][:, :, 0]
 
    """
    plot State index VS Time
    """
    
    plt.figure()
    for i in range(len(time_state[match])):
        plt.plot(time_state[match][i, :, 0], time_state[match][i, :, 1], "--")
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlabel("Time (fs)", size=20)
    plt.ylabel("Occupational state index (int)", size=20)
    plt.title("300MO_QM13", size=20)
    plt.show()

    # """
    # Hist VS Time
    # """  
    # plt.figure()
    # for i in range(len(time_state[match])):
        # plt.plot(time_state[match][i, :, 0], time_state[match][i, :, 1])
    # plt.xticks(size=20)
    # plt.yticks(size=20)
    # plt.xlabel("Time (fs)",size=20)
    # plt.ylabel("Occupational state index (int)",size=20)
    # plt.title("200MO_QM13")
    # plt.show()
    
    index = np.asarray([[51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39],
                        [38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26],
                        [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13],
                        [12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0]], dtype = int)

    # indx = np.asarray(["homo-3", "homo-2", "homo-1", "homo"], dtype = str)
    indx = np.asarray(["States 51-39", "States 38-26", "States 25-13", "States 12-0"], dtype = str)
    
    num_hist = []; sum_hist = []; sSum_hist = []
    
    """
    This (range(np.shape(time_stat['occ_stat'])[1])) will give you occupational 
    state at each time step. You may think it is a clock
    range(len(index)) will do the scan each time step, occupational state will 
    account according to given index state number    
    """
    
    for i in range(np.shape(time_state[match])[1]):
        tmp_=[]; eachorb=[]
        for j in range(len(index)):
            _tmp=[];
            for k in index[j]:
                _tmp.append(np.count_nonzero(time_state[match][:, i, 1] == k))
            tmp_.append(_tmp)
            eachorb.append(sum(_tmp))
        num_hist.append(tmp_)
        sum_hist.append(eachorb)
        sSum_hist.append(sum(eachorb))
    
    sum_hist = np.asarray(sum_hist, dtype = int)
           
    """
    1D plot occuVStime
    """
    
    plt.figure()
    
    for i in range(len(np.transpose(sum_hist))):
        plt.plot(time_state[match][0, :, 0], np.transpose(sum_hist)[i]/len(time_state[match]), "--*")
    plt.xlim([0, t_f])
    plt.ylim([-0.01, 1.01])
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend(indx, loc = "upper right", fontsize = 20)
    plt.ylabel('Population', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)
    
    """
    2D contourf plot
    """
    
    
    # Valence = np.asarray([index[:,0] for i in range(len(sum_hist))], dtype = int)
    
    # plt.figure()
    # plt.contourf(Valence, np.transpose(time_stat['time'])[:,:len(index[:,0])], sum_hist/51 , cmap = 'Blues')
    # plt.title("300MO_QM13", size = 20)
    # plt.tick_params(axis='both', which='major', labelsize=20)
    # plt.xlabel('Occupational state index', fontsize=20)
    # plt.ylabel('Time (fs)', fontsize=20)
    # cbar=plt.colorbar()
    # cbar.set_label("Frequency (number of time visit)", labelpad=0.5, size=20)
          
    return 


def plotKEPE_(time_state, t_f):

    #atomic unit to fs
    time_state[match][:, :, 0] = 0.02418884254 * time_state[match][:, :, 0]
    fig, ax1 = plt.subplots()
    
    ax1.plot(time_state[match][0,:,0],  time_state[match][0,:,1])
    ax1.set_ylabel('Kinetic Energy', color='C0')
    ax1.set_xlabel('Time (fs)', color='C0')
    ax1.tick_params(axis='y', color='C0', labelcolor='C0')
    
    ax1.set_title('KE, PE, Tot_E')
    
    ax2 = ax1.twinx()
    ax2.plot(time_state[match][0,:,0], time_state[match][0,:,2], 'C1')
    ax2.plot(time_state[match][0,:,0], time_state[match][0,:,3], 'C2')
    ax2.set_ylabel('PE', color='C1')
    ax2.tick_params(axis='y', color='C1', labelcolor='C1')
    ax2.spines['right'].set_color('C1')
    ax2.spines['left'].set_color('C0')
    
    plt.show()
    
    return

def plotAD_V_Orbital_E(time_state, t_f):
 
    val_start = 14
    
    plt.hist(time_state[match][:, :, val_start:].flatten(), width = 2, label = "width=2 ev")
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel('Orbital energy (ev)', fontsize=20)
    plt.ylabel('Count', fontsize=20)
        
    return

def hole_pop(time_state, t_f):
    
    atom_sites = np.array([np.arange(len(time_state[match][0, 0, 1:-1])) 
                          for i in range(len(time_state[match][0, :, 0]))])
    time_array = []
    for i in time_state[match][0, :, 0]:
        tmp = []
        for j in range(np.shape(atom_sites)[1]):
            tmp.append(i)
        time_array.append(tmp)
    
    # for each file you may plot     
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, time_state[match][0, :, 1:-1], cmap = 'Blues')
    plt.title("300MO_QM13", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("Hole density", labelpad=0.5, size=20)
    
    # Average of all fileslist_
    ottur_hole = np.zeros(np.shape(time_array))
    for i in range(np.shape(time_array)[0]):
        for j in range(1, np.shape(time_array)[1]):
            ottur_hole[i, j] = np.mean(time_state[match][:, i, j])
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, ottur_hole, cmap = 'Blues')
    plt.title("300MO_QM13_Averaged", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("Hole density", labelpad=0.5, size=20)
    
    return


def partial_charge(time_state, t_f):
    # time_state[match][0, :, 0]
    
    # plt.plot(time_state[match][0, :, 0], time_state[match][0, :, 1])
    
    atom_sites = np.array([np.arange(len(time_state[match][0, 0, 1:-1])) 
                          for i in range(len(time_state[match][0, :, 0]))])
    time_array = []
    for i in time_state[match][0, :, 0]:
        tmp = []
        for j in range(np.shape(atom_sites)[1]):
            tmp.append(i)
        time_array.append(tmp)
    
    # for each file you may plot     
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, time_state[match][0, :, 1:-1], cmap = 'Blues')
    plt.title("300MO_QM13", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("Partial charge density", labelpad=0.5, size=20)
    
    # Average of all files
    ottur_hole = np.zeros(np.shape(time_array))
    for i in range(np.shape(time_array)[0]):
        for j in range(1, np.shape(time_array)[1]):
            ottur_hole[i, j] = np.mean(time_state[match][:, i, j])
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, ottur_hole, cmap = 'Blues')
    plt.title("300MO_QM13_Averaged", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("Partial charge density", labelpad=0.5, size=20)
    
    return

def Plog(time_state, t_f):
    time_state[match][:, :, 0] = 0.02418884254 * time_state[match][:, :, 0]
    
    atom_sites = np.array([np.arange(len(time_state[match][0, 0, 1:-1])) 
                          for i in range(len(time_state[match][0, :, 0]))])
    time_array = []
    for i in time_state[match][0, :, 0]:
        tmp = []
        for j in range(np.shape(atom_sites)[1]):
            tmp.append(i)
        time_array.append(tmp)
    
    # for each file you may plot     
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, time_state[match][0, :, 1:-1], cmap = 'Blues')
    plt.title("300MO_QM13", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("probability to find hole at given state", labelpad=0.5, size=20)
    
    # Average of all files
    ottur_hole = np.zeros(np.shape(time_array))
    for i in range(np.shape(time_array)[0]):
        for j in range(1, np.shape(time_array)[1]):
            ottur_hole[i, j] = np.mean(time_state[match][:, i, j])
    
    plt.figure()  
    plt.contourf(atom_sites, time_array, ottur_hole, cmap = 'Blues')
    plt.title("300MO_QM13_Averaged", size = 20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.xlabel("Atomic site", fontsize=20)
    plt.ylabel('Time (fs)', fontsize=20)
    cbar=plt.colorbar()
    cbar.set_label("probability to find hole at given state", labelpad=0.5, size=20)
    
    return

def hole_center(time_state, t_f):
    
    xyz = np.asarray(time_state[match][:, :, 1:4], dtype=float)
    dist = []
    for i in range(len(xyz)):
        d = []
        for j in range(len(xyz[0, :, :])):
            d.append(np.linalg.norm(xyz[i, j+1, :] - xyz[i, j, :]))
            if j +2 == len(xyz[0, :, :]):
                break
        dist.append(d)
    
    time = np.asarray(time_state[match][0, :, 0][1:], dtype=float)
    
    # plt.plot(time, dist[0], "--")
    # # plt.xlim([0, t_f])
    # # plt.ylim([-0.01, 1.01])
    # plt.xticks(size=20)
    # plt.yticks(size=20)
    # plt.legend("Hole_center_variation", loc = "upper center", fontsize = 20)
    # plt.ylabel('distance (Angs)', fontsize = 20)
    # plt.xlabel('Time (fs)', fontsize = 20)

    """
    Average 
    """
    
    dist = np.asarray(dist, dtype=float)
    plt.plot(time, np.average(dist, axis=0), "-*")
    # plt.xlim([0, t_f])
    # plt.ylim([-0.01, 1.01])
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend("Hole_center_variation", loc = "upper center", fontsize = 20)
    plt.ylabel('Mean distance (Angs)', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)


    return {"hole_center" : xyz}


def charge_center(time_state, t_f):
    
    xyz = np.asarray(time_state[match1][:, :, 1:4], dtype=float)
    dist = []
    for i in range(len(xyz)):
        d = []
        for j in range(len(xyz[0, :, :])):
            d.append(np.linalg.norm(xyz[i, j+1, :] - xyz[i, j, :]))
            if j +2 == len(xyz[0, :, :]):
                break
        dist.append(d)
    
    time = np.asarray(time_state[match1][0, :, 0][1:], dtype=float)
    
    # plt.plot(time, dist[0], "--")
    # # plt.xlim([0, t_f])
    # # plt.ylim([-0.01, 1.01])
    # plt.xticks(size=20)
    # plt.yticks(size=20)
    # plt.legend("Charge_center_variation", loc = "upper center", fontsize = 20)
    # plt.ylabel('distance (Angs)', fontsize = 20)
    # plt.xlabel('Time (fs)', fontsize = 20)
    
    
    """
    Average 
    """
    
    dist = np.asarray(dist, dtype=float)
    plt.plot(time, np.average(dist, axis=0), "-*")
    # plt.xlim([0, t_f])
    # plt.ylim([-0.01, 1.01])
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend("Charge_center_variation", loc = "upper center", fontsize = 20)
    plt.ylabel('Mean distance (Angs)', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)

    return {"charge_center" : xyz}


def distance_hole_charge_center(time_state, hole_cent, charge_cent):
    dist = np.zeros((np.shape(hole_cent["hole_center"])[0], np.shape(hole_cent["hole_center"])[1]))
    
    for i in range(len(dist)):
        dist[i, :] = np.linalg.norm(hole_cent["hole_center"][i, :] - charge_cent['charge_center'][i, :], axis=1)
    
    time = np.asarray(time_state[match][0, :, 0], float)
    
    # plt.plot(time, dist[0, :] , "--")
    # plt.xticks(size=20)
    # plt.yticks(size=20)
    # plt.legend("Charge_center_variation", loc = "upper center", fontsize = 20)
    # plt.ylabel('distance (Angs)', fontsize = 20)
    # plt.xlabel('Time (fs)', fontsize = 20)
    
    """
    Average
    """
    
    # plt.errorbar(time, np.average(dist, axis=0), np.var(dist, axis=0), marker='^')
    plt.errorbar(time, np.average(dist, axis=0), np.std(dist, axis=0), marker='^')
    plt.xticks(size=20)
    plt.yticks(size=20)
    # plt.legend("Charge_center_variation", loc = "upper center", fontsize = 20)
    plt.ylabel('Mean distance between hole and charge (Angs)', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)
    
    
    
    return 


def KE_EN(time_state, t_f):
    
    """
    how many cluster in QM region and MM region 
    """
    
    qm_mol = 13
    QM = time_state[match][:, :, 1:qm_mol*3*3 + 1]
    MM = time_state[match][:, :, qm_mol*3*3 + 1:]
    QM = QM.reshape((np.shape(QM)[0], np.shape(QM)[1], int(np.shape(QM)[2]/3), 3))
    MM = MM.reshape((np.shape(MM)[0], np.shape(MM)[1], int(np.shape(MM)[2]/3), 3))
    
    frame00QM = []; frame00MM = []
    for i in range(np.shape(QM)[1]):
        frame00QM.append(sum((np.linalg.norm(QM[0, i, :], axis = 1))**2))
        frame00MM.append(sum((np.linalg.norm(MM[0, i, :], axis = 1))**2))
    
    x = np.linspace(0,200,401)
    line, = plt.plot(x, frame00QM, "-*")
    line.set_label("QM region water molecules")
    plt.xlim([0, t_f])
    # plt.ylim([-0.01, 1.01])
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend(loc = "upper center", fontsize = 20)
    plt.ylabel('Velocity**2 (A.u)', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)
    
    
    x = np.linspace(0,200,401)
    line, = plt.plot(x, frame00MM, "-*")
    line.set_label("MM region water molecules")
    plt.xlim([0, t_f])
    # plt.ylim([-0.01, 1.01])
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend(loc = "upper center", fontsize = 20)
    plt.ylabel('Velocity**2 (A.u)', fontsize = 20)
    plt.xlabel('Time (fs)', fontsize = 20)
        
    return frame00QM, frame00MM

def main():
    time_state = readlog(list_, t_f)
    time_state1 = readlog1(list1_, t_f)
    
    # plotSlog(time_state, t_f)                  
    # plotKEPE_(time_state, t_f)
    # plotAD_V_Orbital_E(time_state, t_f)
    # hole_pop(time_state, t_f)
    KE_EN(time_state, t_f)
    hole_cent = hole_center(time_state, t_f)
    charge_cent = charge_center(time_state1, t_f)
    distance_hole_charge_center(time_state, hole_cent, charge_cent)
    
    # partial_charge(time_state, t_f)
    # Plog(time_state, t_f)
    
    
    return













"""
hole_populations_w_point_charges
reading whole *.dat files
"""

# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "hole_populations_w_point_charges_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(len(lines)):
#         partial_charge.append(np.asarray(lines[i].split()[1:-1],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(len(partial_charge))]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):
#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='both', which='major', labelsize=20)   
# plt.xlabel('Atom number (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("hole density", size=18)
# cbar.ax.tick_params(labelsize=15)


"""
partial_charges_w_point_charges_
"""


# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "partial_charges_w_point_charges_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(len(lines)):
#         partial_charge.append(np.asarray(lines[i].split()[1:-1],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(len(partial_charge))]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):
#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='both', which='major', labelsize=20)   
# plt.xlabel('Atom number (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("charge density", size=18)
# cbar.ax.tick_params(labelsize=15)

"""
probability to find hole as the given state
"""


# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "P_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(401):
#         partial_charge.append(np.asarray(lines[i].split()[1:-1],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(401)]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):
#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='bothYlGnBu', which='major', labelsize=20)   
# plt.xlabel('Occupation of state (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("charge density", size=18)
# cbar.ax.tick_params(labelsize=15)





"""
hole_density distribution on valence states
reading whole *.dat files
"""

# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "P_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(401):
#         partial_charge.append(np.asarray(lines[i].split()[1:-1],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(401)]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):Wednesday, June 29th:

#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='both', which='major', labelsize=20)   
# plt.xlabel('Occupation of state (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("charge density", size=18)
# cbar.ax.tick_params(labelsize=15)


"""
partial_charges_w_point_charges_
"""


# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "partial_charges_w_point_charges_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(len(lines)):
#         partial_charge.append(np.asarray(lines[i].split()[1:-1],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(len(partial_charge))]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):
#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='both', which='major', labelsize=20)   
# plt.xlabel('Atom number (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("charge density", size=18)
# cbar.ax.tick_params(labelsize=15)


"""
occupied_orbital_energies_w_point_charges
"""


# ls = ["ls","*.dat"] 
# ls_dat = subprocess.Popen(ls, shell=True, stdout=subprocess.PIPE) 
# list_  = [i.decode('ascii') for i in ls_dat.communicate()[0].split()]
# match = "occupied_orbital_energies_w_point_charges_"
# list_ = [i for i in list_ if bool(re.match(match, i)) is True]


# cross=[];En=[]
# for file in list_[:]:
#     f = open(file)
#     lines=f.readlines()
#     partial_charge=[];    
#     for i in range(4,len(lines)):
#         partial_charge.append(np.asarray([float(i) for i in lines[i].split()[14:] if float(i)<100],dtype=float))

# tick_atm = [[i for i in range(1,len(partial_charge[0])+1)] for j in range(len(partial_charge))]

# tick_time=[];time=[i for i in np.arange(0,200.5,0.5)]
# for i in time:
#     tick_=[];
#     for j in range(np.shape(partial_charge)[1]):
#         tick_.append(i)
#     tick_time.append(tick_) 


# plt.contourf(np.asarray(tick_atm,dtype=float),np.asarray(tick_time,dtype=float),partial_charge)
# plt.tick_params(axis='both', which='major', labelsize=20)   
# plt.xlabel('Number of state (Integer)', fontsize=20)
# plt.ylabel('Time (fs)', fontsize=20)
# cbar=plt.colorbar()
# cbar.set_label("Orbital Energy (ev)", size=18)
# cbar.ax.tick_params(labelsize=15)














