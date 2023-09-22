#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Feb 25 12:06:44 2022

@author: rd
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
from optparse import OptionParser
import sys
import warnings
import subprocess
import csv
 
 
"""
basic units in GMX, please look at the following web.
https://manual.gromacs.org/documentation/2019/reference-manual/definitions.html#table-basicunits
 
1 Picoseconds = 41341.374575751 Atomic Unit Of Time
1 nm = 18.8973 Bohr
1 nm/ps = 18.8973/41341.374575751 Bohr/A.u = 0.00045710381413114194
1 nm/ps = 4.571038E-04 Bohr/A.u
"""
 
"""
Reading the input file and covert Angstrom, time to Bohr and atomic unit respectively.
"""
 
"""
Please use python3 "program".py -h to find out
how to use this program properly.
"""
 
class xsample_input(object):
 
    def get_arg(self, Grofile=None, Centerofmolecule=None, **arg): 
        parser=OptionParser() 
        parser = OptionParser(usage="usage: %prog [options] filename")
        parser.add_option('-g','--GRO_file',
                dest='Grofile',
                type='string',
                default="water_prdc_000.gro",
                help='GMX .gro file which is terminated after post production calculation.')
     
        parser.add_option('-n','--cen_mol',
                dest='Centerofmolecule',
                type=int,
                default=13,
                help='​Please enter the number of molecules which you want to center. It is an integer number and not exceeding the number of water molecules in the box​.')
         
        return parser
     
    def start(self, **arg):
        parse=self.get_arg()
        self.opt, self.arg=parse.parse_args()
        return 
    
     
    def readfile(self, filename=str):
        filename = self.opt.Grofile
        f = open(filename, "r")
        _lines = f.readlines()
        bohr = 18.8973
        bohrfs = 4.571038E-04

        global bxabc
         
        xyz=[];vxyz=[];topl=[];elms = [];
        for line in _lines[2:]:
            if len(line.split()) == 3:
                bxabc = np.asarray(line.split(), dtype=float)*bohr
            else:
                topxyzv = line.split()
                xyz.append(np.asarray(topxyzv[3:6], dtype=float)*bohr)
                vxyz.append(np.asarray(topxyzv[6:], dtype=float)*bohrfs)
                topl.append(topxyzv[:2])
                elms.append(topxyzv[1][0])
        return np.asarray(xyz,dtype=float), np.asarray(vxyz,dtype=float), np.asarray(topl,dtype=str), elms
          
    """
    This will read your number of molecules and find out its O atom's coordinate.
    Please remember that we always shift O atoms first than H atoms.
    """
    
    def tallash(self, Ncenter=int):
        Ncenter=self.opt.Centerofmolecule-1
        _xyz=self.readfile()
        newcenmol_xyz=_xyz[0][Ncenter*3]
        if self.readfile()[3][Ncenter*3] == "O":
            print(f'You have chosen {self.opt.Centerofmolecule}th water molecule for center of box Coor:{newcenmol_xyz[0]:9.6f},{newcenmol_xyz[1]:9.6f},{newcenmol_xyz[2]:9.6f}')
        else:
            raise ValueError('Please choose coordinates of O atom.')
        return newcenmol_xyz, _xyz[0], _xyz[1], _xyz[2], _xyz[3]


    """
    Once "prince" molecules have been chosen for center then move
    rest of the molecules which are outside the new box.
    """
     
    def shifting(self, center_new = list, new_xyz = list):
        talla = self.tallash()
        center_new = talla[0]
        new_xyz = talla[1]
        vxyz = np.ndarray.round(talla[2], 4); 
        elems = talla[4]
        
        for j in range(len(new_xyz)):
            """
            choose all O atoms with (mod(j,3)) and calculate distance between center O atom and
            rest of the O atoms
            """
            if np.mod(j, 3) == 0:
                dx = center_new[0] - new_xyz[j][0]
                dy = center_new[1] - new_xyz[j][1]
                dz = center_new[2] - new_xyz[j][2]
                """
                distances between O atoms and their H atoms
                """
                dxOH1 = new_xyz[j + 1, 0] - new_xyz[j, 0]
                dxOH2 = new_xyz[j + 2, 0] - new_xyz[j, 0]
                dyOH1 = new_xyz[j + 1, 1] - new_xyz[j, 1]
                dyOH2 = new_xyz[j + 2, 1] - new_xyz[j, 1]
                dzOH1 = new_xyz[j + 1, 2] - new_xyz[j, 2]
                dzOH2 = new_xyz[j + 2, 2] - new_xyz[j, 2]
                """
                If the distance is larger than half of the box(box size will be the   
                same), The molecule will be shifted right or left according to boundary condition.
                for Y and Z direction should be same.
                """
                if abs(dx) > bxabc[0]/2:
                  org_dx = center_new[0] - bxabc[0]/2
                  if dx > 0:
                    yingi_x = center_new[0] + abs(dx - org_dx)
                  else:
                    yingi_x = center_new[0] - abs(dx - org_dx)
                else:
                    yingi_x = new_xyz[j,0]
                new_xyz[j ,0] = yingi_x
                new_xyz[j + 1, 0] = yingi_x+dxOH1
                new_xyz[j + 2, 0] = yingi_x+dxOH2
                 
                if abs(dy) > bxabc[1]/2:
                  org_dy = center_new[1] - bxabc[1]/2
                  if dy > 0:
                    yingi_y = center_new[1] + abs(dy - org_dy)
                  else:
                    yingi_y = center_new[1] - abs(dy - org_dy)
                else:
                    yingi_y = new_xyz[j,1]
                new_xyz[j, 1] = yingi_y
                new_xyz[j + 1,1] = yingi_y+dyOH1
                new_xyz[j + 2,1] = yingi_y+dyOH2
                 
                if abs(dz) > bxabc[2]/2:
                  org_dz = center_new[2] - bxabc[2]/2
                  if dz > 0:
                    yingi_z = center_new[2] + abs(dz - org_dz)
                  else:
                    yingi_z = center_new[2] - abs(dz - org_dz)
                else:
                    yingi_z = new_xyz[j, 2]
                new_xyz[j, 2] = yingi_z
                new_xyz[j+1, 2] = yingi_z + dzOH1
                new_xyz[j+2, 2] = yingi_z + dzOH2
        return center_new, np.ndarray.round(new_xyz,3), vxyz, elems

     
    """
    New box was created with new center; Than sorted out the distances between
    center of molecules (Oxygen atom) and rest of molecules; than print out coordinates and velocities according to sorted index; Than create input file for "xsamle".
    """
     
    def print_xsample(self, new_xyz = list, center_new = list):
        dist = []; molecules = [];
        velocity = []; elm_mol = []
        shift = self.shifting()
        center_new = shift[0]; new_xyz = shift[1] 
        vxyz = np.ndarray.round(shift[2],4); 
        elems = shift[3]
        
        
        """
        it is bit tricky here, index will take from sorted distances. If there are same distances and return to their
        index would cause potential problem for sorting out velocity. So small decimal would add to corresponding
        coordinates if identical distance exist. This would not change any physics, it is pure mathematical trick. Hence
        result which you get from simulation would be identical.
        """
        for i in range(0,len(new_xyz),3):
            dist.append(np.sqrt((new_xyz[i][0] - center_new[0])**2+\
                                    (new_xyz[i][1] - center_new[1])**2+(new_xyz[i][2] - center_new[2])**2))

     
        """
        According to the sorted index, append coordinates and velocities then printed out.
        """
     
        for idx in np.argsort(dist):
            idx = idx*3
            molecules.append(new_xyz[idx])
            molecules.append(new_xyz[idx + 1])
            molecules.append(new_xyz[idx + 2])
            velocity.append(vxyz[idx])
            velocity.append(vxyz[idx + 1])
            velocity.append(vxyz[idx + 2])
         
        velocity = np.asarray(velocity, dtype=str)
        molecules = np.asarray(molecules, dtype=str)
        elm_mol = np.column_stack([elems, molecules])
         
             
        coor_head = ["$SYSTEM", "$quantum", "$sampling", "$cartPOS", "$cartVEL", "$end"];
        qchem = ["   xunit = bohr", "qchemistry = xmolecule", "type = gs", "   nsample = 1", 
                 "   type = monte-carlo"]
   

         
        f = open(self.opt.Grofile + ".in", "w")
        f.write(coor_head[0] +'\n')
        f.write(qchem[0] +'\n')
        f.write(coor_head[-1] +'\n\n')
        f.write(coor_head[2] +'\n')
        f.write(qchem[-1] +'\n')
        f.write(qchem[-2] +'\n')
        f.write(coor_head[-1] +'\n\n')
        f.write(coor_head[3] +'\n')
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(elm_mol)
        f.write(coor_head[-1] +'\n')
        f.write(coor_head[-2] +'\n')
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(velocity)
        f.write(coor_head[-1] +'\n')
        f.close()
         
         
        f = open(self.opt.Grofile + ".xyz", "w")
        f.write(f"{len(elm_mol)}" +'\n')
        f.write(f"Water_molecules{len(elm_mol)-1:6.0f}" +'\n')
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(elm_mol)
        f.close()
     
        f = open(self.opt.Grofile + ".v", "w")
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(velocity)
        f.close()
         
        return np.asarray(dist)
                 
         
if __name__ == '__main__':
    xs = xsample_input()
    xs.start()
    xs.print_xsample()
