#!/usr/bin/env python
# coding: utf-8

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

from math import *

######################################################################################
## Rotation matrix.

def Rx(theta):
    return np.matrix([[ 1, 0           , 0           ],
                   [ 0, cos(theta),-sin(theta)],
                   [ 0, sin(theta), cos(theta)]])
  
def Ry(theta):
    return np.matrix([[ cos(theta), 0, sin(theta)],
                   [ 0           , 1, 0           ],
                   [-sin(theta), 0, cos(theta)]])
  
def Rz(theta):
    return np.matrix([[ cos(theta), -sin(theta), 0 ],
                   [ sin(theta), cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

######################################################################################
## Reading an initial file where the configuration of a single QD is found.

print('READING THE INITIAL QD')

namefile = 'Single_QD'

if not os.path.isfile(namefile):
    print("Error: The file does not exist.")
    sys.exit(1)

data = pd.read_csv(namefile, header=None, skiprows=1, nrows=1, sep='\s+')
Np   = int(data.values[0][0])

data   = pd.read_csv(namefile, header=None, skiprows=9, nrows=Np, sep='\s+')
Atom   = data.values[:,0]
Type   = data.values[:,1]
Charge = data.values[:,2]
Pos    = data.values[:,3:6]

xcm    = np.mean(Pos[:,0]); Pos[:,0] = Pos[:,0]-xcm
ycm    = np.mean(Pos[:,1]); Pos[:,1] = Pos[:,1]-ycm
zcm    = np.mean(Pos[:,2]); Pos[:,2] = Pos[:,2]-zcm

######################################################################################
## Replicating the QD over the positions of a hexagonal crystal lattice.

### Parameters.
Xsites = 25                               # Size of the grid.
l0     = 85                               # Spacing lattice.


Npsites = Xsites*Xsites
PosLatt = np.zeros([Npsites,3])

x = 0
y = 0
z = 0

### Basis vectors
a1    = np.zeros([3], float)
a1[0] = l0

a2    = np.zeros([3], float)
a2[0] = 0.5*l0
a2[1] = 0.5*sqrt(3.)*l0

print('GENERATING A HEXAGONAL LATTICE')

account = 0
for iy in range(Xsites):
    for ix in range(Xsites):
            
        if (iy%2 == 0):
            PosLatt[account,0] = ix*a1[0]
        else:
            PosLatt[account,0] = ix*a1[0]+a2[0]  

        PosLatt[account,1] = ix*a1[1]+iy*a2[1]
        PosLatt[account,2] = 0

        account = account+1

xcm = np.mean(PosLatt[:,0]); PosLatt[:,0] = PosLatt[:,0]-xcm; Lx = max(PosLatt[:,0])-min(PosLatt[:,0])+0.5*l0
ycm = np.mean(PosLatt[:,1]); PosLatt[:,1] = PosLatt[:,1]-ycm; Ly = max(PosLatt[:,1])-min(PosLatt[:,1])+0.5*sqrt(3)*l0
zcm = np.mean(PosLatt[:,2]); PosLatt[:,2] = PosLatt[:,2]-zcm

######################################################################################
## Writing initial configuration for LAMMPS script.
## Pay attention to the simulation box boundaries. When rotating the QDs, 
## some atoms may end up outside the boundaries.

print('CREATING THE LAMMPS CONFIGURATION')

aux_vec  = np.zeros([3])

fdata = open('Hex_Quantum_Dots', 'w')
fdata.write('LAMMPS data file\n')
fdata.write('\n')
fdata.write('%d atoms\n'  % (Npsites*Np))
fdata.write('2 atom types\n')
fdata.write('\n')
fdata.write('%lf %lf xlo xhi\n' % (-1150.0, 1150.0))
fdata.write('%lf %lf ylo yhi\n' % (-950.0, 950.0))
fdata.write('%lf %lf zlo zhi\n' % (-50, 50))

fdata.write('\n')
fdata.write('Atoms\n')
fdata.write('\n')
id = 0
for n in range(Npsites):
    phi   = pi*np.random.random()    
    theta = pi*np.random.random()    
    psi   = 2.*pi*np.random.random() 
  
    R     = np.eye(3) * np.eye(3) * Rz(psi)            # Rotation along the z-axis.

    # Uncomment these lines to apply a random lateral displacement.
    Deltx = 0. #4*(np.random.random()-0.5)
    Delty = 0. #4*(np.random.random()-0.5)
    Deltz = 0. #4*(np.random.random()-0.5)
    
    for i in range(Np):    
        id = id+1 
                
        aux_vec[0] = Pos[i,0]
        aux_vec[1] = Pos[i,1]
        aux_vec[2] = Pos[i,2]
        aux_rot    = aux_vec*R
        
        aux_rot[0,0] = aux_rot[0,0]+PosLatt[n,0]+Deltx
        aux_rot[0,1] = aux_rot[0,1]+PosLatt[n,1]+Delty
        aux_rot[0,2] = aux_rot[0,2]+PosLatt[n,2]
        
        fdata.write('%d %d %lf %lf %lf %lf\n' % (id, Type[i], Charge[i], aux_rot[0,0], aux_rot[0,1], aux_rot[0,2]))
  
fdata.close()

print('END')

