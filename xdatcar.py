#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 13:55:12 2022

@author: EagleLeslie
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from scipy import stats

# CONSTANTS
# Boltzmann constant eV/K
kb = 8.617332478E-5

# Boltzman constant in J/K
kb_J = 1.380649E-23

# Avogadro's constant
avogadro = 6.0221412927E23

class Xdatcar:
    
    """
    
    Python class for VASP output files
    
    """
    
    def __init__(self):
        
        self.ntimesteps = None
        self.atypes = None
        self.nions = None
        self.nelem = None
        self.totatoms = None
        
        self.xyz = None
        
        self.xdatcar = None
        self.cell = None
        self.volume = None
        self.read_xdatcar()
        
        self.ncomb = None
        
        self.nametype = None
        self.ntype = None
        self.symb = None
        
    def ion_name(self):
        """
        
        Method for linking atom names with XDATCAR xyz positions
        
        """
        
        self.symb = []
        if self.nametype is None:
            self.nametype = [chr(i) for i in range(65,91)[:self.ntype]]

        for i in range(self.ntype):
            self.symb += [np.tile(self.nametype[i].self.nelem[i])]
            
        self.symb = np.concatenate(self.symb)
        
        
    def read_xdatcar(self,File=None):
        """
        
        Borrowed and altered from @QijingZheng
        
        Method for reading in XDATCAR file
        
        """
        if File is None:
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = File
            
        data = [line for line in open(self.xdatcar) if line.strip()]
        #header = str(data[0])
        scale = float(data[1])
        
        self.cell = np.array([line.split() for line in data[2:5]],dtype=float)
        self.cell *= scale
        self.volume = np.linalg.det(self.cell)
        
        atom_types = data[5].split()
        atom_num = data[6].split()
        
        # Check that atom_types is a string in alphabet
        if atom_types[0].isalpha():
            self.atypes = atom_types
            self.nions = len(atom_types)
            self.nelem = np.array(atom_num, dtype=int)
            self.totatoms = self.nelem.sum()
            
    
        # Read in positions from XDATCAR
        # NOTE: READING IN AS DIRECT COORDINATES
        positions = np.array([line.split() for line in data[7:] 
                              if not line.split()[0].isalpha()],dtype=float)
        self.xyz = positions.ravel().reshape((-1,self.totatoms,3))
        self.ntimesteps = self.xyz.shape[0]
        
        
