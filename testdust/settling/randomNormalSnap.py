#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Defines functions to create a pynbody snapshot of gas, normally distributed
along z, with a periodic shape along possibly x and y.  This works in 1D, 2D, 
and 3D.

Created on Thu Feb  9 15:08:05 2017

@author: ibackus
"""
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray


def randomNormal(n, height, baseshape=[]):
    """
    Generate random positions, normally distributed along z.  Base shape can be:
        []          (1D sim)
        [Ly]        (2D sim)
        [Lx, Ly]    (3D sim)
    Where Lx, Ly are lengths along x, y.
    """
    nDim = len(baseshape) + 1
    pos = np.zeros([n, nDim])
    z = np.random.randn(n)
    z *= height
    pos[:,-1] = z
       
    for i in range(nDim - 1):
        
        pos[:, i] = np.random.rand(n) * baseshape[i]
        
    return pos

def normalSnap(n, height, baseshape=[]):
    """
    Generate a snapshot with positions normally distributed along z
    """
    
    snap = pynbody.new(gas=n)
    nDim = len(baseshape) + 1
    pos = randomNormal(n, height, baseshape)
    i0 = 3-nDim
    snap['pos'][:, i0:] = SimArray(pos,'au')
    volume = np.sqrt(2*np.pi) * height
    if nDim > 1:
        volume *= np.prod(baseshape)
    snap['mass'] = volume*SimArray(np.ones(n), 'Msol')/n
    snap['vel'] = SimArray(np.zeros([n,3]), 'km s**-1')
    snap['temp'] = SimArray(np.ones(n),'K')
    snap['eps'] = SimArray(np.ones(n))*height * 5
    snap['rho'] = SimArray(np.ones(n), 'Msol kpc**-3')
    return snap