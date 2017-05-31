# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 13:15:04 2016

@author: ibackus
"""

import numpy as np
import matplotlib.pyplot as plt
import pynbody
SimArray = pynbody.array.SimArray
import os
import diskpy

# Definitions
G = SimArray(1.0, 'G')
kB = SimArray(1.0, 'k')
Msol = SimArray(1.0, 'Msol')


def loadSim(simdir, fprefix = 'settledust'):
    
    paramname = fprefix + '.param'
    paramname = os.path.join(simdir, paramname)
    fnames = diskpy.pychanga.get_fnames(fprefix, simdir)
    fs = [pynbody.load(fname, paramfile=paramname) for fname in fnames]
    parDict = diskpy.utils.configparser(paramname, 'param')
    
    return fs, parDict

def settlingTimes(sim, pars):
    """
    Returns snapshot times for the simulation
    
    Parameters
    ----------
    sim : list 
        A list of SimSnaps for the simulation
    pars : dict
        param dict for the simulation
    
    Returns
    -------
    t : SimArray
        Times of the snapshots.
    """
    units = diskpy.pychanga.units_from_param(pars)
    dtOut = pars['dDelta'] * pars['iOutInterval']
    t = np.arange(len(sim)) * dtOut
    t = SimArray(t, units['t_unit'])
    return t

def settlingPeriod(pars):
    """
    """
    period = 2 * np.pi/settlingOmega(pars)
    return period

def settlingOmega(pars):
    units = diskpy.pychanga.units_from_param(pars)
    M = SimArray(pars['dCentMass'], units['m_unit'])
    R = SimArray(pars['dOrbDist'], units['l_unit'])
    omega = np.sqrt(G * M/ R**3)
    omega.convert_units(1/units['t_unit'])
    return omega

def discreteHistBins(data):
    d = np.diff(np.unique(data)).min()
    left_of_first_bin = data.min() - float(d)/2
    right_of_last_bin = data.max() + float(d)/2
    bins = np.arange(left_of_first_bin, right_of_last_bin + d, d)
    return bins

def snapdefaults(f):
    if not isinstance(f, pynbody.snapshot.SimSnap):
        for snap in f:
            snapdefaults(snap)
        return
    f.g['rho'].convert_units('g cm**-3')
    f['pos'].convert_units('au')