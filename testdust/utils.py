# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:39:07 2016

@author: ibackus
"""
import diskpy
import glob
import os
import shutil
import numpy as np
from types import ModuleType
import pynbody
SimArray = pynbody.array.SimArray

kB = SimArray(1.0, 'k')

def calcGasVals(f, param={}):
    """
    Calculates gas values for a snapshot or subsnap (must be gas particles only)
    for a dusty-gas mixture.  A param dict/path to param can be supplied.
    
    Calculates:
        'rhog': gas density
        'pressure': gas pressure
        'energy': gas internal specific energy
        
    Also sets 'dustFrac' units to 1.
    """
    units = diskpy.pychanga.units_from_param(param)
    molecularWeight = diskpy.pychanga.getpar('dMeanMolWeight', param)
    molecularWeight = SimArray(molecularWeight, 'm_p')
    gamma = diskpy.pychanga.getpar('dConstGamma', param)
    
    # Calculate gas density
    f['dustFrac']
    f['dustFrac'].units = 1
    f['rhog'] = (1.0 - f['dustFrac']) * f['rho']
    
    # Calculate pressure
    pressure = f['rhog'] * kB * f['temp']/molecularWeight
    pressure.convert_units(units['pres_unit'])
    f['pressure'] = pressure
    
    # Calculate internal energy
    energy = f['pressure'] / (f['rhog'] * (gamma - 1))
    specific_energy_unit = (units['v_unit'])**2
    energy.convert_units(specific_energy_unit)
    f['energy'] = energy

def getcs(snap, param):
    """
    From a simulation snapshot and param file (or dict), return the average
    sound speed. (in simulation units)
    """
    
    # Get sound speed
    units = diskpy.pychanga.units_from_param(param)
    mu = units['m_unit']
    tu = units['t_unit']
    lu = units['l_unit']
    m = SimArray(param['dMeanMolWeight'], 'm_p', dtype=float)
    m.convert_units(mu)
    K = pynbody.units.Unit('K')
    kBunit = mu*lu**2/(K*tu**2)
    k = kB.in_units(kBunit)
    T = snap['temp'].mean()
    cs2 = k*T/m
    cs = SimArray(float(np.sqrt(cs2)), lu/tu)
    
    return cs


def setupParamBounds(param, boxShape):
    """
    Sets up the dxPeriod, dyPeriod, and dzPeriod parameters in the .param file
    dict according to boxShape.  boxShape is length 1, 2, or 3 depending on
    how many dimensions the periodic box is in.
    """
    nDim = len(boxShape)
    maxL = max(boxShape)
    if nDim < 3:
        param['dxPeriod'] = 100 * maxL
    else:
        param['dxPeriod'] = float(boxShape[0])
    if nDim < 2:
        param['dyPeriod'] = 100 * maxL
    else:
        param['dyPeriod'] = float(boxShape[-2])
        
    param['dzPeriod'] = float(boxShape[-1])
    
    if 'dPeriod' in param:
        
        param.pop('dPeriod', None)

def gridSize(nParticles, boxshape):
    """
    Makes n-dimensional grid dimensions that will place approximately nParticles
    on an evenly spaced grid of boxshape.  Due to rounding, the number of
    actual particles may be different than nParticles
    """
    volume = float(np.product(boxshape))
    dx = volume**(1./3)
    boxres = []
    for L in boxshape:
        res = int(np.round(L/dx))
        boxres.append(res)
    
    return np.asarray(boxres, dtype=int)
    
def periodicLine(n):
    """
    Places n points on a line length=1, centered at the origin, and periodic
    around x = +/-0.5
    """
    x = np.linspace(-0.5, 0.5, n+1)[0:-1]
    dx = x[1] - x[0]
    x += 0.5 * dx
    return x

def parsePySettings(filename):
    """
    Parses settings in a .py file by running execfile()
    Settings files for this format are just python scripts containing variable
    definitions
    
    Modules and the doc string will be ignored
    
    Returns a dict of the settings
    """
    settings = {}
    globalsDict = {}
    execfile(filename, globalsDict, settings)
    # Ignore the doc string for the file if there is one
    if '__doc__' in settings:
        settings.pop('__doc__', None)
        
    # Ignore any imported modules
    for key in settings.keys():
        
        if isinstance(settings[key], ModuleType):
            
            settings.pop(key, None)
    
    return settings

def copySnapshot(filename, dest):
    """
    Copies a snapshot and auxilliary arrays to dest.  dest can be a directory
    or a new filename
    """
    fnames, destinations = _snapshotSrcDest(filename, dest)
    for fname, destination in zip(fnames, destinations):
        shutil.copyfile(fname, destination)

def moveSnapshot(filename, dest):
    """
    Moves a snapshot and its auxilliary arrays to dest. dest can be a directory
    or a new filename
    """
    fnames, destinations = _snapshotSrcDest(filename, dest)
    for fname, destination in zip(fnames, destinations):
        shutil.move(fname, destination)
    
def _snapshotSrcDest(filename, dest):
    """
    A utility to get all the filenames of a snapshot and its auxilliary arrays
    and set up destination paths for a move or copy operation
    """
    fnames = glob.glob(filename + '*')
    exts = [fname.split(filename)[-1] for fname in fnames]
    
    if os.path.isdir(dest):
        
        destPrefix = os.path.join(dest, filename)
        
    else:
        
        destPrefix = dest
        
    destinations = [destPrefix + ext for ext in exts]
    
    return fnames, destinations

def loadDefaultParam(defaults, userdefaults):
    """
    Loads the parameters in the ChaNGa .param file defaults and overrides 
    them with those in userdefaults
    
    Parameters
    ----------
    defaults, userdefaults : str
        File names to the .param files to load.  userdefaults override defaults
    
    Returns
    -------
    params : dict
        A dictionary of the params
    """
    print 'Loading defaults from:', defaults
    defaults = diskpy.utils.configparser(defaults, 'param')
    print 'Loading user defaults from:', userdefaults
    userdefaults = diskpy.utils.configparser(userdefaults, 'param')
    defaults.update(userdefaults)
    return defaults

def filenamesFromParam(param):
    """
    Retrieves the input filename and output prefix for a .param file (or a 
    loaded .param dict)
    
    Returns (infile, fprefix)
    """
    if isinstance(param, str):
        param = diskpy.utils.configparser(param, 'param')
    
    inFile = param['achInFile']
    prefix = diskpy.pychanga.getpar('achOutName', param)
    return inFile, prefix