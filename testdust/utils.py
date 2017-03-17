# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:39:07 2016

@author: ibackus
"""
import itertools
import diskpy
import glob
import os
import shutil
import numpy as np
from types import ModuleType
import pynbody
SimArray = pynbody.array.SimArray

kB = SimArray(1.0, 'k')


def scaleTimeStep(param, scale=1):
    """
    Rescales the time-step for a simulation.  This is useful if you want to,
    for instance, double the dDelta but force all the outputs to be at the same
    physical time
    
    Parameters
    ----------
    param : dict
        A param dictionary (see diskpy.utils.configparser)
    scale : number
        Factor by which to scale dDelta.  dDelta_new = dDelta * scale, and all
        outputs will be made at the same physical time (or as close to it as
        possible, based on rounding)
    """
    scale = float(scale)
    key = 'dDelta'
    if key in param:        
        param[key] *= scale
        
    key = 'iOutInterval'
    if key in param:
        val = param[key]
        val = np.round(val/scale).astype(int)
        val = max([1, val])
        param[key] = val
        
    key = 'nSteps'
    if key in param:
        val = param[key]
        if val > 0:
            val = np.round(val/scale).astype(int)
            val = max([1, val])
            param[key] = val
        
    key = 'iLogInterval'
    if key in param:
        val = param[key]
        val = np.round(val/scale).astype(int)
        val = max([1, val])
        param[key] = val
        
    key = 'iCheckInterval'
    if key in param:
        val = param[key]
        val = np.round(val/scale).astype(int)
        val = max([1, val])
        param[key] = val
             
    key = 'dDumpFrameStep'
    if key in param:
        val = param[key]
        param[key] = val/scale


def periodicGrid(res, shape=[1,1,1]):
    """
    Generates an ND grid in a periodic box of a given shape, centered on the
    origin.
    
    Parameters
    ----------
    res: float or array-like
        Number of grid points along each axis.  If an integer, the resolution
        will be the same along each axis.
    shape: array-like
        Size of the grid along each axis.  number of dimensions is defined by
        len(shape)
    
    Returns
    -------
    x: array
        Array of positions.  x[i] is a vector of the nD position of point i.
    """
    shape = np.asarray(shape).astype(float)
    ndim = len(shape)
    if ndim == 0:
        # Assume we're trying to do 1D
        shape = shape[None]
        ndim = 1
    res = np.asarray(res).astype(int)
    if res.size == 1:
        # Use same res along each axis
        res = res * np.ones(ndim, dtype=int)
    if len(res) != ndim:
        raise ValueError, "Shape mismatch with res and shape"
    edges = []
    
    for i in range(ndim):
        
        dx = shape[i]/res[i]
        xmax = (shape[i]-dx)/2.
        edge = np.linspace(-xmax, xmax, res[i])
        edges.append(edge)
        
    npoints = np.product(res)
    pos = np.zeros([npoints, ndim])
    
    for i, x in enumerate(itertools.product(*edges)):
        
        pos[i] = x
        
    return pos

def calcGasVals(f, param={}):
    """
    Calculates gas values for a snapshot or subsnap (must be gas particles only)
    for a dusty-gas mixture.  A param dict/path to param can be supplied.
    
    Calculates:
        'rhog': gas density
        'pressure': gas pressure
        'energy': gas internal specific energy
        
    Also sets 'dustFrac' units to 1.
    
    Returns gamma and mean molecular weight
    """
    if isinstance(param, str):
        
        param = diskpy.utils.configparser(param,'param')
        
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
     
    return gamma, molecularWeight

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
    # Sanity checks first
    if not os.path.exists(defaults) or not os.path.isfile(defaults):
        raise IOError, "could not find regular file: {}".format(defaults)
    if not os.path.exists(userdefaults) or not os.path.isfile(userdefaults):
        raise IOError, "could not find regular file: {}".format(userdefaults)    
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