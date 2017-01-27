# -*- coding: utf-8 -*-
"""
This module contains utilities useful for the dust settling test

Created on Wed Oct 19 16:11:02 2016

@author: ibackus
"""
import os
from scipy.interpolate import interp1d
from scipy.special import erf
import numpy as np

import testdust
from testdust.utils import loadDefaultParam as _loadDefaultParam

# Get the path of this file
_directory = os.path.dirname(os.path.realpath(__file__))

# Filenames for the .param files (defaults AND user defaults)
parnameDefaults = {'dust': 'dust_default.param',
                   'gas': 'gas_default.param',
                   'general': 'general_default.param'}
parnameUserDefaults = {'dust': 'dust_user_default.param',
                       'gas': 'gas_user_default.param',
                       'general': 'general_user_default.param'}
# Get full path to all the .param files
for parnames in (parnameDefaults, parnameUserDefaults):
    
    for k, v in parnames.iteritems():
        
        parnames[k] = os.path.join(_directory, v)

def loadDefaultParam(param, fullPars=True):
    """
    Loads the settling .param defaults for param (returned as a dict).  For
    a list of available params do testdust.settling.utils.paramDefaults.keys()
    
    If doing fullPars, the full defaults are generated.  Otherwise, only those
    specific to that step are generated
    """
    if param not in parnameDefaults.keys():
        
        raise RuntimeError, "{0} not a valid param key.  valid options are: {1}"\
        .format(param, parnameDefaults.keys())
    
    parDict = {}
    
    if param is 'gas':
        # First load the general params
        if fullPars:
            parDict = loadDefaultParam('general')
        gasPar = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        parDict.update(gasPar)
        
    elif param is 'dust':
        # Load the gas defaults
        if fullPars:
            parDict = loadDefaultParam('gas')
        dustPar = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        parDict.update(dustPar)
        # Turn off glass damping
        if 'dGlassDamper' in parDict:
            parDict.pop('dGlassDamper', None)
        
    else:
        
        parDict = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        
    return parDict
    
def defaultFilenames(step):
    """
    Get the default filenames for 'gas', 'dust', or 'glass'
    
    step can be 'gas', 'dust', or 'glass'
    
    Returns a tuple of (infile, paramname, fprefix)
    """
    if step is 'glass':
        infile = 'glass.std'
        fprefix = 'glass'
    else:
        parDict = loadDefaultParam(step)
        infile, fprefix = testdust.utils.filenamesFromParam(parDict)
    paramname = fprefix + '.param'
    return infile, paramname, fprefix

def makeInvErf(zmax=30, npts=1e5):
    """
    Generates an inverse error function by making a linear interpolation of
    z vs error.  The error function here is scipy.special.erf
    
    Parameters
    ----------
    zmax : int
        Maximum z to for erf(z) to consider.
    npts : int
        Number of points to use for generating the interpolation spline
    
    Returns
    -------
    invErf(m) : function
        Inverse of the error function, for 0 <= m <= 1
    """
    npts = int(npts)
    z = np.linspace(0, zmax, npts)
    y = erf(z)
    y[-1] = 1
    return interp1d(y, z)

def gaussianGrid(nz=50):
    """
    Generate a set of points centered around zero to create a gaussian
    density profile for SPH.
    
    If you have particles of equal mass, placing the particles at these points
    will create a density profile that decays as exp(-z**2/2).  This gives a 
    standard deviation of 1.
    
    Parameters
    ----------
    nz : int
        Number of points to make
    
    Returns
    -------
    z : array
        Numpy array of points, z >= 0
    """
    erfInv = makeInvErf()
    m = np.linspace(-1, 1, nz+2)[1:-1]
    negs = m < 0
    m = abs(m)
    z = erfInv(m)
    z *= np.sqrt(2)
    z[negs] *= -1
    return z
    
def _estRho(z):
    """
    Estimate density for equal mass particles on a line (up to normalization)
    """
    dz = np.gradient(z, edge_order=2)
    rho = 1.0/dz
    rho /= rho[0]
    return rho
    
def _rhoExpected(z):
    """
    """
    return np.exp(-z**2)