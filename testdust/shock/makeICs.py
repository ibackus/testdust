#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Contains functions for generating the dustyshock tube of Price & Laibe 2015

See shock1d() to make a 1D shock tube.

Price & Laibe shocktube params:
Tube shape:     (2)
Pressure:       (1, 1/8)
Gas Density:    (1, 1/8)
Dust fraction:  (0.5, 0.5)
gamma:          5/3

Created on Fri Feb 24 10:15:21 2017

@author: ibackus
"""

import cPickle as pickle
import numpy as np
import pynbody
import diskpy
import testdust
import os

# Constants
SimArray = pynbody.array.SimArray
kB = pynbody.units.k
defaultParam1d = 'shock_default1d.param'
userDefault1d = 'shock_user_defaults1d.param'
_directory = os.path.dirname(os.path.realpath(__file__))
defaultParam1d = os.path.join(_directory, defaultParam1d)
userDefault1d = os.path.join(_directory, userDefault1d)
argsavename = 'icparams.p'

def shock1d(nParticles=2*569, L=2., P=(1, 1./8), rhog=(1, 1./8), dustFrac=0.5, 
            gamma=5./3, inputParamName=None):
    """
    Creates SPH ICs for a 1D dusty-shock tube a la Price & Laibe 2015
    (see their dustyshock, section 4.2).
    
    Particles are in a line along the z-axis.  Note that since ChaNGa uses 
    periodic boundaries (rather than particles in a tube) the shock tube
    should be twice as long and will actually produce 2 shocks.
    
    Parameters
    ----------
    nParticles: int
        Approximate number of particles to use.  nParticles will be set in the
        end to best approximate the gas density in each region.
    L: float
        Length of the shock tube. (in code units)
    P: array-like
        Pressure in the two regions (left and right)
    rhog: array-like
        Gas density in the two regions
    dustFrac: float
        Dust fraction (between 0 and 1)
    gamma: float
        Adiabatic index
    inputParamName: str
        (optional) Override the default ChaNGa params.  This can also be
        done by changing the appropriate user default param in this directory.
    
    Returns
    -------
    ICs: SimSnap
        Pynbody snapshot of the ICs
    param: dict
        dict of the runtime params
    arguments: dict
        Arguments used when calling this function, including defaults.
    paramsavename: str
        Path the param file is saved to
    """
    # Keep track of arguments (and save at the end)
    arguments = locals()
    
    P = np.asarray(P)
    rhog = np.asarray(P)
    L = float(L)
    gamma = float(gamma)
    
    # --------------------------------------------
    # Initialize
    # --------------------------------------------
    # Load stuff
    if inputParamName is not None:
        param = diskpy.utils.configparser(inputParamName)
    else:
        param = testdust.utils.loadDefaultParam(defaultParam1d, userDefault1d)
#    param = diskpy.utils.configparser(inputParamName, 'param')
    fprefix = param['achOutName']
    savename = param['achInFile']
    molecularWeight = SimArray(param['dMeanMolWeight'], 'm_p')
    units = diskpy.pychanga.units_from_param(param)
    # Initialize quantities
    L = SimArray(L, units['l_unit'])
    P = SimArray(P, units['pres_unit'])
    rhog = SimArray(rhog, units['rho_unit'])
    rho = rhog/(1. - dustFrac)
    # Set-up number of particles to left/right of boundary
    N = np.array(np.round(nParticles * rho/rho.sum()).astype(int))
    nParticles = N.sum()
    # Initialize snapshot
    f = pynbody.new(gas=nParticles)
    # Initialize arrays
    for key in ('temp', 'rho', 'dustFrac'):
        f[key] = 0.
    # Set up left/right slices.  These select the left/right regions
    slices = [slice(0, N[0]), slice(N[0], N[0] + N[1])]
    
    # --------------------------------------------
    # Setup values
    # --------------------------------------------
    # Mass
    totMass = float((rho * (L/2.)).sum())
    f['mass'] = totMass/nParticles
    # Position
    z = [0.5 * L * testdust.utils.periodicLine(n) for n in N]
    f['z'][slices[0]] = z[0] - 0.25 * L
    f['z'][slices[1]] = z[1] + 0.25 * L
    
    # Temperature
    T = molecularWeight * P/(rhog * kB)
    T.convert_units('K')
    
     # Assign arrays, looping over both sides
    for i, s in enumerate(slices):
        f['temp'][s] = T[i]
        f['rho'][s] = rho[i]
        f['dustFrac'][s] = dustFrac
    
    # Other
    f['eps'] = 10*L/nParticles
    
    # Setup params
    testdust.utils.setupParamBounds(param, [float(L)])
    param['dConstGamma'] = gamma
    
    # Save
    f.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
    print 'snapshot saved to:', savename
    paramsavename = fprefix + '.param'
    diskpy.utils.configsave(param, paramsavename)
    print '.param file saved to:', paramsavename
    
    pickle.dump(arguments, open(argsavename, 'w'), 2)
    print 'arguments dict pickled to:', argsavename
    
    return f, param, arguments, paramsavename
