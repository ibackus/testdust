#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Contains functions for generating ICs for the dusty diffusion test of 
Price & Laibe 2015.  Note that to run these simulations, particle positions
must be fixed in place!

Note: This assumes grain size = 1 and grain density = 1 for estimating the
stopping time and thereby estimating the proper SPH particle mass.

UNITS:
lengths are kpc and masses are Msol

Created on Thu Mar  3 12:51:27 2016

@author: ibackus
"""
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import diskpy
import cPickle as pickle
#import json
import os

import testdust
from testdust.utils import periodicGrid, scaleTimeStep

# Settings
kB = SimArray(1.0, 'k')
defaultParam = 'diffusion_defaults.param'
userDefaults = 'diffusion_user_defaults.param'
runParamsName='runparams.p'
# Setup default param paths
_directory = os.path.dirname(os.path.realpath(__file__))
defaultParam = os.path.join(_directory, defaultParam)
userDefaults = os.path.join(_directory, userDefaults)

def makeICs(simdir='difftest', cs=1., boxSize=1., criticalRadius=0.25, dustFrac0=0.1,
            tStop=0.1, rho=1., nSmooth=32, dustFracMin=0., inputSnap=None, 
            dDeltaFactor=1, boxRes=[50, 58, 60]):
    r"""
    Generates ICs for the dustydiffusion test of Price & Laibe 2015.
    
    The ICs are in a periodic cube, centered on the origin, with uniform 
    density and a dust fraction decaying from the origin.  The dust fraction
    profile is set as:
        
    .. math::
        
        \epsilon(r, 0) = \epsilon_0 (1 - (r/r_c)^2)
        
    ICs will be saved to simdir
    
    Parameters
    ----------
    simdir : str
        Directory to save ICs in
    cs : float
        Sound speed (code units)
    boxSize : float
        Size of the periodic cube (code units)
    criticalRadius : float
        :math:`r_c` in the equation above
    dustFrac0 : float
        :math:`\epsilon_0` in the equation above
    tStop : float
        Dust stopping time (code units)
    rho : float
        total density (code units)
    nSmooth : int
        Number of neighbors for smoothing operations
    dustFracMin : float
        Minimum dust fraction
    inputSnap : SimSnap
        (optional) The positions for the simulation can be taken from another
        simulation if one is supplied.
    dDeltaFactor : float
        Factor by which to scale time step.  This will keep outputs at the 
        same simulation time (ie will scale iOutInterval, etc...).  This can 
        be useful e.g. to test multistepping
    boxRes : array-like
        Resolution of the box (number of grid points) along each axis
        
    Returns
    -------
    ICs : SimSnap
        Generated ICs
    param : dict
        params used
    simdir : str
        Directory the ICs are saved in
    paramname : str
        Path to the saved .param file
        
    """
    runParams = locals()
    # -----------------------------------------------
    # Load stuff
    # -----------------------------------------------
    if inputSnap is not None:
        # Load inputSnap file
        if isinstance(inputSnap, str):
            inputSnap = pynbody.load(inputSnap)
        else:
            inputSnap = inputSnap
        ngas = len(inputSnap)
    else:
        # Create grid ICs
        boxRes = np.asarray(boxRes).astype(int)
        ngas = np.product(boxRes)
    # Parse param file
    param = testdust.utils.loadDefaultParam(defaultParam, userDefaults)
    # Get units
    units = diskpy.pychanga.units_from_param(param)
    grainSize = param['dDustSize']
    # filenames
    fname = param['achInFile']
    outparam = param['achOutName'] + '.param'
    molecularWeight = param['dMeanMolWeight']
    
    # -----------------------------------------------
    # set up simulation directory
    # -----------------------------------------------
    if not os.path.exists(simdir):
        
        os.mkdir(simdir)
    
    fRunParams = os.path.join(simdir, runParamsName)
    fname = os.path.join(simdir, fname)
    
    # -----------------------------------------------
    
    snap = pynbody.new(gas=ngas)
    
    # Give everything units
    cs = SimArray(cs, units['v_unit'])
    rho = SimArray(rho, units['rho_unit'])
    boxSize = SimArray(boxSize, units['l_unit'])
    criticalRadius = SimArray(criticalRadius, units['l_unit'])
    tStop = SimArray(tStop, units['t_unit'])
    molecularWeight = SimArray(molecularWeight, 'm_p')
    grainSize = SimArray(grainSize, units['l_unit'])
    # Save/update run parameters
    runParams['criticalRadius'] = criticalRadius
    runParams['cs'] = cs
    runParams['tStop'] = tStop
    pickle.dump(runParams, open(fRunParams, 'w'))
    print "Run parameters saved to:", fRunParams
    
    # Isothermal temperature
    T = (cs**2) * molecularWeight / kB
    T.convert_units('K')
    snap['temp'] = T * np.ones(ngas)
    
    # Mass
    m = (rho * boxSize**3)/ngas
    m.convert_units(units['m_unit'])
    snap['mass'] = m
    
    # Grain density
    grainDensity = rho*cs*tStop*np.sqrt(8/np.pi)/grainSize
    grainDensity.convert_units(units['rho_unit'])
    param['dDustGrainDensity'] = float(grainDensity)
        
    # Setup positions
    if inputSnap is not None:
        pos = inputSnap['pos'].view(float, np.ndarray)
        snap['pos'] = boxSize * pos
        Lfloat = float(boxSize)
        param['dPeriod'] = Lfloat
    else:
        boxShape = boxSize * np.ones(3)
        pos = periodicGrid(boxRes, boxShape)
        snap['pos'] = pos
        Lfloat = float(boxSize)
        param['dPeriod'] = Lfloat
    
    # Define dust fraction
    snap['dustFrac'] = dustFrac0 * (1 - (snap['r']/criticalRadius)**2)
    snap['dustFrac'][snap['dustFrac'] < dustFracMin] = dustFracMin
    
    # Estimate a useable softening length.  Make it big so that it does not
    # constrain the timestepping.  Gravitational forces should be off, so this 
    # is okay
    snap['rho']
    snap['eps'] = 100 * snap['smooth'].mean()
    
    # Save out
    snap.write(fmt=pynbody.tipsy.TipsySnap, filename=fname)
    # Set up the param file
    param['nSmooth'] = nSmooth
    scaleTimeStep(param, dDeltaFactor)
    outparamfull = os.path.join(simdir, outparam)
    diskpy.utils.configsave(param, outparamfull, 'param')
    
    print 'ICs saved to ' + fname    
    return snap, param, simdir, outparam
