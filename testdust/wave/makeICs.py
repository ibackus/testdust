#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 14:26:23 2017

@author: ibackus
"""
import cPickle as pickle
from scipy.interpolate import interp1d
import numpy as np
import pynbody
import os
SimArray = pynbody.array.SimArray
import diskpy
import testdust

kB = SimArray(1.0, 'k')
defaultParam1D = 'wave_defaults.param'
userDefaults1D = 'user_wave_defaults.param'
argname = 'args.p'

_directory = os.path.dirname(os.path.realpath(__file__))
defaultParam1D = os.path.join(_directory, defaultParam1D)
userDefaults1D = os.path.join(_directory, userDefaults1D)

def wavePos(N, drho, L=1., nPeriods=1.):
    """
    Generates particle positions for a density wave, centered on the origin.
    """
    k = 2*np.pi*nPeriods/L
    # Generate a CDF for z from 0 to 1, then scale it
    z = L * np.linspace(-0.5, 0.5, 10000)
    cdf = z - (drho/k) * np.cos(k*z)
    cdf -= cdf[0]
    cdf /= cdf[-1]
    cdfinv = interp1d(cdf, z, kind='linear')
    dm = 1.0/N
    m = np.linspace(dm/2, 1-dm/2, N)
    z = cdfinv(m)
    if pynbody.units.has_units(L):
        z = SimArray(z, L.units)
    return z

def wave1d(rho0=1., dustFrac0=0.5, drho=1e-4, v0=1e-4, gamma=5./3, cs=1., 
            ngas=100, L=1., tStop=5e-4, nSmooth=5, nPeriods=1):
    r"""
    Generate a 1D dusty-gas wave as in Price & Laibe 2015.
    I've only gotten this to work with the M4Kernel and nSmooth = 4 or 5.
    Other kernels or nSmooth create noisy garbage
    
    Default params are stored in 'wave_defaults.param' and can be overridden
    in 'user_wave_defaults.param'
    
    This applies perturbations like:
        
        .. math::
            \rho(z) &= \rho_0 (1 + \delta\rho\sin{(2 \pi N z/L)}) \\
            v(z) &= v_0 \sin{(2 \pi N z/L)} \\
            P(z) &= P_0 + c_s v_0 \rho_0 \sin{(2 \pi N z/L)} 
                [1 + (\delta\rho/2) \sin{(2 \pi N z/L)}]
        
    Which makes a wave moving to the right at the speed of sound :math:`c_s`.
    The pressure perturbation is derived from the density and velocity 
    perturbations using the momentum equation and assuming sinewaves travelling
    to the right at cs
    
    Parameters
    ----------
    All parameters listed are in code units!
    
    rho0: float
        Base total density. :math:`\rho_0`
    dustFrac: float
        Dust fraction of the gas.
    drho: float
        Fractional amplitude of density perturbation
    v0: float
        Velocity perturbation amplitude
    gamma: float
        Adiabatic index
    cs: float
        Sound speed
    ngas: int
        Number of particles
    L: float
        Length of the periodic box
    tStop: float
        Stopping time of the dust
    nSmooth: int
        Number of neighbors for smoothing
    nPeriods: float or int
        Number of wave periods to fit in the box
    
    Returns
    -------
    snapshot: pynbody SimSnap
        ICs for the simulation.  they will also be saved to disk
    param: dict
        Runtime params.  Also saved to a .param file
    paramname: str
        Path to the saved param file
    """
    args = locals()
    nDim = 1
    pickle.dump(args, open(argname, 'w'), 2)
    print 'Saved args to:', argname
    
    param = testdust.utils.loadDefaultParam(defaultParam1D, userDefaults1D)
    
    # -------------------------------------------
    # Parse input param file
    # -------------------------------------------
    # Get units
    units = diskpy.pychanga.units_from_param(param)
    # filenames
    fname = param['achInFile']
    outparam = param['achOutName'] + '.param'
    molecularWeight = param['dMeanMolWeight']
    
    # -------------------------------------------
    # Set-up
    # -------------------------------------------
    mass = rho0/(L * float(ngas))
    # Give everything units
    v0 = SimArray(v0, units['v_unit'])
    cs = SimArray(cs, units['v_unit'])
    rho0 = SimArray(rho0, units['rho_unit'])
    L = SimArray(L, units['l_unit'])
    mass = SimArray(mass, units['m_unit'])
    tStop = SimArray(tStop, units['t_unit'])
    molecularWeight = SimArray(molecularWeight, 'm_p')
    # Initialize snapshot
    snap = pynbody.new(gas=ngas)
    # Initialize arrays
    for key in ('temp', 'rho', 'dustFrac', 'pos', 'mass'):
        
        snap[key] = 0.
    
    k = 2 * np.pi * nPeriods/L
    # -------------------------------------------
    # Positions and mass
    # -------------------------------------------
    z = wavePos(ngas, drho, L, nPeriods)
    snap['z'] = z
    snap['mass'] = mass
    # -------------------------------------------
    # Temperature
    # -------------------------------------------
    P0 = rho0 * (cs**2)/gamma
    P = P0 + cs * v0 * rho0 * np.sin(k * z) * (1 + 0.5 * drho * np.sin(k * z))
    rho = rho0 * (1 + drho * np.sin(k * z))
    T = molecularWeight * P/(kB * (1-dustFrac0) * rho)
    T.convert_units('K')
    print 'Temp:', T
    snap['temp'] = T
    # -------------------------------------------
    # Velocity profile
    # -------------------------------------------
    vz = v0 * np.sin(2*np.pi*snap['z']*nPeriods/L)
    v = SimArray(np.zeros([ngas, 3]), vz.units)
    v[:,2] = vz
    snap['vel'] = v
    # -------------------------------------------
    # Dust stuff
    # -------------------------------------------
    # Dust mass
    grainDensity = rho0*cs*tStop*np.sqrt(8/np.pi/gamma)/units['l_unit']
    grainDensity.convert_units(units['rho_unit'])
    param['dDustSize'] = 1.
    param['dDustGrainDensity'] = float(grainDensity)
    # Define dust fraction
    snap['dustFrac'] = dustFrac0 * np.ones(ngas)
    # Estimate density and eps
    snap['rho'] 
    snap['eps'] = 100 * L/ngas**(nDim/3.)
    
    # -------------------------------------------
    # Write out
    # -------------------------------------------
    testdust.utils.setupParamBounds(param, [float(L)])
    param['nSmooth'] = nSmooth
    param['dConstGamma'] = gamma
    # Write out snapshot and paramfile
    snap.write(fmt=pynbody.tipsy.TipsySnap, filename=fname, double_pos=True, \
    double_vel=True)
    print 'Saved snapshot to:', fname
    diskpy.utils.configsave(param, outparam, 'param')
    print 'Saved .param to:', outparam
    
    return snap, param, outparam