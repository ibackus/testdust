# -*- coding: utf-8 -*-
"""
NOTES:

It is important to time evolve the gas settling enough time.  It can appear to
correct profile spuriously since large oscillation modes exist and need to be
damped out.  This is why I allow the settling to occur several times, setting
the velocties to zero in between.

UNITS NOTES:

Mass unit = Msol
Distance Unit = AU
gas constant R = 1 (in code units)
G = 1 (in code units)
Central Star Mass = 1 Msol (1 in code units)

Also note that the intrinsic dust density and dust size are in somewhat
arbitrary units, since only the quantity (rho_dust * size_dust) matters.  
Therefore (rho_dust * size_dust) is in code units and I just set rho_dust = 1
arbitrarily


Created on Thu Jul 21 14:06:57 2016

@author: ibackus
"""
from pprint import pformat
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import shutil
import json
import pynbody
SimArray = pynbody.array.SimArray
import diskpy
import sphglass
from diskpy.ICgen.ICgen_utils import changa_command, changa_run

import testdust
import utils

class IC():
    """
    IC(self, settingsFileName)
    
    Class for generating tipsy initial conditions for the dust settling test
    of Price & Laibe (2015) (see section 4.4 Dust settling in a 
    protoplanetary disc)
    
    ICs are generated from settings stored in a python settings file.  For an 
    example settings file, look at testdust/settling/settings.py
    
    Examples
    --------
    
    Note, its good to do this a new folder since all files get saved in the 
    current working directory
    
    >>> import testdust
    >>> IC = testdust.settling.makeICs.IC('settings.py')
    
    To use a glass as a starting point:
    
    >>> IC.makeGlass()
    
    To use a snapshot of your own as a starting point, you can just set 
    initialFileName in settings.py . Or (not recommended) you can set
    IC.settings['initialFileName']
    
    Now settle the ICs:
    
    >>> IC.setupGas()
    >>> IC.settleGas()
    
    You can now check to see if the saved snapshot seems to be in 
    equilibrium.  If not, run the settling again
    
    >>> IC.settleGas()
    
    If everything is okay, go ahead and set up the dust settling ICs
    
    >>> IC.setupDust()
    """
    settings = {}
    
    def __init__(self, settingsFileName):
        
        self.readSettings(settingsFileName)
        
    def __repr__(self):
        """
        """
        outString = "< dust settling ICs instance >\n"
        outString += pformat(self.settings)
        return outString
        
    def readSettings(self, settingsFileName=None):
        """
        """
        if settingsFileName is not None:
            
            self.settingsFileName = settingsFileName
            
        self.settings = testdust.utils.parsePySettings(self.settingsFileName)
        self.initCalc()
        
    def initCalc(self):
        """
        """
        x = self.settings
        H, cs, smooth, boxshape, dDelta = initCalc(x['hOverR'], x['R0'], \
            x['boxwidth'], x['nParticles'], x['nSmooth'])
        self.settings.update({'H':H, 'cs': cs, 'smooth': smooth, 'boxshape': boxshape,
                        'dDelta': dDelta})
        
    def makeGlass(self):
        """
        """
        x = self.settings
        output = makeGlass(x['nParticles'], x['boxshape'], x['changaPreset'], 
                         x['initialFileName'], x['verbose'])
        return output
        
    def setupGas(self):
        """
        """
        x = self.settings
        setupGas(x['initialFileName'], x['H'], x['cs'], x['dDelta'], \
            x['boxshape'], x['rho0'], x['smooth'], x['nSoundCross'], x['nSmooth'])
            
    def settleGas(self, changaPreset=None, nRuns=None):
        """
        """
        x = self.settings
        if changaPreset is None:
            changaPreset = x['changaPreset']
        if nRuns is None:
            nRuns = x['nRuns']
            
        output = settleGas(changaPreset, nRuns)
        return output
        
    def setupDust(self):
        """
        """
        x = self.settings
        output = setupDust(x['dustSize'], x['intrinsicDustRho'], x['rho0'], \
            x['R0'], x['numOrbitsRun'], x['dDelta'], x['cs'], x['nSmooth'])
        
        return output

def makeGlass(nParticles, boxshape, changaPreset='default', 
              glassName='glass.std', verbose=True):
    """
    Generate a glass (to be used for the dust settling test)
    This is basically just a wrapper for sphglass, but saves the output to
    glassName
    """
    glass = sphglass.glassBox(nParticles, boxshape, changaPreset, verbose)
    shutil.move('glass.std', glassName)
    return 

def setupGas(filename, height, cs, dDelta, boxshape, rho0, smooth,
             nSoundCross=4, nSmooth=None):
    """
    Sets up the .param and IC snapshot for performing the gas settling portion
    of IC generation (see settleGas).  All quantities are assumed to be in
    simulation units
    
    Parameters
    ----------
    filename : str
        Input snapshot filename.  The input snapshot can be a glass or some
        other ICs
    height : float
        Scale height (usuall defined as cs/omega)
    dDelta : float
        Time step-size
    boxshape : array-like
        Shape of the box to run in (x, y, z)
    rho0 : float
        Density normalization (code units)
    smooth : float
        Approximate mean smoothing length
    nSoundCross : int
        Number of sound crossing times to run for.
    nSmooth : int
        (optional) Number of neighbors for smoothing
    """
    # Set up time stepping
    nSteps = nSoundCross * (2.*height/cs) /float(dDelta)
    nSteps = int(np.round(nSteps))
    # Load gas param defaults
    gasparam = utils.loadDefaultParam('gas')
    if nSmooth is not None:
        gasparam['nSmooth'] = nSmooth
    # Get filenames
    savename, gasParamName, gasPrefix = utils.defaultFilenames('gas')
    # Load the input file (either a glass or something else)
    snap = pynbody.load(filename)
    nParticles = len(snap)
    # Set up the param
    gasparam['dxPeriod'] = boxshape[0]
    gasparam['dyPeriod'] = boxshape[1]
    gasparam['dzPeriod'] = 50 * height
    gasparam['dDelta'] = dDelta
    gasparam['nSteps'] = nSteps
    gasparam['iOutInterval'] = nSteps + 1
    gasparam['iCheckInterval'] = nSteps + 1
    # Set up ICs
    snap['temp'] = cs**2
    boxVolume = boxshape[0] * boxshape[1] * height
    totalMass = rho0 * np.sqrt(2*np.pi) * boxVolume
    snap['mass'] = totalMass/nParticles
    snap['eps'] = smooth*0.5
    # Save snapshot and param file
    diskpy.utils.configsave(gasparam, gasParamName, 'param')
    print ".param file saved to:", gasParamName
    snap.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
    print "snapshot saved to:", savename
    
def settleGas(changaPreset='default', nRuns=2):
    """
    Perform the gas settling portion of the IC generation, using ChaNGa.
    """    
    # Get filenames
    savename, gasParamName, gasPrefix = utils.defaultFilenames('gas')
    # Load the generated param file
    gasparam = diskpy.utils.configparser(gasParamName, 'param')
    
    cmd = changa_command(gasParamName, changaPreset)
    # Now run for more sound crossing times
    outfile = '{0}.{1:06}'.format(gasparam['achOutName'], gasparam['nSteps'])
    for i in range(nRuns):
        # Run
        print "Running ChaNGa with command:", cmd
        changa_run(cmd, require_success=True)
        # Zero out velocities
        snap = pynbody.load(outfile, paramfile=gasParamName)
        snap['vel'] *= 0.
        snap['pos']
        snap['rho']
        # Just to make sure ALL arrays get saved, access each key
        for key in snap.g.loadable_keys():
            snap.g[key]
        snap.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
        
        
def setupDust(dustSize, intrinsicDustRho, rho0, R0, numOrbitsRun, dDelta, cs,
              nSmooth=None):
    """
    Set up the dust .param and IC snapshot files. (after settleGas has 
    been performed)
    
    Parameters
    ----------
    dustSize : float
        Size of the dust (in mm)
    rho0 : float
        Density normalization
    R0 : float
        Distance from star
    numOrbitsRun : float
        Number of orbits to run simulation for
    dDelta : float
        Step size
    cs : float
        Sound speed (in code units)
    nSmooth : int
        (optional) number of neighbors for smoothing
    """
    # Load, using the paramfile for the settled snapshot (gas one)
    dustparam = utils.loadDefaultParam('dust')
    # Get filenames
    filename, inParamName, inPrefix = utils.defaultFilenames('gas')
    savename, dustParamName, dustPrefix = utils.defaultFilenames('dust')
    # Load the output of settleGas
    dustSnap = pynbody.load(filename, paramfile=inParamName)
    # Set up dust params    
    # Set dust fraction such that rho_dust/rho_gas = 1/100
    dustSnap['dustFrac'] = 0.01/(1+0.01)
    # Set up the intrinsic dust size & density
    gamma = dustparam['dMeanMolWeight']
    intrinsicDustRhoCode = 1.
    dustSizeCode = (1.35e-3) * dustSize * (intrinsicDustRho/3.) * rho0 \
    * cs * (R0**1.5) * np.sqrt(8./(np.pi * gamma))/intrinsicDustRhoCode
    dustparam['dDustParam1'] = dustSizeCode
    dustparam['dDustParam2'] = intrinsicDustRhoCode
    # Time stepping
    period = 2*np.pi*R0**1.5
    tRun = period * numOrbitsRun
    dustparam['dDelta'] = dDelta/10.
    dustparam['nSteps'] = np.round(tRun/dustparam['dDelta']).astype(int)
    dustparam['iOutInterval'] = np.round(period/dustparam['dDelta']).astype(int)
    if nSmooth is not None:
        dustparam['nSmooth'] = nSmooth
    
    # ECHO
    tStop = intrinsicDustRhoCode*dustSizeCode*np.sqrt(np.pi*gamma/8.)
    tStop /= rho0 * cs
    print 'stopping time:', tStop
    period = 2*np.pi*R0**1.5
    print '(tStop/tOrbit):', tStop/period
    
    # Save
    diskpy.utils.configsave(dustparam, dustParamName, 'param')
    try:
        # For some reason it helps to access some quantities first to make sure
        # they get loaded
        dustSnap['pos']
        dustSnap['vel']
        dustSnap.write(filename=savename)
    except ValueError:
        dustSnap.write(filename=savename)
    print 'snapshot saved to:', savename
    
def initCalc(hOverR, R0, boxwidth, nParticles, nSmooth=32):
    """
    Initial derivations from settings
    
    Parameters
    ----------
    hOverR  : float
        Ratio of scale height to distance from star
    RO : float
        Distance to star
    boxwidth : float
        Width of the simulation box.  Height will be set as a large number
    nParticles : int
        Number of particles
    nSmooth : int
        Number of neighbors to smooth over for the simulation
    """
    # Initial derivations
    H = hOverR * R0
    cs = hOverR/np.sqrt(R0)
    smooth = (float(nSmooth) * (H*boxwidth**2)/nParticles)**(1./3) # approx smooth length
    boxshape = [boxwidth, boxwidth, 2*H]
    # Time stepping for gas settling
    tCourant = smooth/cs
    dDelta = tCourant
    
    return H, cs, smooth, boxshape, dDelta

if __name__ == '__main__':
    
    pass
#    
#    if do_compare:
#        
#        # Expected values
#        z = np.linspace(-5*H, 5*H, 1000)
#        expected = rho0 * np.exp(-z**2/(2*H**2))
#        
#        fname = 'settledust.std'
#        f = pynbody.load(fname)
#        
#        plt.plot(f['z'], f['rho'],',', z, expected, 'r')
#        plt.legend(['actual', 'expected'])
#        plt.show()
