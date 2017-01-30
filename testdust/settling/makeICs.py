# -*- coding: utf-8 -*-
"""
This module contains functions and classes for generating initial conditions
for the dust settling test of Price & Laibe (2015) (see section 4.4 Dust 
settling in a protoplanetary disc)

The ICs are meant to mimic the vertically-stratified profile of a section of
a protoplanetary disk at some distance R from the star.  

The ICs are in a periodic box which is very extended in z to mimic open
boundary conditions along z.  The box is much longer along z than x-y.  The
vertical density profile is a gaussian of:
    rho(z) = rho0 * exp(-z^2/2 H^2)

IC generation is handled by the IC class.  There are 3 steps in IC generation:
    1)  Generate initial, starting point ICs.  This can a uniform SPH glass,
        or a grid, etc...
    2)  Gas settling.  Before adding dust, the ICs are run to allow them to
        settle into an equilibrium state.  velocities are zeroed out
    3)  Set-up dust parameters

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
import shutil
import itertools
import os
import matplotlib.pyplot as plt

import pynbody
SimArray = pynbody.array.SimArray

import diskpy

from diskpy.ICgen.ICgen_utils import changa_command, changa_run

import testdust
import utils
import ppdgrid

# Setup constants used here
molecularWeight = SimArray(2.0, 'm_p')
kB = SimArray(1.0, 'k')
G = SimArray(1.0, 'G')
Msol = SimArray(1.0, 'Msol')

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
    
    There are 2 ways to generate the initial snapshot.  You can call
    makeInitialSnap (which is wrapper) or call one the makexxx functions:
    
    >>> # Use whatever is specified by settings.py
    >>> IC.makeInitialSnap()
    >>> # OR make a glass
    >>> IC.makeGlass()
    >>> # OR make a 'cubic' grid
    >>> IC.makeGrid()
    
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
        Parses the pythonic settings in settingsFileName 
        See testdust.utils.parsePySettings
        """
        if settingsFileName is not None:
            
            self.settingsFileName = settingsFileName
            
        self.settings = testdust.utils.parsePySettings(self.settingsFileName)
        self.initCalc()
        
    def initCalc(self):
        """
        Perform initial calculations/derivations from the settings
        """
        kind = self.settings['initialSnapKind']
        if kind == 'hexagonal':
            
            self.settings['ndim'] = 2
            
        elif kind == 'line':
            
            self.settings['ndim'] = 1
            
        else:
            
            self.settings['ndim'] = 3
        
        x = self.settings        
        H, cs, smooth, boxshape, dDelta = initCalc(x['hOverR'], x['R0'], \
            x['boxwidth'], x['nParticles'], x['nSmooth'], x['ndim'])
        self.settings.update({'H':H, 'cs': cs, 'smooth': smooth, 'boxshape': boxshape,
                        'dDelta': dDelta})
        
    def makeInitialSnap(self):
        """
        Generate a starting-point, initial snapshot.  A wrapper function
        """
        kind = self.settings['initialSnapKind']
        if kind == 'cubic':
            print 'Making cubic ICs'
            self._makeGrid()
        elif kind == 'glass':
            print 'making glass'
            self._makeGlass()
        elif kind == 'hexagonal':
            print 'Making hexagonal 2D ICs'
            self._makeGrid()
        elif kind == 'line':
            print 'Making line (1D) ICs'
            self._makeGrid()
        else:
            raise ValueError, "Unrecognized initial snapshot kind {0}"\
                .format(kind)
        
    def _makeGlass(self):
        """
        Generate a glass starting point, initial snapshot using 
        makeICs.makeGlass()
        
        Also depends on sphglass
        """
        x = self.settings
        snap = makeGlass(x['nParticles'], x['boxshape'], x['changaPreset'], 
                         x['initialFileName'], x['verbose'])
        return snap
        
    def _makeGrid(self):
        """
        Generate a cubic 'grid' as a starting point, initial snapshot
        """
        x = self.settings
        snap = makeGrid(x['boxres'], x['boxshape'], x['H'], x['R0'], \
            x['initialFileName'], x['dustFrac'], x['rho0'], kind=x['initialSnapKind'])
        return snap
        
    def setupGas(self):
        """
        Set-up the snapshot and .param file used to perform the gas settling.
        See settleGas()
        """
        x = self.settings
        kind = x['initialSnapKind'].lower()
        if kind == 'cubic' or kind == 'hexagonal':
            updateMass = False
        else:
            updateMass = True
            
        setupGas(x['initialFileName'], x['H'], x['cs'], x['dDelta'], \
            x['boxshape'], x['rho0'], x['smooth'], \
            x['nSmooth'], updateMass, x['ndim'], x['changaPreset'], \
             x['tGasSettle'])
            
    def settleGas(self, changaPreset=None, nRuns=None):
        """
        Perform the gas settling portion of the IC generation.  During this
        stage the ICs are time evolved with just gas  (no dust) to allow them
        to settle to an equilibrium position
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
        Set up the dust ICs
        
        After running settleGas(), the dust parameters can be set-up and 
        the final IC generation completed.
        """
        x = self.settings
        output = setupDust(x['dustSize'], x['intrinsicDustRho'], x['rho0'], \
            x['R0'], x['numOrbitsRun'], x['dDelta'], x['cs'], x['nSmooth'])
        
        return output

def hexGrid(boxres, boxshape, H, R0, eps=0.):
    """
    Generate a 2D stretched hexagonal grid in the y-z plane.
    """
    nParticles = boxres[0] * boxres[1]
    y, m = ppdgrid.periodicHexMesh(boxres[0], boxres[1])
    y *= boxshape[1]
    z = ppdgrid.posGen(nParticles, H, R0, eps, m=m)
    pos = np.zeros([nParticles, 3])
    pos[:, 1] = y.flatten()
    pos[:, 2] = z.flatten()
    # Get mass scale required to generate density of rho0 at the origin
    mScales = [ppdgrid.particleMass(z1, rho0=1.) for z1 in z]
    mScale = np.mean(mScales)
    return pos, mScale

def makeGlass(nParticles, boxshape, changaPreset='default', 
              glassName='glass.std', verbose=True):
    """
    Generate a glass (to be used for the dust settling test)
    This is basically just a wrapper for sphglass, but saves the output to
    glassName
    """
    import sphglass
    glass = sphglass.glassBox(nParticles, boxshape, changaPreset, verbose)
    shutil.move('glass.std', glassName)
    return glass
    
def lineGrid(nParticles, H, R0, eps=0.):
    """
    """
    if hasattr(nParticles, '__iter__'):
        nParticles = nParticles[0]
    z = utils.gaussianGrid(nParticles)
    mScale = ppdgrid.posGen(nParticles, H, R0, eps)
    return z, mScale
    
def cubicGrid(boxres, boxshape, H, R0, eps=0.):
    """
    Generates a grid of positions with a gaussian density profile along the 
    z-axis
    """
    nx, ny, nz = boxres
    if len(boxres) != 3:
        raise ValueError, 'boxres must be len(3)'
#    if nz % 2 == 0:
#        raise ValueError, "number of z particles must be odd (for symmetry)"
    x = testdust.utils.periodicLine(nx) * boxshape[0]
    y = testdust.utils.periodicLine(ny) * boxshape[1]
    z = ppdgrid.posGen(nz, H, R0, eps)
    mScale = ppdgrid.particleMass(z, rho0=1.)
    # Now make grid out of them
    grid = np.array([a for a in itertools.product(x, y, z)])
    
    return grid, mScale
    
def makeGrid(boxres, boxshape, H, R0, savename, eps=0., rho0=1., kind='cubic'):
    """
    Generate a tipsy snapshot of gas particles on a grid with an gaussian
    density profile along the z-axis
    """
    nParticles = np.product(boxres)
    # Set up snapshot defaults
    f = pynbody.snapshot.new(gas=nParticles)
    f['temp'] = 1.
    f['vel'] = 0.
    f['eps'] = 1.
    f['rho'] = 0.
    # Generate particle positions
    print 'making {0} grid'.format(kind)
    if kind == 'hexagonal':
        pos, mScale = hexGrid(boxres, boxshape, H, R0, eps)
    elif kind == 'cubic':
        pos, mScale = cubicGrid(boxres, boxshape, H, R0, eps)
    elif kind == 'line':
        z, mScale = lineGrid(boxres, H, R0, eps)
        pos = np.zeros([len(z), 3])
        pos[:,2] = z
    else:
        raise ValueError, 'Unknown grid kind {0}'.format(kind)
    f['pos'] = pos
    f['mass'] = mScale * rho0
    # Now save
    f.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
    print 'snapshot saved to:', savename
    return f

def setupGas(filename, height, cs, dDelta, boxshape, rho0, smooth,
             nSmooth=None, updateMass=True, ndim=3, changaPreset='default', 
             tGasSettle=1000):
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
    nSmooth : int
        (optional) Number of neighbors for smoothing
    updateMass : bool
        (optional) If true, set total mass to be rho0 * box volume
    ndim : int
        Number of dimensions this is meant to be run in.
    tGasSettle : float
        Simulation time (in code units) to run the settling for.  Default is
        the same as in Price & Laibe 2015
    """
    # Set up time stepping
    nSteps = int(np.round(tGasSettle/float(dDelta)))
    assert(nSteps > 0)
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
    gasparam['dzPeriod'] = 50 * height
    gasparam['dxPeriod'] = boxshape[0]
    gasparam['dyPeriod'] = boxshape[1]
    if ndim == 1:
        gasparam['dxPeriod'] = gasparam['dzPeriod'] * 10
        gasparam['dyPeriod'] = gasparam['dzPeriod'] * 10
    elif ndim == 2:
        gasparam['dxPeriod'] = gasparam['dzPeriod'] * 10
    
    gasparam['dDelta'] = dDelta
    gasparam['nSteps'] = nSteps
    gasparam['iOutInterval'] = nSteps + 1
    gasparam['iCheckInterval'] = nSteps + 1
    # Set up ICs
    units = diskpy.pychanga.units_from_param(gasparam)
    v_unit = units['l_unit']/units['t_unit']
    csSimArray = SimArray(cs, v_unit)
    T = (molecularWeight * csSimArray**2)/kB
    T.convert_units('K')
    snap['temp'] = T
    if updateMass:
        print "updating particle masses"
        boxVolume = boxshape[0] * boxshape[1] * height
        totalMass = rho0 * np.sqrt(2*np.pi) * boxVolume
        snap['mass'] = totalMass/nParticles
        
    snap['eps'] = smooth*0.5
    # Save snapshot and param file
    diskpy.utils.configsave(gasparam, gasParamName, 'param')
    print ".param file saved to:", gasParamName
    snap.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
    print "snapshot saved to:", savename
    
    # Run ChaNGa for zero steps to get density
    cmd = changa_command(gasParamName, changaPreset, changa_args='-n 0')
    changa_run(cmd, require_success=True)
    os.system('mv {0}.000000 {1}'.format(gasPrefix, savename))
    snap = pynbody.load(savename, paramfile=gasParamName)
    ind = abs(snap['z']).argsort()[0:30]
    midplaneDensity = snap['rho'][ind].mean()
    denUnit = units['m_unit']/units['l_unit']**3
    rho0physical = SimArray(rho0, denUnit)
    mScale = (rho0physical/midplaneDensity).in_units('1')
    print 'scaling mass by', mScale
    snap['mass'] *= mScale
    snap['rho'] *= mScale
    snap.write()
    
    return gasParamName
    
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
    plt.figure()
    plt.ion()
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
            try:
                snap.g[key]
            except ValueError:
                pass
        snap.write(filename=savename, fmt=pynbody.tipsy.TipsySnap)
        plt.clf()
        plt.plot(snap['z'], snap['rho'], '.')
        plt.xlabel('z')
        plt.ylabel(r'$\rho$')
        plt.title('Gas density after settling')
        plt.show(block=False)
        plt.pause(0.05)
        
        
        
def setupDust(dustSize, intrinsicDustRho, rho0, R0, numOrbitsRun, dDelta, cs,
              nSmooth=None, dustFrac = 0.0099):
    """
    Set up the dust .param and IC snapshot files. (after settleGas has 
    been performed)
    
    Parameters
    ----------
    dustSize : float
        Size of the dust (mm)
    intrinsicDustRho : float
        Intrinsic density of the dust (g/cm^3)
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
    dustFrac : float
        Dust fraction to use
    """
    # Load, using the paramfile for the settled snapshot (gas one)
    gasParamName = utils.defaultFilenames('gas')[1]
    dustparam = diskpy.utils.configparser(gasParamName, 'param')
    dustparam.update(utils.loadDefaultParam('dust', False))
    # Don't do glass damping when running dust!
    if 'dGlassDamper' in dustparam:
        
        dustparam.pop('dGlassDamper', None)
        
    units = diskpy.pychanga.units_from_param(dustparam)
    # Get filenames
    filename, inParamName, inPrefix = utils.defaultFilenames('gas')
    savename, dustParamName, dustPrefix = utils.defaultFilenames('dust')
    # Load the output of settleGas
    dustSnap = pynbody.load(filename, paramfile=inParamName)
    # Set up dust params    
    dustSnap['dustFrac'] = dustFrac
    # Set up the intrinsic dust size & density
    gamma = dustparam['dConstGamma']
    dustSizeCode = SimArray(dustSize, 'mm')
    dustSizeCode.convert_units(units['l_unit'])
    dustparam['dDustSize'] = float(dustSizeCode)
    intrinsicDustRhoCode = SimArray(intrinsicDustRho, 'g cm**-3')
    intrinsicDustRhoCode.convert_units(units['rho_unit'])
    dustparam['dDustGrainDensity'] = float(intrinsicDustRhoCode)
    # Time stepping
    period = 2*np.pi*R0**1.5
    tRun = period * numOrbitsRun
    dustparam['dDelta'] = dDelta
    dustparam['nSteps'] = np.round(tRun/dustparam['dDelta']).astype(int)
    dustparam['iOutInterval'] = np.round(period/dustparam['dDelta']).astype(int)
    if nSmooth is not None:
        dustparam['nSmooth'] = nSmooth
    
    # ECHO
    cs = SimArray(cs, units['v_unit'])
    rho0 = SimArray(rho0, units['rho_unit'])
    tStop = intrinsicDustRhoCode*dustSizeCode*np.sqrt(np.pi*gamma/8.)
    tStop /= rho0 * cs
    R0 = SimArray(R0, units['l_unit'])
    # ASSUME 1-SOLAR MASS STAR
    omega = np.sqrt(G * Msol/R0**3)
    print 't_stop * omega', (omega*tStop).in_units('1')
    
    # Save
    diskpy.utils.configsave(dustparam, dustParamName, 'param')
    try:
        # For some reason it helps to access some quantities first to make sure
        # they get loaded
        dustSnap['pos']
        # Also, just to be safe, set velocities to zero
        dustSnap['vel'] = 0.
        dustSnap.write(filename=savename)
    except ValueError:
        # Also, just to be safe, set velocities to zero
        dustSnap['vel'] = 0.
        dustSnap.write(filename=savename)
    print 'snapshot saved to:', savename
    
def initCalc(hOverR, R0, boxwidth, nParticles, nSmooth=32, ndim=3):
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
    ndim : int
        Number of dimensions running in
    """
    # Initial derivations
    H = hOverR * R0
    cs = hOverR/np.sqrt(R0) # in code units
    volume = H * boxwidth**(ndim-1)
    smooth = (float(nSmooth) * (volume)/nParticles)**(1./ndim) # approx smooth length
    boxshape = [boxwidth, boxwidth, 2*H]
    # Time stepping for gas settling
    tCourant = smooth/cs
    dDelta = tCourant/10.
    
    return H, cs, smooth, boxshape, dDelta
