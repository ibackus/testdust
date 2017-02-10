# -*- coding: utf-8 -*-
"""
This file contains settings required to make ICs for the dust settling test
of Price & Laibe (2015) (see section 4.4 Dust settling in a protoplanetary disc)

Instead of starting with a stretched hexagonal lattice, we start with randomly
placed particles, normally distributed along z.

This file is parsed by testdust.utils.parsePySettings() by running execfile()
on it.  Imports can be included in this file, but all modules will be ignored.

Also, before running the generated dust ICs, you may want to adjust dDelta

THIS IS A 2D TEST.  It requires ChaNGa be compiled with NDSPH=2 (2 dimensions)

Code units are specified in PL15
"""
# --------------------------------------
# Run settings
# --------------------------------------
nParticles = 16 * 56 # Low resolution test
#nParticles = 32 * 111 # Medium resolution test
#nParticles = 64 * 222 # High resolution test

ndim = 2
boxwidth = 0.5
hOverR = 0.05
R0 = 5. # code units
rho0 = 1e-3 # code units
intrinsicDustRho = 3. # in g/cm^3
dustSize = 1. # in mm
numOrbitsRun = 50.
dustFrac = 0.099

# --------------------------------------
# General settings 
# --------------------------------------
# Generate a stretched cubic grid as a starting point
initialSnapKind = 'random' 
initialFileName = 'random.std'
# See diskpy for info on changa presets
changaPreset = 'default' 
verbose = False

# --------------------------------------
# Algorithm settings
# --------------------------------------
#damping = None # automatically set damping force
damping = 0.03 # price and laibe setting

# Sound crossing times to settle gas for
# Price and Laibe use 1000, but you can get away with less, probably
#tGasSettle = 200 # low res
#tGasSettle = 400 # medium res
#tGasSettle = 800 # high res
tGasSettle = 1000 # price and laibe 2015

# For all the runs, how many neighbors to use for smoothing
nSmooth = 28 # PL15 use an "effective" neighbor number of about 28.3
