# -*- coding: utf-8 -*-
"""
This file contains settings required to make ICs for the dust settling test
of Price & Laibe (2015) (see section 4.4 Dust settling in a protoplanetary disc)
This uses a stretched cubic grid with a gaussian vertical density profile as
a starting point (similar to Price & Laibe's hexagonal grid).  This method
for generating ICs the same as in PL15.

This file is parsed by testdust.utils.parsePySettings() by running execfile()
on it.  Imports can be included in this file, but all modules will be ignored.

Also, before running the generated dust ICs, you may want to adjust dDelta

THIS IS A 2D TEST.  It requires ChaNGa be compiled with NDSPH=2 (2 dimensions)

Code units are specified in PL15
"""
# --------------------------------------
# Run settings
# --------------------------------------
boxres = [16, 56] # Low resolution test
#boxres = [32, 111] # Medium resolution test
#boxres = [64, 222] # High resolution test

nParticles = boxres[0] * boxres[1] # Must be specified
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
initialSnapKind = 'hexagonal' 
initialFileName = 'hexgrid.std'
# See diskpy for info on changa presets
changaPreset = 'default' 
verbose = False

# --------------------------------------
# Algorithm settings
# --------------------------------------
# Sound crossing times to settle gas for
# Price and Laibe use 1000, but you can get away with less, probably
tGasSettle = 1000
# how many runs to settle for (zeroing velocity each time)
nRuns = 1
# For all the runs, how many neighbors to use for smoothing
nSmooth = 28 # PL15 use an "effective" neighbor number of about 28.3
