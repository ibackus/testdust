# -*- coding: utf-8 -*-
"""
This file contains settings required to make ICs for the dust settling test
of Price & Laibe (2015) (see section 4.4 Dust settling in a protoplanetary disc)
This uses a glass as a starting point, in contrast to their hexagonal stretched
grid

These mimic the lowest resolution that they use.  We use 57 particles along the
z-direction instead of 56 since our gridding method requires a particle at
the origin

This file is parsed by testdust.utils.parsePySettings() by running execfile()
on it.  Imports can be included in this file, but all modules will be ignored.
"""
# --------------------------------------
# Run settings
# --------------------------------------
nParticles = 16 * 16 * 57
boxwidth = 0.5
hOverR = 0.05
R0 = 5.
rho0 = 1e-3
intrinsicDustRho = 3. # in g/cm^3
dustSize = 1. # in mm
numOrbitsRun = 50.

# --------------------------------------
# General settings 
# --------------------------------------
# Generate an SPH glass as a starting point
initialSnapKind = 'glass' 
initialFileName = 'glass.std'
# See diskpy for info on changa presets
changaPreset = 'default' 
verbose = False

# --------------------------------------
# Algorithm settings
# --------------------------------------
# Sound crossing times to settle gas for
nSoundCross = 4 
# how many runs to settle for (zeroing velocity each time)
nRuns = 2 
# For all the runs, how many neighbors to use for smoothing
nSmooth = 32
