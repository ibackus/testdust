# -*- coding: utf-8 -*-
"""
This file contains settings required to make ICs for the dust settling test
of Price & Laibe (2015) (see section 4.4 Dust settling in a protoplanetary disc)
"""
initialFileName = 'glass.std'
boxres = [8, 8, 29]
#boxres = [64, 228, 64]
#nParticles = 1792
nParticles = boxres[0] * boxres[1] * boxres[2]
boxwidth = 0.5
changaPreset = 'default'
verbose = False
nreglass = 3
hOverR = 0.05
R0 = 5.
rho0 = 1e-3
#nSoundCross = 4 # Sound crossing times to settle gas for
nSoundCross = 1
#nRuns = 2 # how many runs to settle for (zeroing velocity each time)
nRuns = 1
intrinsicDustRho = 3. # in g/cm^3
dustSize = 1. # in mm
numOrbitsRun = 50.
nSmooth = 16
