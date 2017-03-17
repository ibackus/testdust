#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
THIS TEST REQUIRES PARTICLE POSITIONS TO BE FIXED IN PLACE.  Currently,
I have this implemented in ChaNGa as a compile-time flag, which can be
selected by adding
    --enable-testdust=diffusion
to the configure script command line options

Created on Thu Mar 16 16:40:20 2017

@author: ibackus
"""

from testdust import diffusion

import os
import diskpy
import matplotlib.pyplot as plt

# Lets do a low resolution test.  The paper values are [50, 58, 60]
res = [12, 14, 15]
# Set the directory to save the simulatin to
simdir = "pl15_lores"

# --------------------------
# Make the ICs
# --------------------------
print 'Making ICs...'
# Note that I scale dDelta since I'm running at lower resolution and the
# default dDelta is too low.
f, param, simdir, paramname = diffusion.makeICs.makeICs(boxRes=res, \
    simdir=simdir, dDeltaFactor=4.)

# --------------------------
# Run the ICs.
# --------------------------
# NOTE: SLOW! this part takes about 3-4 minutes on a 2014 era laptop.
# You should only really need 1 core to run this
# If you don't have diskpy's ChaNGa presets setup, you can just write your
# own command
# OR if you prefer, run ChaNGa on your own!
print 'Running ICs...'
preset = 'default'
cwd = os.getcwd()
try:
    os.chdir(simdir)
    # First lets run for 0 time steps to plot the ICs properly
    cmd = diskpy.pychanga.changa_command(paramname, preset=preset, 
                                         changa_args='-n 0')
    diskpy.pychanga.changa_run(cmd, require_success=True)
    # Now lets run the simulation
    cmd = diskpy.pychanga.changa_command(paramname, preset=preset)
    diskpy.pychanga.changa_run(cmd, require_success=True)
finally:
    # Go back to the original directory
    os.chdir(cwd)

# --------------------------
# Analyze the results
# --------------------------
print 'Plotting results...'
sim, eps, t, runpars = diffusion.analyze.loadSim(simdir)
plt.clf()
diffusion.plot.dustFracProfile(sim, t, eps)
plt.savefig('diffusion_profile_lores.png')
plt.figure()
diffusion.plot.crossSection(sim, t)
plt.savefig('diffusion_crossection_lores.png')

plt.show()
