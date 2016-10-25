# -*- coding: utf-8 -*-
"""
This example will generate ICs for the dust settling test of Price & Laibe 
(2015) (see section 4.4 Dust settling in a protoplanetary disc)

This is for the lowest resolution version of their test and will use a 
glass OR a cubic grid (instead of a hexagonal grid) as a starting point.  You
can select which one below by commenting/uncommenting the relevant sections.

This will be SLOW!  The long part of the procedure is settling the gas 
(IC.settleGas()).  This takes a long time because the simulation is run for
several dynamical times with ChaNGa in order to reach equilibrium before
adding the dust.

If you stop or things crash, you can load the ICs using
IC = testdust.settling.makeICs.IC(settingsfile)
and restart from where you left off.

Created on Fri Oct 21 17:07:56 2016

@author: ibackus
"""
import testdust
import os

# -----------------------------------
## Generate a glass 
#savedir = 'price_laibe_glass'
#settingsfile = 'price_laibe_glass.py'
# OR use a cubic grid
savedir = 'price_laibe_cubic'
settingsfile = 'price_laibe_cubic_grid.py'
# -----------------------------------

# Initialize the ICs with the settings
IC = testdust.settling.makeICs.IC(settingsfile)
# Its a good idea to generate these ICs in their own folder
if not os.path.exists(savedir):
    os.mkdir(savedir)

os.chdir(savedir)
# Generate initial snapshot
IC.makeInitialSnap()
# Settle the ICs
IC.setupGas()
IC.settleGas()
# Set up the dust
IC.setupDust()

print 'Successfully made dust settling test ICs!'
