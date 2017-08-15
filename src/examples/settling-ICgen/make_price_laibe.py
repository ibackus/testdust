# -*- coding: utf-8 -*-
"""
This example will generate ICs for the dust settling test of Price & Laibe 
(2015) (see section 4.4 Dust settling in a protoplanetary disc)

This will be SLOW!  The long part of the procedure is settling the gas 
(IC.settleGas()).  This takes a long time because the simulation is run for
several dynamical times with ChaNGa in order to reach equilibrium before
adding the dust.

If you stop or things crash, you can load the ICs using
IC = testdust.settling.makeICs.IC(settingsfile)
and restart from where you left off.

To Generate these ICs, ChaNGa MUST be compiled with the glass damping option 
(i.e. the damping force, proportional to the velocity)

There are several example settingsfiles, check them out.  They define the
settling ICs to generate.

Code units are specified in PL15

Created on Fri Oct 21 17:07:56 2016

@author: ibackus
"""
import testdust
import os

# It's a good idea to generate the ICs in a clean folder
savedir = 'dust-settling'

# IC settings are stored in this file
settingsfile = 'settings_hexagonal.py' # PL15 test (stretched hexagonal lattice)
#settingsfile = 'settings_randomnormal_3d.py' # 3D Glass
#settingsfile = 'settings_cubic_grid.py' # Stretched cubic lattice

# Initialize the ICs with the settings
IC = testdust.settling.makeICs.IC(settingsfile)

# Its a good idea to generate these ICs in their own folder
if not os.path.exists(savedir):
    
    os.mkdir(savedir)

cwd = os.getcwd()

try:
    
    os.chdir(savedir)
    # Generate initial snapshot
    IC.makeInitialSnap()
    # Settle the ICs
    IC.setupGas()
    IC.settleGas()
    # Set up the dust
    IC.setupDust()
    
    print 'Successfully made dust settling test ICs!'
    print 'Saved to directory:', savedir

finally:
    
    os.chdir(cwd)
