#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This is a command line utility to  generate ICs for the dust settling 
test of Price & Laibe (2015) (see section 4.4 Dust settling in a 
protoplanetary disc)

This will be SLOW!  The long part of the procedure is settling the gas 
(IC.settleGas()).  This takes a long time because the simulation is run for
several dynamical times with ChaNGa in order to reach equilibrium before
adding the dust.

If you stop or things crash, you can load the ICs using
IC = testdust.settling.makeICs.IC(settingsfile)
and restart from where you left off.

To Generate these ICs, ChaNGa MUST be compiled with the glass damping option 
(i.e. the damping force, proportional to the velocity)

Code units are specified in PL15

Created on Fri Oct 21 17:07:56 2016

@author: ibackus
"""
import sys

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        
        raise RuntimeError, "\n\nNEED 2 ARGUMENTS: <settings file> <output dir>"
    
    import testdust
    import os
    import shutil
    # It's a good idea to generate the ICs in a clean folder
    savedir = sys.argv[2]
    # IC settings are stored in this file
    settingsfile = sys.argv[1]
    # Initialize the ICs with the settings
    IC = testdust.settling.makeICs.IC(settingsfile)
    # Its a good idea to generate these ICs in their own folder
    if not os.path.exists(savedir):
        
        os.mkdir(savedir)
    
    cwd = os.getcwd()
    
    try:
        
        shutil.copyfile(settingsfile, os.path.join(savedir, settingsfile))
        os.chdir(savedir)
        # Generate initial snapshot
        IC.makeInitialSnap()
        print 'Initial snapshot successfully made'
        # Settle the ICs
        IC.setupGas()
        print 'Gas settling snapshot successfully made'
        IC.settleGas()
        print 'Gas settling successfully performed'
        # Set up the dust
        IC.setupDust()
        print 'Successfully made dust settling test ICs!'
        print 'Saved to directory:', savedir
    
    finally:
        
        os.chdir(cwd)
