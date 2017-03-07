#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
An example script to make and run the 1D shock of Price & Laibe then produce
a plot comparing the results to the analytic solution of a dustyshock
with a modified sound speed.

Basically, this is the same as a shock tube but with cs^2 scaled by 1-dustFrac

Created on Tue Mar  7 13:59:45 2017

@author: ibackus
"""
import testdust

# Generate ICs
f, pars, args, paramname = testdust.shock.makeICs.shock1d()

# Run simulation. If diskpy is not installed or ChaNGa presets are not 
# set-up, you can run it separately
import diskpy
cmd = diskpy.pychanga.changa_command(paramname, preset='default')
# You could also make your own command, e.g.:
# cmd = 'charmrun +p1 ~/bin/ChaNGa -p 1 ' + paramname
print 'Running Simulation'
diskpy.pychanga.changa_run(cmd)

# Create plot (and save it)
import matplotlib.pyplot as plt
import pynbody
snapname = '{0}.{1:06}'.format(pars['achOutName'], pars['nSteps'])
results = pynbody.load(snapname, paramfile=paramname)
testdust.shock.analyze.plotPriceLaibe(results, paramname)
plt.suptitle('1D dusty-shock results.  ChaNGa vs analytic')
plt.tight_layout()
plt.savefig('results.png')

# Echo arguments used
print 'Arguments used:'
for k, v in args.iteritems():
    print '\t', k, ':', v

# And show the figure
plt.show()
