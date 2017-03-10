#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
An example script to make and run the 1D wave of Price & Laibe and compare it
to the analytically propagating 1D wave

User defaults for the generated .param file can be placed in:
    testdust/wave/user_wave_defaults.param

Created on Tue Mar  7 13:59:45 2017

@author: ibackus
"""
import testdust
import numpy as np

L = 1.
cs = 1.
drho = 1e-4
v0 = drho
nPeriods = 1
nSmooth = 5

# Generate ICs
f, pars, paramname = testdust.wave.makeICs.wave1d(drho=drho, L=L, cs=cs,
    nSmooth=nSmooth, nPeriods=nPeriods, v0=v0)

# Run simulation. If diskpy is not installed or ChaNGa presets are not 
# set-up, you can run it separately
import diskpy
cmd = diskpy.pychanga.changa_command(paramname, preset='onecore')
# You could also make your own command, e.g.:
# cmd = 'charmrun +p1 ~/bin/ChaNGa -p 1 ' + paramname
print 'Running Simulation'
diskpy.pychanga.changa_run(cmd)

# Create plot (and save it)
import matplotlib.pyplot as plt
import pynbody
snapname = '{0}.{1:06}'.format(pars['achOutName'], pars['nSteps'])
t = pars['nSteps'] * pars['dDelta']
results = pynbody.load(snapname, paramfile=paramname)

period = L/(cs * float(nPeriods))
k = 2 * np.pi * nPeriods / L
z = np.linspace(-0.5*L, 0.5*L, 1000)
ypred = drho * np.sin(k * (z - cs*t))
plt.plot(results['z'], results['rho']/results['rho'].mean() - 1,'o')
plt.plot(z, ypred, 'r')

plt.legend(['ChaNGa', 'analytic'])
plt.ylabel('Fractional density perturbation ')
plt.xlabel('z')
plt.title('Dustywave propagation after t={0} or {1} periods'.format(t, t/period))
plt.tight_layout()
# Save the figure
figname = 'dustywave.png'
plt.savefig(figname)
print 'Saved figure to:', figname
# And show the figure
plt.show()
