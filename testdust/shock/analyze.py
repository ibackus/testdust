# -*- coding: utf-8 -*-
"""
Contains functions for analyzing the results of dusty-shock tests.

See plotPriceLaibe() to generate plots for the test results of the 1D shocktube
test

Created on Mon Feb 22 12:23:31 2016

@author: ibackus
"""
import matplotlib.pyplot as plt
import pynbody
import diskpy
SimArray = pynbody.array.SimArray

from testdust.utils import calcGasVals

def analyticSolution(*args, **kwargs):
    """
    A wrapper for sod.solve().  Requires the sod-shocktube package, available
    at https://github.com/ibackus/sod-shocktube
    See sod.solve() for documentation
    """
    try:
        import sod
    except ImportError:
        print "sod-shocktube package must be installed.  It can be found at"\
        " https://github.com/ibackus/sod-shocktube"
        raise
        
    positions, regions, sol = sod.solve(*args, **kwargs)
    return positions, regions, sol

def plotPriceLaibe(f, paramname, t=None, left_state=(1, 1, 0), 
                   right_state=(0.125, 0.125, 0.), geometry=(0., 1., 0.5)):
    """
    Generates plots for the 1D dustyshock tube, similar to those of 
    Price & Laibe 2015, figure 3.
    
    Parameters
    ----------
    f: SimSnap
        Snapshot to analyze
    paramname: str
        Path to .param file used for the ChaNGa run
    t: float
        (optional) Time of the snapshot (code units).  Used for analytic
        solution.  Defaults to time of the last snapshot.
    left_state, right_state: tuple
        A tuple of the state (pressure, density, velocity) on each side of the
        shocktube barrier for the ICs.  In the case of a dusty-gas, the density
        should be the gas density.
    geometry: tuple
        A tuple of positions for (left boundary, right boundary, barrier)
    """
    # Load stuff
    if isinstance(f, str):
        
        f = pynbody.load(f, paramfile=paramname)
        
    pars = diskpy.utils.configparser(paramname)
    units = diskpy.pychanga.units_from_param(pars)
    
    if t is None:
        
        print "No time supplied.  Assuming final snapshot"
        t = pars['nSteps'] * pars['dDelta']
        print "t =", t
    
    if not pynbody.units.has_units(t):
        t = SimArray(t, units['t_unit'])
    else:
        t.convert_units(units['t_unit'])
    
    gamma, molecularWeight = calcGasVals(f, pars)
    # Analytic solution
    dustFrac = float(f['dustFrac'].mean())
    print 'Mean dust Fraction:', dustFrac
    positions, regions, sol = analyticSolution(left_state, right_state, \
        geometry, t, gamma=gamma, npts=len(f), dustFrac=dustFrac)
    sol['x'] -= geometry[1] * 0.5
    
    # Plot
    plt.clf()
    axs = [plt.subplot(2, 2, 1)]
    axs.extend([plt.subplot(2, 2, i, sharex=axs[0]) for i in range(2,5)])

    ylabels = (r'$v_x$', r'$\rho$', 'u', 'P')
    fkeys = ('vz', 'rho', 'energy', 'pressure')
    solkeys = ('u', 'rho_total', 'energy', 'p')

    for i in range(4):
        ax = axs[i]
        ax.plot(f['z'], f[fkeys[i]], '.k', sol['x'], sol[solkeys[i]],'r')
        ax.set_ylabel(ylabels[i])

    for i in (1, 3):
        axs[i].yaxis.set_major_locator(plt.NullLocator())
        
    for i in (0, 1):
        axs[i].tick_params(axis='x', which='both', bottom='off', 
           labelbottom='off')

    plt.xlim(-0.5, 0.5)
    plt.tight_layout()
    return axs