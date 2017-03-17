# -*- coding: utf-8 -*-
"""
Contains functions for analyzing the output of a dustydiffusion simulation, 
as in Price & Laibe 2015

Created on Wed Aug  3 17:12:21 2016

@author: ibackus
"""
import os
import pynbody
SimArray = pynbody.array.SimArray
import diskpy
import numpy as np
import cPickle as pickle

figdir = 'figures'
fRunParams = 'runparams.json'

# -----------------------------------------
# SETTINGS
# -----------------------------------------
def epsAnalytic(r, t, criticalRadius, tStop, eps0, cs):
    """
    Analytic estimate for the dust fraction profile for the Price and Laibe 
    2015 dustydiffusion test.  Note there seems to be a typo in their paper.
    """
    B = (criticalRadius**2)/eps0
    A = eps0 * abs(B)**(3,5)
    eta = tStop * cs**2
    
    t = diskpy.utils.match_units(t, tStop)[0]
    r = diskpy.utils.match_units(r, criticalRadius)[0]
    
    eps = A * abs(10*eta*t + B)**(-3,5) - (r**2)/(10*eta*t + B)
    
    if (np.size(eps) == 1) and eps < 0:
        eps *= 0
    else:
        eps[eps<0] = 0
    return eps
    
def epsEstimator(runParams):
    """
    A wrapper for epsAnalytic which returns a function to calculate eps
    (analytic) for a set of runParams
    
    returns epsEst(r, t)
    """
    
    def epsEst(r, t):
        """
        """
        return epsAnalytic(r, t, runParams['criticalRadius'], \
        runParams['tStop'], runParams['dustFrac0'], runParams['cs'])
    
    return epsEst

def getTime(f):
    """
    Time for the snapshot(s) f
    """
    if not isinstance(f, pynbody.snapshot.SimSnap):
        
        results = [getTime(x) for x in f]
        units = results[0].units
        output = SimArray(np.zeros(len(results)), units)
        for i, result in enumerate(results):
            
            output[i] = result
            
        return output
        
    t_units = f.infer_original_units('yr')
    t = SimArray(1.0, f.properties['time'] - t_units)
    return t.in_units(t_units)
    #return float(t.in_units(t_units) * len(f)**(-1.0/3))
    
def loadSim(simdir, paramname='snapshot.param'):
    """
    Loads a dustydiffusion simulation from a directory.  In particular, this
    requires the runparameters to be saved in the directory.
    
    Returns
    -------
    fs : list
        list of pynbody snapshots
    eps : fcn
        Estimate the dust fraction profile as a function of (r, t)
    t : SimArray
        Times of the output snapshots
    runparams : dict
        dictionary of the run params
    """
    if not isinstance(simdir, str):
        
        flists = []
        parlists = []
        
        for x in simdir:
            
            f, p = loadSim(x)
            flists.append(f)
            parlists.append(p)
        
        return flists, parlists
        
    
    print simdir
    # Load files
    paramname = os.path.join(simdir, paramname)
    pars = diskpy.utils.configparser(paramname)
    fprefix = diskpy.pychanga.getpar('achOutName', pars)
    fnames = diskpy.pychanga.get_fnames(fprefix, simdir)
    fs = [pynbody.load(fname) for fname in fnames]
    
    runparams = loadRunParams(simdir)
    eps = epsEstimator(runparams)
    t = getTime(fs)
    
    return fs, eps, t, runparams
    
def loadRunParams(simdir):
    
    runParamsPath = os.path.join(simdir, 'runparams.p')
    return pickle.load(open(runParamsPath, 'r'))
        