# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:59:50 2016

@author: ibackus
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pynbody
import diskpy
import os
SimArray = pynbody.array.SimArray
import json

# Derived arrays
@pynbody.derived_array
def dustRho(sim):
    sim['dustFrac'].units = pynbody.units.Unit('1')
    return sim['dustFrac'] * sim['rho']

def loadSim(simdir, fprefix='settledust', settingsFile='runparams.json'):
    
    if hasattr(simdir, '__iter__'):
        
        return [loadSim(x, fprefix) for x in simdir]
        
    paramfile = os.path.join(simdir, fprefix + '.param')
    fnames = diskpy.pychanga.get_fnames(fprefix, simdir)
    fs = [pynbody.load(fname, paramfile=paramfile) for fname in fnames]
    
    settingsFile = os.path.join(simdir, settingsFile)
    settings = json.load(open(settingsFile,'r'))
    
    return fs, settings
    
def getTime(f):
    """
    Normalized by crossing time, for these sims equal to nParticles**1/3
    """
    if not isinstance(f, pynbody.snapshot.SimSnap):
        
        t = [getTime(x) for x in f]
        units = t[0].units
        tnew = SimArray(np.zeros(len(t)), units)
        for i in range(len(t)):
            tnew[i] = t[i]
        return tnew
        
    t_units = f.infer_original_units('yr')
    t = SimArray(1.0, f.properties['time'] - t_units)
    return t.in_units(t_units)
    
def arrayify(fs, key):
    
    array = np.array([f[key] for f in fs])
    return array
    
def loadArrays(fs, keys=['z', 'rho', 'dustFrac']):
    
    if not hasattr(fs, '__iter__'):
        
        raise RuntimeError, 'No __iter__ attribute!'
        
    if not isinstance(fs[0], pynbody.snapshot.SimSnap):
        
        return [loadArrays(x, keys) for x in fs]
        
    results = {}
    
    for key in keys:
        
        results[key] = arrayify(fs, key)
        
    return results
    
def flattenResults(results, keys=None):
    
    vals = []
    if keys is None:
        
        keys = results[0].keys()
    
    for key in keys:
        
        vals.append([results[key]])
        
    return vals, keys
    
def analyze(fs, pars):
    
    keys = ['z', 'rho', 'dustFrac']
    results = loadArrays(fs, keys)
    
    #(z, dustFrac, rho), keys = flattenResults(results, keys)
    results['rhoDust'] = results['rho'] * results['dustFrac']
    results['t'] = getTime(fs)
    results['period'] = 2*np.pi*pars['R0']**1.5
    
    return results
    
def scatter3D(f, qty):
    """
    """
    plt.clf()
    fig = plt.gcf()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('y')
    ax.set_ylabel('z')
    ax.set_zlabel(qty)
    
    ax.scatter(f['y'], f['z'], f[qty], 'o')
    
    fig.canvas.draw()
    
def plotDustRho(results, figname='figures/plot.pdf', **kwargs):
    
    r = results
    xmax = abs(r['z']).max()
    xlim = [-xmax, xmax]
    ylim = [0, r['rhoDust'].max()]
    torb = r['t']/r['period']
    
    with PdfPages(figname) as pdf:
        
        for i in range(9):
            
            figtitle = 't = {:.3} orbits'.format(float(torb[i]))
            print figtitle
            plt.clf()
            plt.plot(r['z'][i], r['rhoDust'][i], ',', rasterized=True, **kwargs)
            plt.ylim(ylim)
            plt.xlim(xlim)
            plt.xlabel('z')
            plt.ylabel('dust density')
            plt.title(figtitle)
            pdf.savefig(bbox_inches='tight')
            
        print 'figure saved to', figname
