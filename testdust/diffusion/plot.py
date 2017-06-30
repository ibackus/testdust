# -*- coding: utf-8 -*-
"""
Contains functions to plot the results of the dustydiffusion test.

@author: ibackus
"""
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import diskpy

#sim, epsEstimator, ts, runpars = analyze.loadSim(simdir)

def crossSection(sim, ts, crossSectionTimes=[0, 1, 10]):
    """
    Reproduces the cross-section plot of dust density of 
    Price & Laibe 2015, fig. 5
    
    Note, sim and ts can be loaded with analyze.loadSim(...)
    
    Parameters
    ----------
    sim : list
        List of SimSnaps for the simulation
    ts : array-like
        Snapshot times
    crossSectionTimes : array-like
        (optional) Sim times to plot (approximate)
    """
    # Select times to plot at
    crossSectionTimes = np.asarray(crossSectionTimes)
    crossSectionTimes = crossSectionTimes.reshape(crossSectionTimes.size)
    if np.ndim(crossSectionTimes) == 0:
        crossSectionTimes = crossSectionTimes[None]
        
    nPlots = len(crossSectionTimes)
    
    # Plot
    axs = diskpy.plot.gridplot(1, nPlots, square=True)
    fig = plt.gcf()
    
    for iPlot in range(nPlots):
        
        ax = axs[iPlot]
        iTime = abs(ts - crossSectionTimes[iPlot]).argmin()
        t = ts[iTime]
        f = sim[iTime]
        
        im=pynbody.plot.sph.image(f, 'dustFrac', width=1, log=False, vmin=0,
                               vmax = 0.11, cmap='cubehelix_r', 
                               show_cbar=False, subplot=ax, ret_im=True)
        
        ax.set_xlabel('t={:.2g}'.format(float(t)))
                               
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.05, 0.6])
    fig.colorbar(im, cax=cbar_ax)
    fig.set_size_inches(8.5, 3.4, forward=True)
    plt.suptitle('Cross section of dust fraction in z=0 plane\n'\
    'See Price & Laibe (2015)')

def dustFracProfile(sim, ts, epsEstimator, 
                    epsPlotTimes=[0., 0.1, 0.3, 1, 3, 10], nr=200,
                    colorcode=True, legend=True):
    """
    Note, sim and ts and epsEstimator can be loaded with analyze.loadSim(...)
    
    Parameters
    ----------
    sim : list
        List of SimSnaps for the simulation
    ts : array-like
        Snapshot times
    epsEstimator : function
        A function of (r, t) that returns the analytic estimate of the dust
        fraction density profile of P&L15 dustydiffusion
    epsPlotTimes : array-like
        Approximate times to plot at
    nr : int
        Number of radial bins
    """
    # Make plot times an array
    epsPlotTimes = np.asarray(epsPlotTimes)
    epsPlotTimes = epsPlotTimes.reshape(epsPlotTimes.size)
    
    nt = len(epsPlotTimes)
    actualPlotTimes = np.zeros(nt)
    title = 'plot times: '
    
    if colorcode:
        markercolor = None
    else:
        markercolor = 'k'
    
    for iPlot in range(nt):
        
        iTime = abs(ts - epsPlotTimes[iPlot]).argmin()
        
        # Calculate stuff
        f = sim[iTime]
        t = ts[[iTime]]
        actualPlotTimes[iPlot] = t
        print t
        r = np.linspace(0, f['r'].max(), nr)
        epsAnalytic = epsEstimator(r, t)
        
        # Plot
        scatter=plt.plot(f['r'], f['dustFrac'], 'o', markersize=3, 
                 markeredgecolor='none', label='t={:.2g}'.format(float(t)),
                 color=markercolor)
        line=plt.plot(r, epsAnalytic, 'r')
        if colorcode:
            # Make lines and points the same color
            line[0].set_color(scatter[0].get_color())
        title += '{:.2g}, '.format(float(t))
        
    # Set-up plot
    plt.ylim(0, 0.11)
    plt.xlim(0, 0.5)
    plt.ylabel('Dust fraction')
    plt.xlabel('r')
    if legend:
        plt.legend(loc='best', markerscale=2)
    plt.title(title)
    
