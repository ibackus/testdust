# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:17:42 2016

@author: ibackus
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz

def _hexline(nx, firstSpacing=1):
    """
    """
    x = np.zeros(nx)
    nx0 = int((nx + 1)/2)
    nx1 = nx - nx0
    x0 = 3 * np.arange(nx0, dtype=float)
    x[0::2] = x0
    x[1::2] = x0[0:nx1] + firstSpacing
    return x
    
def periodicHexMesh(nx, ny):
    """
    Return a hexagonal mesh lattice that is periodic along x and centered
    on the origin
    """
    if (nx % 2) != 0:
        raise ValueError, 'nx must be even'
        
    x, y = hexmesh(nx, ny)
    L = 1 + x.max()
    # Add half a grid spacing
    x += 0.5
    # Shift and scale
    x -= 0.5 * L
    x /= L
    y /= L
    yL = (y.max() - y.min())
    y += 0.5*yL - y.max()
    return x, y
    

def hexmesh(nx, ny):
    """
    Generate a hexagonal mesh lattice.  The lattice spacing is 1
    
    Returns
    -------
    xmat, ymat : arrays
        x, y positions of the lattice sites, shape (nx x ny)
    """
    xeven = _hexline(nx, firstSpacing=2)
    xodd = _hexline(nx, firstSpacing=1)
    xshift = np.zeros(ny)
    xshift[1::2] = 0.5
    
    y = (np.sqrt(3)/2.) * np.arange(ny, dtype=float)
    
    xmat = np.ones([nx, ny])
    ymat = np.dot(np.ones([nx, 1]), y[None, :])
    
    for i in range(0, ny, 2):
        # Even loop
        xmat[:,i] = xshift[i] + xeven
        
    for i in range(1, ny, 2):
        # odd loop
        xmat[:,i] = xshift[i] + xodd
        
    return xmat, ymat
    

def rho(w, h, eps=0., rho0=1.):
    """
    Vertical density profile as a function of dimensionless z
    
    For an isothermal PPD with no self gravity.
    
    Parameters
    ----------
    w : array-like
        Dimensionless z-coordinate (z/R0)
    h : array-like
        Dimensionless scale height (H/R0)
    eps : float
        Dust fraction (between 0 and 1)
    rho0 : float
        Density at z=0
    
    Returns
    -------
    density : array-like
        Analytic solution, \
        :math:`\\rho = \\rho_0 \\exp(- (\\frac{\\sqrt{1+w^2} - 1} \
        {h^2 (1 - \\epsilon) \\sqrt{1+w^2}}  ))`
    """
    a = np.sqrt(1. + w**2)
    den = np.exp( -(a - 1)/((1-eps) * a * h**2) )
    return den*rho0
    
def CDF(w, h, eps=0.):
    """
    Normalized integral of vertical density profile as function of dimensionless
    height.  (integral of rho, see rho()).  Note, this density profile is 
    not normalizable and therefore requires a maximum-height cutoff
    
    For an isothermal PPD with no self gravity.
    
    Parameters
    ----------
    w : array-like
        Dimensionless z-coordinate (z/R0)
    h : array-like
        Dimensionless scale height (H/R0)
    eps : float
        Dust fraction (between 0 and 1)
        
    Returns
    -------
    cdf : array
    """
    rho0 = 1
    den = rho(w, h, eps, rho0)
    # Force density to zero at bound
    den[-1] = 0
    denInt = cumtrapz(den, w, initial=0.)
    # normalize
    denInt /= denInt.max()
    return denInt
    
def CDFInv(h, eps=0., nPoints=1e5, nScaleHeight=5.):
    """
    Generates an interpolating function that returns the CDF inverse for
    a vertical density profile as a function of dimensionless z.
    For an isothermal PPD with no self gravity.  Note, this density profile is 
    not normalizable and therefore requires a maximum-height cutoff
    
    Parameters
    ----------
    h : array-like
        Dimensionless scale height (H/R0)
    eps : float
        Dust fraction (between 0 and 1)
    nPoints : int
        Number of spline points to generate
    nScaleHeight : int
        Number of scale heights to look out to (wmax = nScaleHeight=h)
    
    Returns
    -------
    cdfInvSpl : function
        Function  cdfInvSpl(m) = w. w is z/R0 and 0 < m < 1
    """
    nPoints = int(nPoints)
    wmax = nScaleHeight*h
    w = np.linspace(0, wmax, nPoints)
    # Calculate CDF
    cdf = CDF(w, h, eps)
    # Create the spline interpolation
    cdfInvSpl = interp1d(cdf, w)
    return cdfInvSpl

def rhoEst(z, m=1):
    """
    Estimate (linear) density from sampled points, assuming constant mass, 
    using rho = dm/dz
    """
    dz = np.gradient(z, edge_order=2)
    dm = float(m)/len(z)
    rho = dm/dz
    
    return rho
    
def posGen(n, H, R0, eps, nScaleHeight=5, nSplinePts=1e5, m=None):
    """
    Generate z positions for an isothermal PPD with no self gravity.
    
    Parameters
    ----------
    n : int
        Number of positions to make
    H : float-like
        Disk scale height.  :math:`H = c_s / \\Omega`
    R0 : float-like
        Distance from star (used as length normalization for H and z)
    eps : float
        Dust fraction.  :math:`0 \\leq \\epsilon < 1`
    nScaleHeight : int
        Defines largest z value to calculate.
    nSplinePts : int
        Number of spline points used in CDF inverse calculation
    
    Returns
    -------
    z : array-like
        z - points that will create a correct vertical density profile
    """
    # Dimensionless scale height
    h = H/R0
    # Generate positions
    cdfInvSpl = CDFInv(h, eps, nPoints=1e5)
    if m is None:
        m = np.linspace(-1, 1, n)
    negatives = (m < 0)
    m = abs(m)
    # Generate dimensionless z positions
    w = cdfInvSpl(m)
    w[negatives] *= -1
    # Position
    z = R0 * w
    return z
    
def particleMass(z, rho0):
    """
    Calculate particle mass required to make the one-dimensional density profile
    defined by z have a max value equal to rho0.
    """
    # Density estimate if m=1
    den = rhoEst(z, m=1.)
    m = rho0/den.max()
    return m
    