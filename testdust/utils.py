# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:39:07 2016

@author: ibackus
"""
import diskpy
import glob
import os
import shutil
from types import ModuleType

def parsePySettings(filename):
    """
    Parses settings in a .py file by running execfile()
    Settings files for this format are just python scripts containing variable
    definitions
    
    Modules and the doc string will be ignored
    
    Returns a dict of the settings
    """
    settings = {}
    globalsDict = {}
    execfile(filename, globalsDict, settings)
    # Ignore the doc string for the file if there is one
    if '__doc__' in settings:
        settings.pop('__doc__', None)
        
    # Ignore any imported modules
    for key in settings.keys():
        
        if isinstance(settings[key], ModuleType):
            
            settings.pop(key, None)
    
    return settings

def copySnapshot(filename, dest):
    """
    Copies a snapshot and auxilliary arrays to dest.  dest can be a directory
    or a new filename
    """
    fnames, destinations = _snapshotSrcDest(filename, dest)
    for fname, destination in zip(fnames, destinations):
        shutil.copyfile(fname, destination)

def moveSnapshot(filename, dest):
    """
    Moves a snapshot and its auxilliary arrays to dest. dest can be a directory
    or a new filename
    """
    fnames, destinations = _snapshotSrcDest(filename, dest)
    for fname, destination in zip(fnames, destinations):
        shutil.move(fname, destination)
    
def _snapshotSrcDest(filename, dest):
    """
    A utility to get all the filenames of a snapshot and its auxilliary arrays
    and set up destination paths for a move or copy operation
    """
    fnames = glob.glob(filename + '*')
    exts = [fname.split(filename)[-1] for fname in fnames]
    
    if os.path.isdir(dest):
        
        destPrefix = os.path.join(dest, filename)
        
    else:
        
        destPrefix = dest
        
    destinations = [destPrefix + ext for ext in exts]
    
    return fnames, destinations

def loadDefaultParam(defaults, userdefaults):
    """
    Loads the parameters in the ChaNGa .param file defaults and overrides 
    them with those in userdefaults
    
    Parameters
    ----------
    defaults, userdefaults : str
        File names to the .param files to load.  userdefaults override defaults
    
    Returns
    -------
    params : dict
        A dictionary of the params
    """
    defaults = diskpy.utils.configparser(defaults, 'param')
    userdefaults = diskpy.utils.configparser(userdefaults, 'param')
    defaults.update(userdefaults)
    return defaults

def filenamesFromParam(param):
    """
    Retrieves the input filename and output prefix for a .param file (or a 
    loaded .param dict)
    
    Returns (infile, fprefix)
    """
    if isinstance(param, str):
        param = diskpy.utils.configparser(param, 'param')
    
    inFile = param['achInFile']
    prefix = diskpy.pychanga.getpar('achOutName', param)
    return inFile, prefix