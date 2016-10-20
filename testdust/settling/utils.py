# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:11:02 2016

@author: ibackus
"""
import os

import testdust
from testdust.utils import loadDefaultParam as _loadDefaultParam

# Get the path of this file
_directory = os.path.dirname(os.path.realpath(__file__))

# Filenames for the .param files (defaults AND user defaults)
parnameDefaults = {'dust': 'dust_default.param',
                   'gas': 'gas_default.param',
                   'general': 'general_default.param'}
parnameUserDefaults = {'dust': 'dust_user_default.param',
                       'gas': 'gas_user_default.param',
                       'general': 'general_user_default.param'}
# Get full path to all the .param files
for parnames in (parnameDefaults, parnameUserDefaults):
    
    for k, v in parnames.iteritems():
        
        parnames[k] = os.path.join(_directory, v)

def loadDefaultParam(param):
    """
    Loads the settling .param defaults for param (returned as a dict).  For
    a list of available params do testdust.settling.utils.paramDefaults.keys()
    
    """
    if param not in parnameDefaults.keys():
        
        raise RuntimeError, "{0} not a valid param key.  valid options are: {1}"\
        .format(param, parnameDefaults.keys())
        
    if param is 'gas':
        # First load the general params
        parDict = loadDefaultParam('general')
        gasPar = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        parDict.update(gasPar)
        
    elif param is 'dust':
        # Load the gas defaults
        parDict = loadDefaultParam('gas')
        dustPar = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        parDict.update(dustPar)
        
    else:
        
        parDict = _loadDefaultParam(parnameDefaults[param], parnameUserDefaults[param])
        
    return parDict
    
def defaultFilenames(step):
    """
    Get the default filenames for 'gas', 'dust', or 'glass'
    
    step can be 'gas', 'dust', or 'glass'
    
    Returns a tuple of (infile, paramname, fprefix)
    """
    if step is 'glass':
        infile = 'glass.std'
        fprefix = 'glass'
    else:
        parDict = loadDefaultParam(step)
        infile, fprefix = testdust.utils.filenamesFromParam(parDict)
    paramname = fprefix + '.param'
    return infile, paramname, fprefix