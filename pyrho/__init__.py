"""A Python module for fitting, characterising and simulating rhodopsin photocurrents"""
#__doc__ = 
# Main module file for PyRhO

#if sys.version_info < (3,0) #IPython.__version__ < "3":
from __future__ import division # a/b -> float
from __future__ import absolute_import, print_function, unicode_literals 

import platform
IPY = False
try:
    import IPython
    IPY = True
except ImportError:
    IPY = False

from pkg_resources import get_distribution, DistributionNotFound

# Necessary?
import matplotlib as mpl   
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import lmfit


# Place all submodule functions and variables into namespace
from pyrho.config import *
from pyrho.parameters import *
from pyrho.utilities import *
from pyrho.expdata import * #import loadData
from pyrho.models import * #import models
from pyrho.simulators import *
from pyrho.protocols import * #import modProtocols
from pyrho.fitting import * #import fitStates
from pyrho.jupytergui import * # IPythonGUI

#TODO
#__all__ = ['config', 'utilities', 'parameters', 'expdata', 'models', 'protocols', 'simulators', 'fitting', 'IPythonGUI']

__project__ = 'pyrho'
# http://stackoverflow.com/questions/17583443/what-is-the-correct-way-to-share-package-version-with-setup-py-and-the-package
#__version__ = get_distribution(__project__).version #'0.8.0'
try:
    import os
    _dist = get_distribution(__project__)
    distLoc = os.path.normcase(_dist.location) # Normalise case for Windows systems
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(distLoc, __project__)): # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = __project__ + '-' + '(local)'
else:
    __version__ = _dist.version

if __name__ == '__main__': 
    try:
        __IPYTHON__
    except NameError:
        pass
    else: # and IPython. See also get_ipython()
        print('Loading IPython GUI!')
        loadGUI()

# TODO Move everything except imports elsewhere
        
def runAll(listOfModels=[6]):
    """Run all protocols with the Python simulator on a list of models with default parameters!
    listOfModels    := individual or list of integers or strings 
                        e.g. [3,4,6], 3, '4', ['4','6'], modelList
    """
    
    if not isinstance(listOfModels, (list, tuple)):
        listOfModels = [listOfModels] # ints or strs
    listOfModels = [str(m) for m in listOfModels]
    
    for model in listOfModels:
        ### Select generative model
        RhO = models[model]()
        for prot in protocols:
            ### Select simulation protocol
            Prot = protocols[prot]()
            for sim in ['Python']:#simulators:
                Sim = simulators[sim](Prot, RhO)            
                print("\nRunning Protocol '{}' on the {}-state model...".format(prot, model))
                print('--------------------------------------------------------------------------------\n')
                Sim.run()                
                Sim.plot()
                print("\nFinished!")
                print('================================================================================\n\n')

                
def printVersions():
    """Display version information for PyRhO and its dependencies"""

    #import platform
    print("Python version: ", platform.python_version())
    if IPY:    
        try:
            #import IPython
            #__IPYTHON__
            print("IPython version: ", IPython.__version__)
        except ImportError: # IPython not found
            pass
    #deps = [numpy, scipy, matplotlib, lmfit, warnings, os, pickle, collections, platform]
    #depsGUI = [IPython, ast, base64]
    #for mod in dependencies:
    #    print("{} version: {}".format(mod, mod.__version__))
    #import numpy
    print("NumPy version: ", np.__version__)
    #import scipy
    print("SciPy version: ", sp.__version__)
    #import matplotlib
    print("Matplotlib version: ", mpl.__version__)
    #import lmfit
    print("Lmfit version: ", lmfit.__version__)
    
    print("PyRhO version: ", __version__)

