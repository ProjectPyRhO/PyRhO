
# Main module file for PyRhO

#if sys.version_info < (3,0) #IPython.__version__ < "3":
from __future__ import division # a/b -> float
from __future__ import absolute_import, print_function, unicode_literals 


#__all__ = ['config', 'utilities', 'parameters', 'loadData', 'models', 'protocols', 'simulators', 'fitting', 'IPythonGUI']

__doc__ = """A Python module for fitting, characterising and simulating rhodopsin photocurrents"""
from pkg_resources import get_distribution, DistributionNotFound
__project__ = 'pyrho'
# http://stackoverflow.com/questions/17583443/what-is-the-correct-way-to-share-package-version-with-setup-py-and-the-package
#__version__ = get_distribution(__project__).version #'0.8.0'
try:
    import os
    _dist = get_distribution(__project__)
    distLoc = os.path.normcase(_dist.location) # Normalize case for Windows systems
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(distLoc, __project__)): # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = __project__ + '-' + '(local)'
else:
    __version__ = _dist.version


from pyrho.config import *

import matplotlib as mp    
import matplotlib.pyplot as plt
import numpy as np

# Place all submodule functions and variables into namespace
from pyrho.parameters import *
from pyrho.utilities import *
from pyrho.loadData import * #import loadData
from pyrho.models import * #import models
from pyrho.simulators import *
from pyrho.protocols import * #import modProtocols
from pyrho.fitting import * #import fitStates
from pyrho.IPythonGUI import *


if __name__ == '__main__': 
    try:
        __IPYTHON__
    except NameError:
        pass
    else: # and IPython. See also get_ipython()
        print('Loading IPython GUI!')
        loadGUI()

        
def runAll():
    """Run all protocols with all simulators on all models!"""
    
    for model in modelList: #models: #nStates in [3,4,6]:
        ### Select generative model
        RhO = models[model]()
        for sim in ['Python']:#simulators:
            Sim = simulators[sim](RhO)            
            for prot in protocols: #protocol, init in protParams.items():
                print("\nRunning Protocol '{}' on the {}-state model...".format(prot, model))
                print('--------------------------------------------------------------------------------\n')
                ### Select simulation protocol
                Prot = protocols[prot]()
                Prot.run(Sim, RhO)                
                Prot.plot(Sim, RhO)
                print("\nFinished!")
                print('================================================================================\n\n')

                
def printVersions():
    ### Display version information

    import platform
    print("Python version: ", platform.python_version())
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
    import numpy
    print("NumPy version: ", numpy.__version__)
    import scipy
    print("SciPy version: ", scipy.__version__)
    import matplotlib
    print("Matplotlib version: ", matplotlib.__version__)
    import lmfit
    print("Lmfit version: ", lmfit.__version__)
    
    print("PyRhO version: ", __version__)

