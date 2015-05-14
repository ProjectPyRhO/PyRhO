
# Main module file for pyRhO

#if IPython.__version__ < "3":
from __future__ import division # a/b -> float
from __future__ import absolute_import, print_function, unicode_literals 


#__all__ = ['parameters', 'loadData', 'models', 'protocols', 'fitting', 'IPythonGUI'] # 'config'

__doc__ = """A Python module for fitting, characterising and simulating rhodopsin photocurrents"""
__version__ = '0.5.0'

from .config import *

#deps = [numpy, scipy, matplotlib, lmfit, warnings, os, pickle, collections, platform]
#depsGUI = [IPython, ast, base64]

# global verbose
# if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
    # verbose = 1 # Text notification output level [0,1,2]

# global addTitles
# if 'addTitles' not in vars() or 'addTitles' not in globals() or addTitles is None:
    # print('Plotting Titles')
    # addTitles = True
    


import matplotlib as mp    
#from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
#import IPython

#from scipy.signal import argrelextrema
#from scipy.optimize import curve_fit
#from lmfit import minimize, Parameters, fit_report
#import sys, os, errno

#import pickle

#import numpy as np
#import matplotlib.pyplot as plt
#from scipy.signal import * 
#import warnings
#from lmfit import Parameters

# Place all submodule functions and variables into namespace
from .parameters import *
from .loadData import * #import loadData
from .models import * #import models
from .simulators import *
from .protocols import * #import modProtocols
from .fitting import * #import fitStates
from .IPythonGUI import *


# Plot configuration
#colours = ['b','g','r','c','m','y','k']
#styles = ['-', '--', '-.', ':']
#golden_ratio  = (1.0 + np.sqrt(5)) / 2.0
#figWidth = 16
#figHeight = figWidth / golden_ratio
#mp.rcParams['figure.figsize'] = (figWidth, figHeight) #pylab.rcParams
#mp.rcParams['axes.labelsize'] = 16
#mp.rcParams['ytick.labelsize'] = 14
#mp.rcParams['xtick.labelsize'] = 14
#mp.rcParams['legend.fontsize'] = 12
#mp.rcParams['axes.titlesize'] = 'x-large'
#eqSize = 18
#latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
#if latexInstalled:
#    mp.rcParams['text.usetex'] = True
#plt.rcParams.update(rcdef) # Set back to defaults

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
    
    for model in models: #nStates in [3,4,6]:
        ### Select generative model
        RhO = models[model]()
        #nStates = int(model)
        #RhO = selectModel(nStates)
        
        for sim in simulators:
            
            Sim = simulators[sim](RhO)
            
            for prot in protocols: #protocol, init in protParams.items():
                ### Select simulation protocol
                #protocols = ['custom', 'saturate', 'rectifier', 'shortPulse', 'recovery']
                #protocol = protocols[2] #'recovery'#'shortPulse' # Set this interactively with radio buttons?
                
                print("\nRunning Protocol '{}' on the {}-state model...".format(protocol,nStates))
                print('--------------------------------------------------------------------------------\n')
                #Prot = selectProtocol(protocol)
                Prot = protocols[prot]()
                Prot.run(Sim,RhO)
                
                # Plot settings
                #plotPeakRecovery = False
                #Prot.plotStateVars = False
                #plotKinetics = False
                Prot.plot(Sim,RhO)
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
    #for mod in dependencies:
    #    print("{} version: {}".format(mod, mod.__version__))
    #import numpy
    print("NumPy version: ", numpy.__version__)
    import scipy
    print("SciPy version: ", scipy.__version__)
    #import matplotlib
    print("Matplotlib version: ", matplotlib.__version__)
    import lmfit
    print("lmfit version: ", lmfit.__version__)
    
    print("PyRhO version: ", pyrho.__version__)

