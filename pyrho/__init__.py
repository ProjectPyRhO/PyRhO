
# Main module file for pyRhO

#if IPython.__version__ < "3":
from __future__ import division # a/b -> float
from __future__ import absolute_import, print_function, unicode_literals 


#__all__ = ['parameters', 'loadData', 'models', 'protocols', 'fitStates', 'IPythonGUI'] # 'config'

__doc__ = """A Python module to fit and characterise rhodopsin photocurrents"""

from .config import *

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
from .protocols import * #import modProtocols
from .fitStates import * #import fitStates
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