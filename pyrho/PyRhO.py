# Main module file for pyRhO

#if IPython.__version__ < "3":
from __future__ import division # a/b -> float
from __future__ import absolute_import, print_function, unicode_literals 

global verbose
if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
    verbose = 1 # Text notification output level [0,1,2]


    
from matplotlib import pyplot as plt
import matplotlib as mp
#import matplotlib.pyplot as plt
import numpy as np
#import IPython


from scipy.signal import argrelextrema
#from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, fit_report
import sys, os, errno

import pickle

#import numpy as np
#import matplotlib.pyplot as plt
from scipy.signal import * 
import warnings
#from lmfit import Parameters

# Place all submodule functions and variables into namespace
from parameters import *
from loadData import * #import loadData
from models import * #import models
from modProtocols import * #import modProtocols
from fitStates import * #import fitStates
