# config.py

import matplotlib as mp
import numpy as np
import warnings
import os


##### Output level settings #####
global verbose
if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
    verbose = 1 # Text notification output level [0,1,2]


### Set data and figure directories to defaults
if 'dDir' not in vars() or 'dDir' not in globals() or dDir is None:
    dDir = 'data/'
if not os.path.isdir(dDir):
    createDir(dDir)
    
if 'fDir' not in vars() or 'fDir' not in globals() or fDir is None:
    fDir = 'figs/'
if not os.path.isdir(fDir):
    createDir(fDir)

def createDir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise    
    
    
##### Plot settings #####

display = 'paper' #'screen' #'paper'

colours = ['b','g','r','c','m','y','k']
styles = ['-', '--', '-.', ':']

latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mp.rcParams['text.usetex'] = True
#plt.rcParams.update(rcdef) # Set back to defaults

golden_ratio  = (1.0 + np.sqrt(5)) / 2.0

# http://matplotlib.org/users/customizing.html
if display == 'screen':
    dpi = 300
    #figFormat = 'retina'
    saveFigFormat = 'png'    
    figWidth = 12
    figHeight = figWidth / golden_ratio
    tickSize = 14
    labelSize = 16
    legendSize = 12
    titleSize = 'x-large'
    eqSize = 18
    linewidth = 1.5
    addTitles = True
    addStimulus = True

elif display == 'paper':
    #http://www.nature.com/nature/authors/gta/3c_Final_artwork.pdf
    #http://www.plosone.org/static/figureSpecifications
    dpi = 300 # Set DPI 300 - 600
    saveFigFormat = 'png'#'eps'           # Supported formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz.
    #figFormat = 'eps'
    figWidth = 5#7
    figHeight = figWidth / golden_ratio
    tickSize = 6
    labelSize = 7
    legendSize = 7
    titleSize = 'medium'
    eqSize = 7
    linewidth = 1
    addTitles = False
    addStimulus = True

else:
    warnings.warn('Warning: Unknown display type selected - using matplotlib defaults!')

mp.rcParams['figure.figsize'] = (figWidth, figHeight) #pylab.rcParams
mp.rcParams['figure.dpi'] = dpi
mp.rcParams['savefig.dpi'] = dpi
mp.rcParams['axes.labelsize'] = labelSize
mp.rcParams['ytick.labelsize'] = tickSize
mp.rcParams['xtick.labelsize'] = tickSize
mp.rcParams['legend.fontsize'] = legendSize
mp.rcParams['axes.titlesize'] = titleSize
mp.rcParams['lines.linewidth'] = linewidth


plotPeakRecovery = False
plotStateVars = False
plotKinetics = False



