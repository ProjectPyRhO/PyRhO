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
    if os.path.isdir(path):
        return
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise    
    
#pyrhoPath = os.path.abspath(pyrho.__file__) # http://stackoverflow.com/questions/247770/retrieving-python-module-path
NMODLfiles = ['RhO3.mod', 'RhO4.mod', 'RhO6.mod']
#pyrhoNEURONpath = os.path.join(pyrhoPath, 'NEURON')
#NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
    
##### Plot settings #####

display = 'screen' #'paper'

colours = ['b','g','r','c','m','y','k']
styles = ['-', '--', '-.', ':']

latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mp.rcParams['text.usetex'] = True
#mp.rcParams.update(rcdef) # Set back to defaults

def setFigOutput(display='screen'):
    golden_ratio  = (1.0 + np.sqrt(5)) / 2.0

    # http://matplotlib.org/users/customizing.html
    if display == 'screen':
        dpi = 300 #mp.rcParams['savefig.dpi'] #300
        #figFormat = 'retina'
        saveFigFormat = 'png'    
        figWidth = 12
        figHeight = figWidth / golden_ratio
        tickSize = 14
        labelSize = 16
        legendSize = 12
        titleSize = 'small' #'x-large'
        eqSize = 18
        linewidth = 1.5
        addTitles = True
        addStimulus = True

    elif display == 'paper':
        #http://www.nature.com/nature/authors/gta/3c_Final_artwork.pdf
        #http://www.plosone.org/static/figureSpecifications
        #http://www.ploscompbiol.org/static/figureSpecifications
        #http://www.frontiersin.org/about/AuthorGuidelines#_FigureandTableGuidelines
        dpi = 1200 # Set DPI 300 - 600
        saveFigFormat = 'png'#'eps'             # Supported formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz.
        #figFormat = 'eps'
        figWidth = 18/2.54 #5#7
        figHeight = figWidth / golden_ratio
        tickSize = 8
        labelSize = 10
        legendSize = 10
        titleSize = 'medium'
        eqSize = 12
        linewidth = 2
        addTitles = False
        addStimulus = True

        mp.rcParams['savefig.dpi'] = dpi        # http://nbviewer.ipython.org/gist/minrk/3301035
                                                # http://wiki.scipy.org/Cookbook/Matplotlib/AdjustingImageSize
        
    else:
        warnings.warn('Warning: Unknown display type selected - using matplotlib defaults!')
    
    # http://matplotlib.org/users/customizing.html
    mp.rcParams['figure.figsize'] = (figWidth, figHeight) #pylab.rcParams
    mp.rcParams['figure.dpi'] = dpi
    #mp.rcParams['savefig.dpi'] = dpi           # http://nbviewer.ipython.org/gist/minrk/3301035
    mp.rcParams['axes.labelsize'] = labelSize
    mp.rcParams['ytick.labelsize'] = tickSize
    mp.rcParams['xtick.labelsize'] = tickSize
    mp.rcParams['legend.fontsize'] = legendSize
    mp.rcParams['axes.titlesize'] = titleSize
    mp.rcParams['lines.linewidth'] = linewidth
    mp.rcParams['font.size'] = eqSize
    
    
    return saveFigFormat, addTitles, addStimulus, eqSize

# plotPeakRecovery = False
# plotStateVars = False
# plotKinetics = False

saveFigFormat, addTitles, addStimulus, eqSize = setFigOutput(display)

try:
    __IPYTHON__
    import IPython
except NameError:
    pass
else: # and IPython. See also get_ipython()
    figFormat = 'retina' # 'svg'
    IPython.display.set_matplotlib_formats(figFormat)
    print("Default display figure format set: "+figFormat)


#IPython.core.magic.Magics.basic.pprint()
#%pprint
#%precision %.6g
#np.set_printoptions(precision=6) # Number of decimal places for floats
# Long    Scientific
# 1e-4    1e-5
# 1e15    1e16

