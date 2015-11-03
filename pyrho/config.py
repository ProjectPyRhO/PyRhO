# config.py

import matplotlib as mp
import numpy as np
import warnings
import os
#import glob
import platform
import shutil
import subprocess
import sys
import time

pyVer = sys.version_info

# Set timing function
if pyVer < (3,3):
    if sys.platform == 'win32':
        # On Windows, the best timer is time.clock
        wallTime = time.clock
    else: #sys.platform.startswith('linux') 'cygwin' 'darwin'
        # On most other platforms the best timer is time.time
        wallTime = time.time
else:
    wallTime = time.perf_counter

    
#pyrhoPath = os.path.abspath(pyrho.__file__) # http://stackoverflow.com/questions/247770/retrieving-python-module-path
pyrhoPath = os.path.dirname(os.path.abspath(__file__)) # modulePath
#print(pyrhoPath)


pyrhoNEURONpath = os.path.join(pyrhoPath, 'NEURON')
NMODLfiles = ['RhO3.mod', 'RhO4.mod', 'RhO6.mod']
#HOCfiles = glob.glob("*.hoc")
HOCfiles = [h for h in os.listdir(pyrhoNEURONpath) if h.endswith('.hoc')]
#NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
#print(NMODLfilesIncPath)

##### Output level settings #####
global verbose
if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
    verbose = 1 # Text notification output level [0,1,2]


def createDir(path):
    if os.path.isdir(path):
        return
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise    
    
    
### Set data and figure directories to defaults
if 'dDir' not in vars() or 'dDir' not in globals() or dDir is None:
    dDir = 'data/'
#if not os.path.isdir(dDir):
    createDir(dDir)
    
if 'fDir' not in vars() or 'fDir' not in globals() or fDir is None:
    fDir = 'figs/'
#if not os.path.isdir(fDir):
    createDir(fDir)

def setupGUI(path=None):
    if path is None: # Copy image folders to home directory
        path = os.getcwd()
    createDir('gui')
    # ...
    return
    
def simAvailable(sim):
    sim = sim.lower()
    if sim is 'python': #sim is 'Python' or sim is 'python':
        return True
    elif sim is 'neuron': #sim is 'NEURON' or sim is 'Neuron' or sim is 'neuron':
        return checkNEURON()
    elif sim is 'brian' or sim is 'brian2': #sim is 'Brian' or sim is 'brian' or sim is 'Brian2' or sim is 'brian2':
        return checkBrian()
    else:
        return False
    
    
def checkNEURON(test=False):
    try:
        import neuron as nrn
        found = True
        if verbose > 1:
            print('NEURON module found!')
        if test:
            print('Running suite of tests...')
            nrn.test()
    except ImportError:
        found = False
        if verbose > 1:
            print('NEURON module not found!')
    return found
    
def setupNEURON(path=None):
    # Check for a working NEURON installation...
    #checkNEURON()
    cwd = os.getcwd()
    if not checkNEURON():
        print('Please install NEURON and rerun!')
        return
    
    ### To load mod files:
    # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
    if path is None:
        try:
            path = os.environ['NRN_NMODL_PATH']
            print("Using path from 'NRN_NMODL_PATH' for hoc & nmodl files: ")
        except KeyError:
            path = cwd
            os.environ['NRN_NMODL_PATH'] = path
            #os.putenv(key, value) # Direct assignment is preferable
            print('Using current directory for hoc & nmodl files: ')
    else:
        print("Using suplied path for hoc & nmodl files: ")
    print(path)
        
    # mod files should be compiled in the working directory for NEURON to find them
    # Alternatively, change os.environ['NRN_NMODL_PATH'] ?
    
    #print('NEURON module path: ', path)
    print('Copying mod {} and hoc {} files from {}...'.format(NMODLfiles, HOCfiles, pyrhoNEURONpath))
    if os.path.isdir(path): # Check path
        NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
        for f in NMODLfilesIncPath:
            shutil.copy2(f, path)
        #for f in NMODLfiles:
        #    shutil.copy2(pyrhoNEURONpath, path, f) # >= 3.3
        # for h in os.listdir("pyrhoNEURONpath"): #HOCfiles:
            # if h.endswith(".hoc"):
                # shutil.copy2(pyrhoNEURONpath, path, h)
        HOCfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in HOCfiles]
        for f in HOCfilesIncPath:
            shutil.copy2(f, path)
        #for f in HOCfiles:
        #    shutil.copy2(pyrhoNEURONpath, path, f)
    
    nrnivmodl = shutil.which('nrnivmodl')
    if nrnivmodl is None:
        arch = platform.machine()
        nrnivmodl = os.path.join(path, "nrn", arch, "bin", "nrnivmodl") # Should this be another variable?
    print('Compiling mod files with ', nrnivmodl)
    try:
        cwd = os.getcwd()
        os.chdir(path)
        retcode = subprocess.call(nrnivmodl) #nrn/x86_64/bin/nrnivmodl
        os.chdir(cwd)
        # Use subprocess.run() for Python >= 3.5
        if retcode < 0:
            print('NMODL compilation was terminated by signal', -retcode, file=sys.stderr)
        else:
            print("NMODL compilation returned", retcode) #, file=sys.stderr)
    except OSError as e:
        print('NMODL Compilation failed: ', e, file=sys.stderr)
    return
    
def checkBrian(test=False):
    try:
        import brian2 as br
        found = True
        if verbose > 1:
            print('Brian module found!')
        if test:
            print('Running suite of tests...')
            br.prefs.codegen.target = 'cython' #'numpy #'auto' #'weave' - Python 2 only
            out = br.codegen.cpp_prefs.get_compiler_and_args()
            print(out)
            br.test()
    except ImportError:
        found = False
        if verbose > 1:
            print('Brian module not found!')
        
    return found
    
def setupBrian():
    return
    


    
##### Plot settings #####

figDisplay = 'screen' #'paper'

colours = ['b','g','r','c','m','y','k']
styles = ['-', '--', '-.', ':']

latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mp.rcParams['text.usetex'] = True
#mp.rcParams.update(rcdef) # Set back to defaults

xLabelPos = 0.98

def setFigOutput(figDisplay='screen', width=None):
    golden_ratio  = (1.0 + np.sqrt(5)) / 2.0
    global saveFigFormat
    global addTitles
    global addStimulus
    global eqSize
    
    # http://matplotlib.org/users/customizing.html
    if figDisplay == 'screen':
        dpi = 300 #mp.rcParams['savefig.dpi'] #300
        #figFormat = 'retina'
        saveFigFormat = 'png'    
        if width is None:
            width = 12
        figWidth = width
        figHeight = figWidth / golden_ratio
        tickSize = 14
        labelSize = 16
        legendSize = 12
        titleSize = 'small' #'x-large'
        eqSize = 18
        linewidth = 1.5
        addTitles = True
        addStimulus = True

    elif figDisplay == 'paper':
        #http://www.nature.com/nature/authors/gta/3c_Final_artwork.pdf
        #http://www.plosone.org/static/figureSpecifications
        #http://www.ploscompbiol.org/static/figureSpecifications
        #http://www.frontiersin.org/about/AuthorGuidelines#_FigureandTableGuidelines
        dpi = 300 #1200 # Set DPI 300 - 600
        saveFigFormat = 'png'#'eps'             # Supported formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz.
        #figFormat = 'eps'
        if width is None:
            width = 18/2.54 #5#7
        figWidth = width
        figHeight = figWidth / golden_ratio
        tickSize = 8
        labelSize = 10
        legendSize = 10
        titleSize = 'medium'
        eqSize = 10
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
    #mp.rcParams['figure.autolayout'] = True
    
    
    #return saveFigFormat, addTitles, addStimulus, eqSize

# plotPeakRecovery = False
# plotStateVars = False
# plotKinetics = False

#saveFigFormat, addTitles, addStimulus, eqSize = setFigOutput(display)
setFigOutput(figDisplay)

fancyPlots = False #True

def setFigStyle(fancyPlots=False): # Merge this with setFigOutput
    if fancyPlots:
        try:
            import seaborn
        except ImportError:
            warnings.warn('Seaborn not found - using default plotting scheme.')
        else: # Try another package?
            pass
        finally: # Always do this last
            pass
    else:
        if 'seaborn' in sys.modules:
            seaborn.reset_orig()
        #setFigOutput(figDisplay)
        
def resetPlot():
    if 'seaborn' in sys.modules:
        seaborn.reset_orig()
    setFigOutput(figDisplay)
        
try:
    __IPYTHON__
    import IPython
    #from IPython import display
except NameError:
    pass
else: # and IPython. See also get_ipython()
    figFormat = 'retina' # 'svg'
    IPython.display.set_matplotlib_formats(figFormat)
    #display.set_matplotlib_formats(figFormat)
    if verbose > 1:
        print("Default display figure format set: "+figFormat)


#IPython.core.magic.Magics.basic.pprint()
#%pprint
#%precision %.6g
#np.set_printoptions(precision=6) # Number of decimal places for floats
# Long    Scientific
# 1e-4    1e-5
# 1e15    1e16

