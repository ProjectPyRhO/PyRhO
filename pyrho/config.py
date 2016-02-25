# config.py

import matplotlib as mpl
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


home = os.path.expanduser('~')
pyrhoPath = os.path.dirname(os.path.abspath(__file__)) # modulePath

pyrhoNEURONpath = os.path.join(pyrhoPath, 'NEURON')
NMODLfiles = ['RhO3c.mod', 'RhO4c.mod', 'RhO6c.mod'] #['RhO3.mod', 'RhO4.mod', 'RhO6.mod']
#HOCfiles = glob.glob("*.hoc")
HOCfiles = [h for h in os.listdir(pyrhoNEURONpath) if h.endswith('.hoc')]
#NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
NEURONinstallScript = 'install_neuron.sh'

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
    dDir = 'data' + os.sep
    createDir(dDir)
    
if 'fDir' not in vars() or 'fDir' not in globals() or fDir is None:
    fDir = 'figs' + os.sep
    createDir(fDir)

GUIdir = 'gui'
def setupGUI(path=None):
    if path is None: # Copy image folders to home directory
        path = os.path.join(os.getcwd(), GUIdir)
        
    createDir(path)
    pyrhoGUIpath = os.path.join(pyrhoPath, GUIdir)
    pngFiles = [f for f in os.listdir(pyrhoGUIpath) if f.endswith('.png')]
    for f in pngFiles:
        shutil.copy2(os.path.join(pyrhoGUIpath, f), path)
    
    return
    
def simAvailable(sim):
    sim = sim.lower()
    if sim is 'python':
        return True
    elif sim is 'neuron':
        return checkNEURON()
    elif sim is 'brian' or sim is 'brian2':
        return checkBrian()
    else:
        return False
    
    
def checkNEURON(test=False):
    """Check for NEURON installation and optionally run tests"""
    
    if 'NRN_NMODL_PATH' in os.environ:
        nmodlPath = os.environ['NRN_NMODL_PATH']
        if verbose > 1:
            print("'NRN_NMODL_PATH' environment variable is set: {}".format(nmodlPath))
            modList = [file for file in os.listdir(nmodlPath) if file.endswith(".mod")]
            hocList = [file for file in os.listdir(nmodlPath) if file.endswith(".hoc")]
            #for file in os.listdir(nmodlPath):
            #    if file.endswith(".hoc"):
            #        print(file)
            print("MOD files found: ", modList)
            print("HOC files found: ", hocList)
            if not set(NMODLfiles).issubset(set(modList)):
                warnings.warn('Missing mod files in {}: {}'.format(nmodlPath, set(NMODLfiles) - set(modList)))
            if not set(HOCfiles).issubset(set(hocList)):
                warnings.warn('Missing hoc files in {}: {}'.format(nmodlPath, set(HOCfiles) - set(hocList)))
            
    else:
        warnings.warn("'NRN_NMODL_PATH' is not set - add e.g. 'export NRN_NMODL_PATH=/path/to/NEURON' to your bash profile.")
    
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
    
    
def setupNEURON(path=None, NEURONpath=None):
    """Setup the NEURON simulator to work with PyRhO
        path        := Path to PyRhO's working directory containing hoc and mod files (default=pwd)
        NEURONpath  := Path to NEURON installation directory containing nrn (and iv)"""
        
    cwd = os.getcwd()
    # Boiler plates to expand '~'
    if path is not None:
        path = os.path.expanduser(path)
    if NEURONpath is not None:
        NEURONpath = os.path.expanduser(NEURONpath)
        
    if not checkNEURON():   # Check for a working NEURON installation...
        if False: # TODO: Make NEURONinstallScript work with subprocess
            if sys.platform == 'win32':
                # TODO: Create cmd script for compiling on Windows
                warnings.warn('Compilation on Windows is not yet supported - please install NEURON manually and rerun!')
                return
            else:
                try:
                    if NEURONpath is None:
                        print("Please specify an installation path for NEURON with 'NEURONpath' and run again")
                        return
                    #NEURONscriptPath = os.path.join(home, 'NEURON')
                    #NEURONscriptPath = pyrhoNEURONpath
                    NEURONscriptIncPath = os.path.join(pyrhoNEURONpath, NEURONinstallScript)
                    exitcode = subprocess.call([NEURONscriptIncPath, NEURONpath], shell=True) # .check_call
                except:
                    shutil.copy2(os.path.join(NEURONscriptPath, NEURONinstallScript), cwd)
                    print('Unable to install NEURON - please install manually with the script copied to {}.'.format(cwd))
        else:    
            shutil.copy2(os.path.join(NEURONscriptPath, NEURONinstallScript), cwd)
            print("NEURON must be compiled from source to work with Python 3. Please use the script '{}' copied to '{}' and then rerun setupNEURON.".format(NEURONinstallScript, cwd))
            print("E.g.: ./{} /abs/path/to/install/NEURON/".format(NEURONinstallScript))
    
    ### To load mod files:
    # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
    if path is None:
        try:
            path = os.environ['NRN_NMODL_PATH']
            print("Using path from 'NRN_NMODL_PATH' for hoc & nmodl files: ", end='')
        except KeyError:
            path = cwd
            os.environ['NRN_NMODL_PATH'] = path
            #os.putenv(key, value) # Direct assignment is preferable
            print('Using current directory for hoc & nmodl files: ', end='')
    else:
        print("Using suplied path for hoc & nmodl files: ", end='')
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
        if NEURONpath is not None:
            nrnivmodl = os.path.join(NEURONpath, "nrn", arch, "bin", "nrnivmodl")
        else:
            nrnivmodl = os.path.join(path, "nrn", arch, "bin", "nrnivmodl") # Try the hoc/mod path
    
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
        print('Try setting the NEURON directory as an environment variable or passing it as the NEURONpath argument.')
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
    if not checkBrian():
        import pip
        pip.main(['install', 'brian2']) # pip install brian2
    return
    


    
##### Plot settings #####

figDisplay = 'screen' #'paper'

colours = ['b','g','r','c','m','y','k']
styles = ['-', '--', '-.', ':']

latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mpl.rcParams['text.usetex'] = True
#mpl.rcParams.update(rcdef) # Set back to defaults

xLabelPos = 0.98

def setFigOutput(figDisplay='screen', width=None):
    golden_ratio  = (1.0 + np.sqrt(5)) / 2.0
    global saveFigFormat
    global addTitles
    global addStimulus
    global eqSize
    
    # http://matplotlib.org/users/customizing.html
    if figDisplay == 'screen':
        dpi = 300 #mpl.rcParams['savefig.dpi'] #300
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
        markerSize = 6 # Default: 6
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
        markerSize = 6
        addTitles = False
        addStimulus = True

        mpl.rcParams['savefig.dpi'] = dpi        # http://nbviewer.ipython.org/gist/minrk/3301035
                                                # http://wiki.scipy.org/Cookbook/Matplotlib/AdjustingImageSize
        
    else:
        warnings.warn('Warning: Unknown display type selected - using matplotlib defaults!')
    
    # http://matplotlib.org/users/customizing.html
    mpl.rcParams['figure.figsize'] = (figWidth, figHeight) #pylab.rcParams
    mpl.rcParams['figure.dpi'] = dpi
    #mpl.rcParams['savefig.dpi'] = dpi           # http://nbviewer.ipython.org/gist/minrk/3301035
    mpl.rcParams['axes.labelsize'] = labelSize
    mpl.rcParams['ytick.labelsize'] = tickSize
    mpl.rcParams['xtick.labelsize'] = tickSize
    mpl.rcParams['legend.fontsize'] = legendSize
    mpl.rcParams['axes.titlesize'] = titleSize
    mpl.rcParams['lines.linewidth'] = linewidth
    mpl.rcParams['lines.markersize'] = markerSize
    mpl.rcParams['font.size'] = eqSize
    #mpl.rcParams['figure.autolayout'] = True
    
    #return saveFigFormat, addTitles, addStimulus, eqSize

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

