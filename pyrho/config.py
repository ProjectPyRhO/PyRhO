"""General configuration variables and functions"""

from __future__ import print_function, division
import warnings
import logging
import os
#import glob
import platform
import shutil
import subprocess
import sys
import time
import importlib
import pkgutil

import matplotlib as mpl
import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import numpy as np
#import scipy as sp
#import lmfit

__all__ = ['setupGUI', 'simAvailable', 'setupNEURON', 'setupBrian', 'check_package',
           'setFigOutput', 'setFigStyle', 'resetPlot', 'setOutput', 'logger']
# 'wallTime',
# TODO: Place in dict i.e. CONFIG_PARAMS['dDir'] or class with setter methods e.g. to call setOutput
#, 'colours', 'styles', 'verbose', 'dDir', 'fDir', 'DASH_LINE', 'DOUB_DASH_LINE'

pyVer = sys.version_info

# Set timing function
if pyVer < (3, 3):
    if sys.platform == 'win32':
        # On Windows, the best timer is time.clock
        wallTime = time.clock
    else:  # sys.platform.startswith('linux') 'cygwin' 'darwin'
        # On most other platforms the best timer is time.time
        wallTime = time.time
else:
    wallTime = time.perf_counter


if pyVer[0] == 3:
    if pyVer <= (3, 3):
        def check_package(pkg):
            """Test if 'pkg' is available"""
            return importlib.find_loader(pkg) is not None
    elif pyVer >= (3, 4):
        def check_package(pkg):
            """Test if 'pkg' is available"""
            return importlib.util.find_spec(pkg) is not None
elif pyVer[0] == 2:
    def check_package(pkg):
        """Test if 'pkg' is available"""
        return pkgutil.find_loader(pkg) is not None

_DASH_LINE = '-' * 80
_DOUB_DASH_LINE = '=' * 80


# Output level settings #
#global verbose
#if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
verbose = 1  # Text notification output level [0,1,2]
logLevels = [logging.CRITICAL, logging.ERROR, logging.WARNING,
             logging.INFO, logging.DEBUG, logging.NOTSET]

logger = logging.getLogger(__name__)
logging.basicConfig(filename='PyRhO.log', level=logLevels[verbose])

def setOutput(level):
    verbose = level
    logging.setLevel(level)

logger.info('Starting PyRhO')
logger.debug('Initialised Logger')
home = os.path.expanduser('~')
pyrhoPath = os.path.dirname(os.path.abspath(__file__)) # modulePath

pyrhoNEURONpath = os.path.join(pyrhoPath, 'NEURON')
NMODLfiles = ('RhO3c.mod', 'RhO4c.mod', 'RhO6c.mod')  # ['RhO3.mod', 'RhO4.mod', 'RhO6.mod']
#HOCfiles = glob.glob("*.hoc")
HOCfiles = [h for h in os.listdir(pyrhoNEURONpath) if h.endswith('.hoc')]
#NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
NEURONinstallScript = 'install_neuron.sh'


def createDir(path):
    """Create directory"""
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

    """
    Setup the Jupyter notebook GUI.

    Parameters
    ----------
    path : str, optional
        Specify the path for where to set up the GUI folders.
        Defaults to the current working directory

    """

    if path is None: # Copy image folders to home directory
        path = os.path.join(os.getcwd(), GUIdir)

    createDir(path)
    pyrhoGUIpath = os.path.join(pyrhoPath, GUIdir)
    pngFiles = [f for f in os.listdir(pyrhoGUIpath) if f.endswith('.png')]
    for f in pngFiles:
        shutil.copy2(os.path.join(pyrhoGUIpath, f), path)

    return

def simAvailable(simName, test=False):

    """
    Check if a simulator is available.

    Parameters
    ----------
    simName : str {'python', 'neuron', 'brian', 'brian2'}
        Specify the simulator to check
    test : bool, optional
        Specify whether to run the simulator's test suite (the default is False)

    Returns
    -------
    bool
        True if the simulator is available, otherwise False
    """

    simName = simName.lower()
    if simName is 'python':
        return True
    elif simName is 'neuron':
        return checkNEURON(test)
    elif simName is 'brian' or simName is 'brian2':
        return checkBrian(test)
    else:
        return False


def checkNEURON(test=False):
    """Check for NEURON installation and optionally run tests."""

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
        warnings.warn("'NRN_NMODL_PATH' is not set - add 'export NRN_NMODL_PATH=/path/to/NEURON' to your bash profile.")

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


def setupNEURON(path=None): # , NEURONpath=None):

    """
    Setup the NEURON simulator to work with PyRhO.

    Parameters
    ----------
    path : str, optional
        Path to NEURON installation directory containing nrn (and iv) containing hoc and mod files (default=pwd)
    """

    #path : str, optional
    #    Path to PyRhO's working directory containing hoc and mod files (default=pwd)
    #NEURONpath : str, optional
    #    Path to NEURON installation directory containing nrn (and iv)

    cwd = os.getcwd()
    # Boiler plates to expand '~'
    if path is not None:
        path = os.path.expanduser(path)
    #if NEURONpath is not None:
    #    NEURONpath = os.path.expanduser(NEURONpath)

    # TODO: Add instructions to compile/install for Python 2
    if not checkNEURON():   # Check for a working NEURON installation...
        if False: # TODO: Make NEURONinstallScript work with subprocess
            if sys.platform == 'win32':
                # TODO: Create bat/cmd script for compiling on Windows
                warnings.warn('Compilation on Windows is not yet supported - please install NEURON manually and rerun!')
                return
            else:
                try:
                    if path is None: #NEURONpath
                        ans = input('The `path` argument was not specified - please enter one or press `Enter` to use the current working directory ({}): '.format(cwd))
                        if ans == '' or ans is None:
                            path = cwd
                        else:
                            path = os.path.expanduser(ans)
                        # Set NRN_NMODL_PATH environment variable to path
                        #print('Please specify an installation path for NEURON with `NEURONpath` and run again')
                        #return
                    #NEURONscriptPath = os.path.join(home, 'NEURON')
                    #NEURONscriptPath = pyrhoNEURONpath

                    NEURONscriptIncPath = os.path.join(pyrhoNEURONpath, NEURONinstallScript)
                    exitcode = subprocess.call([NEURONscriptIncPath, path], shell=True) #NEURONpath # .check_call
                except:
                    createDir(path)
                    shutil.copy2(os.path.join(pyrhoNEURONpath, NEURONinstallScript), path) #cwd
                    print('Unable to install NEURON - please install manually with the script copied to {}.'.format(path)) #cwd
        else:
            if path is None:
                path = cwd
            createDir(path)
            shutil.copy2(os.path.join(pyrhoNEURONpath, NEURONinstallScript), path) #cwd
            print('NEURON must be compiled from source to work with Python 3. ')
            print('Please use the script `{}` in `{}` and then rerun setupNEURON.'.format(NEURONinstallScript, path)) #cwd
            print("E.g.: ./{} /abs/path/to/install/NEURON/".format(NEURONinstallScript))
            return

    ### To load mod files:
    # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
    if path is None:
        try:
            path = os.environ['NRN_NMODL_PATH']
            print('Using path from `NRN_NMODL_PATH` for hoc & nmodl files: ', end='')
        except KeyError:
            path = cwd
            os.environ['NRN_NMODL_PATH'] = path
            #os.putenv(key, value) # Direct assignment is preferable
            print('Using current directory for hoc & nmodl files: ', end='')
    else:
        print('Using suplied path for hoc & nmodl files: ', end='')
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


    if pyVer >= (3, 3):
        nrnivmodl = shutil.which('nrnivmodl')
    else:
        from distutils.spawn import find_executable
        nrnivmodl = find_executable('nrnivmodl')

    if nrnivmodl is None:
        arch = platform.machine()
        #if NEURONpath is not None:
        #    nrnivmodl = os.path.join(NEURONpath, "nrn", arch, "bin", "nrnivmodl")
        #else:
        nrnivmodl = os.path.join(path, 'nrn', arch, 'bin', 'nrnivmodl') # Try the hoc/mod path

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
        print('Try setting the NEURON directory as an environment variable or passing it as the `path` argument.') #NEURONpath
    return

def checkBrian(test=False):
    """Check for the Brian2 simulator"""
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
    """Setup the Brian2 simulator."""
    if not checkBrian():
        try:
            import pip
            pip.main(['install', 'brian2']) # pip install brian2
        except:
            warnings.warn('Unable to install Brian2 - please run `pip install brian2` from the terminal')
    return




##### Plot settings #####

figDisplay = 'screen' #'paper'

colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
styles = ['-', '--', '-.', ':']

latexInstalled = False # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mpl.rcParams['text.usetex'] = True
#mpl.rcParams.update(rcdef) # Set back to defaults

xLabelPos = 0.98

def setFigOutput(figDisplay='screen', width=None):

    """
    Set figure plotting options.

    Parameters
    ----------
    figDisplay : str, {'screen', 'paper'}, optional
        Specify the output medium to set figure properties for
    width : int or float, optional
        Figure width in inches. Default depends on figDisplay.
    """

    golden_ratio = (1.0 + np.sqrt(5)) / 2.0
    global saveFigFormat
    global addTitles
    global addStimulus
    global eqSize

    # http://matplotlib.org/users/customizing.html
    if figDisplay == 'screen':
        dpi = 80 #300 #mpl.rcParams['savefig.dpi'] #300
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

    elif figDisplay == 'nb' or figDisplay == 'notebook':
        mpl.use('nbagg') # Switch backend - TODO: Move this to __init__ before import matplotlib.pyplot
        dpi = 80
        saveFigFormat = 'png'
        if width is None:
            width = 12
        figWidth = width
        figHeight = figWidth / golden_ratio
        tickSize = 14
        labelSize = 12
        legendSize = 8
        titleSize = 'small' #'x-large'
        eqSize = 12
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
                                                # https://github.com/matplotlib/matplotlib/issues/5945

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

#global fancyPlots
fancyPlots = False #True

def setFigStyle(): #fancyPlots=False): # Merge this with setFigOutput
    """Set figure style"""
    global fancyPlots
    global colours
    if fancyPlots:
        try:
            import seaborn as sns
            #cp = sns.color_palette()
            colours = sns.color_palette()
        except ImportError:
            warnings.warn('Seaborn not found - using default plotting scheme.')
        else: # Try another package?
            pass
        finally: # Always do this last
            fancyPlots = True
    else:
        if 'seaborn' in sys.modules:
            sns.reset_orig()
        fancyPlots = False
        #setFigOutput(figDisplay)

def resetPlot():
    """Reset figure style"""
    global fancyPlots
    if 'seaborn' in sys.modules and fancyPlots:
        sns.reset_orig()
        fancyPlots = False
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
