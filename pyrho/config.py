"""General configuration variables and functions"""

import importlib
import logging
import os
import pkgutil
#import glob
import platform
import shutil
import subprocess
import sys
import time
import warnings

import matplotlib as mpl
#import matplotlib.pyplot as plt
import numpy as np

#from pyrho.__init__ import printVersions

__all__ = ['setupGUI', 'simAvailable', 'setupNEURON', 'setupBrian', 'check_package',
           'setFigOutput', 'setFigStyle', 'resetPlot']
# 'wall_time',
# TODO: Place in dict i.e. CONFIG_PARAMS['dDir'] or class with setter methods e.g. to call set_output
#, 'colours', 'styles', 'verbose', 'dDir', 'fDir', 'DASH_LINE', 'DOUB_DASH_LINE'

pyVer = sys.version_info

if pyVer < (2, 7):
    raise RuntimeError('Only Python versions >= 2.7 are supported')

wall_time = time.perf_counter

def check_package(pkg):
    """Test if 'pkg' is available"""
    return importlib.util.find_spec(pkg) is not None

_LINE_LENGTH = 80
_DASH_LINE = '-' * _LINE_LENGTH
_DOUB_DASH_LINE = '=' * _LINE_LENGTH


# Output level settings #
#global verbose
#if 'verbose' not in vars() or 'verbose' not in globals() or verbose is None:
verbose = 1  # Text notification output level [0,1,2]


home = os.path.expanduser('~')
pyrhoPath = os.path.dirname(os.path.abspath(__file__))  # modulePath

pyrhoNEURONpath = os.path.join(pyrhoPath, 'NEURON')
NMODLfiles = ('RhO3c.mod', 'RhO4c.mod', 'RhO6c.mod')  # ['RhO3.mod', 'RhO4.mod', 'RhO6.mod']
#HOCfiles = glob.glob("*.hoc")
HOCfiles = [h for h in os.listdir(pyrhoNEURONpath) if h.endswith('.hoc')]
#NMODLfilesIncPath = [os.path.join(pyrhoNEURONpath, f) for f in NMODLfiles]
NEURONinstallScript = 'install_neuron.sh'


# Set data and figure directories to defaults
if 'dDir' not in vars() or 'dDir' not in globals():  # or dDir is None:
    dDir = 'data' + os.sep
    # createDir(dDir)
    os.makedirs(dDir, exist_ok=True)

if 'fDir' not in vars() or 'fDir' not in globals():  # or fDir is None:
    fDir = 'figs' + os.sep
    # createDir(fDir)
    os.makedirs(fDir, exist_ok=True)

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

    if path is None:  # Copy image folders to home directory
        path = os.path.join(os.getcwd(), GUIdir)

    # createDir(path)
    os.makedirs(path, exist_ok=True)
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
        Specify whether to run the simulator's test suite (default: False)

    Returns
    -------
    bool
        True if the simulator is available, otherwise False
    """

    simName = simName.lower()
    if simName == 'python':
        return True
    elif simName == 'neuron':
        return checkNEURON(test)
    elif simName == 'brian' or simName == 'brian2':
        return checkBrian(test)
    else:
        return False


def checkNEURON(test=False):
    """Check for NEURON installation and optionally run tests."""

    if 'NRN_NMODL_PATH' in os.environ:
        nmodlPath = os.environ['NRN_NMODL_PATH']
        if verbose > 1:
            print(f"'NRN_NMODL_PATH' environment variable is set: {nmodlPath}")
            modList = [file for file in os.listdir(nmodlPath) if file.endswith(".mod")]
            hocList = [file for file in os.listdir(nmodlPath) if file.endswith(".hoc")]
            print("MOD files found: ", modList)
            print("HOC files found: ", hocList)
            if not set(NMODLfiles).issubset(set(modList)):
                warnings.warn(f'Missing mod files in {nmodlPath}: {set(NMODLfiles) - set(modList)}')
            if not set(HOCfiles).issubset(set(hocList)):
                warnings.warn(f'Missing hoc files in {nmodlPath}: {set(HOCfiles) - set(hocList)}')

    else:
        warnings.warn("'NRN_NMODL_PATH' is not set - add it to environment "
                      "e.g. add 'export NRN_NMODL_PATH=/path/to/NEURON' "
                      "to your bash profile.")

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


def setupNEURON(path=None):  # , NEURONpath=None):
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
        if path is None:
            path = cwd
        os.makedirs(path, exist_ok=True)
        shutil.copy2(os.path.join(pyrhoNEURONpath, NEURONinstallScript), path) #cwd
        print('NEURON must be compiled from source to work with Python 3. ')
        print(f'Please use the script `{NEURONinstallScript}` in `{path}` and then rerun setupNEURON.') #cwd
        print(f"E.g.: ./{NEURONinstallScript} /abs/path/to/install/NEURON/")
        return

    # To load mod files:
    # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
    if path is None:
        try:
            path = os.environ['NRN_NMODL_PATH']
            print('Using path from `NRN_NMODL_PATH` for hoc & nmodl files: ', end='')
        except KeyError:
            path = cwd
            os.environ['NRN_NMODL_PATH'] = path
            # os.putenv(key, value) # Direct assignment is preferable
            print('Using current directory for hoc & nmodl files: ', end='')
    else:
        print('Using suplied path for hoc & nmodl files: ', end='')
    print(path)

    # mod files should be compiled in the working directory for NEURON to find them
    # Alternatively, change os.environ['NRN_NMODL_PATH'] ?

    #print('NEURON module path: ', path)
    print(f'Copying mod {NMODLfiles} and hoc {HOCfiles} files from {pyrhoNEURONpath}...')
    if os.path.isdir(path):  # Check path
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
        #if NEURONpath is not None:
        #    nrnivmodl = os.path.join(NEURONpath, "nrn", arch, "bin", "nrnivmodl")
        #else:
        nrnivmodl = os.path.join(path, 'nrn', arch, 'bin', 'nrnivmodl')  # Try the hoc/mod path

    print('Compiling mod files with ', nrnivmodl)
    try:
        cwd = os.getcwd()
        os.chdir(path)
        retcode = subprocess.call(nrnivmodl)  # nrn/x86_64/bin/nrnivmodl
        os.chdir(cwd)
        # Use subprocess.run() for Python >= 3.5
        if retcode < 0:
            print('NMODL compilation was terminated by signal', -retcode, file=sys.stderr)
        else:
            print("NMODL compilation returned", retcode)  # , file=sys.stderr)
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
            br.prefs.codegen.target = 'cython'  # 'numpy #'auto' #'weave' - Python 2 only
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
            pip.main(['install', 'brian2'])  # pip install brian2
        except:
            warnings.warn('Unable to install Brian2 - please run `pip install brian2` from the terminal')
    return


##### Plot settings #####
# TODO: Simplify style/colour cycling with axess.prop_cycle
# http://matplotlib.org/users/whats_new.html#added-axes-prop-cycle-key-to-rcparams

figDisplay = 'screen'  # 'paper'
colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
styles = ['-', '--', '-.', ':']
xLabelPos = 0.98
latexInstalled = False  # Change this!
# http://nipunbatra.github.io/2014/08/latexify/
# http://www.tex.ac.uk/ctan/macros/latex/contrib/pythontex/pythontex_install.py
if latexInstalled:
    mpl.rcParams['text.usetex'] = True
#mpl.rcParams.update(rcdef) # Set back to defaults


def setFigOutput(figDisplay='screen', width=None, height=None):

    """
    Set figure plotting options.

    Parameters
    ----------
    figDisplay : str, {'screen', 'paper'}, optional
        Specify the output medium to set figure properties for
    width : int or float, optional
        Figure width in inches. Default depends on figDisplay.
    """

    global saveFigFormat
    global addTitles
    global addStimulus
    global eqSize

    # http://matplotlib.org/users/customizing.html
    if figDisplay == 'screen':
        dpi = 80  # 300 #mpl.rcParams['savefig.dpi'] #300
        #figFormat = 'retina'
        saveFigFormat = 'png'
        if width is None:
            width = 12
        tickSize = 14
        labelSize = 16
        legendSize = 12
        titleSize = 'small'  # 'x-large'
        eqSize = 18  # General annotations
        linewidth = 1.5
        markerSize = 6  # Default: 6
        addTitles = True
        addStimulus = True

    elif figDisplay == 'nb' or figDisplay == 'notebook':
        mpl.use('nbagg')  # Switch backend - TODO: Move this to __init__ before import matplotlib.pyplot
        dpi = 80
        saveFigFormat = 'png'
        if width is None:
            width = 12
        tickSize = 14
        labelSize = 12
        legendSize = 8
        titleSize = 'small'  # 'x-large'
        eqSize = 12  # General annotations
        linewidth = 1.5
        markerSize = 6  # Default: 6
        addTitles = True
        addStimulus = True

    elif figDisplay == 'paper':
        #http://www.nature.com/nature/authors/gta/3c_Final_artwork.pdf
        #http://www.plosone.org/static/figureSpecifications
        #http://www.ploscompbiol.org/static/figureSpecifications
        #http://www.frontiersin.org/about/AuthorGuidelines#_FigureandTableGuidelines
        dpi = 300 #1200 # Set DPI 300 - 600
        saveFigFormat = 'eps'  # Supported formats: eps, pdf, pgf, png, ps, raw, rgba, svg, svgz.
        #figFormat = 'eps'
        if width is None:
            width = 18/2.54
        tickSize = 8
        labelSize = 10
        legendSize = 10
        titleSize = 'medium'
        eqSize = 12  # General annotations
        linewidth = 2
        markerSize = 6
        addTitles = False
        addStimulus = True

        mpl.rcParams['savefig.dpi'] = dpi       # http://nbviewer.ipython.org/gist/minrk/3301035
                                                # http://wiki.scipy.org/Cookbook/Matplotlib/AdjustingImageSize
                                                # https://github.com/matplotlib/matplotlib/issues/5945

    else:
        warnings.warn('Warning: Unknown display type selected'
                      ' - using matplotlib defaults!')

    if height is None:
        golden_ratio = (1.0 + np.sqrt(5)) / 2.0
        height = width / golden_ratio

    # http://matplotlib.org/users/customizing.html
    mpl.rcParams['figure.figsize'] = (width, height)  # pylab.rcParams
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
fancyPlots = False


def setFigStyle():  # fancyPlots=False): # TODO: Merge this with setFigOutput
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
        else:  # Try another package?
            pass
        finally:  # Always do this last
            fancyPlots = True
    else:
        if 'seaborn' in sys.modules:
            sns.reset_orig()
        fancyPlots = False
        #setFigOutput(figDisplay)


def resetPlot():
    """Reset figure style."""
    global fancyPlots
    if 'seaborn' in sys.modules and fancyPlots:
        try:
            import seaborn as sns
            sns.reset_orig()
        except ImportError:
            sns = None
        fancyPlots = False
    setFigOutput(figDisplay)

try:
    __IPYTHON__
    import IPython

    #from IPython import display
except NameError:
    pass
else:  # and IPython. See also get_ipython()
    figFormat = 'retina'  # 'svg'
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
