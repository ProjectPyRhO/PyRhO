"""A Python module for fitting, characterising and simulating rhodopsin photocurrents."""
#__doc__ =
# Main module file for PyRhO


#from pkg_resources import get_distribution, DistributionNotFound
import logging
import os
import platform

import lmfit
# Necessary?
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
import scipy as sp

# Place all submodule functions and variables into namespace
from pyrho import config
from pyrho.config import *
from pyrho.config import _DASH_LINE, _DOUB_DASH_LINE
from pyrho.expdata import *
from pyrho.fitting import *
from pyrho.models import *
from pyrho.parameters import *
from pyrho.protocols import *
from pyrho.simulators import *
from pyrho.utilities import *

# from pyrho.jupytergui import *

# TODO
#__all__ = ['config', 'utilities', 'parameters', 'expdata',
#           'models', 'protocols', 'simulators', 'fitting', 'jupytergui']

__project__ = 'pyrho'
# http://stackoverflow.com/questions/17583443/what-is-the-correct-way-to-share-package-version-with-setup-py-and-the-package
#__version__ = pkg_resources.get_distribution(__project__).version #'0.8.0'
try:
    _dist = pkg_resources.get_distribution(__project__)
    _distLoc = os.path.normcase(_dist.location)  # Normalise case for Windows
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(_distLoc, __project__)):
        # not installed, but there is another version that *is*
        raise pkg_resources.DistributionNotFound
except pkg_resources.DistributionNotFound:
    __version__ = __project__ + '-' + '(local)'
else:
    __version__ = _dist.version

# TODO: Refactor environment checks.
# In notebooks, __IPYTHON__ and get_ipython() are defined.
#IPY = False
#try:
#    import IPython
#    IPY = True
#except ImportError:
#    IPY = False

# IPY = check_package('IPython')
try:
    __IPYTHON__
    IPY = True
except NameError:
    IPY = False
# https://github.com/ipython/ipython/issues/9791
#try:
#    from IPython import get_ipython
#    ip = get_ipython()
#    IPY = True
#except ImportError, NameError:
#    IPY = False

if IPY:
    from pyrho.jupytergui import *

if __name__ == '__main__':
    try:
        __IPYTHON__
    except NameError:
        pass
    else:  # and IPython. See also get_ipython()
        print('Loading IPython GUI!')
        loadGUI()

# TODO Move everything except imports elsewhere

# TODO: Rename global models, protocols and simulators dictionaries

def run(mods=6, prots='step', sims='Python', plot=True):
    """
    Run all protocols on a list of models with default parameters.

    Parameters
    ----------
    mods : int, str, list
        Individual or list of integers or strings specifying the models to run
        e.g. [3, 4, 6], 3, '4', ['4', '6'].

    prots : str, list
        String or list of strings of the names of protocols to run (default: 'step').
        e.g. 'ramp', ['chirp', 'sinusoid'].

    sims : str, list
        List of strings of the names of simulators to use (default: 'Python').
        e.g. 'Brian', ['Python', 'NEURON'].
    """

    if not isinstance(mods, (list, tuple)):
        assert isinstance(mods, (int, str))
        mods = [mods]  # ints or strs
    mods = [str(m) for m in mods]

    if not isinstance(prots, (list, tuple)):
        assert isinstance(prots, str)  # TODO: or Protocol class
        prots = [prots]

    if not isinstance(sims, (list, tuple)):
        assert isinstance(sims, str)
        sims = [sims]

    results = {simulator: {protocol: {} for protocol in prots} for simulator in sims}
    for model in mods:
        rho = models[model]()  # Select generative model
        for protocol in prots:
            prot = protocols[protocol]()  # Select simulation protocol
            for simulator in sims:
                sim = simulators[simulator](prot, rho)
                print(f"\nUsing {simulator} to run Protocol '{protocol}' on the {model}-state model...")
                print(_DASH_LINE, '\n')
                results[simulator][protocol][model] = sim.run()
                if plot:
                    sim.plot()
                print("\nFinished!")
                print(_DOUB_DASH_LINE, '\n\n')
    return results


def print_versions():
    """Display version information for PyRhO and its dependencies."""
    # import platform
    print("Python version:        ", platform.python_version())
    if IPY:
        try:
            import IPython

            # __IPYTHON__
            print("IPython version:       ", IPython.__version__)
        except ImportError:  # IPython not found
            pass
    #deps = [numpy, scipy, matplotlib, lmfit, warnings, os, pickle, collections, platform]
    #depsGUI = [IPython, ast, base64]
    #for mod in dependencies:
    #    print("{} version: {}".format(mod, mod.__version__))
    #import numpy
    print("NumPy version:         ", np.__version__)
    #import scipy
    print("SciPy version:         ", sp.__version__)
    #import matplotlib
    print("Matplotlib version:    ", mpl.__version__)
    #import lmfit
    print("Lmfit version:         ", lmfit.__version__)

    print("PyRhO version:         ", __version__)

    try:
        import neuron
        print(f"NEURON version:         {neuron.__version__}")
    except ImportError:
        pass

    try:
        import brian2
        print(f"Brian version:          {brian2.__version__}")
    except ImportError:
        pass


def get_versions_table():
    """Display version information for PyRhO and its dependencies."""

    table = ["     Module | Version    "]
    table.append("============|============")
    table.append(f"{'Python':>11} | {platform.python_version()}")
    if IPY:
        try:
            import IPython

            # __IPYTHON__
            table.append(f"{'IPython':>11} | {IPython.__version__}")
        except ImportError:  # IPython not found
            pass
    table.append(f"{'NumPy':>11} | {np.__version__}")
    table.append(f"{'SciPy':>11} | {sp.__version__}")
    table.append(f"{'Matplotlib':>11} | {mpl.__version__}")
    table.append(f"{'Lmfit':>11} | {lmfit.__version__}")
    table.append(f"{'PyRhO':>11} | {__version__}")
    try:
        import neuron
        table.append(f"{'NEURON':>11} | {neuron.__version__}")
    except ImportError:
        pass
    try:
        import brian2
        table.append(f"{'Brian':>11} | {brian2.__version__}")
    except ImportError:
        pass
    table.append("-------------------------\n")

    return "\n".join(table)


logLevels = [logging.CRITICAL, logging.ERROR, logging.WARNING,
             logging.INFO, logging.DEBUG, logging.NOTSET]


def setOutput(logger, level):
    verbose = level
    logger.setLevel(level)

logging.basicConfig(filename='PyRhO.log', level=logLevels[config.verbose],
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%y-%m-%d %H:%M:%S', filemode='w')

logger = logging.getLogger(__name__)
setOutput(logger, config.verbose)

logger.info('Starting PyRhO')
logger.debug('Initialised Logger')
logger.info('Module versions table\n' + get_versions_table())
