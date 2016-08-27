"""General utility functions used throughout PyRhO"""

from __future__ import print_function, division

import os
import copy
import warnings
import pickle
from string import Template

import numpy as np
import matplotlib.pyplot as plt

from pyrho import config
#from pyrho.utilities import wallTime

__all__ = ['Timer', 'saveData', 'loadData', 'getExt', 'getIndex', 'calcV1',
           'lam2rgb', 'irrad2flux', 'flux2irrad', 'times2cycles', 'cycles2times',
           'plotLight', 'setCrossAxes', 'round_sig']
# 'printParams', 'compareParams', 'texIt', 'expDecay', 'biExpDecay', 'biExpSum', 'calcgbar'


class Timer:
    """
    Class for timing blocks of code.
    http://preshing.com/20110924/timing-your-code-using-pythons-with-statement/
    
    Examples
    --------
    >>> with Timer() as t:
    >>>    run code...
    Execution took <t>s)
    """
    #interval = 0
    def __init__(self):
        self.start = 0
        self.end = 0
        self.interval = 0

    def __enter__(self):
        self.start = config.wallTime() #time.clock()
        return self

    def __exit__(self, *args):
        self.end = config.wallTime() #time.clock()
        self.interval = self.end - self.start
        print('{:.3g}s'.format(self.interval))

    def __str__(self):
        return '{:.3g}s'.format(self.interval)

    def reset(self):
        """Reset timer to 0"""
        self.start = 0
        self.end = 0
        self.interval = 0


# The following two functions are only used in fitModel
def printParams(params):
    """
    Print an LMFIT Parameters object
    """
    vd = params.valuesdict()
    report = '------------------------\n'
    report += '       Parameters\n'
    report += '------------------------\n'
    for k, v in vd.items():
        if isinstance(v, (int, float, complex)):
            report += '{:>7} = {:8.3g}\n'.format(k, v)
        else: # Check for bool?
            report += '{:>7} = {:8}\n'.format(k, str(v))
    report += '========================\n'
    print(report)


def compareParams(origParams, newParams):
    """
    Print two sets of LMFIT Parameters with the percentage (or absolute) difference
    """
    ovd = origParams.valuesdict()
    nvd = newParams.valuesdict()
    report  = '--------------------------------------------\n'
    report += '          Original        New    Change     \n'
    report += '--------------------------------------------\n'
    for k, nv in nvd.items():
        ov = ovd[k]
        if origParams[k].vary:
            if isinstance(nv, (int, float, complex)):
                if ov > 1e-4: #ov != 0:
                    report += '{:>7} = {:8.3g} --> {:8.3g} ({:+.3g}%)\n'.format(k, ov, nv, (nv-ov)*100/ov)
                else:
                    report += '{:>7} = {:8.3g} --> {:8.3g} (Abs: {:+.3g})\n'.format(k, ov, nv, nv-ov)
            else: # Check for bool?
                report += '{:>7} = {:8}\n'.format(k, str(nv))
        else:
            report += '{:>7} = {:8.3g} --> {:8.3g}   ~ Fixed ~\n'.format(k, ov, nv)
    report += '============================================\n'
    print(report)


def texIt(texString):
    """Function to add '$' signs around a (La)TeX string"""
    texTemplate = Template('${content}$')
    return texTemplate.substitute(content=texString)


def saveData(data, pkl, path=None):
    """
    Pickle data in ``dDir`` or as specified in optional ``path`` argument. 
    
    Parameters
    ----------
    data : object
        Variable to pickle. 
    pkl : str
        Filename to save (without extension). 
    path : str, optional
        Optionally specify a path for where to save the file (default=``config.dDir``). 
    
    Returns
    -------
    str
        Pickle filename (including extension). 
    """
    
    # if pkl is None:
        # pkl = data.__name__
    if path is None:
        path = config.dDir
    pklFile = os.path.join(path, pkl+".pkl")
    with open(pklFile, 'wb') as fh:
        pickle.dump(data, fh)
    if config.verbose > 0:
        print("Data saved to disk: {}".format(pklFile))
    return pklFile


def loadData(pkl, path=None):
    """
    Load a pickled dataSet. 

    The function searches paths in the following order:
        1) Optional ``path`` argument
        2) Present working directory
        3) ``config.dDir``
    
    Parameters
    ----------
    pkl : str
        Pickle filename (optionally including the extension)
    path : str, optional
        Optionally specify a path for where to find the file. 
    
    Returns
    -------
    dataSet
        The unpickled contents of ``pkl``
    """
    
    if pkl.lower().endswith('.pkl'):
        pklFile = pkl
    else:
        pklFile = pkl + '.pkl'
    if path is None:
        if os.path.isfile(pklFile):
            pass
        else:
            pklFile = os.path.join(config.dDir, pklFile)
    else:
        pklFile = os.path.join(path, pklFile)
    with open(pklFile, 'rb') as fh:
        dataSet = pickle.load(fh)
    return dataSet


def getExt(vector, ext='max'):
    """
    Get the most extreme value from a vector. 
    
    Parameters
    ----------
    vector : ndarray
        Array of values. 
    ext : {'max', 'min'}, optional
        Specify whether to find the maximum or minimum value. 
    
    Returns
    -------
    tuple
        (mVal := the extreme value, mInd := the index of mVal)
    """
    
    if ext == 'max':
        mVal = max(vector)
    elif ext == 'min':
        mVal = min(vector)
    mInd = np.searchsorted(vector, mVal)
    return mVal, mInd


def getIndex(valList, val):
    """
    Return the index of val in valList.
    
    This handles lists containing ``None``.
    
    Parameters
    ----------
    valList : list
        List of values to search through. 
    val : int or float
        Value to index. 
    
    Returns
    -------
    int
        The index of ``val`` in ``valList`` or ``None`` if not found. 
    """

    # valList types: list, array, number
    #   +/- None
    # val types: list, array, number
    #   +/- None

    # Vs=[-100,-70,-40,-10]
    # V = -70
    # np.where(np.isclose(Vs,V))[0] # Array of indices

    # locList = copy.copy(valList)
    # if isinstance(valList, (list, tuple)):
        # try:
            # ind = valList.index(val)
        # except:
            # pass
    # elif isinstance(valList, (np.ndarray, np.generic)):
        # try:
            # cl = np.isclose(valList, val)
            # ind = np.searchsorted(cl, True)
            ###ind = np.searchsorted(cl, True, equal_nan=True)
        # except:
            # pass
        # else:
            # locList = valList.tolist()
    # elif isinstance(valList, (int, long, float, complex)):
        # locList = list([copy.copy(valList)])
    # else:
        # raise TypeError("Value list must be a list, array or number")

    locList = list(copy.copy(valList))
    if val is None:
        try:
            ind = locList.index(None)
        except ValueError:
            raise
    else:
        try:
            iNone = locList.index(None)
            locList[iNone] = np.nan
        except:
            pass
        cl = list(np.isclose(locList, val))
        try:
            ind = cl.index(True)
        except ValueError:
            ind = None
    return ind


def calcV1(E, v0):
    """
    Calculate :math:`v_1` from :math:`v_0` and :math:`E` to satisfy the definition: :math:`f_v(-70):= 1`. 
    
    Parameters
    ----------
    E : float
        The opsin reversal potential in millivolts [mV]. 
    v0 : float
        The model rectifier parameter :math:`v_0`. 
    
    Returns
    -------
    float
        The calculated value of the rectifier scaling parameter :math:`v_1`. 
    """
    return (70+E)/(np.exp((70+E)/v0)-1)
    #return (-70-E)/(1-np.exp(-(-70-E)/v0))


def lam2rgb(wav, gamma=0.8, output='norm'):
    """
    Convert a wavelength [nm] to an RGB colour tuple in the visible spectrum. 
    
    This converts a given wavelength of light to an
    approximate RGB colour value with edge attenuation.
    The wavelength must be given in nanometres in the
    range from 380 nm - 750 nm (789 THz - 400 THz).

    Adapted from: http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    
    Parameters
    ----------
    wav : int or float
        Wavelength of light in nanometres [nm]. 
    gamma : float, optional
        Gamma correction exponent (default=0.8). 
    output : {'norm', 'hex'}, optional
        Specify whether to return an RGB tuple or hexadecimal string and RGB tuple (default=(RGB)). 
    
    Returns
    -------
    tuple
        If 'norm' output an RGB tuple with ints in [0, 255] or if 'hex', output a hexadecimal string followed by an RGB tuple. 
    """
    
    # == A few notes about colour ==

    # Color   Wavelength(nm) Frequency(THz)
    # Red     620-750        484-400
    # Orange  590-620        508-484
    # Yellow  570-590        526-508
    # Green   495-570        606-526
    # Blue    450-495        668-606
    # Violet  380-450        789-668

    # f is frequency (cycles per second)
    # l (lambda) is wavelength (meters per cycle)
    # e is energy (Joules)
    # h (Plank's constant) = 6.6260695729 x 10^-34 Joule*seconds
    #                      = 6.6260695729 x 10^-34 m^2*kg/seconds
    # c = 299792458 meters per second
    # f = c/l
    # l = c/f
    # e = h*f
    # e = c*h/l

    # List of peak frequency responses for each type of
    # photoreceptor cell in the human eye:
    #     S cone: 437 nm
    #     M cone: 533 nm
    #     L cone: 564 nm
    #     rod:    550 nm in bright daylight, 498 nm when dark adapted.
    #             Rods adapt to low light conditions by becoming more sensitive.
    #             Peak frequency response shifts to 498 nm.

    wav = float(wav)
    if wav >= 380 and wav < 440:
        attenuation = 0.3 + 0.7 * (wav - 380) / (440 - 380)
        R = ((-(wav - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wav >= 440 and wav < 490:
        R = 0.0
        G = ((wav - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wav >= 490 and wav < 510:
        R = 0.0
        G = 1.0
        B = (-(wav - 510) / (510 - 490)) ** gamma
    elif wav >= 510 and wav < 580:
        R = ((wav - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wav >= 580 and wav < 645:
        R = 1.0
        G = (-(wav - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wav >= 645 and wav <= 750:
        attenuation = 0.3 + 0.7 * (750 - wav) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else: # Outside the visible spectrum
        R = 0.0
        G = 0.0
        B = 0.0

    if output == 'norm':
        return (int(R), int(G), int(B))
    elif output == 'hex':
        R *= 255
        R = int(max(0, min(round(R), 255)))
        G *= 255
        G = int(max(0, min(round(G), 255)))
        B *= 255
        B = int(max(0, min(round(B), 255)))
        #return (int(R), int(G), int(B)) # int() truncates towards 0
        return "#{0:02x}{1:02x}{2:02x}".format(R, G, B), (R, G, B)


### Model functions ###

### Physical constants
_h = 6.6260695729e-34    # Planck's constant (Js)
_c = 2.99792458e8        # Speed of light    (m*s^-1)
#NA = 6.0221413e23      # Avogadro's Number (mol^-1)
    
def irrad2flux(E, lam=470):   # E2phi
    """
    Converts irradiance [mW * mm^-2] and wavelength (default: 470) [nm]
    to flux [photons * mm^-2 * s^-1]. 
    
    Parameters
    ----------
    E : float
        Irradiance [mW * mm^-2]. 
    lam : int or float, optional
        Wavelength (default=470) [nm]. 
    
    Returns
    -------
    float
        Flux [photons * mm^-2 * s^-1]. 
    """

    Ep = 1e12 * _h * _c / lam # Energy per photon [mJ] (using lambda in [nm])
    return E / Ep             # Photon flux (phi) scaled to [photons * s^-1 * mm^-2]


def flux2irrad(phi, lam=470):
    """
    Converts flux [photons * mm^-2 * s^-1] and wavelength (default: 470) [nm]
    to irradiance [mW * mm^-2]. 

    Parameters
    ----------
    phi : float
        Flux [photons * mm^-2 * s^-1]. 
    lam : int or float, optional
        Wavelength (default=470) [nm]. 
    
    Returns
    -------
    float
        Irradiance [mW * mm^-2]. 
    """

    Ep = 1e12 * _h * _c / lam # Energy per photon [mJ] (using lambda in [nm])
    return phi * Ep           # Irradiance (E) scaled to [mW * mm^-2]


def _calcgbar(Ip, Vclamp, A=1):
    # Unused
    """Estimate (lower bound) the cell's maximum conductance from its peak current
    Ip      :=  Peak current [nA]
    Vclamp  :=  Clamp Voltage [mV]
    A       :=  Cell surface area [um^2]
    return gbar [pS/um^2]"""
    Gmax = Ip/Vclamp  # Maximum conductance for the whole cell
    gbar = Gmax/A     # Maximum conductance pS / um^2
    return gbar * (1e6) # 1e-12 S / (1e-6 m)^2 = (1e-6)*(1e-9 A / 1e-3 V)/(1e-6 m)^2


def times2cycles(times, Dt_tot):       # TODO revise to handle negative delay times c.f. PhotoCurrent
    r"""
    Convert times (absolute events) to pulse cycles (durations). 
    
    Parameters
    ----------
    times : list or array
        List of pulse times e.g. :math:`[[t_{on,0}, t_{off,0}], ..., [t_{on,N-1}, t_{off,N-1}]]`. 
    Dt_tot : float
        Total protocol time, :math:`t_{tot}`
    
    Returns
    -------
    tuple
        Array of cycles and delay duration e.g. :math:`[[\Delta t_{on,0}, \Delta t_{off,0}], ..., [\Delta t_{on,N-1}, \Delta t_{off,N-1}]], \Delta t_{del}`. 
    """
    times = np.array(times, copy=True)
    #nPulses = times.shape[0]
    assert(times.shape[1] <= 2)
    delD = times[0, 0] # This assumes that the times have not been shifted
    onDs = [row[1]-row[0] for row in times] # pulses[:,1] - pulses[:,0]   # Pulse Durations
    offDs = np.append(times[1:, 0], Dt_tot) - times[:, 1]
    cycles = np.vstack((onDs, offDs)).transpose()
    return (cycles, delD)


def cycles2times(cycles, delD):
    r"""
    Convert pulse cycles (durations) to times (absolute events)
    
    Parameters
    ----------
    cycles : list or array
        Array of cycles and delay duration e.g. :math:`[[\Delta t_{on,0}, \Delta t_{off,0}], ..., [\Delta t_{on,N-1}, \Delta t_{off,N-1}]], \Delta t_{del}`. 
    delD : float
        Delay duration, :math:`\Delta t_{del}`
    
    Returns
    -------
    tuple
        List of pulse times and total protocol time 
        e.g. :math:`[[t_{on,0}, t_{off,0}], ..., [t_{on,N-1}, t_{off,N-1}]], \Delta t_{tot}`. 
    """

    # TODO: Generalise to delDs c.f. recovery
    cycles = np.array(cycles)
    nPulses = cycles.shape[0]
    assert(cycles.shape[1] <= 2)
    times = np.zeros((nPulses, 2)) #[delD,delD+cycles[row,0] for row in pulses]
    lapsed = delD
    for p in range(nPulses):
        times[p, 0] = lapsed
        times[p, 1] = lapsed+cycles[p, 0]
        lapsed += sum(cycles[p, :])
    return (times, lapsed)


def plotLight(times, ax=None, light='shade', dark=None, lam=470, alpha=0.2):
    """
    Plot light pulse(s)
    
    Parameters
    ----------
    times : list or array
        List of pulse times e.g. :math:`[[t_{on,0}, t_{off,0}], ..., [t_{on,N-1}, t_{off,N-1}]]`. 
    ax : axes, optional
        Figure axes to plot on [default: gca()]
    light : {'shade', 'borders', 'greyscale', 'hatch', 'spectral'}, optional
        Representation type for plotting the light (default: 'shade'). 
    dark : float, optional
        Lightness of the background [0 (black), 1 (white)] (default:``None``). 
    lam : float, optional
        Wavelength [nm] (default=470). 
    alpha :float, optional
        Light region's transparency level [0, 1] (default=0.2). 
    """

    ### Change plt.axvspan to ax.axvspan etc.
    if ax is None:
        ax = plt.gca()
    else:
        plt.sca(ax)
    nPulses = times.shape[0]

    if dark is None:
        pass
    else:
        ax.set_axis_bgcolor(str(dark))
        for p in range(nPulses):
            ax.axvspan(times[p][0], times[p][1], facecolor='w')

    if light == 'shade':
        for p in range(nPulses):
            ax.axvspan(times[p][0], times[p][1], facecolor='y', alpha=alpha)
    elif light == 'borders':
        for p in range(0, nPulses):
            ax.axvline(x=times[p][0], linestyle='--', color='k')
            ax.axvline(x=times[p][1], linestyle='--', color='k')
    elif light == 'greyscale':
        # Set background to grey and illumination to white
        ax.set_axis_bgcolor('0.6') #'0.3'
        for p in range(nPulses):
            ax.axvspan(times[p][0], times[p][1], facecolor='w')
    elif light == 'hatch':
        for p in range(nPulses):
            ax.axvspan(times[p][0], times[p][1], hatch='/') #'*'
    elif light == 'spectral':
        # Plot the colour corresponding to the wavelength
        if 380 <= lam <= 750:
            rgb = lam2rgb(lam)
            for p in range(0, nPulses):
                ax.axvspan(times[p][0], times[p][1], facecolor=rgb, alpha=alpha)
        else: # Light is not in the visible spectrum - plot borders instead
            for p in range(0, nPulses):
                ax.axvline(x=times[p][0], linestyle='--', color='k')
                ax.axvline(x=times[p][1], linestyle='--', color='k')
                ax.axvspan(times[p][0], times[p][1], hatch='/') #'*'
    elif light == 'None' or light is None:
        pass
    else:
        warnings.warn('Warning: Unrecognised light representation: {}!'.format(light))
    return


def setCrossAxes(ax, zeroX=True, zeroY=False):
    """
    Remove box and set axes to run through zero
    
    Parameters
    ----------
    zeroX : bool
        Set the x-axis to run through y=0 (default=True). 
    zeroY : bool
        Set the y-axis to run through x=0 (default=False). 
    """
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    if zeroY:
        ax.spines['left'].set_position('zero') # y-axis
    if zeroX:
        ax.spines['bottom'].set_position('zero') # x-axis
    # 'center' -> ('axes', 0.5)
    # 'zero'   -> ('data', 0.0)
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
    #ax.yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter(useMathText=True))


def round_sig(x, n=3):
    """
    Round a number to ``n`` significant digits (default=3). 
    
    Parameters
    ----------
    x : float
        Number to round. 
    n : int, optional
        Number of significant digits to round to (default=3). 
    """
    
    if abs(x) == 0 or np.isinf(x) or np.isnan(x):
        return x
    else:
        return round(x, n-int(np.floor(np.log10(abs(x))))-1)


### Used in analysing kinetics (protocols) and steady-state (loadData)
def expDecay(t, a, b, c):
    r"""
    Calculate single exponential decay function.
    
    Parameters
    ----------
    t : ndarray
        Time array. 
    a, b, c : float
        Exponential function coefficients. 
    
    Returns
    -------
    ndarray
        :math:`a \cdot e^{-t/b} + c`
    """
    return a * np.exp(-t/b) + c

def biExpDecay(t, a1, tau1, a2, tau2, I_ss):
    """Calculate double exponential decay function"""
    return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2) + I_ss

def biExpSum(t, a_act, tau_act, a_deact, tau_deact, a0):
    """Calculate the sum of two opposite exponential functions"""
    return a0 + a_act*(1-np.exp(-t/tau_act)) + a_deact*np.exp(-t/tau_deact)

