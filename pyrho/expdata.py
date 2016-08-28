"""Classes for storing and processing experimental photocurrent data"""

from __future__ import print_function, division
import warnings
import logging
import copy

import numpy as np
# import scipy.io as sio # Use for Matlab files < v7.3
# import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyrho.utilities import getIndex, times2cycles, setCrossAxes, round_sig, plotLight
from pyrho.config import check_package
from pyrho import config

__all__ = ['PhotoCurrent', 'ProtocolData']


# TODO: Import/Export from/to python electrophysiology modules
# Neo: http://neuralensemble.org/neo/
# G-Node: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3942789/ (uses Neo and odML)
# Stimfit: http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00016/full
# fit_neuron: http://pythonhosted.org/fit_neuron/tutorial_easy.html


# TODO: Move to utilities.py
# import h5py
# f = h5py.File('myfile.hdf5','r')

# from StringIO import StringIO   # StringIO behaves like a file object
# c = StringIO("0 1\n2 3")
# np.loadtxt(c)
# array([[ 0.,  1.],
       # [ 2.,  3.]])

# This should be left to the user
#def loadMatFile(filename):
    ### Extend to load pkl files too
    #try:
    #import scipy.io as sio # Use for Matlab files < v7.3
    #sio.whosmat(filename)
#    data = sio.loadmat(filename)
    #except:
    #    import h5py
    #    fh = h5py.File(filename,'r')
    #    data = fh.get("var")
    #    fh.close()
#    return data


class PhotoCurrent(object):
    """
    Data storage class for an individual Photocurrent and associated properties

    Attributes
    ----------
    I : ndarray(float)
        Array of photocurrent values in nanoamps [nA].
    t : ndarray(float)
        Array of time values in milliseconds [ms] corresponding to ``I``.
    pulses : array, shape = [n_pulses, 2]
        Pairs of time points describing the beginning and end of stimulation
        e.g. [[t_on0, t_off0], [t_on1, t_off1], ..., [t_onN-1, t_offN-1]].
    phi : float
        Stimulating flux (max) [ph*mm^-2*s^-1].
    V : float or ``None``
        Voltage clamp potential [mV] or ``None`` if no clamp was used.
    stimuli : ndarray(float) or ``None``, optional
        Optional array of flux values corresponding to ``I``. Expect nStimuli x nSamples row vectors.
    states : ndarray(float) or ``None``, optional
        Optional array of model state variables if the data are synthetic (shape=nStates x nSamples).
    stateLabels : list(str) or ``None``, optional
        Optional list of LaTeX strings labelling each of the state variables.
    """
    #pulses : list(list(float))

    '''
    Dt_phi      t_phi
    cycles      times
    periods     pulses
    durations
    '''
    # TODO: Prefer "reverse notation". e.g. t_end, t_beg
    '''

    t_start        # t[0]          HIDE _t_start_
    t_end        # t[-1]         HIDE _t_end_
    t_peak_      # Time of biggest current peak     Replace with t_peaks_[0]?
    t_peaks_     # Times of current peaks in each pulse
    Dt_total        # t[-1] - t[0]
    Dt_delay        # Delay duration before the first pulse
    Dt_delays       # Total delay before each pulse
    Dt_ons        # On-phase durations
    Dt_IPIs        # Inter-pulse-intervals t_off <-> t_on
    Dt_offs       # Off-phase durations
    nPulses     # Number of pulses
    pulseCycles # Pulse durations # TODO: Rename cycles to durations/periods/Dt_phis
    nStimuli    # Number of stimuli
    nStates     # Number of model states
    synthetic   # Modelling data
    clamped     # Membrane potential was clamped
    lam         # Stimulus wavelength       RETHINK c.f. stimuli --> phi_lambda

    I_peak_       # Biggest current peak                        (.peak_ published)
    I_peaks_      # Current peaks in each pulse

    Dt_lag_        # Lag of first pulse      REMOVE
    Dt_lags_       # t_lag = Dt_act := t_peak - t_on
    I_sss_        # Steady-state currents for each pulse
    I_ss_         # Steady-state current of first pulse   REMOVE (.ss_ published)
    I_range_      # [Imin, Imax]
    I_span_       # Imax - Imin
    type        # Polarity of current     RENAME

    pulseAligned# Aligned to (a) pulse    RETHINK...
    alignPoint  # {0:=t_on, 1:=t_peak, 2:=t_off}
    p0          # Time shift
    overlap     # Periods are up to *and including* the start of the next phase
    dt          # Sampling time step                dt_
    sr          # Sampling rate [Hz] [samples/s]    sr_
    _idx_peak_    # Index of biggest current peak                     HIDE
    _idx_peaks_   # Indexes of current peaks in each pulse            HIDE
    _idx_pulses_   # Indexes for the start of each on- and off-phases  HIDE
    isFiltered  # Data hase been filtered               HIDE
    _I_orig       # Original photocurrent (unfiltered)    HIDE __I_orig_
    _I_prev       # Previous photocurrent                 HIDE __I_prev_
    _offset_     # Current offset calculated to zero dark current    HIDE
    on_         # Current at t_on[:]        REMOVE
    off_        # Current at t_off[:]       REMOVE
    '''

    # Move states (and Vm, stimuli) into .extras['states']...

    # TODO: Make this a setter which calls _find_idx_peaks and findSteadyState when changed
    overlap = True  # Periods are up to *and including* the start of the next e.g. onPhase := t[onInd] <= t <? t[offInd]

    def __init__(self, I, t, pulses, phi, V, stimuli=None, states=None, stateLabels=None, label=None):
        """
        I       := Photocurrent [nA]
        t       := Time series (or time step) [ms]
        phi     := Stimulating flux (max) [ph*mm^-2*s^-1]
        V       := Voltage clamp potential [mV]
        pulses  := Pairs of time points describing the beginning and end of stimulation
                    e.g. [[t_on1,t_off1],[t_on2,t_off2],...]

        # I = [0., 0., 0., ..., -0.945, ..., 0.]
        # t = [0, 0.1, 0.2, ..., tmax]
        # pulses = [[100, 250], [300, 500]]
        # photocurrent(I, t, V, phi, pulses)
        """

        ### Load data
        self.I = np.copy(I)                     # Array of photocurrent values np.copy(I) == np.array(I, copy=True) == np.array(I)
        self.nSamples = len(self.I)             # Number of samples

        if len(t) == len(I):
            assert(len(t) > 1)
            self.t = np.copy(t)                 # Corresponding array of time points [ms] #np.array copies by default
            tdiff = self.t[1:] - self.t[:-1]
            self.dt = tdiff.sum()/len(tdiff)    # (Average) step size
            self.sr = 1000/(self.dt)            # Sampling rate [samples/s]
        elif len(t) == 1:                       # Assume time step is passed rather than time array
            assert(t > 0)
            self.t = np.array(t*range(len(I)))
            self.dt = t                         # Step size
            self.sr = 1000/t                    # Sampling rate [samples/s]
        else:
            raise ValueError("Dimension mismatch: |t|={}; |I|={}. t must be either an array of the same length as I or a scalar defining the timestep!".format(len(t), len(I)))


        self.t_start = self.t[0]                   # Beginning trial time
        self.t_end = self.t[-1]                  # Last trial time point
        self.Dt_total = self.t_end - self.t_start       # Total trial time #max(self.t) # Handles negative delays

        if stimuli is not None:
            # TODO: Remove stimuli and make it equivalent to t i.e. float or array
            # Expect nStimuli x nSamples row vectors
            self.stimuli = np.copy(stimuli)
            ndim = self.stimuli.ndim
            shape = self.stimuli.shape
            if ndim == 1:
                self.nStimuli = 1
                assert(shape[0] == self.nSamples)
            elif ndim == 2:
                self.nStimuli = shape[0]
                assert(shape[1] == self.nSamples)
            else:
                raise ValueError('Dimension mismatch with stimuli: {}; shape: {}!'.format(ndim, shape))
        else:
            self.stimuli = None

        if states is not None:
            self.states = np.copy(states)
            self.nStates = self.states.shape[1] # len(stateLabels)
            self.stateLabels = copy.copy(stateLabels)
            assert(len(self.stateLabels) == self.nStates)
            self.synthetic = True
            assert(self.states.shape[0] == self.nSamples)
        else:
            self.synthetic = False

        '''
        if Vm is not None:
            self.Vm = np.copy(Vm)               # Array of membrane voltage values
            assert(len(self.Vm) == self.nSamples)

        if spikes is not None:
            self.spikes = np.copy(spikes)       # Array of spike times
        '''

        ### Load metadata
        self.pulses = np.array(pulses)          # nPulses x 2 array [t_on, t_off] # assumes deepcopy(pulses)
        self.nPulses = self.pulses.shape[0]
        self.pulseCycles, _ = times2cycles(self.pulses, self.t_end) #self.Dt_total)

        self.Dt_delay = self.pulses[0, 0] - self.t_start                   # Handles negative delays
        self.Dt_delays = np.array(self.pulses[:, 0] - self.t_start)        # Delay Durations
        self.Dt_ons = np.array(self.pulses[:, 1] - self.pulses[:, 0]) # Pulse Durations
        #if self.nPulses > 1:
        #    self.Dt_IPIs = np.zeros(self.nPulses - 1)
        #    for p in range(0, self.nPulses-1):
        #        self.Dt_IPIs[p] = self.pulses[p+1,0] - self.pulses[p,1]
        self.Dt_IPIs = np.array([self.pulses[p+1, 0] - self.pulses[p, 1] for p in range(self.nPulses-1)]) # end <-> start
        self.Dt_offs = np.append(self.Dt_IPIs, self.t_end-self.pulses[-1, 1])
        # self.Dt_offs = [self.Dt_total-((Dt_on+pOff)*nPulses)-Dt_delay for pOff in pulseCycles[:,1]]

        #for p in self.nPulses: # List comprehension instead?
        #   self._idx_pulses_[p,0] = np.searchsorted(self.t, pulses[p,0], side="left")  # CHECK: last index where value <= t_on
        #   self._idx_pulses_[p,1] = np.searchsorted(self.t, pulses[p,1], side="left")  # CHECK: last index where value <= t_off
        #self._idx_pulses_ = np.array([[np.searchsorted(self.t, pulses[p,time]) for time in range(2)] for p in range(self.nPulses)])
        self._idx_pulses_ = np.array([np.searchsorted(self.t, self.pulses[p, :]) for p in range(self.nPulses)], dtype=np.int)

        ### Record Experimental constants
        self.V = copy.copy(V)           # Clamp Voltage [mV]: None if no clamp was used
        self.clamped = bool(V != None)  # Flag for voltage-clamped recording
        self.phi = copy.copy(phi)       # Light intensity
        # Future inclusions
        self.lam = 470 #Lambda          # Wavelength [nm]
        # self.pH                       # pH
        # self.Temp                     # Temperature
        self.label = copy.copy(label)   # Optional trial label e.g. "saturate"

        self.isFiltered = False
        self._I_orig = None
        self._I_prev = None
        #self.filterData()              # Smooth the data with a moving average

        ### Calibrate - correct any current offset in experimental recordings
        Idel, _ = self.getDelayPhase()
        I_offset = np.mean(Idel[:int(round(0.9*len(Idel)))+1]) # Calculate the mean over the first 90% to avoid edge effects
        if abs(I_offset) > 0.01 * abs(max(self.I) - min(self.I)): # Recalibrate if the offset is more than 1% of the span
            self.I -= I_offset
            if config.verbose > 0:
                print("Photocurrent recalibrated by {} [nA]".format(I_offset))
            self._offset_ = I_offset

        #if pulses[0][0] > 0: # Check for an initial delay period
        #    onInd = self._idx_pulses_[0,0]
        #    trim = int(round(0.1*onInd)) # Discount the first and last 10% of the delay period to remove edge effects
        #    offset = np.mean(self.I[trim:onInd-trim+1]) # self.I[:onInd]
        #    if 0.01*abs(offset) > abs(max(self.I) - min(self.I)):
        #        self.I -= offset
        #        if verbose > 0:
        #            print("Photocurrent recalibrated by {} [nA]".format(offset))

            # Subtract delay from time vector to start at 0 with the first on period
            #self.t -= pulses[0][0]
            #self.t_end = max(self.t)



        ### Derive properties from the data
        self.on_ = np.array([self.I[pInd[0]] for pInd in self._idx_pulses_])     # Current at t_on[:]
        self.off_ = np.array([self.I[pInd[1]] for pInd in self._idx_pulses_])    # Current at t_off[:]

        # Add this to findPeaks
        self.I_range_ = [min(self.I), max(self.I)]
        self.I_span_ = self.I_range_[1] - self.I_range_[0]
        #if abs(self.Irange[0]) > abs(self.Irange[1]):
        #    self.Ipeak = self.Irange[0] # Min
        #    self.Ipeaks = np.asarray([min(self.getCycle(p)) for p in range(self.nPulses)]) # Peak may occur after stimulation #np.asarray([min(self.I[self._idx_pulses_[p,0]:self._idx_pulses_[p,1]]) for p in range(self.nPulses)])
        #else:
        #    self.Ipeak = self.Irange[1] # Max
        #    self.Ipeaks = np.asarray([max(self.getCycle(p)) for p in range(self.nPulses)])
            #self.Ipeaks = np.asarray([max(self.I[self._idx_pulses_[p,0]:self._idx_pulses_[p,1]]) for p in range(self.nPulses)])
        #np.asarray([max(abs(self.I[self._idx_pulses_[p,0]:self._idx_pulses_[p,1]])) for p in range(self.nPulses)])

        self._idx_peak_ = np.argmax(abs(self.I)) #np.searchsorted(self.I, self.Ipeak)
        self.t_peak_ = self.t[self._idx_peak_]
        self.I_peak_ = self.I[self._idx_peak_]

        #self._idx_peaks_ = np.array([np.argmax(abs(self.getCycle(p)[0])) for p in range(self.nPulses)]) #np.searchsorted(self.I, self.Ipeaks)
        self._idx_peaks_ = self._find_idx_peaks()
        self.t_peaks_ = self.t[self._idx_peaks_]
        self.I_peaks_ = self.I[self._idx_peaks_]

        self.Dt_lags_ = np.array([self.t_peaks_[p] - self.pulses[p, 0] for p in range(self.nPulses)]) # t_lag = t_peak - t_on
        self.Dt_lag_ = self.Dt_lags_[0]
        # For Go: t[peakInds[0]]-self.pulses[0,1]

        self.I_sss_ = np.array([self.findSteadyState(p) for p in range(self.nPulses)])
        self.I_ss_ = self.I_sss_[0]

        if self.I_peak_ < 0 and self.I_ss_ < 0:
            self.type = 'excitatory'  # Depolarising
        else:
            self.type = 'inhibitory'  # Hyperpolarising

        # Align t_0 to the start of the first pulse
        self.pulseAligned = False
        self.alignPoint = 0
        self.p0 = None
        self.alignToPulse()

        #self.findKinetics()

        if config.verbose > 1:
            print("Photocurrent data loaded! nPulses={}; Total time={}ms; Range={}nA".format(self.nPulses, self.Dt_total, str(self.I_range_)))


    def __len__(self):
        return self.Dt_total

    def __str__(self):
        """Print out summary details of the photocurrent"""
        if self.nPulses > 1:
            plural = 's'
        else:
            plural = ''
        if self.clamped: # self.V is not None
            clStr = '@ {:.3g} mV'.format(self.V)
        else:
            clStr = '(unclamped)'
        str = 'Photocurrent with {} pulse{} {} sampled at {:.3g} samples/s over {:.3g} ms {}; {:.3g} ph/s/mm^2'.format(self.nPulses, plural, self.pulses, self.sr, self.Dt_total, clStr, self.phi)
        return str

    def __call__(self, incTime=False):
        if incTime:
            return self.I, self.t
        else:
            return self.I

    #TODO: Finish this!
    def toDF(self):
        """Export to pandas dictionary"""
        if not check_package('pandas'):
            warnings.warn('Pandas not found!')
            return
        else:
            import pandas as pd
        df = pd.DataFrame({
        't'    : self.t,
        'I'    : self.I
                          })
        if self.synthetic:
            #for si, st in enumerate(self.stateLabels):
            #    df[st] = self.states[si, :]
            df[self.stateLabels] = self.states #.T?
        if self.stimuli is not None:
            df['stimuli'] = self.stimuli # TODO: Check this works with matrices
        if self.isFiltered:
            df['_I_orig'] = self._I_orig
        return df


    #TODO: Finish this - avoid circular imports or move to fitting.py!
    def fitKinetics(self, p=0, method='powell'): # trim=0.1, # defMethod
        r"""
        Fit exponentials to a photocurrent to find time constants of kinetics.

        Plot the time-constants along with the photocurrent:
            * :math:`\tau_{act} :=` The activation time-constant of :math:`[I_{on}:I_{peak}]`
            * :math:`\tau_{inact} :=` The inactivation time-constant of :math:`[I_{peak}:I_{off}]`
            * :math:`\tau_{deact} :=` The deactivation time-constant(s) of :math:`[I_{off}:]`. A single and double exponential function are fit the the off-curve.

        Parameters
        ----------
        p : int
            Specify which pulse to use (default=0) ``0 <= p < nPulses``.
        method : str
            Optimisation method (default=defMethod)
        """

        def calcOn(p, t):
            """Fit a biexponential curve to the on-phase to find lambdas"""
            v = p.valuesdict()
            return v['a0'] + v['a_act']*(1-np.exp(-t/v['tau_act'])) + v['a_deact']*np.exp(-t/v['tau_deact'])

        #def jacOn(p,t):
        #    v = p.valuesdict()
        #    return [(v['a1']/v['tau_act'])*np.exp(-t/v['tau_act']) - (v['a2']/v['tau_deact'])*np.exp(-t/v['tau_deact'])]

        def residOn(p, I, t):
            return I - calcOn(p, t)

        def calcOff(p, t):
            v = p.valuesdict()
            return v['a0'] + v['a1']*np.exp(-v['Gd1']*t) + v['a2']*np.exp(-v['Gd2']*t)

        def residOff(p, I, t):
            return I - calcOff(p, t)


        plt.figure()
        self.plot()

        # These are used only in fitKinetics
        from lmfit import Parameters, minimize
        #from pyrho.fitting import methods, defMethod
        from pyrho.fitting import reportFit

        ### On phase ###

        # t=0 :         I = a0 + a_deact = 0    ==> a0 = -a_deact
        # t=t_off :     I = a0 + a_act = Iss    ==> a_act = Iss - a0
        #                                           a_act = Iss + a_deact
        # Iss = a_act - a_deact

        Ion, ton = self.getOnPhase(p)

        pOn = Parameters()

        Iss = self.I_ss_
        #Ipeak = self.I_peak_

        if Iss < 0: # Excitatory
            pOn.add('a_act', value=Iss, min=-1e3, max=1e-9) #1
            pOn.add('a_deact', value=Iss*0.1, min=-1e3, max=1e-9, expr='a_act - {}'.format(Iss)) #0.1
        else: # Inhibitory
            pOn.add('a_act', value=Iss, min=1e-9, max=1e3) # peak_?
            pOn.add('a_deact', value=Iss*0.1, min=1e-9, max=1e3, expr='a_act - {}'.format(Iss))
        pOn.add('a0', value=0, min=-1e3, max=1e3, expr='-a_deact') # redundant
        pOn.add('tau_act', value=5, min=1e-9, max=1e3)
        pOn.add('tau_deact', value=50, min=1e-9, max=1e3)

        # Dictionary unpacking also works if preferred
        # from pyrho.utilities import biExpSum
        # def residBiExpSum(p, I, t):
             #v = p.valuesdict()
            # return I - biExpSum(t, **p.valuesdict())#v['a_act'], v['tau_act'], v['a_deact'], v['tau_deact'], v['a0'])
        # minRes = minimize(residBiExpSum, pOn, args=(Ion,ton), method=method)

        minRes = minimize(residOn, pOn, args=(Ion, ton), method=method)

        fpOn = minRes.params #pOn
        v = fpOn.valuesdict()
        print('tau_{{act}} = {:.3g}, tau_{{deact}} = {:.3g}'.format(v['tau_act'], v['tau_deact']))
        print('a_{{act}} = {:.3g}, a_{{deact}} = {:.3g}, a_0 = {:.3g}'.format(v['a_act'], v['a_deact'], v['a0']))
        if config.verbose > 1:
            #print(fit_report(minRes))
            reportFit(minRes, 'On-phase sum of exponentials', method)

        plt.plot(ton, calcOn(fpOn, ton), label=r'On-Fit $\tau_{{act}}={:.3g}, \tau_{{deact}}={:.3g}$'.format(v['tau_act'], v['tau_deact']))
        #plt.plot(ton, biExpSum(ton, **fpOn.valuesdict()), label=r'On-Fit $\tau_{{act}}={:.3g}, \tau_{{deact}}={:.3g}$'.format(v['tau_act'], v['tau_deact']))


        ### Add a check for steady-state before fitting the off-curve

        ### Off phase ###

        # t0 = t_off
        # t=0 :     I = a0 + a1 + a2 = Iss

        Iss = self.I_ss_ #fpOn['a0'].value + fpOn['a1'].value
        Ioff, toff = self.getOffPhase(p)

        # Single exponential
        pOffs = Parameters()
        if Iss < 0: # Excitatory
            pOffs.add('a0', value=0, min=Iss*.001, max=-Iss*.001, vary=True) # Add some tolerance # expr='{}-a1-a2'.format(Iss))
            pOffs.add('a1', value=0, min=-1e3, max=-1e-9, vary=True, expr='{}-a0'.format(Iss))
            #pOffs.add('a2', value=0, min=-1e3, max=0, vary=False)
        else: # Inhibitory
            pOffs.add('a0', value=0, min=-Iss*.001, max=Iss*.001, vary=True) # expr='{}-a1-a2'.format(Iss))
            pOffs.add('a1', value=0, min=1e-9, max=1e3, vary=True, expr='{}-a0'.format(Iss))
            #pOffs.add('a2', value=0, min=0, max=1e3, vary=False)

        pOffs.add('Gd1', value=10, min=1e-3, max=1e3)
        pOffs.add('Gd2', value=0, min=0, max=1e3, vary=False) #, expr='Gd1')#, min=1e-9)
        pOffs.add('a2', value=0, min=-1e-9, max=1e-9, vary=False)

        minRes = minimize(residOff, pOffs, args=(Ioff, toff-toff[0]), method=method)
        fpOffs = minRes.params #pOff
        print('tau_{{off}} = {:.3g}'.format(1/fpOffs['Gd1'].value))
        if config.verbose > 1:
            #print(fit_report(minRes))
            reportFit(minRes, 'Off-phase mono-exponential decay', method)

        plt.plot(toff, calcOff(fpOffs, toff-toff[0]), label=r'Off-Fit (Mono-Exp) $\tau_{{off}}={:.3g}$'.format(1/fpOffs['Gd1'].value))

        # Double exponential
        pOffd = Parameters()
        if Iss < 0: # Excitatory
            pOffd.add('a0', value=0, min=-1e3, max=1e3, vary=False)
            pOffd.add('a1', value=0.8*Iss, min=-1e3, max=-1e-9)
            pOffd.add('a2', value=0.2*Iss, min=-1e3, max=-1e-9, expr='{}-a0-a1'.format(Iss))
        else: # Inhibitory
            pOffd.add('a0', value=0, min=-1e3, max=1e3, vary=False)
            pOffd.add('a1', value=0.8*Iss, min=1e-9, max=1e3)
            pOffd.add('a2', value=0.2*Iss, min=-1e3, max=-1e-9, expr='{}-a0-a1'.format(Iss))
        pOffd.add('Gd1', value=0.1, min=1e-9, max=1e3)
        pOffd.add('Gd2', value=0.01, min=1e-9, max=1e3)#, vary=True) #, expr='Gd1')#, min=1e-9)

        minRes = minimize(residOff, pOffd, args=(Ioff, toff-toff[0]), method=method)
        fpOffd = minRes.params #pOff
        print('tau_{{off1}} = {:.3g}, tau_{{off2}} = {:.3g}'.format(1/fpOffd['Gd1'].value, 1/fpOffd['Gd2'].value))
        if config.verbose > 1:
            #print(fit_report(minRes))
            reportFit(minRes, 'Off-phase bi-exponential decay', method)

        plt.plot(toff, calcOff(fpOffd, toff-toff[0]), label=r'Off-Fit (Bi-Exp) $\tau_{{off1}}={:.3g}, \tau_{{off2}}={:.3g}$'.format(1/fpOffd['Gd1'].value, 1/fpOffd['Gd2'].value))
        #plt.show(block=False)
        plt.legend(loc='best')


        # TODO: Move this to fitting subpackage
        if config.verbose > 1:
            def solveGo(tlag, Gd, Go0=1000, tol=1e-9):
                Go, Go_m1 = Go0, 0
                #print(tlag, Gd, Go, Go_m1)
                while abs(Go_m1 - Go) > tol:
                    Go_m1 = Go
                    Go = ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
                    #Go_m1, Go = Go, ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
                    #print(Go, Go_m1)
                return Go

            E = 0 ### Find this from fitting fV first!!!

            GoA = solveGo(self.Dt_lag_, Gd=1/fpOn['tau_deact'].value)
            GoB = solveGo(self.Dt_lag_, Gd=max(fpOffd['Gd1'].value, fpOffd['Gd2'].value))

            corrFac = lambda Gact, Gdeact: 1 + Gdeact / Gact
            Gd = max(fpOffd['Gd1'].value, fpOffd['Gd2'].value)

            print('Lag method (tau_deact): Go = {}, cf={} --> g0 = {}'.format(GoA, corrFac(GoA, 1/fpOn['tau_deact'].value), 1e6 * self.I_peak_ * corrFac(GoA, 1/fpOn['tau_deact'].value) / (self.V - E))) #(1 + 1 / (GoA * pOn['tau_deact'].value))
            print('Lag method (max(Gd1,Gd2)): Go = {}, cf={} --> g0 = {}'.format(GoB, corrFac(GoB, Gd), 1e6 * self.I_peak_ * corrFac(GoB, Gd) / (self.V - E))) #(1 + max(pOff['Gd1'].value, pOff['Gd2'].value)/GoB)
            print('Exp method (tau_deact): Gact = {}, cf={} --> g0 = {}'.format(1/fpOn['tau_act'].value, corrFac(1/fpOn['tau_act'].value, 1/fpOn['tau_deact'].value), 1e6 * self.I_peak_ * corrFac(1/fpOn['tau_act'].value, 1/fpOn['tau_deact'].value) / (self.V - E))) #(1 + pOn['tau_act'].value / pOn['tau_deact'].value)
            print('Exp method (max(Gd1,Gd2)): Gact = {}, cf={} --> g0 = {}'.format(1/fpOn['tau_act'].value, corrFac(1/fpOn['tau_act'].value, Gd), 1e6 * self.I_peak_ * corrFac(1/fpOn['tau_act'].value, Gd) / (self.V - E))) #(1 + pOn['tau_act'].value * max(pOff['Gd1'].value, pOff['Gd2'].value))



        ### Segment the photocurrent into ON, INACT and OFF phases (Williams et al., 2013)
        # I_p := maximum (absolute) current
        # I_ss := mean(I[400ms:450ms])
        # ON := 10ms before I_p to I_p ?!
        # INACT := 10:110ms after I_p
        # OFF := 500:600ms after I_p



        # from scipy.optimize import curve_fit
        # from .parameters import p0on, p0inact, p0off

        # def monoExp(t, r, Imax):
            # return Imax * np.exp(-r*t) - Imax

        # def biExp(t, a1, tau1, a2, tau2, I_ss):
            # return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2) + I_ss

        ### Fit curve for tau_on
        #Iact, tact = pc.getActivation(p)
        #popt, pcov = curve_fit(monoExp, tact, Iact, p0=(-1, -0.2, -1)) #Needs ball-park guesses (0.3, 125, 0.5)
        #print("Activation: ", popt)

        ### Fit curve for tau_inact
        #Iinact, tinact = pc.getDeactivation(p)
        #popt, pcov = curve_fit(monoExp, tinact, Iinact, p0=(-1, 0.02, -1)) #Needs ball-park guesses (0.3, 125, 0.5)
        #print("Inactivation: ", popt)

        ### Fit curve for tau_off (bi-exponential)
        #Ioff, toff = pc.getOffPhase(p)
        #popt, pcov = curve_fit(monoExp, toff, Ioff, p0=(-0.1, 0.1, -0.1)) #Needs ball-park guesses (0.3, 125, 0.5)
        #print("Off (Mono-Exp): ", popt)

        #popt, pcov = curve_fit(biExp, toff, Ioff, p0=(-1, 7.5, -1, 35, -1)) #Needs ball-park guesses (0.3, 125, 0.5)
        #print("Off (Bi-Exp): ", popt)


        # Taken from protocols.py
        '''
        # TODO: Incorporate into Photocurrent class
        def _plotKinetics(self):
            ### Segment the photocurrent into ON, INACT and OFF phases (Williams et al., 2013)
            # I_p := maximum (absolute) current
            # I_ss := mean(I[400ms:450ms])
            # ON := 10ms before I_p to I_p ?!
            # INACT := 10:110ms after I_p
            # OFF := 500:600ms after I_p

            if not peakInds: # Prevent indexing problems when no peak was found
                peakInds = [0]
            else:
                ### Analyse kinetics for the first pulse
                ### Fit curve for tau_on
                if verbose > 1:
                    print('Analysing on-phase decay...')
                onBegInd = np.searchsorted(t,Dt_delay,side="left")
                self.fitPeaks(t[onBegInd:peakInds[0]], I_RhO[onBegInd:peakInds[0]], expDecay, p0on, '$I_{{on}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
                ### Plot tau_on vs Irrad (for curves of V)
                ### Plot tau_on vs V (for curves of Irrad)

            ### Fit curve for tau_inact
            if verbose > 1:
                print('Analysing inactivation-phase decay...')
            onEndInd = np.searchsorted(t,Dt_on+Dt_delay,side="left") # Add one since upper bound is not included in slice
            popt, _, _ = self.fitPeaks(t[peakInds[0]:onEndInd + 1], I_RhO[peakInds[0]:onEndInd + 1], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
            if verbose > 1:
                print("$\tau_{{inact}} = {}$; $I_{{ss}} = {}$".format(popt[1],popt[2]))
            Iss=popt[2]
            IssVals[run][phiInd][vInd] = Iss
            ### Plot tau_inact vs Irrad (for curves of V)
            ### Plot tau_inact vs V (for curves of Irrad)

            ### Fit curve for tau_off (bi-exponential)
            if verbose > 1:
                print('Analysing off-phase decay...')
    #                 endInd = -1 #np.searchsorted(t,Dt_off+Dt_on+Dt_delay,side="right") #Dt_total
            popt, _, _ = self.fitPeaks(t[onEndInd:], I_RhO[onEndInd:], biExpDecay, p0off, '$I_{{off}} = {:.3}e^{{-t/{:g}}} {:+.3}e^{{-t/{:g}}} {:+.3}$','')
            ### Plot tau_off vs Irrad (for curves of V)
            ### Plot tau_off vs V (for curves of Irrad)

            # Draw boundary between ON and INACT phases
            for p in peakInds:
                plt.axvline(x=t[p],linestyle=':',color='m')
                plt.axhline(y=I_RhO[peakInds[0]],linestyle=':',color='r')
                plt.axhline(y=Iss,linestyle=':',color='b')

            plt.legend(loc='best')
            return
        '''

        plt.show()

        return


    def alignToPulse(self, pulse=0, alignPoint=0):
        """Set time array so that the first pulse occurs at t=0 (with negative delay period)"""
        if not self.pulseAligned or alignPoint != self.alignPoint: #and abs(self.pulses[pulse,0]) > 1e-12:
            if alignPoint == 0:         # Start
                self.p0 = self.pulses[pulse, 0]
            elif alignPoint == 1:       # Peak
                self.p0 = self.t_peaks_[pulse]
            elif alignPoint == 2:       # End
                self.p0 = self.pulses[pulse, 1]
            else:
                raise NotImplementedError
            self.t -= self.p0           # Time array
            self.pulses -= self.p0      # Pulse times
            self.t_start = self.t[0]       # Beginning Time of Trial
            self.t_end = self.t[-1]      # End Time of Trial
            self.t_peak_ = self.t[self._idx_peak_]
            self.t_peaks_ = self.t[self._idx_peaks_]
            self.pulseAligned = True
            self.alignPoint = alignPoint


    def alignToTime(self, t=None):
        """Set time array so that it begins at t=0 [default] (with the first pulse at t>0)"""
        if t is None:
            if self.pulseAligned:
                self.p0 = self.t[0]
            else:
                self.p0 = 0
        else:
            self.p0 = t

        #if self.pulseAligned: # and abs(self.pulses[0,0]) < 1e-12:
        #    self.p0 = self.t[0]
        self.t -= self.p0           # Time array
        self.pulses -= self.p0      # Pulse times
        self.t_start = self.t[0]       # Beginning Time of Trial
        self.t_end = self.t[-1]      # End Time of Trial
        self.t_peak_ = self.t[self._idx_peak_]
        self.t_peaks_ = self.t[self._idx_peaks_]
        self.pulseAligned = False


    def _find_idx_peaks(self): #, pulse=0):
        """
        Find the indicies of the photocurrent peaks for each pulse

        Returns
        -------
        ndarry(int)
            Array of peak indexes for each pulse (shape=nPulses)
            np.array([_idx_peak_0, _idx_peak_1, ..., _idx_peak_n-1])
        """

        offsetInd = len(self.getDelayPhase()[0]) - int(self.overlap) #1
        peakInds = np.zeros((self.nPulses,), dtype=np.int)
        for p in range(self.nPulses):
            peakInds[p] = np.argmax(abs(self.getCycle(p)[0])) + offsetInd
            offsetInd += len(self.getCycle(p)[0]) - int(self.overlap) #1

        return peakInds

    ### Move findPeaks from models.py to here?
    # def findPeaks(self): ### See findPeaks in models.py
        # self.peakInds = findPeaks(self.I) ### This needs some careful tweaking for real data...
        # self.t_peaks = self.t[self.peakInds]
        # self.I_peaks = self.I[self.peakInds]

    def findSteadyState(self, pulse=0, tail=0.05, method=0): #, window=tFromOff): ### c.f. findPlateauCurrent() in models.py
        """Find the steady-state current either as the last ``tail`` proportion of the on-phase or by fitting a decay function"""
        assert(0 <= pulse < self.nPulses)

        #offInd = self._idx_pulses_[pulse][1] #np.searchsorted(t,Dt_on+Dt_delay,side="left")

        #if self.Dt_ons[pulse] < window:
        #    raise ValueError('Error: The plateau buffer must be shorter than the on phase!')
            #windowInd = int(round(p*len(I_phi))) #np.searchsorted(t,t[onEndInd]-100,side="left") # Generalise
            #I_ss = np.mean(I_phi[-windowInd:])
        #    return None

        Ion, ton = self.getOnPhase(pulse)

        # Calculate step change (or gradient with t[1:] - t[:-1])
        cutInd = max(2, int(round(tail*len(Ion)))) # Need at least 2 points to calculate the difference
        if cutInd < 5: # On-phase is too short
            warnings.warn('Duration Warning: The on-phase is too short for steady-state convergence!')
            #return None
            method = 0

        dI = Ion[-cutInd+1:] - Ion[-cutInd:-1]
        if abs(np.mean(dI)) > 0.01 * self.I_span_:
            logging.warn('Steady-state Convergence Warning: The average step size is larger than 1% of the current span!')
            #return None
            method = 0

        if method == 0: # Empirical: Calculate Steady-state as the mean of the last 5% of the On phase

            Iss = np.mean(Ion[-cutInd:])

            # Calculate Steady-state as the mean of the last 50ms of the On phase
            #tFromOffInd = np.searchsorted(t,t[offInd]-window,side="left")
            #self.Iss = np.mean(self.I[tFromOffInd:offInd+1])
            #windowInd = int(round(p*len(self.I))) #np.searchsorted(t,t[onEndInd]-100,side="left") # Generalise
            #self.Iss = np.mean(self.I[-windowInd:])

            # Calculate step change (or gradient with t[1:] - t[:-1])
            #Idiff = self.I[tFromOffInd+1:offInd+1] - self.I[tFromOffInd:offInd]
            #if abs(np.mean(Idiff)) > 0.01 * self.Ispan:
            #    warnings.warn('Steady-state Convergence Warning: The average step size is larger than 1% of the current span!')

        elif method == 1: # Theoretical: Fit curve from peak to end of on phase

            postPeak = slice(self._idx_peaks_[pulse], self._idx_pulses_[pulse, 1]+int(self.overlap)) #1 # t_peak : t_off+1
            #popt = fitPeaks(self.t[postPeak], self.I[postPeak], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
            #Iss = popt[2]

            from pyrho.parameters import p0inact
            from pyrho.utilities import expDecay
            from scipy.optimize import curve_fit

            t = self.t[postPeak]
            I = self.I[postPeak]
            shift = t[0]
            popt, pcov = curve_fit(expDecay, t-shift, I, p0=p0inact)
            Iss = popt[2]
            if config.verbose > 1:
                peakEq = '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$'.format(*[round_sig(p, 3) for p in popt])
                print(peakEq)

        return Iss


    def getdIdt(self, offset=1):
        """Calculate the first derivative of the photocurrent"""
        dI = self.I[offset:] - self.I[:-offset]
        dt = self.t[offset:] - self.t[:-offset]
        #return (dI/dt, np.cumsum(dt) - dt/2)
        #return (dI/dt, self.t[:-offset] + dt/2)
        return (dI/dt, self.t[offset:] - dt/2)


    def getd2Idt2(self, offset=1):
        """Calculate the second derivative of the photocurrent"""
        dI = self.I[offset:] - self.I[:-offset]
        dt = self.t[offset:] - self.t[:-offset]
        d2I = dI[offset:] - dI[:-offset]
        tp = self.t[offset:] - dt/2
        dt2 = tp[offset:] - tp[:-offset]
        #dt2 = dt[offset:] - dt[:-offset]
        return (d2I/dt2, tp[offset:] - dt2/2)
        #dIdt, t = self.getdIdt(offset)
        #dt = (t - self.t[offset:]) * -2
        #dI = dIdt * dt
        #d2I = dI[offset:] - dI[:-offset]
        #dt2 = dt[offset:] - dt[:-offset]
        #return (d2I/dt2, t[offset:] - dt2/2)

    #d3I = d2I[offset:] - d2I[:-offset]
    #dt3 = dt2[offset:] - dt2[:-offset]
    #plt.plot(dt2[offset:] - dt3/2, d3I/dt3)


    def getDelayPhase(self):
        """Return Idel, tdel"""
        delSlice = slice(0, self._idx_pulses_[0, 0]+int(self.overlap)) # [:_idx_pulses_[0, 0]+overlap]
        return (self.I[delSlice], self.t[delSlice])


    def getOnPhase(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the on-phase (Ion, ton) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        onSlice = slice(self._idx_pulses_[pulse, 0], self._idx_pulses_[pulse, 1]+int(self.overlap)) # [_idx_pulses_[pulse,0]:_idx_pulses_[pulse,1]+overlap]
        #return (self.I[self._idx_pulses_[pulse,0]:self._idx_pulses_[pulse,1]+1], self.t[self._idx_pulses_[pulse,0]:self._idx_pulses_[pulse,1]+1])
        return (self.I[onSlice], self.t[onSlice])


    def getOffPhase(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the off-phase (Ioff, toff) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        if 0 <= pulse < self.nPulses-1:
            offSlice = slice(self._idx_pulses_[pulse, 1], self._idx_pulses_[pulse+1, 0]+int(self.overlap))
            #return self.I[self._idx_pulses_[pulse,1]:self._idx_pulses_[pulse+1,0]+1]
        elif pulse == self.nPulses-1:   # Last Pulse
            offSlice = slice(self._idx_pulses_[pulse, 1], None) # [_idx_pulses_[pulse, 1]:]
            #return self.I[self._idx_pulses_[pulse,1]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
        return (self.I[offSlice], self.t[offSlice])


    def getCycle(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the on- and off-phase (Ip, tp) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        if 0 <= pulse < self.nPulses-1:
            cycleSlice = slice(self._idx_pulses_[pulse, 0], self._idx_pulses_[pulse+1, 0]+int(self.overlap))
            #return self.I[self._idx_pulses_[pulse,0]:self._idx_pulses_[pulse+1,0]+1] # Consider removing the +1
        elif pulse == self.nPulses-1:   # Last Pulse
            cycleSlice = slice(self._idx_pulses_[pulse, 0], None)
            #return self.I[self._idx_pulses_[pulse,0]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
        return (self.I[cycleSlice], self.t[cycleSlice])

    def getActivation(self, pulse=0):
        """Return I [nA] and t [ms] arrays from beginning of the on-phase to the peak for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        actSlice = slice(self._idx_pulses_[pulse, 0], self._idx_peaks_[pulse]+int(self.overlap))
        return (self.I[actSlice], self.t[actSlice])


    def getDeactivation(self, pulse=0): # Inactivation, Deactivation, Desensitisation???
        """Return I [nA] and t [ms] arrays from the peak to the end of the on-phase for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        deactSlice = slice(self._idx_peaks_[pulse], self._idx_pulses_[pulse, 1]+int(self.overlap))
        return (self.I[deactSlice], self.t[deactSlice])


    def plot(self, ax=None, light='shade', dark=None, addFeatures=True, colour=None, linestyle=None):

        """
        Plot the photocurrent
        Optional arguments:
        ax          :=  Specify axes on which to plot
        light       :=  Specify style of plotting for the stimulus. See plotLight()
        addFeatures :=  Plot additional features including peak, steady-state and stimulation times
        colour      :=  Specify colour for the photocurrent data
        linestyle   :=  Specify linestyle for the photocurrent data
        """

        if ax == None:
            ax = plt.gca()
        else:
            plt.sca(ax)

        fig = plt.gcf()

        # TODO: Implement optional stimulus plotting
        '''
        if addFeatures and self.stimuli is not None:
            inner_grid = mpl.gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=ax)
            axStim = plt.Subplot(fig, inner_grid[0,0])
            axStim.plot(self.t, self.stimuli)
            plt.setp(axStim.get_xticklabels(), visible=False)
            #plt.setp(axLine.get_xticklabels(), visible=False)

            axPC = plt.Subplot(fig, inner_grid[1:,0], sharex=True)
            ax = axPC
            #fig, (ax1, ax2) = plt.subplots(1,2, sharex=True, sharey=True)
            #gs = plt.GridSpec(4, 1)
            #axStim = fig.add_subplot(gs[1,1]) #, sharex=axLine)

            #axStim.plot(self.t, self.stimuli)

            #axPC = fig.add_subplot(gs[2:,1])
            #plt.sca(axPC)
        '''

        if colour is None or linestyle is None:
            plt.plot(self.t, self.I)
        else:
            plt.plot(self.t, self.I, color=colour, linestyle=linestyle)

        plotLight(self.pulses, ax=ax, light=light, dark=dark, lam=470, alpha=0.2)

        plt.xlabel(r'$\mathrm{Time\ [ms]}$', position=(config.xLabelPos, 0), ha='right')
        plt.xlim((self.t_start, self.t_end))
        plt.ylabel(r'$\mathrm{Photocurrent\ [nA]}$')

        setCrossAxes(ax)

        if addFeatures:

            #p = 0
            #plt.axvline(x=self.t_peaks_[p], linestyle=':', color='k')
            #plt.axhline(y=self.I_peaks_[p], linestyle=':', color='k')

            toffset = round(0.1 * self.t_end)

            for p in range(self.nPulses):
                # Add Pointer to peak currents
                #ax.arrow(self.t_peaks_[p], 0.8*self.I_peaks_[p], 0, 0.05*self.I_peaks_[p], head_width=0.05, head_length=0.1, fc='k', ec='k')
                # ax.annotate("", xy=(self.t_peaks_[p], self.I_peaks_[p]), xycoords='data',
                    # xytext=(self.t_peaks_[p], 0.9*self.I_peaks_[p]), textcoords='data', #textcoords='axes fraction',
                    # arrowprops=dict(arrowstyle="wedge,tail_width=1.", facecolor='red', shrinkB=10), #, shrinkB=5 , shrink=0.05
                    # horizontalalignment='center', verticalalignment='top')

                # plt.text(self.t_peaks_[p], 1.02*self.I_peaks_[p], '$I_{{peak}} = {:.3g}\mathrm{{nA}};\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(self.I_peaks_[p], self.Dt_lags_[0]), ha='left', va='top', fontsize=eqSize)

                if self.I_peaks_[p] is not None:
                    ax.annotate(r'$I_{{peak}} = {:.3g}\mathrm{{nA}};\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(self.I_peaks_[p], self.Dt_lags_[0]),
                                xy=(self.t_peaks_[p], self.I_peaks_[p]),
                                xytext=(toffset+self.t_peaks_[p], self.I_peaks_[p]),
                                arrowprops=dict(arrowstyle="wedge,tail_width=0.6", shrinkA=5, shrinkB=5, facecolor='red'),
                                horizontalalignment='left', verticalalignment='center', fontsize=config.eqSize)

                # Add pointer to steady-state currents
                if self.I_sss_[p] is not None:
                    #plt.text(1.1*self.pulses[p,1], self.I_ss_, '$I_{{ss}} = {:.3g}\mathrm{{nA}}$'.format(self.I_ss_), ha='left', va='center', fontsize=eqSize)
                    #if toffset+self.pulses[p,1] > self.t_end:
                    #xPos = 0
                    ax.annotate(r'$I_{{ss}} = {:.3g}\mathrm{{nA}}$'.format(self.I_sss_[p]), xy=(self.pulses[p, 1], self.I_sss_[p]),
                                xytext=(toffset+self.pulses[p, 1], self.I_sss_[p]),
                                arrowprops=dict(arrowstyle="wedge,tail_width=0.6", shrinkA=5, shrinkB=5),
                                horizontalalignment='left', verticalalignment='center', fontsize=config.eqSize)

                # Add labels for on and off phases
                #ymin, ymax = plt.ylim()
                #plt.ylim(round_sig(ymin,3), round_sig(ymax,3))
                #pos = 0.95 * abs(ymax-ymin)
                #arrowy = 0.085 #0.075
                #texty = 0.05

                # TODO: Fix positioning - data or axes proportions?
                texty = -round(0.1 * self.I_peak_)
                arrowy = 1.5 * texty
                #awidth=10
                ax.annotate('', xy=(self.pulses[p,0], arrowy),
                            xytext=(self.pulses[p,1], arrowy),
                            arrowprops=dict(arrowstyle='<->', color='blue', shrinkA=0, shrinkB=0))
                plt.text(self.pulses[p,0]+self.Dt_ons[p]/2, texty,
                         r'$\Delta t_{{on_{}}}={:.3g}\mathrm{{ms}}$'.format(p, self.Dt_ons[p]),
                        ha='center', va='bottom', fontsize=config.eqSize)
                if p < self.nPulses-1:
                    end = self.pulses[p+1,0]
                else:
                    end = self.t_end
                ax.annotate('', xy=(self.pulses[p,1], arrowy), xytext=(end, arrowy),
                            arrowprops=dict(arrowstyle='<->', color='green', shrinkA=0, shrinkB=0))
                plt.text(self.pulses[p,1]+self.Dt_offs[p]/2, texty,
                         r'$\Delta t_{{off_{}}}={:.3g}\mathrm{{ms}}$'.format(p, self.Dt_offs[p]),
                        ha='center', va='bottom', fontsize=config.eqSize)

        #plt.show()

        return # ax


    def plotStates(self, plotPies=True, pulse=None, name=None): #, verbose=config.verbose):

        #phi = self.phi # Use the value at t_off if the stimulus if a function of time
        t = self.t
        states = self.states
        pulses = self.pulses
        peakInds = self._idx_peaks_
        labels = self.stateLabels

        plotSum = False


        if pulse is None:
            piePulses = list(range(self.nPulses))
        else:
            if isinstance(pulse, (list, tuple)):
                piePulses = pulse
            else:
                piePulses = [pulse]

        if plotPies:
            plotInit = bool(len(piePulses) > 1)
            plotPeaks = bool(peakInds is not None)
            plotSS = True
            plotSSinf = hasattr(self, 'ssInf') # not plotInit #
            count = sum([plotInit, plotPeaks, plotSS, plotSSinf])
        else:
            count = 1
            piePulses = []
        nPulses = len(piePulses)

        figWidth, figHeight = mpl.rcParams['figure.figsize']
        fig = plt.figure(figsize=(figWidth, (1+nPulses/2)*figHeight)) # 1.5*
        gs = plt.GridSpec(2+nPulses, count)

        t_start, t_end = t[0], t[-1] # self.t_start, self.t_end

        # Plot line graph of states
        axLine = fig.add_subplot(gs[0,:])
        plt.plot(t, states)
        plt.setp(axLine.get_xticklabels(), visible=False)

        if plotSum:
            sig, = plt.plot(t, np.sum(states,axis=1), color='k', linestyle='--')
            labelsIncSum = np.append(labels, r'$\Sigma s_i$')
            plt.legend(labelsIncSum, loc=6)
        else:
            plt.legend(labels, loc=6)

        plt.ylabel(r'$\mathrm{State\ occupancy}$')
        plt.xlim((t_start, t_end))
        #plt.ylim((-0.1,1.1))
        plt.ylim((0, 1))
        if config.addTitles:
            plt.title(r'$\mathrm{State\ variables\ through\ time}$')
            #plt.title('State variables through time: $v={} \mathrm{{mV}},\ \phi={:.3g} \mathrm{{photons}} \cdot \mathrm{{s}}^{{-1}} \cdot \mathrm{{cm}}^{{-2}}$'.format(V,phiOn))
        plotLight(pulses, axLine)
        ### New plot format (plus change in ylims)
        #axLine.spines['left'].set_position('zero') # y-axis
        #axLine.spines['right'].set_color('none')
        #axLine.spines['bottom'].set_position('zero') # x-axis
        #axLine.spines['bottom'].set_color('none')
        #axLine.spines['top'].set_color('none')
        axLine.spines['left'].set_smart_bounds(True)
        axLine.spines['bottom'].set_smart_bounds(True)
        #axLine.xaxis.set_ticks_position('bottom')
        #axLine.yaxis.set_ticks_position('left')


        ### Plot stack plot of state variables
        axStack = fig.add_subplot(gs[1,:], sharex=axLine)
        plt.stackplot(t, states.T)
        plt.ylim((0, 1))
        plt.xlim((t_start, t_end))
        plotLight(pulses, axStack, 'borders')
        if config.addTitles:
            axStack.title.set_visible(False)
        plt.xlabel(r'$\mathrm{Time\ [ms]}$')
        plt.ylabel(r'$\mathrm{State\ occupancy}$')


        if plotPies:
            #TODO: Remove this again and fix bug with piechart sns colours
            #if config.fancyPlots and check_package('seaborn'):
            #    import seaborn as sns
            #    cp = sns.color_palette()
            #else:
            cp = config.colours

            for p in piePulses:
                pieInd = 0
                if plotInit:
                    axS0 = fig.add_subplot(gs[p+2, pieInd])
                    #initialStates = states[0,:] * 100 #self.s0 * 100
                    initialStates = states[self._idx_pulses_[p,0],:] * 100
                    if config.verbose > 1:
                        pct = {l:s for l,s in zip(labels, initialStates)}
                        print('Initial state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(initialStates, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title(r'$\mathrm{Initial\ state\ occupancies}$')
                    else:
                        #axS0.text(-1, 1, '$t_{0}$')
                        axS0.annotate('$t_{0}$', xycoords='axes fraction', xy=(0, 1))
                    if pieInd == 0:
                        axS0.annotate('$pulse={}$'.format(p), xycoords='axes fraction', xy=(0, 0))
                    pieInd += 1

                if plotPeaks: #peakInds is not None: ### Plot peak state proportions
                    pInd = peakInds[p] # Plot the first peak
                    axLine.axvline(x=t[pInd], linestyle=':', color='k')
                    axStack.axvline(x=t[pInd], linestyle=':', color='k')
                    axPeak = fig.add_subplot(gs[p+2, pieInd])
                    sizes = states[pInd,:] * 100
                    #sizes = [s*100 for s in sizes]
                    #explode = (0,0,0.1,0.1,0,0)
                    if config.verbose > 1:
                        pct = {l:s for l,s in zip(labels,sizes)}
                        print('Peak state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)#, explode=explode)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title(r'$\mathrm{Simulated\ peak\ state\ occupancies}$')
                    else:
                        #axPeak.text(-1, 1, '$t_{peak}$')
                        axPeak.annotate('$t_{peak}$', xycoords='axes fraction', xy=(0, 1))
                    if pieInd == 0:
                        axPeak.annotate('$pulse={}$'.format(p), xycoords='axes fraction', xy=(0, 0))
                    pieInd += 1

                if plotSS: #not plotInit: # Plot steady-state proportions
                    axSS = fig.add_subplot(gs[p+2, pieInd])
                    offInd = self._idx_pulses_[p, 1] # TODO: Revise
                    ss = states[offInd,:] * 100
                    if config.verbose > 1:
                        pct = {l:s for l,s in zip(labels, ss)}
                        print('Steady-state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(ss, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title(r'$\mathrm{Simulated\ steady-state\ occupancies}$')
                    else:
                        #axSS.text(-1, 1, '$t_{peak}$')
                        axSS.annotate('$t_{ss}$', xycoords='axes fraction', xy=(0, 1))
                    pieInd += 1

                # TODO: Generalise to use phi(t=t_off)
                if plotSSinf: # hasattr(self, 'ssInf'): #phi > 0: ### Plot steady state proportions
                    axInf = fig.add_subplot(gs[p+2, pieInd])
                    #ssInf = self.calcSteadyState(phi) * 100 # Convert array of proportions to %
                    ssInf = self.ssInf[p, :] * 100
                    if config.verbose > 1:
                        pct = {l:s for l,s in zip(labels, ssInf)}
                        print('Analytic steady-state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(ssInf, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp) #, explode=explode
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title(r'$\mathrm{Analytic\ steady-state\ occupancies}$')
                    else:
                        #axInf.text(-1, 1, r'$t_{\inf}$')#, fontsize=mpl.rcParams['legend.fontsize'])
                        axInf.annotate(r'$t_{\infty}$', xycoords='axes fraction', xy=(0, 1))

        plt.tight_layout()
        plt.show()

        if name is not None:
            from os import path
            figName = path.join(config.fDir, name+'.'+config.saveFigFormat)
            plt.savefig(figName, format=config.saveFigFormat)

        return

    '''
    def genPhiArray(self,phiOn,t_ons,t_offs,tstep):
        # t_ons and t_offs are the *start* of the on and off periods
        self.nPulses = len(t_ons)
        assert(self.nPulses == len(t_offs))
        phi = np.zeros(t_ons[0]/tstep - 1)
        for p in range(nPulses):
            phi = phi.append(phi,phiOn*np.ones((t_offs[p]-t_ons[p])/tstep - 1))
        # Consider the final off
        # Alternative...
        phi = np.zeros(len(self.t))
        for p in range(nPulses):
            phi[t_ons[p]/tstep:t_offs[p]/tstep] = phiOn
    '''

    def filterData(self, t_window=1):
        """
        Pass frequency bands to filter out
        t_window    := Time window [ms] over which to calculate the moving average
        """

        if not self.isFiltered:
            self._I_orig = np.copy(self.I) # TODO: Put in __init__
            self.isFiltered = True
            I = self._I_orig
        else:
            self._I_prev = np.copy(self.I)
            I = self._I_prev

        # Moving average
        nPoints = int(round(t_window/self.dt))
        self.I = np.convolve(I, np.ones(nPoints)/nPoints, mode='same')



class ProtocolData(object):
    """
    Container for PhotoCurrent data from parameter variations in the same protocol

    Attributes
    ----------
    protocol : str
        Label corresponding to the stimulation protocol used to collect the data.
    nRuns : int >= 1
        The number of `runs` recorded in the protocol (default=1).
    phis : list[float]
        The flux values used in the protocol.
    nPhis : int >= 1
        The number of flux values used == len(phis).
    Vs : list[float, None]
        The membrane clamp voltages used in the protocol.
    nVs : int >= 1 or None
        The number of membrane clamp voltages == len(Vs) or `None` if no voltage clamps were used.
    trials : list[list[list[PhotoCurrent]]]
        Nested lists of the PhotoCurrents for each protocol condition.
        nRuns x nPhis x nVs
    runLabels : list[str]
        List of series labels specifying the independent variable for each run.
    """

    # TODO: Replace lists with dictionaries or pandas data structures

    def __init__(self, protocol, nRuns, phis, Vs):

        self.protocol = protocol

        self.nRuns = nRuns
        self.phis = phis
        self.nPhis = len(phis)
        self.Vs = Vs
        self.nVs = len(Vs)

        self.trials = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)] # Array of PhotoCurrent objects
        self.runLabels = None

        #PD = [{'PC':PC, 'run':run, 'phi':phi, 'V':V, ...},{},... ]
        self.metaData = {'nRuns':self.nRuns, 'nPhis':self.nPhis, 'nVs':self.nVs}

        ### Extract the parameters from the pcs and file them accordingly
        ### Alternatively, put them all in a flat list/set and use getTrials to retrieve relevant pcs.

    ''' # This may be overkill since the metadata should be set upon initialisation and not changed
    @property
    def phis(self):
        return self._phis

    @phis.setter
    def phis(self, phis):
        self._phis = phis
        self.nPhis = len(phis)
    '''

    def __str__(self):
        return 'Protocol data set: [nRuns={}, nPhis={}, nVs={}]'.format(self.nRuns, self.nPhis, self.nVs)

    #def __info__:
    #    """Report data set features"""
    #    pass

    def __iter__(self):
        """Iterator to return the pulse sequence for the next trial"""
        self.run = 0
        self.phiInd = 0
        self.vInd = 0
        return self

    def __next__(self):

        if self.run >= self.nRuns:
            raise StopIteration
        pc = self.trials[self.run][self.phiInd][self.vInd]
        #self.vInd = (self.vInd + 1) % self.nVs
        self.vInd += 1
        if self.vInd >= self.nVs:
            self.phiInd += 1
            self.vInd = 0
        if self.phiInd >= self.nPhis:
            self.run += 1
            self.phiInd = 0

        return pc


    def addTrials(self, photocurrents, run=0):
        """
        Add a photocurrent to the ProtocolData object without the need to specify the phi or V index.

        Parameters
        ----------
        photocurrents : list[PhotoCurrent] or PhotoCurrent
            PhotoCurrent or list of PhotoCurrent data to add to the ProtocolData object.
        run : int, optional
            Specify the run index, as this is not apparent from the PhotoCurrent (default=0).
        """

        assert(0 <= run < self.nRuns)

        if not isinstance(photocurrents, (list, tuple)):
            if not isinstance(photocurrents, (PhotoCurrent)):
                raise TypeError("Trials must be either a PhotoCurrent or a list of PhotoCurrent objects")
            photocurrents = list([photocurrents])

        indices = []

        Vs = [pc.V for pc in photocurrents]
        phis = [pc.phi for pc in photocurrents]
        if np.allclose(Vs, np.ones_like(Vs)*Vs[0]) and np.allclose(phis, np.ones_like(phis)*phis[0]):
            V = Vs[0]
            iV = getIndex(self.Vs, V) #self._getVindex(V)
            if iV is None:
                iV = 0

            phi = phis[0]
            iPhi = getIndex(self.phis, phi)
            if iPhi is None:
                iPhi = 0

            for run, pc in enumerate(photocurrents):
                self.trials[run][iPhi][iV] = pc
            return

        for pc in photocurrents:
            V = pc.V
            iV = getIndex(self.Vs, V) #self._getVindex(V)
            if iV is None:
                iV = 0

            phi = pc.phi
            iPhi = getIndex(self.phis, phi)
            if iPhi is None:
                iPhi = 0

            self.trials[run][iPhi][iV] = pc

            if config.verbose > 1:
                print("PhotoCurrent added to run={}, iPhi={}, iV={}".format(run, iPhi, iV))
            indices.append((run, iPhi, iV))

        return indices


    def addTrial(self, photocurrent, run=0):

        assert(0 <= run < self.nRuns)

        V = photocurrent.V
        iV = getIndex(self.Vs, V) #self._getVindex(V)
        if iV is None:
            iV = 0

        phi = photocurrent.phi
        iPhi = getIndex(self.phis, phi)
        if iPhi is None:
            iPhi = 0

        self.trials[run][iPhi][iV] = photocurrent

        if config.verbose > 1:
            print("PhotoCurrent added to run={}, iPhi={}, iV={}".format(run, iPhi, iV))

        return (run, iPhi, iV)


    def getTrials(self, **kwargs):
        """Pass arrays of values for any/all/none of 'runs', 'phis' and 'Vs'
            Return a flattened list of matching photocurrents"""

        for k, v in kwargs.items():
            if not isinstance(v, (list, tuple)):
                if isinstance(v, (np.ndarray, np.generic)):
                    kwargs[k] = v.tolist()
                else:
                    kwargs[k] = list([v])

        runs = kwargs["runs"] if "runs" in kwargs else list(range(self.nRuns))
        #assert(0 <= runs < self.nRuns)
        assert(0 < len(runs) <= self.nRuns)

        if config.verbose > 0:
            print('runs = ', runs)

        if 'phis' in kwargs:
            phis = kwargs['phis']
            cl = list(np.isclose(self.phis, phis))
            iPhis = [i for i,el in enumerate(cl) if el]
        else:
            iPhis = list(range(self.nPhis))

        # phis = kwargs["phis"] if "phis" in kwargs else self.phis
        # cl = list(np.isclose(self.phis, phis))
        # iPhis = [i for i,el in enumerate(cl) if el]
        assert(0 < len(iPhis) <= self.nPhis)

        if config.verbose > 0:
            print('phis = ', iPhis)

        #cl = list(np.isclose(Vs, Vs))
        #iVs = [i for i,el in enumerate(cl) if el]
        if 'Vs' in kwargs:
            iVs = [getIndex(self.Vs, V) for V in kwargs['Vs']]
        else:
            iVs = [getIndex(self.Vs, V) for V in self.Vs]

        if config.verbose > 0:
            print('Vs = ', iVs)

        # Vs = kwargs["Vs"] if "Vs" in kwargs else self.nVs
        # assert(0 <= len(Vs) < self.nVs)

        trials = []
        for run in runs:
            for iPhi in iPhis:
                for iV in iVs:
                    if self.trials[run][iPhi][iV]: # Skip missing PhotoCurrents
                        trials.append(self.trials[run][iPhi][iV])

        return trials

    # def _getVindex(self, V): ### Generalise!!!
        # Vs = list(copy.copy(self.Vs))
        # if V is None:
            # try:
                # iV = Vs.index(None)
            # except ValueError:
                # raise
        # else:
            # try:
                # iNone = Vs.index(None)
                # Vs[iNone] = np.nan
            # except:
                # pass
            # cl = list(np.isclose(Vs, V))
            # try:
                # iV = cl.index(True)
            # except ValueError:
                # iV = None
        # return iV


    def getLineProps(self, run, phiInd, vInd):
        colours = config.colours
        styles = config.styles

        if config.verbose > 1 and (self.nRuns > len(colours) or len(self.phis) > len(colours) or len(self.Vs) > len(colours)):
            warnings.warn("Warning: only {} line colours are available!".format(len(colours)))
        if config.verbose > 0 and self.nRuns > 1 and len(self.phis) > 1 and len(self.Vs) > 1:
            warnings.warn("Warning: Too many changing variables for one plot!")
        if config.verbose > 2:
            print("Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}".format(run, self.nRuns, phiInd, len(self.phis), vInd, len(self.Vs)))
        if self.nRuns > 1:
            col = colours[run % len(colours)]
            if len(self.phis) > 1:
                style = styles[phiInd % len(styles)]
            elif len(self.Vs) > 1:
                style = styles[vInd % len(styles)]
            else:
                style = '-'
        else:
            if len(self.Vs) > 1:
                col = colours[vInd % len(colours)]
                if len(self.phis) > 1:
                    style = styles[phiInd % len(styles)]
                else:
                    style = '-'
            else:
                if len(self.phis) > 1:
                    col = colours[phiInd % len(colours)]
                    style = '-'
                else:
                    col = 'b'   ### colours[0]
                    style = '-' ### styles[0]

        return col, style


    def plot(self, ax=None, light='shade', addFeatures=True):
        """
        Convenience method to plot the data set.

        Parameters
        ----------
        ax : axis, optional
            Existing axis on which to plot the PhotoCurrent data (default=None).
            If `None`, create a new axis.
        light : str {'shade', ...}, optional
            Specify the light style to plot.
        addFeatures : bool, optional
            Flag to pass to PhotoCurrent.plot for adding extra features e.g. peak and steady-state

        See Also
        --------
        utilities.plotLight
        """

        if ax == None:
            ax = plt.gca()
        else:
            plt.sca(ax)
        self.legLabels = []

        #Dt_ons = []
        t_starts, t_ends = [], []
        #pulseSet = [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
        self.nPulses = self.trials[0][0][0].nPulses # Assume all trials have the same number
        pulseSet = np.zeros((self.nPulses, 2, self.nRuns))
        for run in range(self.nRuns):
            for phiInd, phi in enumerate(self.phis):
                for vInd, V in enumerate(self.Vs):
                    pc = self.trials[run][phiInd][vInd]
                    #pc.alignToPulse()
                    t_starts.append(pc.t_start)
                    t_ends.append(pc.t_end)
                    #pulseSet[run][phiInd][vInd] = pc.pulses
                    pulseSet[:, :, run] = pc.pulses
                    #Dt_ons.append(pc.Dt_ons)
                    col, style = self.getLineProps(run, phiInd, vInd)
                    pc.plot(ax=ax, light=None, addFeatures=False, colour=col, linestyle=style)
                    label = ''
                    if self.nRuns > 1 and self.runLabels is not None: #hasattr(self, 'runLabels'):
                        label += self.runLabels[run]
                    if self.nPhis > 1:
                        label += r'$\phi = {:.3g}\ \mathrm{{[ph. \cdot mm^{{-2}} \cdot s^{{-1}}]}}$ '.format(phi)
                    if self.nVs > 1:
                        label += r'$\mathrm{{V}} = {:+}\ \mathrm{{[mV]}}$ '.format(V)
                    self.legLabels.append(label)

                    # if run==0 and phiInd==0 and vInd==0:
                        # #self.t_start, self.t_end = pc.t_start, pc.t_end
                        # pulses = pc.pulses
                        # plotLight(pulses, ax=ax, light=light, lam=470, alpha=0.2)
                    # else:
                        # for p in range(pc.nPulses):
                            # if np.allclose(pulses[p], pc.pulses[p]):
                                # pass
                            # elif np.allclose(pulses[p,0], pc.pulses[p,0]) or np.allclose(pulses[p,1], pc.pulses[p,1]):
                                # pass
                            # else:
                                # plotLight(np.asarray([pc.pulses[p]]), ax=ax, light=light, lam=470, alpha=0.2)
                    # if pc.t_start < self.t_start:
                        # self.t_start = pc.t_start
                    # if pc.t_end > self.t_start:
                        # self.t_end = pc.t_end
                    #plotLight(pc.pulses, ax=ax, light=light, lam=470, alpha=0.2)

        self.t_start, self.t_end = min(t_starts), max(t_ends)

        #for run in range(self.nRuns):
        #run = 0 # Arbitrary choice TODO: Reconsider!
        # Add stimuli
        for p in range(self.nPulses):
            sameStart, sameEnd = False, False
            if np.allclose(pulseSet[p, 0, :], np.tile(pulseSet[p, 0, 0], (1, 1, self.nRuns))): #pth t_on are the same
                sameStart = True
            if np.allclose(pulseSet[p, 1, :], np.tile(pulseSet[p, 1, 0], (1, 1, self.nRuns))): #pth t_off are the same
                sameEnd = True

            if sameStart and sameEnd: #np.allclose(pulseSet[p,:,run], np.tile(pulseSet[p,:,0], (1,1,self.nRuns))): #pth pulses are the same
                plotLight(np.asarray([pulseSet[p, :, 0]]), ax=ax, light=light, lam=470, alpha=0.2)

            elif not sameStart and not sameEnd: # No overlap
                for run in range(self.nRuns):
                    plotLight(np.asarray([pulseSet[p, :, run]]), ax=ax, light=light, lam=470, alpha=0.2)

            else: #not (sameStart and sameEnd): # One or the other - xor
                pass # This applies to shortPulse only at present - do not shade!
                #for run in range(self.nRuns):
                    # Plot bars
                #for run in range(self.nRuns):
                #    plotLight(np.asarray([pulseSet[p, :, run]]), ax=ax, light=light, lam=470, alpha=0.2)


        ### Move to protocols...
        if len(self.Vs) == 1:
            ncol = 1
        else:
            ncol = len(self.phis)
        if label != '':
            ax.legend(self.legLabels, loc='best', borderaxespad=0, ncol=ncol, fancybox=True)

        # Freeze y-limits
        #ax.set_ylim(ax.get_ylim())
        ax.set_ybound(ax.get_ylim())

        #tickLabels = [item.get_text() for item in ax.get_yticklabels(which='both')]
        ax.set_xlabel(r'$\mathrm{Time\ [ms]}$', position=(config.xLabelPos, 0), ha='right')
        plt.xlim((self.t_start, self.t_end))
        ax.set_ylabel(r'$\mathrm{Photocurrent\ [nA]}$')

        setCrossAxes(ax)

        #ax.set_yticklabels(tickLabels)
        #if np.all([Dt_on == Dt_ons[0] for Dt_on in Dt_ons]):
        #if np.allclose(Dt_ons, Dt_ons[0] * np.ones(len(Dt_ons))):
        #plotLight(pc.pulses, ax=ax, light=light, lam=470, alpha=0.2)
        #else:
            # Plot bars for on periods
        #    pass
        #plotLight(self.getProtPulses(), ax=ax, light=light, lam=470, alpha=0.2)
        #ax.tight_layout()

        #plt.show()
        return # ax


    def getIpmax(self, vInd=None):
        """
        Find the maximum peak current for the whole data set.

        This is useful when the 'delta' protocol is absent.

        Parameters
        ----------
        vInd : int, optional
            Optionally restrict the search to a particular membrane potential (default=None).

        Returns
        -------
        Ipmax_ : float
            Maximum absolute (most extreme) peak current value found in data set.
        tuple
            Indexes of the most extreme value found (rmax, pmax, vmax)
        """

        self.Ipmax_ = 0
        if vInd is None:
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    for vInd in range(self.nVs):
                        if abs(self.trials[run][phiInd][vInd].I_peak_) > abs(self.Ipmax_):
                            self.Ipmax = self.trials[run][phiInd][vInd].I_peak_
                            rmax = run
                            pmax = phiInd
                            vmax = vInd
        else:
            assert(vInd < self.nVs)
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    if abs(self.trials[run][phiInd][vInd].I_peak_) > abs(self.Ipmax_):
                        self.Ipmax_ = self.trials[run][phiInd][vInd].I_peak_
                        rmax = run
                        pmax = phiInd
                        vmax = vInd
        return self.Ipmax_, (rmax, pmax, vmax)

    # reduce(lambda a,b: a if (a > b) else b, list)

    def getProtPeaks(self):
        """
        Return the set of maximum absolute (most extreme) peak currents across a whole set of photocurrents

        Returns
        -------
        Ipeaks : list[list[list[float]]]
            Nested lists of peak values: nRuns x nPhis x nVs.
        tpeaks : list[list[list[float]]]
            Nested lists of peak value times: nRuns x nPhis x nVs.
        """
        if self.nRuns > 1:
            phiInd = 0
            vInd = 0
            self.IrunPeaks = [self.trials[run][phiInd][vInd].I_peak_ for run in range(self.nRuns)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].t_peak_ for run in range(self.nRuns)]
            Ipeaks = self.IrunPeaks
            tpeaks = self.trunPeaks
        if self.nPhis > 1:
            run = 0
            vInd = 0
            self.IphiPeaks = [self.trials[run][phiInd][vInd].I_peak_ for phiInd in range(self.nPhis)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].t_peak_ for phiInd in range(self.nPhis)]
            Ipeaks = self.IphiPeaks
            tpeaks = self.trunPeaks
        if self.nVs > 1:
            run = 0
            phiInd = 0
            self.IVPeaks = [self.trials[run][phiInd][vInd].I_peak_ for vInd in range(self.nVs)]
            self.tVPeaks = [self.trials[run][phiInd][vInd].t_peak_ for vInd in range(self.nVs)]
            Ipeaks = self.IVPeaks
            tpeaks = self.tVPeaks
        return Ipeaks, tpeaks


    def getSteadyStates(self, run=0, phiInd=None):
        """
        Return Iss, Vss.

        Find the steady-state currents and the corresponding voltage clamp potentials.

        Parameters
        ----------
        run : int, optional
            Specify the run index to collect the data from (default=0).
        phiInd : int, optional
            Specify the phi index to collect the data from or collect from all flux values (default=None).

        Returns
        -------
        Iss : list[list[float]] or list[float]
            Nested lists of peak values: len == nPhis x nVs or nVs.
        Vss : list[list[float]] or list[float]
            Nested lists of steady-state membrane potentials: len == nPhis x nVs or nVs
        """

        assert(self.nVs > 1)
        if phiInd is None: # Return 2D array
            self.Isss_ = np.zeros((self.nPhis, self.nVs))
            self.Vss_ = np.zeros((self.nPhis, self.nVs))
            for phiInd, phi in enumerate(self.phis):
                for vInd, V in enumerate(self.Vs):
                    self.Isss_[phiInd, vInd] = self.trials[run][phiInd][vInd].I_ss_ # Variations along runs are not useful here
                    self.Vss_[phiInd, vInd] = self.trials[run][phiInd][vInd].V
        else:
            self.Isss_ = np.zeros(self.nVs)
            self.Vss_ = np.zeros(self.nVs)
            for vInd, V in enumerate(self.Vs):
                self.Isss_[vInd] = self.trials[run][phiInd][vInd].I_ss_
                self.Vss_[vInd] = self.trials[run][phiInd][vInd].V
        return self.Isss_, self.Vss_



'''
from collections import defaultdict

class DataSet():
    """Container for photocurrent data used to produce arrays for parameter extraction"""

    def __init__(self, fluxSet, saturate=None, recovery=None, rectifier=None, shortPulses=None):
        self.data = defaultdict(list)

        #self.protocol = protocol

        # # if protocol == shortPulse:
            # # self.data = data
        # # elif protocol == recovery:
            # # self.Is = Is
            # # self.ts = ts
            # # pulseCycles=np.column_stack((Dt_on*np.ones(len(Dt_IPIs)),[IPI-Dt_on for IPI in Dt_IPIs])) # [:,0] = on phase duration; [:,1] = off phase duration

    def addData(self, photoData, protocol):
        self.data[protocol].append(photoData)
        self.phis.append(photoData.phi)
'''

# ProtocolData.append(TrialData(I,t,V,phi,pulses))

# print "Sort the dataSet in place by V ..."
# import operator
# dataSet.sort(key=operator.attrgetter('V'))
