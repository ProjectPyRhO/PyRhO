"""Classes for storing and processing experimental photocurrent data"""

import warnings
import copy

import numpy as np
import scipy.io as sio # Use for Matlab files < v7.3
#import h5py
from lmfit import Parameters, minimize, fit_report #*
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyrho.fitting import methods, defMethod
from pyrho.utilities import * # For times2cycles and cycles2times, expDecay, findPlateauCurrent
#from pyrho.parameters import tFromOff
from pyrho.config import * #verbose, colours, styles
from pyrho.config import xLabelPos
from pyrho import config

# See also python electrophysiology modules
# Neo: http://neuralensemble.org/neo/
# G-Node: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3942789/ (uses Neo and odML)
# Stimfit: http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00016/full
# fit_neuron: http://pythonhosted.org/fit_neuron/tutorial_easy.html

__all__ = ['PhotoCurrent', 'ProtocolData']

# TODO: Move to utilities.py
# import h5py
# f = h5py.File('myfile.hdf5','r')

# from StringIO import StringIO   # StringIO behaves like a file object
# c = StringIO("0 1\n2 3")
# np.loadtxt(c)
# array([[ 0.,  1.],
       # [ 2.,  3.]])

# I = [0.,0.,0.,...-0.945,...,0.]
# t = [0,0.1,0.2,...,tmax]
# pulses = [[100,250],[300,500]]
# photocurrent(I,t,V,phi,pulses)


def loadMatFile(filename):
    ### Extend to load pkl files too
    #try:
    #import scipy.io as sio # Use for Matlab files < v7.3
    #sio.whosmat(filename)
    data = sio.loadmat(filename)
    #except: 
    #    import h5py
    #    fh = h5py.File(filename,'r')
    #    data = fh.get("var")
    #    fh.close()
    return data

'''
class RhodopsinStates():
    """Data storage class for models states and their associated properties"""
    
    def __init__(self, states, t, varLabels):
        ### Load data
        self.states = np.copy(states)                     # Array of state values
        self.nStates = len(varLabels)
        self.nPoints = states.size/self.nStates
        
        if len(t) == len(self.nPoints):
            assert(len(t) > 1)
            self.t = np.copy(t)                 # Corresponding array of time points [ms] #np.array copies by default
            tdiff = t[1:] - t[:-1]
            self.dt = tdiff.sum()/len(tdiff)    # (Average) step size
            self.sr = 1000/(self.dt)            # Sampling rate [samples/s]
        elif len(t) == 1:                       # Assume time step is passed rather than time array
            assert(t > 0)
            self.t = np.array(t*range(self.nPoints))
            self.dt = t                         # Step size
            self.sr = 1000/t                    # Sampling rate [samples/s]
        else:
            raise ValueError("Dimension mismatch: t must be either an array of the same length as I or a scalar defining the timestep!")
        
        #...
'''

class PhotoCurrent():
    """Data storage class for an individual Photocurrent and its associated properties"""
    # TODO: Make this a setter which calls findPeakInds and findSteadyState when changed
    overlap = True  # Periods are up to *and including* the start of the next e.g. onPhase := t[onInd] <= t <? t[offInd]
    
    def __init__(self, I, t, pulses, phi, V, states=None, stateLabels=None, label=None):
        """ I       := Photocurrent [nA]
            t       := Time series (or time step) [ms]
            phi     := Stimulating flux (max) [ph*mm^-2*s^-1]
            V       := Voltage clamp potential [mV]
            pulses  := Pairs of time points describing the beginning and end of stimulation 
                        e.g. [[t_on1,t_off1],[t_on2,t_off2],...]"""
                        
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
        
        self.begT = self.t[0]                   # Beginning trial time
        self.endT = self.t[-1]                  # Last trial time point
        self.totT = self.endT - self.begT       # Total trial time #max(self.t) # Handles negative delays
        
        
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
        self.pulseCycles, _ = times2cycles(self.pulses, self.endT) #self.totT)
                
        self.delD = self.pulses[0,0] - self.begT                    # Handles negative delays
        self.delDs = np.array(self.pulses[:,0] - self.begT)         # Delay Durations
        self.onDs = np.array(self.pulses[:,1] - self.pulses[:,0])   # Pulse Durations
        #if self.nPulses > 1:
        #    self.IPIs = np.zeros(self.nPulses - 1)
        #    for p in range(0, self.nPulses-1):
        #        self.IPIs[p] = self.pulses[p+1,0] - self.pulses[p,1]
        self.IPIs = np.array([self.pulses[p+1,0] - self.pulses[p,1] for p in range(self.nPulses-1)]) # end <-> start
        self.offDs = np.append(self.IPIs, self.endT-self.pulses[-1,1])
        # self.offDs = [self.totT-((onD+pOff)*nPulses)-delD for pOff in pulseCycles[:,1]]    
        
        #for p in self.nPulses: # List comprehension instead?
        #   self.pulseInds[p,0] = np.searchsorted(self.t, pulses[p,0], side="left")  # CHECK: last index where value <= t_on
        #   self.pulseInds[p,1] = np.searchsorted(self.t, pulses[p,1], side="left")  # CHECK: last index where value <= t_off
        #self.pulseInds = np.array([[np.searchsorted(self.t, pulses[p,time]) for time in range(2)] for p in range(self.nPulses)])
        self.pulseInds = np.array([np.searchsorted(self.t, self.pulses[p,:]) for p in range(self.nPulses)], dtype=np.int)
        
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
        #self.filterData()              # Smooth the data with a moving average
        
        ### Calibrate - correct any current offset in experimental recordings
        Idel, _ = self.getDelayPhase()
        I_offset = np.mean(Idel[:int(round(0.9*len(Idel)))+1]) # Calculate the mean over the first 90% to avoid edge effects
        if abs(I_offset) > 0.01 * abs(max(self.I) - min(self.I)): # Recalibrate if the offset is more than 1% of the span
            self.I -= I_offset
            if verbose > 0:
                print("Photocurrent recalibrated by {} [nA]".format(I_offset))
            self.offset_ = I_offset
        
        #if pulses[0][0] > 0: # Check for an initial delay period
        #    onInd = self.pulseInds[0,0]
        #    trim = int(round(0.1*onInd)) # Discount the first and last 10% of the delay period to remove edge effects
        #    offset = np.mean(self.I[trim:onInd-trim+1]) # self.I[:onInd]
        #    if 0.01*abs(offset) > abs(max(self.I) - min(self.I)):
        #        self.I -= offset
        #        if verbose > 0:
        #            print("Photocurrent recalibrated by {} [nA]".format(offset))
                
            # Subtract delay from time vector to start at 0 with the first on period
            #self.t -= pulses[0][0]
            #self.endT = max(self.t)
        

        
        ### Derive properties from the data
        self.on_ = np.array([self.I[pInd[0]] for pInd in self.pulseInds])     # Current at t_on[:]
        self.off_ = np.array([self.I[pInd[1]] for pInd in self.pulseInds])    # Current at t_off[:]
        
        # Add this to findPeaks
        self.range_ = [min(self.I), max(self.I)]
        self.span_ = self.range_[1] - self.range_[0]
        #if abs(self.Irange[0]) > abs(self.Irange[1]):
        #    self.Ipeak = self.Irange[0] # Min
        #    self.Ipeaks = np.asarray([min(self.getCycle(p)) for p in range(self.nPulses)]) # Peak may occur after stimulation #np.asarray([min(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]]) for p in range(self.nPulses)])
        #else:
        #    self.Ipeak = self.Irange[1] # Max
        #    self.Ipeaks = np.asarray([max(self.getCycle(p)) for p in range(self.nPulses)])
            #self.Ipeaks = np.asarray([max(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]]) for p in range(self.nPulses)])
        #np.asarray([max(abs(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]])) for p in range(self.nPulses)])
        
        self.peakInd_ = np.argmax(abs(self.I)) #np.searchsorted(self.I, self.Ipeak)
        self.tpeak_ = self.t[self.peakInd_]
        self.peak_ = self.I[self.peakInd_]
        
        #self.peakInds_ = np.array([np.argmax(abs(self.getCycle(p)[0])) for p in range(self.nPulses)]) #np.searchsorted(self.I, self.Ipeaks)
        self.peakInds_ = self.findPeakInds()
        self.tpeaks_ = self.t[self.peakInds_]
        self.peaks_ = self.I[self.peakInds_]
        
        self.lags_ = np.array([self.tpeaks_[p] - self.pulses[p,0] for p in range(self.nPulses)]) # t_lag = t_peak - t_on
        self.lag_ = self.lags_[0]
        # For Go: t[peakInds[0]]-self.pulses[0,1]
                
        self.sss_ = np.array([self.findSteadyState(p) for p in range(self.nPulses)])
        self.ss_ = self.sss_[0]
        
        if self.peak_ < 0 and self.ss_ < 0:
            self.type = 'excitatory'
        else:
            self.type = 'inhibitory'
        
        # Align t_0 to the start of the first pulse
        self.pulseAligned = False
        self.alignPoint = 0
        self.alignToPulse()
        
        #self.findKinetics()
        
        if verbose > 1:
            print("Photocurrent data loaded! nPulses={}; Total time={}ms; Range={}nA".format(self.nPulses, self.totT, str(self.range_)))
    
    
    def __len__(self):
        return self.totT
    
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
        str = 'Photocurrent with {} pulse{} {} sampled at {:.3g} samples/s over {:.3g} ms {}; {:.3g} ph/s/mm^2'.format(self.nPulses, plural, self.pulses, self.sr, self.totT, clStr, self.phi)
        return str
        
    def __call__(self, incTime=False):
        if incTime:
            return self.I, self.t
        else:
            return self.I
    
    #TODO: Finish this!
    def toDF():
        """Export to pandas dictionary"""
        df = DataFrame({
                        't'    : self.t,
                        'I'    : self.I
                        })
        if self.synthetic:
            #for si, st in enumerate(self.stateLabels):
            #    df[st] = self.states[si, :]
            df[self.stateLabels] = self.states #.T?
        if self.isFiltered:
            df['Iorig'] = self.Iorig
        return df
    
    
    #TODO: Finish this!
    def fitKinetics(self, p=0, trim=0.1, method=defMethod):
        """
        Fit exponentials to a photocurrent to find time constants of kinetics
        p       : specify which pulse to use (default=0)
        
        method  : optimisation method (default=defMethod)
        """
        
        def calcOn(p,t):
            """Fit a biexponential curve to the on-phase to find lambdas"""
            v = p.valuesdict()
            return v['a0'] + v['a_act']*(1-np.exp(-t/v['tau_act'])) + v['a_deact']*np.exp(-t/v['tau_deact'])
        
        #def jacOn(p,t):
        #    v = p.valuesdict()
        #    return [(v['a1']/v['tau_act'])*np.exp(-t/v['tau_act']) - (v['a2']/v['tau_deact'])*np.exp(-t/v['tau_deact'])]
        
        def residOn(p,I,t):
            return I - calcOn(p,t)
        
        def calcOff(p,t):
            v = p.valuesdict()
            return v['a0'] + v['a1']*np.exp(-v['Gd1']*t) + v['a2']*np.exp(-v['Gd2']*t)
        
        def residOff(p,I,t):
            return I - calcOff(p,t)
        
        
        plt.figure()
        self.plot()
        
        from pyrho.fitting import reportFit
        
        ### On phase ###
        
        # t=0 :         I = a0 + a_deact = 0    ==> a0 = -a_deact
        # t=t_off :     I = a0 + a_act = Iss    ==> a_act = Iss - a0
        #                                           a_act = Iss + a_deact
        # Iss = a_act - a_deact        
        
        Ion, ton = self.getOnPhase(p)
        
        pOn = Parameters()
        
        Iss = self.ss_
        #Ipeak = self.peak_
        
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
#        from pyrho.utilities import biExpSum
#        def residBiExpSum(p, I, t):
#            #v = p.valuesdict()
#            return I - biExpSum(t, **p.valuesdict())#v['a_act'], v['tau_act'], v['a_deact'], v['tau_deact'], v['a0'])
#        minRes = minimize(residBiExpSum, pOn, args=(Ion,ton), method=method)

        minRes = minimize(residOn, pOn, args=(Ion,ton), method=method)

        fpOn = minRes.params #pOn
        v = fpOn.valuesdict()
        print('tau_{{act}} = {:.3g}, tau_{{deact}} = {:.3g}'.format(v['tau_act'], v['tau_deact']))
        print('a_{{act}} = {:.3g}, a_{{deact}} = {:.3g}, a_0 = {:.3g}'.format(v['a_act'], v['a_deact'], v['a0']))
        if config.verbose > 1:
            #print(fit_report(minRes))
            reportFit(minRes, 'On-phase sum of exponentials', method)
        
        plt.plot(ton, calcOn(fpOn,ton), label=r'On-Fit $\tau_{{act}}={:.3g}, \tau_{{deact}}={:.3g}$'.format(v['tau_act'], v['tau_deact']))
        #plt.plot(ton, biExpSum(ton, **fpOn.valuesdict()), label=r'On-Fit $\tau_{{act}}={:.3g}, \tau_{{deact}}={:.3g}$'.format(v['tau_act'], v['tau_deact']))
        
        
        ### Add a check for steady-state before fitting the off-curve
        
        ### Off phase ###
        
        # t0 = t_off
        # t=0 :     I = a0 + a1 + a2 = Iss
        
        Iss = self.ss_ #fpOn['a0'].value + fpOn['a1'].value        
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
        
        minRes = minimize(residOff, pOffs, args=(Ioff,toff-toff[0]), method=method)
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
        
        minRes = minimize(residOff, pOffd, args=(Ioff,toff-toff[0]), method=method)
        fpOffd = minRes.params #pOff
        print('tau_{{off1}} = {:.3g}, tau_{{off2}} = {:.3g}'.format(1/fpOffd['Gd1'].value, 1/fpOffd['Gd2'].value))
        if config.verbose > 1:
            #print(fit_report(minRes))
            reportFit(minRes, 'Off-phase bi-exponential decay', method)
        
        plt.plot(toff, calcOff(fpOffd,toff-toff[0]), label=r'Off-Fit (Bi-Exp) $\tau_{{off1}}={:.3g}, \tau_{{off2}}={:.3g}$'.format(1/fpOffd['Gd1'].value, 1/fpOffd['Gd2'].value))        
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
            
            GoA = solveGo(self.lag_, Gd=1/fpOn['tau_deact'].value)
            GoB = solveGo(self.lag_, Gd=max(fpOffd['Gd1'].value, fpOffd['Gd2'].value))
            
            corrFac = lambda Gact, Gdeact: 1 + Gdeact / Gact
            Gd = max(fpOffd['Gd1'].value, fpOffd['Gd2'].value)
            
            print('Lag method (tau_deact): Go = {}, cf={} --> g0 = {}'.format(GoA, corrFac(GoA, 1/fpOn['tau_deact'].value), 1e6 * self.peak_ * corrFac(GoA, 1/fpOn['tau_deact'].value) / (self.V - E))) #(1 + 1 / (GoA * pOn['tau_deact'].value))
            print('Lag method (max(Gd1,Gd2)): Go = {}, cf={} --> g0 = {}'.format(GoB, corrFac(GoB, Gd), 1e6 * self.peak_ * corrFac(GoB, Gd) / (self.V - E) )) #(1 + max(pOff['Gd1'].value, pOff['Gd2'].value)/GoB)
            print('Exp method (tau_deact): Gact = {}, cf={} --> g0 = {}'.format(1/fpOn['tau_act'].value, corrFac(1/fpOn['tau_act'].value, 1/fpOn['tau_deact'].value), 1e6 * self.peak_ * corrFac(1/fpOn['tau_act'].value, 1/fpOn['tau_deact'].value) / (self.V - E) )) #(1 + pOn['tau_act'].value / pOn['tau_deact'].value)
            print('Exp method (max(Gd1,Gd2)): Gact = {}, cf={} --> g0 = {}'.format(1/fpOn['tau_act'].value, corrFac(1/fpOn['tau_act'].value, Gd), 1e6 * self.peak_ * corrFac(1/fpOn['tau_act'].value, Gd) / (self.V - E) )) #(1 + pOn['tau_act'].value * max(pOff['Gd1'].value, pOff['Gd2'].value))
        


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
                onBegInd = np.searchsorted(t,delD,side="left")
                self.fitPeaks(t[onBegInd:peakInds[0]], I_RhO[onBegInd:peakInds[0]], expDecay, p0on, '$I_{{on}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
                ### Plot tau_on vs Irrad (for curves of V)
                ### Plot tau_on vs V (for curves of Irrad)
            
            ### Fit curve for tau_inact
            if verbose > 1:
                print('Analysing inactivation-phase decay...')
            onEndInd = np.searchsorted(t,onD+delD,side="left") # Add one since upper bound is not included in slice
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
    #                 endInd = -1 #np.searchsorted(t,offD+onD+delD,side="right") #totT
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



        return
    
    
    def alignToPulse(self, pulse=0, alignPoint=0):
        """Set time array so that the first pulse occurs at t=0 (with negative delay period)"""
        if not self.pulseAligned or alignPoint != self.alignPoint: #and abs(self.pulses[pulse,0]) > 1e-12:
            if alignPoint == 0:         # Start
                self.p0 = self.pulses[pulse,0]
            elif alignPoint == 1:       # Peak
                self.p0 = self.tpeaks_[pulse]
            elif alignPoint == 2:       # End
                self.p0 = self.pulses[pulse,1]
            else:
                raise NotImplementedError
            self.t -= self.p0           # Time array
            self.pulses -= self.p0      # Pulse times
            self.begT = self.t[0]       # Beginning Time of Trial
            self.endT = self.t[-1]      # End Time of Trial
            self.tpeak_ = self.t[self.peakInd_]
            self.tpeaks_ = self.t[self.peakInds_] 
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
        self.begT = self.t[0]       # Beginning Time of Trial
        self.endT = self.t[-1]      # End Time of Trial
        self.tpeak_ = self.t[self.peakInd_]
        self.tpeaks_ = self.t[self.peakInds_]
        self.pulseAligned = False
    
    
    def findPeakInds(self): #, pulse=0):
        """Find the indicies of the photocurrent peaks for each pulse
            OUT:    np.array([peakInd_0, peakInd_1, ..., peakInd_n])"""
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
        """Find the steady-state current either as the last tail % of the on-phase or by fitting a decay function"""
        assert(0 <= pulse < self.nPulses)
        
        #offInd = self.pulseInds[pulse][1] #np.searchsorted(t,onD+delD,side="left")
        
        #if self.onDs[pulse] < window:
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
        if abs(np.mean(dI)) > 0.01 * self.span_:
            warnings.warn('Steady-state Convergence Warning: The average step size is larger than 1% of the current span!')
            return None
            #method = 0
        
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
            
            postPeak = slice(self.peakInds_[pulse], self.pulseInds[pulse, 1]+int(self.overlap)) #1 # t_peak : t_off+1
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
                peakEq = '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$'.format(*[round_sig(p,3) for p in popt])
                print(peakEq)
        
        return Iss
    
    
    def getdIdt(self, offset=1):
        dI = self.I[offset:] - self.I[:-offset]
        dt = self.t[offset:] - self.t[:-offset]
        #return (dI/dt, np.cumsum(dt) - dt/2)
        #return (dI/dt, self.t[:-offset] + dt/2)
        return (dI/dt, self.t[offset:] - dt/2)
    
    
    def getd2Idt2(self, offset=1):
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
        delSlice = slice(0, self.pulseInds[0,0]+int(self.overlap))
        #return self.I[:self.pulseInds[0,0]+1]
        return (self.I[delSlice], self.t[delSlice])
    
    
    def getOnPhase(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the on-phase (Ion, ton) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        onSlice = slice(self.pulseInds[pulse,0], self.pulseInds[pulse,1]+int(self.overlap))
        #return (self.I[self.pulseInds[pulse,0]:self.pulseInds[pulse,1]+1], self.t[self.pulseInds[pulse,0]:self.pulseInds[pulse,1]+1])
        return (self.I[onSlice], self.t[onSlice])
    
    
    def getOffPhase(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the off-phase (Ioff, toff) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        if 0 <= pulse < self.nPulses-1:
            offSlice = slice(self.pulseInds[pulse,1], self.pulseInds[pulse+1,0]+int(self.overlap))
            #return self.I[self.pulseInds[pulse,1]:self.pulseInds[pulse+1,0]+1]
        elif pulse == self.nPulses-1:   # Last Pulse
            offSlice = slice(self.pulseInds[pulse,1], None)
            #return self.I[self.pulseInds[pulse,1]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
        return (self.I[offSlice], self.t[offSlice])
    
    
    def getCycle(self, pulse=0):
        """Return I [nA] and t [ms] arrays from the on- and off-phase (Ip, tp) for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        if 0 <= pulse < self.nPulses-1:
            cycleSlice = slice(self.pulseInds[pulse,0], self.pulseInds[pulse+1,0]+int(self.overlap))
            #return self.I[self.pulseInds[pulse,0]:self.pulseInds[pulse+1,0]+1] # Consider removing the +1
        elif pulse == self.nPulses-1:   # Last Pulse
            cycleSlice = slice(self.pulseInds[pulse,0], None)
            #return self.I[self.pulseInds[pulse,0]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
        return (self.I[cycleSlice], self.t[cycleSlice])
        
    def getActivation(self, pulse=0):
        """Return I [nA] and t [ms] arrays from beginning of the on-phase to the peak for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        actSlice = slice(self.pulseInds[pulse,0], self.peakInds_[pulse]+int(self.overlap))
        return (self.I[actSlice], self.t[actSlice])
    
    
    def getDeactivation(self, pulse=0): # Inactivation, Deactivation, Desensitisation???
        """Return I [nA] and t [ms] arrays from the peak to the end of the on-phase for a given pulse"""
        assert(0 <= pulse < self.nPulses)
        deactSlice = slice(self.peakInds_[pulse], self.pulseInds[pulse,1]+int(self.overlap))
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
        
        if colour is None or linestyle is None:
            plt.plot(self.t, self.I)
        else:
            plt.plot(self.t, self.I, color=colour, linestyle=linestyle)
        
        plotLight(self.pulses, ax=ax, light=light, dark=dark, lam=470, alpha=0.2)
        
        plt.xlabel(r'$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
        plt.xlim((self.begT, self.endT))
        plt.ylabel(r'$\mathrm{Photocurrent\ [nA]}$')
        
        setCrossAxes(ax)
        
        if addFeatures:
            #p = 0
            #plt.axvline(x=self.tpeaks_[p], linestyle=':', color='k')
            #plt.axhline(y=self.peaks_[p], linestyle=':', color='k')
            
            toffset = round(0.1 * self.endT)
            
            for p in range(self.nPulses):
                # Add Pointer to peak currents
                #ax.arrow(self.tpeaks_[p], 0.8*self.peaks_[p], 0, 0.05*self.peaks_[p], head_width=0.05, head_length=0.1, fc='k', ec='k')
                # ax.annotate("", xy=(self.tpeaks_[p], self.peaks_[p]), xycoords='data',
                    # xytext=(self.tpeaks_[p], 0.9*self.peaks_[p]), textcoords='data', #textcoords='axes fraction',
                    # arrowprops=dict(arrowstyle="wedge,tail_width=1.", facecolor='red', shrinkB=10), #, shrinkB=5 , shrink=0.05
                    # horizontalalignment='center', verticalalignment='top')
                
                # plt.text(self.tpeaks_[p], 1.02*self.peaks_[p], '$I_{{peak}} = {:.3g}\mathrm{{nA}};\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(self.peaks_[p], self.lags_[0]), ha='left', va='top', fontsize=eqSize)
                
                if self.peaks_[p] is not None:
                    ax.annotate('$I_{{peak}} = {:.3g}\mathrm{{nA}};\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(self.peaks_[p], self.lags_[0]), xy=(self.tpeaks_[p], self.peaks_[p]), xytext=(toffset+self.tpeaks_[p], self.peaks_[p]), arrowprops=dict(arrowstyle="wedge,tail_width=0.6", shrinkA=5, shrinkB=5, facecolor='red'), horizontalalignment='left', verticalalignment='center', fontsize=config.eqSize)
                
                # Add pointer to steady-state currents
                if self.sss_[p] is not None:
                    #plt.text(1.1*self.pulses[p,1], self.ss_, '$I_{{ss}} = {:.3g}\mathrm{{nA}}$'.format(self.ss_), ha='left', va='center', fontsize=eqSize)
                    #if toffset+self.pulses[p,1] > self.endT:
                    #xPos = 0
                    ax.annotate('$I_{{ss}} = {:.3g}\mathrm{{nA}}$'.format(self.sss_[p]), xy=(self.pulses[p,1], self.sss_[p]), xytext=(toffset+self.pulses[p,1], self.sss_[p]), arrowprops=dict(arrowstyle="wedge,tail_width=0.6", shrinkA=5, shrinkB=5), horizontalalignment='left', verticalalignment='center', fontsize=config.eqSize)
                
                # Add labels for on and off phases
                #ymin, ymax = plt.ylim()
                #plt.ylim(round_sig(ymin,3), round_sig(ymax,3))
                #pos = 0.95 * abs(ymax-ymin)
                #arrowy = 0.085 #0.075
                #texty = 0.05
                
                texty = -round(0.1 * self.peak_)
                arrowy = 1.5 * texty
                #awidth=10
                ax.annotate('', xy=(self.pulses[p,0], arrowy), xytext=(self.pulses[p,1], arrowy), arrowprops=dict(arrowstyle='<->',color='blue',shrinkA=0,shrinkB=0))
                plt.text(self.pulses[p,0]+self.onDs[p]/2, texty, '$\Delta t_{{on_{}}}={:.3g}\mathrm{{ms}}$'.format(p, self.onDs[p]), ha='center', va='bottom', fontsize=config.eqSize)
                if p < self.nPulses-1:
                    end = self.pulses[p+1,0]
                else:
                    end = self.endT
                ax.annotate('', xy=(self.pulses[p,1], arrowy), xytext=(end, arrowy), arrowprops=dict(arrowstyle='<->',color='green',shrinkA=0,shrinkB=0))
                plt.text(self.pulses[p,1]+self.offDs[p]/2, texty, '$\Delta t_{{off_{}}}={:.3g}\mathrm{{ms}}$'.format(p, self.offDs[p]), ha='center', va='bottom', fontsize=config.eqSize)
        
        
        return # ax
    
    
    def plotStates(self, plotPies=True, pulse=None, name=None, verbose=config.verbose):
        
        phi = self.phi # Use the value at t_off if the stimulus if a function of time
        t = self.t
        states = self.states
        pulses = self.pulses
        peakInds = self.peakInds_
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
        
        begT, endT = t[0], t[-1] # self.begT, self.endT
        
        # Plot line graph of states
        axLine = fig.add_subplot(gs[0,:])
        plt.plot(t, states)
        plt.setp(axLine.get_xticklabels(), visible=False)
        
        if plotSum:
            sig, = plt.plot(t, np.sum(states,axis=1), color='k', linestyle='--')
            labelsIncSum = np.append(labels, '$\Sigma s_i$')
            plt.legend(labelsIncSum, loc=6)
        else:
            plt.legend(labels, loc=6)
        
        plt.ylabel('$\mathrm{State\ occupancy}$')
        plt.xlim((begT, endT))
        #plt.ylim((-0.1,1.1))
        plt.ylim((0, 1))
        if config.addTitles:
            plt.title('$\mathrm{State\ variables\ through\ time}$') 
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
        plt.xlim((begT, endT))
        plotLight(pulses, axStack, 'borders')
        if config.addTitles:
            axStack.title.set_visible(False)
        plt.xlabel('$\mathrm{Time\ [ms]}$')
        plt.ylabel('$\mathrm{State\ occupancy}$')
        
        
        if plotPies:
            if config.fancyPlots:
                import seaborn as sns
                cp = sns.color_palette()
            else:
                cp = config.colours
            
            for p in piePulses:
                pieInd = 0
                if plotInit:
                    axS0 = fig.add_subplot(gs[p+2, pieInd])
                    #initialStates = states[0,:] * 100 #self.s0 * 100
                    initialStates = states[self.pulseInds[p,0],:] * 100
                    if verbose > 1:
                        pct = {l:s for l,s in zip(labels, initialStates)}
                        print('Initial state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(initialStates, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title('$\mathrm{Initial\ state\ occupancies}$')
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
                    if verbose > 1:
                        pct = {l:s for l,s in zip(labels,sizes)}
                        print('Peak state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)#, explode=explode)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title('$\mathrm{Simulated\ peak\ state\ occupancies}$')
                    else:
                        #axPeak.text(-1, 1, '$t_{peak}$')
                        axPeak.annotate('$t_{peak}$', xycoords='axes fraction', xy=(0, 1))
                    if pieInd == 0:
                        axPeak.annotate('$pulse={}$'.format(p), xycoords='axes fraction', xy=(0, 0))
                    pieInd += 1
                
                if plotSS: #not plotInit: # Plot steady-state proportions
                    axSS = fig.add_subplot(gs[p+2, pieInd])
                    offInd = self.pulseInds[p, 1] # TODO: Revise
                    ss = states[offInd,:] * 100
                    if verbose > 1:
                        pct = {l:s for l,s in zip(labels, ss)}
                        print('Steady-state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(ss, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp)
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title('$\mathrm{Simulated\ steady-state\ occupancies}$')
                    else:
                        #axSS.text(-1, 1, '$t_{peak}$')
                        axSS.annotate('$t_{ss}$', xycoords='axes fraction', xy=(0, 1))
                    pieInd += 1
                
                # TODO: Generalise to use phi(t=t_off)
                if plotSSinf: # hasattr(self, 'ssInf'): #phi > 0: ### Plot steady state proportions
                    axInf = fig.add_subplot(gs[p+2, pieInd])
                    #ssInf = self.calcSteadyState(phi) * 100 # Convert array of proportions to %
                    ssInf = self.ssInf[p, :] * 100
                    if verbose > 1:
                        pct = {l:s for l,s in zip(labels, ssInf)}
                        print('Analytic steady-state occupancies (%):', sorted(pct.items(), key=lambda x: labels.index(x[0])))
                    patches, texts, autotexts = plt.pie(ssInf, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False, colors=cp) #, explode=explode
                    for lab in range(len(labels)):
                        texts[lab].set_fontsize(mpl.rcParams['ytick.labelsize'])
                        autotexts[lab].set_fontsize(mpl.rcParams['axes.labelsize'])
                    plt.axis('equal')
                    if config.addTitles:
                        plt.title('$\mathrm{Analytic\ steady-state\ occupancies}$')
                    else:
                        #axInf.text(-1, 1, r'$t_{\inf}$')#, fontsize=mpl.rcParams['legend.fontsize'])
                        axInf.annotate(r'$t_{\infty}$', xycoords='axes fraction', xy=(0, 1))

        plt.tight_layout()
        
        if name is not None:
            from os import path
            figName = path.join(fDir, name+'.'+config.saveFigFormat)
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
            self.Iorig = np.copy(self.I) # TODO: Put in __init__
            self.isFiltered = True
            I = self.Iorig
        else:
            self.Iprev = np.copy(self.I)
            I = self.Iprev
        
        # Moving average
        nPoints = int(round(t_window/self.dt))
        self.I = np.convolve(I, np.ones(nPoints)/nPoints, mode='same')
        


class ProtocolData():
    """Container for PhotoCurrent data from parameter variations in the same protocol"""
    
    # TODO: Replace lists with dictionaries or pandas data structures    
    
    def __init__(self, protocol, nRuns, phis, Vs):
        
        self.protocol = protocol
        
        self.nRuns = nRuns
        self.phis = phis
        self.nPhis = len(phis)
        self.Vs = Vs
        self.nVs = len(Vs)
        
        self.trials = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)] # Array of PhotoCurrent objects
        
        #PD = [{'PC':PC, 'run':run, 'phi':phi, 'V':V, ...},{},... ]
        self.metaData = {'nRuns':self.nRuns, 'nPhis':self.nPhis, 'nVs':self.nVs}
        
        ### Extract the parameters from the pcs and file them accordingly
        ### Alternatively, put them all in a flat list/set and use getTrials to retrieve relevant pcs. 

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
            
            if verbose > 1:
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
        
        if verbose > 1:
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
        
        if verbose > 0:
            print("runs = ", runs)
        
        if "phis" in kwargs:
            phis = kwargs["phis"]
            cl = list(np.isclose(self.phis, phis))
            iPhis = [i for i,el in enumerate(cl) if el]
        else:
            iPhis = list(range(self.nPhis))
        
        # phis = kwargs["phis"] if "phis" in kwargs else self.phis
        # cl = list(np.isclose(self.phis, phis))
        # iPhis = [i for i,el in enumerate(cl) if el]
        assert(0 < len(iPhis) <= self.nPhis)
        
        if verbose > 0:
            print("phis = ", iPhis)
        
        #cl = list(np.isclose(Vs, Vs))
        #iVs = [i for i,el in enumerate(cl) if el]
        if "Vs" in kwargs:
            iVs = [getIndex(self.Vs, V) for V in kwargs["Vs"]]
        else:
            iVs = [getIndex(self.Vs, V) for V in self.Vs]
        
        if verbose > 0:
            print("Vs = ", iVs)
        
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
        #global colours
        #global styles
        if verbose > 1 and (self.nRuns>len(colours) or len(self.phis)>len(colours) or len(self.Vs)>len(colours)):
            warnings.warn("Warning: only {} line colours are available!".format(len(colours)))
        if verbose > 0 and self.nRuns>1 and len(self.phis)>1 and len(self.Vs)>1:
            warnings.warn("Warning: Too many changing variables for one plot!")
        if verbose > 2:
            print("Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}".format(run,self.nRuns,phiInd,len(self.phis),vInd,len(self.Vs)))
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
        if ax == None:
            ax = plt.gca()
        else:
            plt.sca(ax)
        self.legLabels = []
        
        #onDs = []
        begTs, endTs = [], []
        #pulseSet = [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
        self.nPulses = self.trials[0][0][0].nPulses # Assume all trials have the same number
        pulseSet = np.zeros((self.nPulses, 2, self.nRuns))
        for run in range(self.nRuns):
            for phiInd, phi in enumerate(self.phis):
                for vInd, V in enumerate(self.Vs):
                    pc = self.trials[run][phiInd][vInd]
                    #pc.alignToPulse()
                    begTs.append(pc.begT)
                    endTs.append(pc.endT)
                    #pulseSet[run][phiInd][vInd] = pc.pulses
                    pulseSet[:,:,run] = pc.pulses
                    #onDs.append(pc.onDs)
                    col, style = self.getLineProps(run, phiInd, vInd)
                    pc.plot(ax=ax, light=None, addFeatures=False, colour=col, linestyle=style)
                    label = ""
                    if self.nRuns > 1 and hasattr(self, 'runLabels'):
                        label += self.runLabels[run]
                    if self.nPhis > 1:
                        label += "$\phi = {:.3g}\ \mathrm{{[ph. \cdot mm^{{-2}} \cdot s^{{-1}}]}}$ ".format(phi)
                    if self.nVs > 1:
                        label += "$\mathrm{{V}} = {:+}\ \mathrm{{[mV]}}$ ".format(V)
                    self.legLabels.append(label)

                    # if run==0 and phiInd==0 and vInd==0:
                        # #self.begT, self.endT = pc.begT, pc.endT
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
                    # if pc.begT < self.begT:
                        # self.begT = pc.begT
                    # if pc.endT > self.begT:
                        # self.endT = pc.endT
                    #plotLight(pc.pulses, ax=ax, light=light, lam=470, alpha=0.2)
        
        self.begT, self.endT = min(begTs), max(endTs)
        
        # Add stimuli
        for p in range(self.nPulses):
            sameStart, sameEnd = False, False
            if np.allclose(pulseSet[p,0,run], np.tile(pulseSet[p,0,0], (1,1,self.nRuns))): #pth t_on are the same
                sameStart = True
            if np.allclose(pulseSet[p,1,run], np.tile(pulseSet[p,1,0], (1,1,self.nRuns))): #pth t_off are the same
                sameEnd = True
            
            if sameStart and sameEnd: #np.allclose(pulseSet[p,:,run], np.tile(pulseSet[p,:,0], (1,1,self.nRuns))): #pth pulses are the same
                plotLight(np.asarray([pulseSet[p,:,0]]), ax=ax, light=light, lam=470, alpha=0.2)
            
            elif not sameStart and not sameEnd: # No overlap
                for run in range(self.nRuns):
                    plotLight(np.asarray([pulseSet[p,:,run]]), ax=ax, light=light, lam=470, alpha=0.2)
            
            else: #not (sameStart and sameEnd): # One or the other - xor
                pass
                #for run in range(self.nRuns):
                    # Plot bars

        ### Move to protocols...
        if len(self.Vs) == 1:
            ncol = 1
        else:
            ncol = len(self.phis)
        if label != "":
            lgd = ax.legend(self.legLabels, loc='best', borderaxespad=0, ncol=ncol, fancybox=True)
        
        # Freeze y-limits
        #ax.set_ylim(ax.get_ylim())
        ax.set_ybound(ax.get_ylim())
        
        #tickLabels = [item.get_text() for item in ax.get_yticklabels(which='both')]
        ax.set_xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
        plt.xlim((self.begT, self.endT))
        ax.set_ylabel('$\mathrm{Photocurrent\ [nA]}$')
        
        setCrossAxes(ax)
        
        #ax.set_yticklabels(tickLabels)
        #if np.all([onD == onDs[0] for onD in onDs]):
        #if np.allclose(onDs, onDs[0] * np.ones(len(onDs))):
        #plotLight(pc.pulses, ax=ax, light=light, lam=470, alpha=0.2)
        #else:
            # Plot bars for on periods
        #    pass
        #plotLight(self.getProtPulses(), ax=ax, light=light, lam=470, alpha=0.2)
        #ax.tight_layout()
        
        return # ax
        
        
    def getIpmax(self, vInd=None): 
        """Find the maximum peak current for the whole data set. This is useful when the 'delta' protocol is absent"""
        self.Ipmax_ = 0
        if vInd is None:
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    for vInd in range(self.nVs):
                        if abs(self.trials[run][phiInd][vInd].peak_) > abs(self.Ipmax_):
                            self.Ipmax = self.trials[run][phiInd][vInd].peak_
                            rmax = run
                            pmax = phiInd
                            vmax = vInd
        else:
            assert(vInd < self.nVs)
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    if abs(self.trials[run][phiInd][vInd].peak_) > abs(self.Ipmax_):
                        self.Ipmax_ = self.trials[run][phiInd][vInd].peak_
                        rmax = run
                        pmax = phiInd
                        vmax = vInd
        return self.Ipmax_, (rmax, pmax, vmax)
    
    # reduce(lambda a,b: a if (a > b) else b, list)
    
    def getProtPeaks(self):
        """Return the set of maximum (absolute) peak currents across a whole set of photocurrents"""
        if self.nRuns > 1:
            phiInd = 0
            vInd = 0
            self.IrunPeaks = [self.trials[run][phiInd][vInd].peak_ for run in range(self.nRuns)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].tpeak_ for run in range(self.nRuns)]
            Ipeaks = self.IrunPeaks
            tpeaks = self.trunPeaks
        if self.nPhis > 1:
            run = 0
            vInd = 0
            self.IphiPeaks = [self.trials[run][phiInd][vInd].peak_ for phiInd in range(self.nPhis)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].tpeak_ for phiInd in range(self.nPhis)]
            Ipeaks = self.IphiPeaks
            tpeaks = self.trunPeaks
        if self.nVs > 1:
            run = 0
            phiInd = 0
            self.IVPeaks = [self.trials[run][phiInd][vInd].peak_ for vInd in range(self.nVs)]
            self.tVPeaks = [self.trials[run][phiInd][vInd].tpeak_ for vInd in range(self.nVs)]
            Ipeaks = self.IVPeaks
            tpeaks = self.tVPeaks
        return Ipeaks, tpeaks
        
    
    def getSteadyStates(self, run=0, phiInd=None):
        """Return Iss, Vss"""
        assert(self.nVs > 1)
        if phiInd is None: # Return 2D array
            self.Isss_ = np.zeros((self.nPhis,self.nVs))
            self.Vss_ = np.zeros((self.nPhis,self.nVs))
            for phiInd, phi in enumerate(self.phis): 
                for vInd, V in enumerate(self.Vs): 
                    self.Isss_[phiInd,vInd] = self.trials[run][phiInd][vInd].ss_ # Variations along runs are not useful here
                    self.Vss_[phiInd,vInd] = self.trials[run][phiInd][vInd].V
        else:
            self.Isss_ = np.zeros(self.nVs)
            self.Vss_ = np.zeros(self.nVs)
            for vInd, V in enumerate(self.Vs): 
                self.Isss_[vInd] = self.trials[run][phiInd][vInd].ss_
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
            # # pulseCycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI-onD for IPI in IPIs])) # [:,0] = on phase duration; [:,1] = off phase duration

    def addData(self, photoData, protocol):
        self.data[protocol].append(photoData)
        self.phis.append(photoData.phi)
'''

# ProtocolData.append(TrialData(I,t,V,phi,pulses))

# print "Sort the dataSet in place by V ..."
# import operator
# dataSet.sort(key=operator.attrgetter('V'))
