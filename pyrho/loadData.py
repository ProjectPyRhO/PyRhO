import numpy as np
import scipy.io as sio # Use for Matlab files < v7.3
from lmfit import Parameters, minimize, fit_report #*
from pyrho.utilities import * # For times2cycles and cycles2times, expDecay, findPlateauCurrent
from pyrho.parameters import tFromOff
from pyrho.config import * #verbose, colours, styles
from pyrho import config
import warnings
import copy


# See also python electrophysiology modules
# Neo: http://neuralensemble.org/neo/
# G-Node: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3942789/ (uses Neo and odML)
# Stimfit: http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00016/full
# fit_neuron: http://pythonhosted.org/fit_neuron/tutorial_easy.html


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
    #sio.whosmat(matfile)
    data = sio.loadmat(matfile)
    #except: 
    #    import h5py
    #    fh = h5py.File(matfile,'r')
    #    data = fh.get("var")
    #    fh.close()
    return data


class RhodopsinStates():
    """Data storage class for models states and their associated properties"""
    
    def __init__(self, states, t, varLabels):
        ### Load data
        self.states = np.copy(states)                     # Array of state values
        
        if len(t) == len(I):
            assert(len(t) > 1)
            self.t = np.copy(t)                 # Corresponding array of time points [ms] #np.array copies by default
            tdiff = t[1:] - t[:-1]
            self.dt = tdiff.sum()/len(tdiff)    # (Average) step size
            self.sr = 1000/(self.dt)            # Sampling rate [samples/s]
        elif len(t) == 1:                       # Assume time step is passed rather than time array
            assert(t > 0)
            self.t = np.array(t*range(len(I)))
            self.dt = t                         # Step size
            self.sr = 1000/t                    # Sampling rate [samples/s]
        else:
            raise ValueError("Dimension mismatch: t must be either an array of the same length as I or a scalar defining the timestep!")
        
        #...


class PhotoCurrent():
    """Data storage class for an individual Photocurrent and its associated properties"""
    overlap = True  # Periods are up to *and including* the start of the next e.g. onPhase := t[onInd] <= t <? t[offInd]
    
    def __init__(self, I, t, pulses, phi, V, label=None):
        """ I       := Photocurrent [nA]
            t       := Time series (or time step) [ms]
            phi     := Stimulating flux (max) [ph*mm^-2*s^-1]
            V       := Voltage clamp potential [mV]
            pulses  := Pairs of time points describing the beginning and end of stimulation 
                        e.g. [[t_on1,t_off1],[t_on2,t_off2],...]"""
                        
        ### Load data
        self.I = np.copy(I)                     # Array of photocurrent values np.copy(I) == np.array(I, copy=True) == np.array(I)
        
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
        
        self.nSamples = len(self.I)             # Number of samples
        
        self.begT = self.t[0]                   # Beginning trial time
        self.endT = self.t[-1]                  # Last trial time point
        self.totT = self.endT - self.begT       # Total trial time #max(self.t) # Handles negative delays
        
        '''
        if states is not None:
            self.states = np.copy(states)
            self.nStates = self.states.shape[0]
            self.stateLabels = copy.copy(stateLabels)
            assert(len(self.stateLabels) == self.nStates)
            self.synthetic = True
            assert(self.states.shape[1] == self.nSamples)
        
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
        
        # Align t_0 to the start of the first pulse
        self.pulseAligned = False
        self.alignPoint = 0
        self.alignToPulse()
        
        #self.findKinetics()
        
        if verbose > 1:
            print("Photocurrent data loaded! nPulses={}; Total time={}ms; Range={}nA".format(self.nPulses, self.totT, str(self.Irange)))
    
    
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
    def findKinetics(self, p=0, trim=0.1):
        ### Segment the photocurrent into ON, INACT and OFF phases (Williams et al., 2013)
        # I_p := maximum (absolute) current
        # I_ss := mean(I[400ms:450ms])
        # ON := 10ms before I_p to I_p ?!
        # INACT := 10:110ms after I_p
        # OFF := 500:600ms after I_p
        
        method = 'powell'
        
        def calcOn(p,t):
            """Fit a biexponential curve to the on-phase to find lambdas"""
            v = p.valuesdict()
            return -(v['a0'] + v['a1']*(1-np.exp(-t/v['tau_act'])) + v['a2']*np.exp(-t/v['tau_deact']))
        
        #def jacOn(p,t):
        #    v = p.valuesdict()
        #    return [(v['a1']/v['tau_act'])*np.exp(-t/v['tau_act']) - (v['a2']/v['tau_deact'])*np.exp(-t/v['tau_deact'])]
        
        def residOn(p,I,t):
            return I - calcOn(p,t)
        
        def calcOff(p,t):
            v = p.valuesdict()
            return -(v['a0'] + v['a1']*np.exp(-v['Gd1']*t) + v['a2']*np.exp(-v['Gd2']*t))
        
        def residOff(p,I,t):
            return I - calcOff(p,t)
        
        def monoExp(t, A, B, C):
            #C = -A
            return A * np.exp(-B*t) + C

        def biExp(t, a1, tau1, a2, tau2, I_ss):
            return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2) + I_ss
        
        plt.figure()
        self.plot()
        
        ### On phase
        Ion, ton = self.getOnPhase(p)
        
        pOn = Parameters()
        pOn.add('a0', value=0, expr='-a2')
        pOn.add('a1', value=1, min=1e-9)
        pOn.add('a2', value=0.1, min=1e-9)
        pOn.add('tau_act', value=5, min=1e-9)
        pOn.add('tau_deact', value=50, min=1e-9)
        minRes = minimize(residOn, pOn, args=(Ion,ton), method=method)
        fp = pOn
        print('tau_{{act}} = {:.3g}, tau_{{deact}} = {:.3g}'.format(fp['tau_act'].value, fp['tau_deact'].value))
        if verbose > 1:
            print(fit_report(minRes))
        
        plt.plot(ton, calcOn(fp,ton), label='On-Fit $\\tau_{{act}}={:.3g}, \\tau_{{deact}}={:.3g}$'.format(fp['tau_act'].value, fp['tau_deact'].value))

        
        ### Add a check for steady-state before fitting the off-curve
        
        ### Off phase
        Iss = fp['a0'].value + fp['a1'].value
        pOff = Parameters() # copy.deepcopy(pOn)
        Ioff, toff = self.getOffPhase(p)
        
        # Single exponential
        pOff.add('a0', value=0, expr='{}-a1-a2'.format(Iss))
        pOff.add('a1', value=0, vary=True)
        pOff.add('a2', value=-0, vary=False)
        pOff.add('Gd1', value=0.1)#, min=1e-9)
        pOff.add('Gd2', value=0, vary=False) #, expr='Gd1')#, min=1e-9)
        minRes = minimize(residOff, pOff, args=(Ioff,toff-toff[0]), method=method)
        fp = pOff
        print('tau_{{off}} = {:.3g}'.format(1/fp['Gd1'].value))
        if verbose > 1:
            print(fit_report(minRes))
        
        plt.plot(toff, calcOff(fp, toff-toff[0]), label='Off-Fit (Mono-Exp) $\\tau_{{off}}={:.3g}$'.format(1/fp['Gd1'].value))
        
        # Double exponential
        pOff = Parameters()
        pOff.add('a0', value=0, vary=False)
        pOff.add('a1', value=0.1)
        pOff.add('a2', value=-0.1, expr='{}-a0-a1'.format(Iss))
        pOff.add('Gd1', value=0.1)#, min=1e-9)
        pOff.add('Gd2', value=0.01)#, vary=True) #, expr='Gd1')#, min=1e-9)
        minRes = minimize(residOff, pOff, args=(Ioff,toff-toff[0]), method=method)
        fp = pOff
        print('tau_{{off1}} = {:.3g}, tau_{{off2}} = {:.3g}'.format(1/fp['Gd1'].value, 1/fp['Gd2'].value))
        if verbose > 1:
            print(fit_report(minRes))
        
        
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
        
        GoA = solveGo(self.lag_, Gd=1/pOn['tau_deact'].value)
        GoB = solveGo(self.lag_, Gd=max(pOff['Gd1'].value, pOff['Gd2'].value))
        
        corrFac = lambda Gact, Gdeact: 1 + Gdeact / Gact
        Gd = max(pOff['Gd1'].value, pOff['Gd2'].value)
        
        print('Lag method (tau_deact): Go = {}, cf={} --> g0 = {}'.format(GoA, corrFac(GoA, 1/pOn['tau_deact'].value), 1e6 * self.peak_ * corrFac(GoA, 1/pOn['tau_deact'].value) / (self.V - E))) #(1 + 1 / (GoA * pOn['tau_deact'].value))
        print('Lag method (max(Gd1,Gd2)): Go = {}, cf={} --> g0 = {}'.format(GoB, corrFac(GoB, Gd), 1e6 * self.peak_ * corrFac(GoB, Gd) / (self.V - E) )) #(1 + max(pOff['Gd1'].value, pOff['Gd2'].value)/GoB)
        print('Exp method (tau_deact): Gact = {}, cf={} --> g0 = {}'.format(1/pOn['tau_act'].value, corrFac(1/pOn['tau_act'].value, 1/pOn['tau_deact'].value), 1e6 * self.peak_ * corrFac(1/pOn['tau_act'].value, 1/pOn['tau_deact'].value) / (self.V - E) )) #(1 + pOn['tau_act'].value / pOn['tau_deact'].value)
        print('Exp method (max(Gd1,Gd2)): Gact = {}, cf={} --> g0 = {}'.format(1/pOn['tau_act'].value, corrFac(1/pOn['tau_act'].value, Gd), 1e6 * self.peak_ * corrFac(1/pOn['tau_act'].value, Gd) / (self.V - E) )) #(1 + pOn['tau_act'].value * max(pOff['Gd1'].value, pOff['Gd2'].value))
        
        plt.plot(toff, calcOff(pOff,toff-toff[0]), label='Off-Fit (Bi-Exp) $\\tau_{{off1}}={:.3g}, \\tau_{{off2}}={:.3g}$'.format(1/pOff['Gd1'].value, 1/pOff['Gd2'].value))  
        
        #plt.show(block=False)
        plt.legend(loc='best')
        
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
            popt = fitPeaks(self.t[postPeak], self.I[postPeak], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
            Iss = popt[2]
        
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
    
    
    def plot(self, ax=None, light='shade', dark=None, addFeatures=True, colour=None, linestyle=None): #colour, linestyle
        """Plot the photocurrent
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
        
        plt.xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
        plt.xlim((self.begT, self.endT))
        plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
        
        #ax.spines['left'].set_position('zero') # y-axis
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero') # x-axis
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
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
                    ax.annotate('$I_{{ss}} = {:.3g}\mathrm{{nA}}$'.format(self.sss_[p]), xy=(self.pulses[p,1], self.sss_[p]), xytext=(toffset+self.pulses[p,1], self.sss_[p]), arrowprops=dict(arrowstyle="wedge,tail_width=0.6", shrinkA=5, shrinkB=5), horizontalalignment='left', verticalalignment='center', fontsize=config.eqSize)
                
                # Add labels for on and off phases
                #ymin, ymax = plt.ylim()
                #plt.ylim(round_sig(ymin,3), round_sig(ymax,3))
                #pos = 0.95 * abs(ymax-ymin)
                arrowy = 0.085 #0.075
                texty = 0.05
                #awidth=10
                ax.annotate('', xy=(self.pulses[p,0], arrowy), xytext=(self.pulses[p,1], arrowy), arrowprops=dict(arrowstyle='<->',color='blue'))
                plt.text(self.pulses[p,0]+self.onDs[p]/2, texty, '$\Delta on_{}={:.3g}\mathrm{{ms}}$'.format(p, self.onDs[p]), ha='center', va='bottom', fontsize=config.eqSize)
                if p < self.nPulses-1:
                    end = self.pulses[p+1,0]
                else:
                    end = self.endT
                ax.annotate('', xy=(self.pulses[p,1], arrowy), xytext=(end, arrowy), arrowprops=dict(arrowstyle='<->',color='green'))
                plt.text(self.pulses[p,1]+self.offDs[p]/2, texty, '$\Delta off_{}={:.3g}\mathrm{{ms}}$'.format(p, self.offDs[p]), ha='center', va='bottom', fontsize=config.eqSize)
        
        
        return # ax
    
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
            self.Iorig = copy(self.I)
            self.isFiltered = True
            I = self.Iorig
        else:
            self.Iprev = copy(self.I)
            I = self.Iprev
        
        # Moving average
        nPoints = int(round(t_window/self.dt))
        self.I = np.convolve(I, np.ones(nPoints)/nPoints, mode='same')
        


class ProtocolData():
    """Container for PhotoCurrent data from parameter variations in the same protocol"""
    
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
        
        #ax.spines['left'].set_position('zero') # y-axis
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero') # x-axis
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
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
