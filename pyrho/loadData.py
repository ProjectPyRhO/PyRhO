import numpy as np
import scipy.io as sio # Use for Matlab files < v7.3
from lmfit import Parameters
from .models import * # For times2cycles and cycles2times
#from config import verbose
import warnings


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

#def loadPhotoCurrent(key,dict):
#    I_phi = dict[key]
# Check for sign of current, subtract any shift
#    return


class PhotoCurrent():
    """Data storage class for an individual Photocurrent and its associated properties"""
    
    def __init__(self, I, t, phi, V, pulses, label=None):
        """ I       := Photocurrent [nA]
            t       := Time series (or time step) [ms]
            phi     := Stimulating flux (max) [ph*mm^-2*s^-1]
            V       := Voltage clamp potential [mV]
            pulses  := Pairs of time points describing the beginning and end of stimulation 
                        e.g. [[t_on1,t_off1],[t_on2,t_off2],...]"""
                        
        ### Load data
        self.I = I                  # Array of photocurrent values
        
        if len(t) == len(I):
            assert(len(t) > 1)
            self.t = t              # Corresponding array of time points [ms]
            tdiff = t[1:] - t[:-1]
            self.dt = tdiff.sum()/len(tdiff)    # (Average) step size
            self.sr = 1000/(self.dt)# Sampling rate [samples/s]
        elif len(t) == 1:           # Assume time step is passed rather than time array
            assert(t > 0)
            self.t = t*range(len(I))
            self.dt = t             # Step size
            self.sr = 1000/t        # Sampling rate [samples/s]
        else:
            raise ValueError("Dimension mismatch: t must be either an array of the same length as I or a scalar defining the timestep!")
        
        self.begT = self.t[0]       # Beginning trial time
        self.endT = self.t[-1]      # Last trial time point
        self.totT = self.endT - self.begT      # Total trial time #max(self.t) # Handles negative delays
        
        ### Load metadata
        pulses = np.asarray(pulses)
        
        #if isinstance(pulses[0], int) or isinstance(pulses[0], int):
        #    pulses = [pulses]
        self.pulses = pulses # nPulses x 2 array [t_on, t_off]      
        self.nPulses = self.pulses.shape[0]
        self.pulseCycles = times2cycles(self.pulses, self.endT) #self.totT)
        
        #if self.nPulses > 1: ### Revise this!!!
        #self.delDs = [row[0] for row in pulses] # pulses[:,0]    # Delay Durations
        self.delD = self.pulses[0,0] - self.begT # Handles negative delays
        self.onDs = [row[1]-row[0] for row in pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        #if self.nPulses > 1:
        #    self.IPIs = np.zeros(self.nPulses - 1)
        #    for p in range(0, self.nPulses-1):
        #        self.IPIs[p] = self.pulses[p+1,0] - self.pulses[p,1]
        self.IPIs = np.asarray([self.pulses[p+1,0] - self.pulses[p,1] for p in range(self.nPulses-1)]) # end <-> start
        self.offDs = np.append(self.IPIs, self.endT-self.pulses[-1,1])
        # self.offDs = [self.totT-((onD+pOff)*nPulses)-delD for pOff in pulseCycles[:,1]]    
        
        #for p in self.nPulses: # List comprehension instead?
        #   self.pulseInds[p,0] = np.searchsorted(self.t, pulses[p,0], side="left")  # CHECK: last index where value <= t_on
        #   self.pulseInds[p,1] = np.searchsorted(self.t, pulses[p,1], side="left")  # CHECK: last index where value <= t_off
        #self.pulseInds = np.array([[np.searchsorted(self.t, pulses[p,time]) for time in range(2)] for p in range(self.nPulses)])
        self.pulseInds = np.array([np.searchsorted(self.t, pulses[p,:]) for p in range(self.nPulses)])
        
        #if self.nPulses > 1:
        #    self.IPIs = [self.pulses[p,0]-self.pulses[p-1,1] for p in range(1,self.nPulses)]  # end <-> start
            # self.IPIs = [self.pulses[p,0]-self.pulses[p-1,0] for p in range(1,self.nPulses)]    # start <-> start
        
        ### Record Experimental constants
        self.V = V                  # Clamp Voltage [mV]: None if no clamp was used
        self.phi = phi              # Light intensity
        # Future inclusions
        # self.Lambda            # Wavelength
        # self.pH                   # pH
        # self.Temp                 # Temperature
        self.label = label          # Optional trial label e.g. "saturate"
        
        ### Calibrate - correct any current offset in experimental recordings
        Idel = self.getDelayPhase()
        offset = np.mean(Idel[:int(round(0.9*len(Idel)))+1]) # Calculate the mean over the first 90% to avoid edge effects
        if abs(offset) > 0.01 * abs(max(self.I) - min(self.I)): # Recalibrate if the offset is more than 1% of the span
            self.I -= offset
            if verbose > 0:
                print("Photocurrent recalibrated by {} [nA]".format(offset))
        
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
        
        
        
        
        ### Extract properties from the data
        # Add this to findPeaks
        self.Irange = [min(self.I), max(self.I)]
        self.Ispan = self.Irange[1] - self.Irange[0]
        #if abs(self.Irange[0]) > abs(self.Irange[1]):
        #    self.Ipeak = self.Irange[0] # Min
        #    self.Ipeaks = np.asarray([min(self.getCycle(p)) for p in range(self.nPulses)]) # Peak may occur after stimulation #np.asarray([min(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]]) for p in range(self.nPulses)])
        #else:
        #    self.Ipeak = self.Irange[1] # Max
        #    self.Ipeaks = np.asarray([max(self.getCycle(p)) for p in range(self.nPulses)])
            #self.Ipeaks = np.asarray([max(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]]) for p in range(self.nPulses)])
        #np.asarray([max(abs(self.I[self.pulseInds[p,0]:self.pulseInds[p,1]])) for p in range(self.nPulses)])
        
        self.peakInd = np.argmax(abs(self.I)) #np.searchsorted(self.I, self.Ipeak)
        self.tpeak = self.t[self.peakInd]
        self.Ipeak = self.I[self.peakInd]
        
        self.peakInds = np.asarray([np.argmax(abs(self.getCycle(p))) for p in range(self.nPulses)]) #np.searchsorted(self.I, self.Ipeaks)
        self.tpeaks = self.t[self.peakInds]
        self.Ipeaks = self.I[self.peakInds]

        if label == 'saturate':     ##### Remove...
            self.Isat = self.Ipeak
            self.Ipmax = self.Ipeak # Deprecate
        
        
        self.Iss = findPlateauCurrent(self.I) # self.findPlateaus???
        
        if verbose > 1:
            print("Photocurrent data loaded! nPulses={}; Total time={}ms; Range={}nA".format(self.nPulses, self.totT, str(self.Irange)))
    
    def alignToPulse(self):
        """Set time array so that the first pulse occurs at t=0 (with negative delay period)"""
        if abs(self.pulses[0,0]) > 1e-12: # HACK!!! Write compare floats function
            p0 = self.delD #pulses[0,0]
            self.t -= p0            # Time array
            self.pulses -= p0       # Pulse times
            self.begT = self.t[0]   # Beginning Time of Trial
            self.endT = self.t[-1]  # End Time of Trial
            self.tpeak = self.t[self.peakInd]
            self.tpeaks = self.t[self.peakInds]
            
        
    def alignToTime(self):
        """Set time array so that it begins at t=0 (with the first pulse at t>0)"""
        if abs(self.pulses[0,0]) < 1e-12:
            p0 = self.delD
            self.t += p0            # Time array
            self.pulses += p0       # Pulse times
            self.begT = self.t[0]   # Beginning Time of Trial
            self.endT = self.t[-1]  # End Time of Trial
            self.tpeak = self.t[self.peakInd]
            self.tpeaks = self.t[self.peakInds]
    
    
    #def printSummary(self):
    def __str__(self):
        """Print out summary details of the photocurrent"""
        str = 'Photocurrent with {} stimulation periods sampled at {} samples/s over {} ms'.format(self.nPulses, self.sr, self.totT)
        return str

    
    ### Move findPeaks from models.py to here?
    # def findPeaks(self): ### See findPeaks in models.py
        # self.peakInds = findPeaks(self.I) ### This needs some careful tweaking for real data...
        # self.t_peaks = self.t[self.peakInds]
        # self.I_peaks = self.I[self.peakInds]
        
    def findPlateaus(self, pulse=0, method=0, window=tFromOff): ### c.f. findPlateauCurrent() in models.py
        # Find plateau
        
        assert(0 <= pulse < self.nPulses)
        
        offInd = self.pulseInds[pulse][1] #np.searchsorted(t,onD+delD,side="left")
        
        if self.onDs[pulse] < window:
            raise ValueError('Error: The plateau buffer must be shorter than the on phase!')
            #windowInd = int(round(p*len(I_phi))) #np.searchsorted(t,t[onEndInd]-100,side="left") # Generalise
            #I_ss = np.mean(I_phi[-windowInd:])
        
        if method == 0: # Empirical
            # Calculate Steady-state as the mean of the last 50ms of the On phase
            tFromOffInd = np.searchsorted(t,t[offInd]-window,side="left")
            self.Iss = np.mean(self.I[tFromOffInd:offInd+1])
            #windowInd = int(round(p*len(self.I))) #np.searchsorted(t,t[onEndInd]-100,side="left") # Generalise
            #self.Iss = np.mean(self.I[-windowInd:])
            
            # Calculate step change (or gradient with t[1:] - t[:-1])
            Idiff = self.I[tFromOffInd+1:offInd+1] - self.I[tFromOffInd:offInd]
            if abs(np.mean(Idiff)) > 0.01 * self.Ispan:
                warnings.warn('Steady-state Convergence Warning: The average step size is larger than 1% of the current span!')

        elif method == 1: # Analytical # Try this first, then resort to empirical method?
            # Fit curve from peak to end of on phase
            popt=fitPeaks(t[peakInds[0]:offInd+1], self.I[peakInds[0]:offInd+1], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
            self.Iss = popt[2]
        
        return self.Iss

    def getDelayPhase(self):
        return self.I[:self.pulseInds[0,0]+1]
        
    def getOnPhase(self, pulse=0):
        return self.I[self.pulseInds[pulse,0]:self.pulseInds[pulse,1]+1]
        
    def getOffPhase(self, pulse=0):
        if 0 <= pulse < self.nPulses-1:
            return self.I[self.pulseInds[pulse,1]:self.pulseInds[pulse+1,0]+1]
        elif pulse == self.nPulses-1:   # Last Pulse
            return self.I[self.pulseInds[pulse,1]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
    
    def getCycle(self, pulse=0):
        if 0 <= pulse < self.nPulses-1:
            return self.I[self.pulseInds[pulse,0]:self.pulseInds[pulse+1,0]+1] # Consider removing the +1
        elif pulse == self.nPulses-1:   # Last Pulse
            return self.I[self.pulseInds[pulse,0]:]
        else:
            raise IndexError("Error: Selected pulse out of range!")
    
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
        
    def filterData(self):
        """Pass frequency bands to filter out"""
        # http://wiki.scipy.org/Cookbook/SignalSmooth
        # http://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
        # http://www.nehalemlabs.net/prototype/blog/2013/04/05/an-introduction-to-smoothing-time-series-in-python-part-i-filtering-theory/
        pass


        
# class Empirical():
    # """Container for empirically derived parameters"""
    # def __init__(self, E=0.0, phi0=None, gam=None, A=None):
        # #self.protocol = protocol
        # self.E = E          
        # self.phi0 = phi0
        # self.gam = gam
        # self.A = A
        
        #Name  #Value   #Units                   #Description
        #E    = 0      # [mV]                     Reversal potential for opsin

        # 4 & 6 State models only
        #phi0 = 1e14   # [photons * s^-1 * mm^-2] Normalising photon flux density ~ 1/10 of threshold
        #gam  = 0.05   # [Dimensionless]          Ratio of single channel conductances: O2/O1
        #A    = 31192  # [mu m^2]                 Effective area of the cell / compartment


        
        
        
        
### Add functions to bundle Protocol data with Empirical parameters into a dataSet dictionary, check for minimum requirements and report which models/features can be fit from the dataSet. ###




# bundle['varyIPI'][IPI].I


# def loadDataSet(self):
    # filename = '%s-in-%s.hoc'%(opsin,loc[0])
    # print("filename =",filename)
    # try:
        # h.load_file(filename)
        # print("loaded %s %s successful"%(opsin,loc[0]))
    # except:
        # print("Error: could not load hoc file for opsin/location: ", filename)
    # ss = '%s_in_%s(%s)'%(opsin,loc[0],','.join([str(x) for x in loc[1]]))
    # print(ss)
    # h(ss)



class ProtocolData():
    """Container for PhotoCurrent data from parameter variations in the same protocol"""
    
    def __init__(self,protocol,nRuns,phis,Vs):
        #self.trials = [] # Array of PhotoCurrent objects
        # self.trials.append(pc)
        self.protocol = protocol
        
        self.nRuns = nRuns
        self.phis = phis
        self.nPhis = len(phis)
        self.Vs = Vs
        self.nVs = len(Vs)
        
        self.trials = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)] # Array of PhotoCurrent objects
        
        #PD = [{'PC':PC, 'run':run, 'phi':phi, 'V':V, ...},{},... ]
        self.metaData = {'nRuns':self.nRuns, 'nPhis':self.nPhis, 'nVs':self.nVs}
        
        
        
    def __str__(self):
        return 'Protocol data set: [nRuns={}, nPhis={}, nVs={}]'.format(self.nRuns, self.nPhis, self.nVs)
        
    #def __info__:
    #    """Report data set features"""
    #    pass
    
    # def getSlice(self):
        
        # for run in range(self.nRuns):
            # for phiInd in range(self.nPhis):
                # for vInd in range(self.nVs):
                    
    
    
    # Function to get array based on independent variable???
    # e.g. IPIs, onDs, phis etc. 
    
        
    def getIpmax(self): 
        """Find the maximum peak current for the whole data set. This is useful when the 'saturate' protocol is absent"""
        self.Ipmax = 0
        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
                    if abs(self.trials[run][phiInd][vInd].Ipeak) > abs(self.Ipmax):
                        self.Ipmax = self.trials[run][phiInd][vInd].Ipeak
                        rmax = run
                        pmax = phiInd
                        vmax = vInd
        return self.Ipmax, (rmax, pmax, vmax)
    
    
    def getProtPeaks(self):
        """Return the set of maximum (absolute) peak currents across a whole set of photocurrents"""
        if self.nRuns > 1:
            phiInd = 0
            vInd = 0
            self.IrunPeaks = [self.trials[run][phiInd][vInd].Ipmax for run in range(self.nRuns)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].tpeak for run in range(self.nRuns)]
            Ipeaks = self.IrunPeaks
            tpeaks = self.trunPeaks
        if self.nPhis > 1:
            run = 0
            vInd = 0
            self.IphiPeaks = [self.trials[run][phiInd][vInd].Ipmax for phiInd in range(self.nPhis)]
            self.trunPeaks = [self.trials[run][phiInd][vInd].tpeak for phiInd in range(self.nPhis)]
            Ipeaks = self.IphiPeaks
            tpeaks = self.trunPeaks
        if self.nVs > 1:
            run = 0
            phiInd = 0
            self.IVPeaks = [self.trials[run][phiInd][vInd].Ipmax for vInd in range(self.nVs)]
            self.tVPeaks = [self.trials[run][phiInd][vInd].tpeak for vInd in range(self.nVs)]
            Ipeaks = self.IVPeaks
            tpeaks = self.tVPeaks
        return Ipeaks, tpeaks
        
    
    def getSteadyStates(self,run=0):
        
        assert(self.nVs > 1)
        self.Iplats = np.zeros((self.nPhis,self.nVs))
        self.Vplats = np.zeros((self.nPhis,self.nVs))
        for phiInd, phi in enumerate(self.phis): 
            for vInd, V in enumerate(self.Vs): 
                self.Iplats[phiInd,vInd] = self.trials[run][phiInd][vInd].Iss # Variations along runs are not useful here
                self.Vplats[phiInd,vInd] = self.trials[run][phiInd][vInd].V
        return self.Iplats, self.Vplats
        
        
        # self.Is = Is  # [run][phiInd][vInd]
        # self.phis = phis # List of light intensities during stimulation
        # self.Vs = Vs
        
        # self.nPulses = pulses.shape[0]
        # self.pulseInds = [[np.searchsorted(self.t, pulses[p,time]) for time in range(2)] for p in range(self.nPulses)]

        # self.totT = max(self.t)     # Total trial time
        
        # IPI = 75 # Inter-Pulse-Interval
        # pDs = [1,2,3,5,8,10,20]
        # nRuns = len(pDs)
        # pulseCycles=np.column_stack((pDs,[IPI-pD for pD in pDs])) # [:,0] = on phase duration; [:,1] = off phase duration
        
        # IPIs = [500,1000,1500,2500,5000,7500,10000]
        # nRuns = len(IPIs)

        
        # ts = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        # Is = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        # IpInds = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        # IssVals = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        # PulseInds = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]

        
from collections import defaultdict        
        
class DataSet():
    """Container for photocurrent data used to produce arrays for parameter extraction"""
    
    def __init__(self,protocol):
        self.data = defaultdict(list)
        
        self.protocol = protocol
        
        # # if protocol == varyPL:
            # # self.data = data
        # # elif protocol == varyIPI:
            # # self.Is = Is
            # # self.ts = ts
            # # pulseCycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI-onD for IPI in IPIs])) # [:,0] = on phase duration; [:,1] = off phase duration

    def addData(self, photoData, protocol):
        self.data[protocol].append(photoData)
        self.phis.append(photoData.phi)



        
# def MyClass(object):
    # def __init__(number):
        # self.number=number

# my_objects = []

# for i in range(100) :
    # my_objects.append(MyClass(i))

# #later

# for obj in my_objects :
    # print obj.number

# # Or:
# ProtocolData.append(TrialData(I,t,V,phi,pulses))
    
# print "Sort the dataSet in place by V ..."
# import operator
# dataSet.sort(key=operator.attrgetter('V'))

# print "Print all phi == 1e14"
# look = "Sanosi"
# for person in personList:
    # if look in person.name:
        # print "%s: \"%s\"" % (person.name, person.quote)

# DataSet.data[pulse[1,0]].I
