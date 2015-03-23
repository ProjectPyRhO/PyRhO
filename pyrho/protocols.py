import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import *
from scipy.optimize import curve_fit
from .parameters import *
from .loadData import * #import loadData
from .models import *
from .config import verbose, saveFigFormat, eqSize, addTitles, addStimulus, colours, styles, dDir, fDir
import pickle
import warnings
import os

#from config import verbose

### Select simulation protocol
#protocols = ['custom', 'saturate', 'inwardRect', 'varyPL', 'varyIPI']
#protocol = protocols[4] #'varyIPI'#'varyPL' # Set this interactively with radio buttons?

# if 'eqSize' not in vars() or 'eqSize' not in globals() or eqSize is None:
    # eqSize = 18 # Move to config

# #plt.tight_layout()
# if 'addTitles' not in vars() or 'addTitles' not in globals() or addTitles is None:
    # addTitles = True

# ### Set default plotting colour and style cycles
# if 'colours' not in vars() or 'colours' not in globals() or colours is None:
    # colours = ['b','g','r','c','m','y','k']
# if 'styles' not in vars() or 'styles' not in globals() or styles is None:
    # styles = ['-', '--', '-.', ':']

# if 'saveFigFormat' not in vars() or 'saveFigFormat' not in globals() or saveFigFormat is None:
    # saveFigFormat = 'png'



### Deprecated!!!
def selectProtocol(protocol):
    """Protocol selection function"""
    if protocol == 'custom':
        return protCustom()
    elif protocol == 'step':
        return protStep()
    elif protocol == 'sinusoid':
        return protSinusoid()
    elif protocol == 'chirp':
        return protChirp()
    elif protocol == 'ramp':
        return protRamp()
    elif protocol == 'saturate':
        return protSaturate()
    elif protocol == 'inwardRect':
        return protInwardRect()
    elif protocol == 'varyPL':
        return protVaryPL()
    elif protocol == 'varyIPI':
        return protVaryIPI()
    else:
        raise NotImplementedError(protocol)
        #print("Error in selecting protocol - please choose from 'custom', 'saturate', 'inwardRect', 'varyPL' or 'varyIPI'")
        
### Protocols to be included in the next version:
### - Temperature (Q10)
### - pH (intracellular and extracellular)
### - Wavelength (lambda)


def characterise(RhO):
    """Run small signal analysis on Rhodopsin"""
    for protocol in smallSignalAnalysis: # .keys()
        RhO.setLight(0.0)
        Prot = selectProtocol(protocol)
        #Prot = protocols[protocol]()
        Prot.runProtocol(RhO)
        Prot.plotProtocol(RhO)
    return

class Protocol(object):
    """Common base class for all protocols"""
    
    def __init__(self, protocol='custom', saveData=True):
        self.protocol = protocol
        #self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        #self.plotPeakRecovery = plotPeakRecovery
        #self.plotStateVars = plotStateVars
        #self.plotKinetics = plotKinetics
    
    def __str__(self):
        return "Protocol type: "+self.protocol
    
    def setParams(self, params):
        for p in params.keys():
            #if p in self.__dict__:
            self.__dict__[p] = params[p].value
            #else:
            #    warnings.warn('Warning: "{p}" not found in {self}'.format(p,self))
        self.prepare()
        self.lam = 470 # Default wavelength [nm]

    def exportParams(self, params):
        """Export parameters to lmfit dictionary"""
        for p in self.__dict__.keys():
            params[p].value = self.__dict__[p]
        return params
            
    def printParams(self):
        for p in self.__dict__.keys():
            print(p,' = ',self.__dict__[p])
            
    def runProtocol(self, Sim, RhO, verbose=verbose): ### rename to run()
        
        self.prepare()
        Sim.prepare(self.dt)
        
        ### Setup containers and run variables
        label = ""
        

        ### Temporary variables to avoid renaming everything...!!!
        protocol = self.protocol
        phis=self.phis
        Vs=self.Vs
        if self.protocol == 'varyPL' or self.protocol == 'varyIPI' or self.protocol == 'sinusoid':
            if self.protocol == 'sinusoid':
                delDs=self.delDs
            pulseCycles = self.pulseCycles
            padDs = self.padDs
        else:
            pulses=self.pulses
        nPulses=self.nPulses
        delDs=self.delDs
        onDs=self.onDs
        offDs=self.offDs
        totT=self.totT
        nRuns=self.nRuns
        dt=self.dt
        
        ### HACKS!!! - revise runTrial()
        delD = delDs[0] 
        onD = onDs[0]
        offD = offDs[0]
        padD = 0.0
        


        # Create container lists for storing all arrays of solution variables
        self.multisoln = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)] # multisoln[runs][phis][Vs][t,s]
        self.ts = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.Is = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.IpInds = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.IpVals = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.IssVals = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.PulseInds = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.labels = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.data = [[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.phiFuncs = [[None for p in range(len(phis))] for r in range(nRuns)]

        if self.saveData:
            # if 'dataTag' not in globals():
                # dataTag = str(RhO.nStates)+"s"
                # self.dataTag = dataTag
                # print('Saving data with tag: '+dataTag)
            self.dataTag = str(RhO.nStates)+"s"
            if verbose > 0:
                print('Saving data with tag: '+self.dataTag)
            dataTag = self.dataTag
            self.PD = ProtocolData(protocol,nRuns,phis,Vs)
            # Add empirical parameters here or to dataSet dictionary?
        
        if verbose > 1: 
            self.printParams()
        
        ### Solve the system of ODEs

        # Loop over the number of runs...             ### Place within V & phi loops to test protocols at different V & phi?
        for run in range(nRuns):
            
            if nRuns > 1:
                if nPulses > 1 or protocol == 'varyPL':
                    onD = pulseCycles[run,0]
                    offD = pulseCycles[run,1]
                    padD = padDs[run]
                else: #if protocol == 'sinusoid':
                    onD = pulseCycles[0,0]
                    offD = pulseCycles[0,1]
                    padD = padDs[0]
                    
            # runExperiment()
            
            # Loop over light intensity...
            for phiInd, phiOn in enumerate(phis): 
                
                if verbose > 1 and (len(phis)>1 or (run==0 and phiInd==0)): # len(phis)>0
                    #print(RhO.phi); print(type(RhO.phi))
                    RhO.dispRates()
                
                # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                for vInd, V in enumerate(Vs): 
                    
                    if protocol == "varyPL":
                        label = "$\mathrm{{Pulse}}={}\mathrm{{ms}}$ ".format(onD)
                        #figTitle += "for varying pulse length "
                    elif protocol == "varyIPI":
                        label = "$\mathrm{{IPI}}={}\mathrm{{ms}}$ ".format(self.IPIs[run])
                        #figTitle += "for varying inter-pulse-interval "
                    elif protocol == 'sinusoid':
                        label = "$f={}\mathrm{{Hz}}$ ".format(round_sig(self.fs[run],3))
                    #else:
                        #figTitle += "\n "
                    
                    if len(phis)>1:
                        label += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    #else:
                        #figTitle += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    
                    if len(Vs)>1:
                        label += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V)
                    #else:
                        #figTitle += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V)
                    


                    ### Run the trial ###
                    start, end = 0.0, totT #0.00, stimD
                    pStart, pEnd = delD, (delD+onD)
                    
                    if self.squarePulse: #protocol in squarePulses: ##### Change after changing 'custom'
                        phi_t = InterpolatedUnivariateSpline([pStart,pEnd],[phiOn,phiOn], k=1, ext=1) 
                        I_RhO,t,soln = Sim.runTrial(RhO, nPulses, V,phiOn,delD,onD,offD,padD,dt)
                        
                    else: # Arbitrary functions of time: phi(t)
                        
                        if protocol == 'sinusoid':
                            t = np.linspace(0.0,onD,(onD*self.sr/1000)+1, endpoint=True)
                            phi_t = InterpolatedUnivariateSpline(pStart + t, self.A0[0] + 0.5*phiOn*(1-np.cos(self.ws[run]*t)), ext=1) # Extrapolate=0 # A0[r]
                        elif protocol == 'chirp':
                            t = np.linspace(0.0,onD,(onD*self.sr/1000)+1, endpoint=True)
                            ft = self.f0*(self.fT/self.f0)**(t/pEnd)
                            phi_t = InterpolatedUnivariateSpline(pStart + t, self.A0[0] + 0.5*phiOn*(np.cos((ft/1000)*t)+1), ext=1) # 0.5*phiOn*(1-np.cos((ft/1000)*t))
                        elif protocol == 'ramp':
                            phi_t = InterpolatedUnivariateSpline([pStart,pEnd], [self.phi_ton,phiOn], k=1, ext=1) #[start,delD,end,totT], [0,self.phi_ton,phiOn,0] 
                        elif protocol == 'custom':
                            t = np.linspace(start,end,((end-start)/self.dt)+1, endpoint=True) # sr or dt... Check!!!
                            #phi_t = InterpolatedUnivariateSpline(t, phi_ft)
                            phi_t = InterpolatedUnivariateSpline([pStart,pEnd],[phiOn,phiOn], k=1, ext=1) ##### Hack!!!!! Remove when custom is generalised to arbitrary functions
                        #elif protocol in squarePulses: #self.squarePulse: ##### Change after changing 'custom'
                            #pulses = np.array([[delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD] for p in range(nPulses)])
                            #phi_t = InterpolatedUnivariateSpline([pStart,pEnd], [phiOn,phiOn], k=1, ext=1) ### Generalise for nPulses!!!
                        
                        I_RhO,t,soln = Sim.runTrialPhi_t(RhO,V,phi_t,delD,onD,totT,dt)
                        
                    self.phiFuncs[run][phiInd] = phi_t
                    
                    if verbose > 1:
                        print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]; label=<{}>'.format(run,nRuns,phiInd,len(phis),vInd,len(Vs),I_RhO.min(),I_RhO.max(),label))
                    #####################
                    
                    
                    # Save simulation results
                    self.multisoln[run][phiInd][vInd] = soln
                    self.ts[run][phiInd][vInd] = t
                    self.Is[run][phiInd][vInd] = I_RhO
                    self.PulseInds[run][phiInd][vInd] = RhO.pulseInd
                    self.labels[run][phiInd][vInd] = label
                    label = ""
                    
                    ### Find peaks in photocurrent - just take min() or max() within each pulse?
                    # if (V < RhO.E): #if abs(min(I_RhO)) > abs(max(I_RhO)):
                        # minmax = np.less
                    # else: 
                        # minmax = np.greater
                    
                    if protocol == 'custom': ### This should be within the loops...!!!
                        #PC.findPeaks()
                        #PC.peakInds = findPeaks(PC.I,minmax)
                        #PC.peakInds = findPeaks(PC.I)
                        #PC.t_peaks = PC.t[peakInds]
                        #PC.I_peaks = PC.I[peakInds]
                        peakInds = findPeaks(PC.I)
                        t_peaks = t[peakInds]
                        I_peaks = I[peakInds]
                        
                        ### Generalise for multiple pulses
                        #PC.findPlateaus(0,1)
                        onEndInd = np.searchsorted(t,onD+delD,side="left")
                        #offInd = PC.pulseInds[0][1]
                        tFromOffInd = np.searchsorted(t,t[onEndInd]-tFromOff,side="left") ### Move into loops # Generalise!!!
                        Iss = np.mean(PC.I[tFromOffInd:onEndInd+1])
                        #self.IssVals[run][phiInd][vInd] = Iss
                
                    elif protocol == "saturate":
        #                 startInd = PulseInds[run][phiInd][vInd][0,0]
        #                 print(startInd)
        #                 extOrder = int(onD/dt)
        #                 peakInds = findPeaks(I_RhO,startInd,extOrder)
        #                 Ipeak = I_RhO[peakInds]
                        if (V < RhO.E): #abs(min(I_RhO)) > abs(max(I_RhO)): #(V < RhO.E): # Find Minima
                            #peakInds = argrelextrema(I_RhO, np.less, order=extOrder)
                            Ipmax = min(I_RhO)
                            peakInds = [np.argmin(I_RhO)]
                        else:       # Find Maxima
                            #peakInds = argrelextrema(I_RhO, np.greater, order=extOrder)
                            Ipmax = max(I_RhO)
                            peakInds = [np.argmax(I_RhO)]# + startInd
                        if verbose > 1:
                            print("I_peak = {}nA; t_peak={}ms; peakInds={}".format(Ipmax,t[peakInds],peakInds))
                    #     Ip = I_RhO[IpInds[run][phiInd][vInd]]
                    
                    elif protocol == 'sinusoid' or protocol == 'chirp':
                        if (V < RhO.E): # Find Minima
                            #peakInds = argrelextrema(I_RhO, np.less, order=extOrder)
                            Ipeak = min(I_RhO)
                            peakInds = [np.argmin(I_RhO)]
                        else:       # Find Maxima
                            #peakInds = argrelextrema(I_RhO, np.greater, order=extOrder)
                            Ipeak = max(I_RhO)
                            peakInds = [np.argmax(I_RhO)]# + startInd
                        if verbose > 1:
                            print("I_peak = {}nA; t_peak={}ms; peakInds={}".format(Ipeak,t[peakInds],peakInds))
                        onEndInd = np.searchsorted(t,onD+delD,side="left")
                        tFromOffInd = np.searchsorted(t,t[onEndInd]-tFromOff,side="left") # Generalise!!!
                        Iss = np.mean(I_RhO[tFromOffInd:onEndInd+1])
                        self.IssVals[run][phiInd][vInd] = Iss
                        
                    elif protocol == 'inwardRect':
                        #### Calculate Steady-state as the mean of the last tFromOff ms of the On phase
                        
                        ### Add test to check that steady state has been reached
                        #if abs(dIdt) < tol:
                        
                        onEndInd = np.searchsorted(t,onD+delD,side="left")
                        tFromOffInd = np.searchsorted(t,t[onEndInd]-tFromOff,side="left") # Generalise!!!
                        Iss = np.mean(I_RhO[tFromOffInd:onEndInd+1])
                        self.IssVals[run][phiInd][vInd] = Iss
                        startInd = self.PulseInds[run][phiInd][vInd][0,0] # Start of first pulse
                        extOrder = int(round(onD/dt))
                        #peakInds = findPeaks(I_RhO,minmax,startInd,extOrder)
                        peakInds = findPeaks(I_RhO,startInd,extOrder)

                    elif protocol == 'varyIPI':
                        ### Search only within the on phase of the second pulse
                        startInd = self.PulseInds[run][phiInd][vInd][1,0]
                        endInd = self.PulseInds[run][phiInd][vInd][1,1]
                        extOrder = int(1+endInd-startInd) #100#int(round(len(I_RhO)/5))
                        #peakInds = findPeaks(I_RhO[:endInd+extOrder+1],minmax,startInd,extOrder)
                        peakInds = findPeaks(I_RhO[:endInd+extOrder+1],startInd,extOrder)
                        if len(peakInds) > 0: # Collect data at the (second) peak
                            self.IpIPI[run] = I_RhO[peakInds[0]] #-1 peaks
                            self.tpIPI[run] = t[peakInds[0]] # tPeaks
                    
                    else: ### Check what this should apply to...
                        startInd = self.PulseInds[run][phiInd][vInd][0,0] # Start of first pulse
                        extOrder = int(round(onD/dt))
                        #peakInds = findPeaks(I_RhO,minmax,startInd,extOrder)
                        peakInds = findPeaks(I_RhO,startInd,extOrder)
                        ### PLOT
                        ### if plotPeakRecovery and peakInds and len(peakInds) > 0: ### Fit curves for recovery kinetics #(nPulses > 1): # and (nRuns == 1): #run+1 == nRuns:
                        ###     plt.figure(Ifig.number)
                        ###     fitPeaks(t[peakInds], I_RhO[peakInds], expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$') #(1,100,1)
                        ###     plt.legend(loc='best')
                        
                    self.IpInds[run][phiInd][vInd] = peakInds
                    self.IpVals[run][phiInd][vInd] = I_RhO[peakInds]
                    
                    
                    
                    if self.saveData:
                        ### Save data to pickle files
                        pulses = np.array([[delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD] for p in range(nPulses)]) #np.array([[delD,onD]])
                        PC=PhotoCurrent(I_RhO,t,phiOn,V,pulses,protocol)
                        self.data[run][phiInd][vInd] = PC
                        #PD.trials.append(PC)
                        self.PD.trials[run][phiInd][vInd] = PC

                        
                    
        
        if protocol == "saturate": # Move into loops?
            # Take the max over all runs, phis and Vs?
            # Ipmax = minmax(self.IpVals[run][phiInd][vInd][:])# = I_RhO[peakInds]
            try: #if V != RhO.E:
                Gmax = Ipmax / (V - RhO.E) # Assuming [O_p] = 1 ##### Should fV also be used?
            except ZeroDivisionError: #else:
                print("The clamp voltage must be different to the reversal potential!")
            
            # if RhO.nStates == 6:
                # try:
                    # gbar_est = (Gmax/RhO.A) * 1e6
                # except ZeroDivisionError:
                    # print("An area > 0 must be specified!")
            # else:
            gbar_est = Gmax * 1e6
            
            if verbose > 0:
                print("Estimated maximum conductance (gmax) = {} uS".format(round_sig(gbar_est,3)))

        
        

        
        if self.saveData:           ##### This should all be PD!!!
            if protocol == 'custom': ### This should be within the loops...!!!
                #PC.findPeaks()
                #PC.peakInds = findPeaks(PC.I,minmax)
                PC.peakInds = findPeaks(PC.I)
                PC.t_peaks = PC.t[peakInds]
                PC.I_peaks = PC.I[peakInds]
                #PC.findPlateaus(0,1)
                #offInd = PC.pulseInds[0][1]
                #tFromOffInd = np.searchsorted(t,t[offInd]-tFromOff,side="left") ### Move into loops # Generalise!!!
                #PC.Iss = np.mean(PC.I[tFromOffInd:offInd+1])              ### Move into loops
        #         popt = fitPeaks(t[peakInds[0]:offInd+1], PC.I[peakInds[0]:offInd+1], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
        #         PC.Iss = popt[2]
                #PC.IssVals = self.IssVals
                
            elif protocol == 'saturate':
                #PC.Ipeak = Ipeak
                PC.Ipmax = Ipmax
                PC.gbar_est = gbar_est
            elif protocol == 'inwardRect':
                PC.IssVals = self.IssVals # Change name here to avoid confusion with scalar Iss?
                PC.Vs = Vs
            elif protocol == 'varyPL':
                pass
            elif protocol == 'varyIPI':
                PC.Ip = self.IpIPI #peaks
                PC.tp = self.tpIPI #tPeaks
                #PC.tau_r = popt[1]
            
            if verbose > 0:
                print("Saving data to disk")
            fh = open(dDir+protocol+dataTag+".pkl","wb")
            pickle.dump(PC,fh)
            fh.close()
            
        return self.PD
    
    

    def plotStimulus(self,phi_t,pulses,totT,ax=None,light='shade'):
        t = np.linspace(0,totT,10*int(round(totT/self.dt))+1) #10001) # 
        
        # if self.protocol == 'chirp':
            # fig, ax1 = plt.subplot()
            # ax2 = twinx()
            # ax1.plot(t,phi_t(t),'b')
            # ax1.set_ylabel('Amplitude')
            # ax2.set_yscale('log')
            # ft = self.f0*(self.fT/self.f0)**(t/self.onD)
            # ax2.plot(t,ft,'g')
            # ax2.sset_ylabel('Instantaneous frequency')
        # else:
        if ax == None:
            #fig = plt.figure()    
            ax = plt.gca()
        else:
            #plt.figure(fig.number)
            plt.sca(ax)
            
        plt.plot(t,phi_t(t))
        ### Finish function to plot the shape of the light stimuli

        if light == 'spectral':
            plotLight(pulses, ax=ax, light='spectral', lam=self.lam)
        else:
            plotLight(pulses, ax=ax, light=light)
        
        plt.xlabel('$\mathrm{Time\ [ms]}$') #(r'\textbf{Time} [ms]')
        plt.xlim((0,totT))
        #plt.ylabel('$\mathrm{\phi\ [photons \cdot s^{-1} \cdot mm^{-2}]}$')
        plt.ylabel('$\mathrm{\phi\ [ph. / s / mm^{2}]}$')
        
        return ax
    
    
    def getLineProps(self, run, vInd, phiInd):
        #global colours
        #global styles
        if verbose > 1 and (self.nRuns>len(colours) or len(self.phis)>len(colours) or len(self.Vs)>len(colours)):
            warnings.warn("Warning: only {} line colours are available!".format(len(colours)))
        if verbose > 0 and self.nRuns>1 and len(self.phis)>1 and len(self.Vs)>1:
            warnings.warn("Warning: Too many changing variables for one plot!")
        if verbose > 2:
            print("Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}".format(run,self.nRuns,phiInd,len(self.phis),vInd,len(self.Vs)))
        if self.nRuns>1:
            col = colours[run%len(colours)]
            if len(self.phis)>1:
                style=styles[phiInd%len(styles)]
            elif len(self.Vs)>1:
                style=styles[vInd%len(styles)]
            else:
                style = '-'
        else:
            if len(self.Vs)>1:
                col = colours[vInd%len(colours)]
                if len(self.phis)>1:
                    style = styles[phiInd%len(styles)]
                else:
                    style = '-'
            else:
                if len(self.phis)>1:
                    col = colours[phiInd%len(colours)]
                    style = '-'
                else:
                    col = 'b'   ### colours[0]
                    style = '-' ### styles[0]
        return col, style
    
    
    ### Move fitPeaks and fitfV to fitting.py ###
    def fitPeaks(self, t_peaks, I_peaks, curveFunc, p0, eqString, fig=None):
        #print(p0)
        shift = t_peaks[0] # ~ delD
    #     if protocol == 'varyIPI':
    #         plt.ylim(ax.get_ylim()) # Prevent automatic rescaling of y-axis
        popt, pcov = curve_fit(curveFunc, t_peaks-shift, I_peaks, p0=p0) #Needs ball-park guesses (0.3, 125, 0.5)
        peakEq = eqString.format(*[round_sig(p,3) for p in popt]) # *popt rounded to 3s.f.
        
        if fig:
            plt.figure(fig.number) # Select figure
    #     ext = 10 # Extend for ext ms either side
    #     xspan = t_peaks[-1] - t_peaks[0] + 2*ext 
    #     xfit=np.linspace(t_peaks[0]-ext-shift,t_peaks[-1]+ext-shift,xspan/dt)
            plt.plot(t_peaks, I_peaks, linestyle='', color='r', marker='*')
            xfit=np.linspace(-shift,self.totT-shift,self.totT/self.dt) #totT
            yfit=curveFunc(xfit,*popt)
            
            plt.plot(xfit+shift,yfit,linestyle=':',color='#aaaaaa',linewidth=1.5*mp.rcParams['lines.linewidth'])#,label="$v={:+} \mathrm{{mV}}$, $\phi={:.3g}$".format(V,phiOn)) # color='#aaaaaa' 
            #ylower = copysign(1.0,I_peaks.min())*ceil(abs((I_peaks.min()*10**ceil(abs(log10(abs(I_peaks.min())))))))/10**ceil(abs(log10(abs(I_peaks.min()))))
            #yupper = copysign(1.0,I_peaks.max())*ceil(abs((I_peaks.max()*10**ceil(abs(log10(abs(I_peaks.max())))))))/10**ceil(abs(log10(abs(I_peaks.max()))))
        #     if (len(Vs) == 1) and (len(phis) == 1) and (nRuns == 1):
        #         x, y = 0.8, 0.9
        #     else:
            x = 0.8
            y = yfit[-1] #popt[2]
            
            plt.text(x*self.totT,y,peakEq,ha='center',va='bottom',fontsize=eqSize) #, transform=ax.transAxes)
        
        print(peakEq)
        if verbose > 1:
            print("Parameters: {}".format(popt))
            if type(pcov) in (tuple, list):
                print("$\sigma$: {}".format(np.sqrt(pcov.diagonal())))
            else:
                print("Covariance: {}".format(pcov))
        return popt, pcov, peakEq
    
    
    # def calcIssfromfV(V,v0,v1,E):#,G): # Added E as another parameter to fit
        # ##[s1s, s2s, s3s, s4s, s5s, s6s] = RhO.calcSteadyState(RhO.phiOn)
        # ##psi = s3s + (RhO.gam * s4s) # Dimensionless
        
        # #E = RhO.E
        # if type(V) != np.ndarray:
            # V = np.array(V)
        # fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
        # fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
        # ##psi = RhO.calcPsi(RhO.steadyStates) ### This is not necessary for fitting!!!
        # ##g_RhO = RhO.gbar * psi * fV # Conductance (pS * mu m^-2)
        # ##I_ss = RhO.A * g_RhO * (V - E) # Photocurrent: (pS * mV)
        # #I_ss = G * fV * (V-E)
        # ##return I_ss * (1e-6) # 10^-12 * 10^-3 * 10^-6 (nA)
        # return fV * (V - E)
    
    def fitfV(self, Vs, Iss, curveFunc, p0, RhO, fig=None):#, eqString): =plt.gcf()
        if fig==None:
            fig=plt.gcf()
        markerSize=40
        eqString = r'$f(V) = \frac{{{v1:.3}}}{{V-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{V-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'
        psi = RhO.calcPsi(RhO.steadyStates)
        #sf = RhO.A * RhO.gbar * psi * 1e-6 # Six-state only
        sf = RhO.g * psi * 1e-6 
        fVs = np.asarray(Iss)/sf # np.asarray is not needed for the six-state model!!!
        popt, pcov = curve_fit(curveFunc, Vs, fVs, p0=p0) # (curveFunc, Vs, Iss, p0=p0)
        pFit = [round_sig(p,3) for p in popt]
        #peakEq = eqString.format(pFit[0],pFit[2],pFit[2],pFit[1])
        peakEq = eqString.format(v1=pFit[0],E=pFit[2],v0=pFit[1])
        
        Vrange = max(Vs)-min(Vs)
        xfit=np.linspace(min(Vs),max(Vs),Vrange/.1) #Prot.dt
        yfit=curveFunc(xfit,*popt)*sf
        
        #peakEq = eqString.format(*[round_sig(p,3) for p in popt])
        
        fig.plot(xfit,yfit)#,label=peakEq)#,linestyle=':', color='#aaaaaa')
        #col, = getLineProps(Prot, 0, 0, 0) #Prot, run, vInd, phiInd
        #plt.plot(Vs,Iss,linestyle='',marker='x',color=col)
        fig.scatter(Vs,Iss,marker='x',color=colours,s=markerSize)#,linestyle=''
        
        # x = 1 #0.8*max(Vs)
        # y = 1.2*yfit[-1]#max(IssVals[run][phiInd][:])
        # plt.text(-0.8*min(Vs),y,peakEq,ha='right',va='bottom',fontsize=eqSize)#,transform=ax.transAxes)
        
        if verbose > 1:
            print(peakEq)
        return popt, pcov, peakEq
    
    
    
    
    
    def plotProtocol(self, Sim, RhO, verbose=verbose): # Remove RhO as an argument? # Rename to plot()
        
        protocol = self.protocol
        phis=self.phis
        Vs=self.Vs
        if self.protocol == 'varyPL':
            pulseCycles = self.pulseCycles
            padDs = self.padDs
            pDs = self.pDs
        elif self.protocol == 'varyIPI' or self.protocol == 'sinusoid':
            pulseCycles = self.pulseCycles
            padDs = self.padDs
        else:
            pass#pulses=self.pulses
        nPulses=self.nPulses
        delDs=self.delDs
        onDs=self.onDs
        offDs=self.offDs
        totT=self.totT
        nRuns=self.nRuns
        #dt=self.dt
        ts=self.ts
        Is=self.Is
        multisoln=self.multisoln
        labels=self.labels
        IpInds=self.IpInds
        IssVals=self.IssVals
        
        ### HACKS!!! - revise runTrial()
        delD = delDs[0] 
        onD = onDs[0]
        offD = offDs[0]
        # padD = 0.0
        
        ### Move all plotting in here?
        # Ifig = self.Ifig
        # ax = self.ax
        
        if self.saveData:# and 'dataTag' not in globals():
            dataTag = self.dataTag
        
        figTitle = "Photocurrent through time "
        
        
        ### Additional plots
        # plt.figure(Ifig.number)
        # if protocol == "varyPL":
            # ymin, ymax = plt.ylim()
            # pos = 0.02 * abs(ymax-ymin)
            # print("[{},{}] --> {}".format(ymin,ymax,pos))
            
        # elif protocol == "varyIPI":
            # ymin, ymax = plt.ylim()
            # pos = 0.02 * abs(ymax-ymin)
            # plt.ylim(ax.get_ylim())
            # #plt.ylim(ymin*1.1,ymax)
        # else:
            # plt.ylim(ax.get_ylim()) # Prevent automatic rescaling of y-axis
        # print("y-axis limits: {}".format(plt.ylim()))

        # Loop over the number of runs ### Place within V & phi loops to test protocols at different V & phi?
        for run in range(nRuns): 
            if nRuns > 1:
                if nPulses > 1 or protocol == 'varyPL':
                    onD = pulseCycles[run,0]
                    offD = pulseCycles[run,1]
                    padD = padDs[run]
                else: #if protocol == 'sinusoid':
                    onD = pulseCycles[0,0]
                    offD = pulseCycles[0,1]
                    padD = padDs[0]
            
            # Loop over light intensity...
            for phiInd, phiOn in enumerate(phis): #for phiInd in range(0, len(phis)):
                #print('phis[{}]={}'.format(phiInd,phiOn))
        #         phiOn = phis[phiInd]
                
                # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                for vInd, V in enumerate(Vs): #range(0, len(Vs)):
        #             V = Vs[vInd]
                    
                    if verbose > 1:
                        print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}'.format(run,nRuns,phiInd,len(phis),vInd,len(Vs)))
                    
                    t = ts[run][phiInd][vInd]
                    I_RhO = Is[run][phiInd][vInd]
                    label = labels[run][phiInd][vInd]
                    
                    peakInds = IpInds[run][phiInd][vInd]
                    
                    if protocol == "varyPL":
                        #label = "$\mathrm{{Pulse}}={}\mathrm{{ms}}$ ".format(onD)
                        figTitle += "for varying pulse length "
                    elif protocol == "varyIPI":
                        #label = "$\mathrm{{IPI}}={}\mathrm{{ms}}$ ".format(IPIs[run])
                        figTitle += "for varying inter-pulse-interval "
                    elif protocol == "sinusoid":
                        figTitle += "$\phi_0 = {:.3g}$ ".format(self.A0[0]) # A0[r]
                    elif protocol == 'chirp':
                        figTitle += "$\phi_0 = {:.3g}; f_0={}, f_T={}$ ".format(self.A0[0],self.f0,self.fT)
                    else:
                        figTitle += "\n "
                    
                    if len(phis)>1:
                        pass
                        #label += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    else:
                        figTitle += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    
                    if len(Vs)>1:
                        pass
                        #label += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V)
                    else:
                        figTitle += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V)
                    
                    ##### PLOT #####
                    
                    
                    ### Plot photocurrent
                    col, style = self.getLineProps(run, vInd, phiInd)
                    
                    if (vInd == 0) and (phiInd == 0) and (run == 0): #or (len(Vs) == 1):
                        figWidth, figHeight = mp.rcParams['figure.figsize']
                        ### Plot Stimulus
                        if 'pulses' in self.__dict__:
                            if not addStimulus: # Create a separate stimulus figure
                                stimFig = self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,ax=None,light='spectral')
                                if max(phis)/min(phis) >= 100:
                                    plt.yscale('log')
                            # else:
                                # if protocol == 'saturate':
                                    # figHeight *= 1.5
                        
                        Ifig = plt.figure(figsize=(figWidth, figHeight))
                        if protocol == 'varyPL':
                            gsPL = plt.GridSpec(2,3)
                            axLag = Ifig.add_subplot(gsPL[0,-1])
                            axPeak = Ifig.add_subplot(gsPL[1,-1],sharex=axLag)
                            axI = Ifig.add_subplot(gsPL[:,:-1])
                            
                        elif  protocol == 'saturate':
                            if addStimulus: 
                                gsStim = plt.GridSpec(4,1)
                                axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
                                axI = Ifig.add_subplot(gsStim[1:,:],sharex=axS) # Photocurrent axes
                                self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,axS,light='spectral')
                                plt.setp(axS.get_xticklabels(), visible=False)
                                plt.xlabel('')
                                if max(phis)/min(phis) >= 100:
                                    plt.yscale('log')
                            else:
                                axI = Ifig.add_subplot(111)
                            
                            plotLight(self.pulses,axI)
                        
                        elif  protocol == 'ramp':
                            if addStimulus: 
                                gsStim = plt.GridSpec(4,1)
                                axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
                                axI = Ifig.add_subplot(gsStim[1:,:],sharex=axS) # Photocurrent axes
                                self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,axS,light='spectral')
                                plt.setp(axS.get_xticklabels(), visible=False)
                                plt.xlabel('')
                                #if phis[-1]/phis[0] >= 100:
                                #    plt.yscale('log')
                            else:
                                axI = Ifig.add_subplot(111)
                                
                            plotLight(self.pulses,axI)
                        
                        elif  protocol == 'chirp':
                            if addStimulus: 
                                gsStim = plt.GridSpec(4,1)
                                axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
                                axI = Ifig.add_subplot(gsStim[1:,:],sharex=axS) # Photocurrent axes
                                self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,axS,light='spectral')
                                plt.setp(axS.get_xticklabels(), visible=False)
                                plt.xlabel('')
                                axS.set_ylim(self.A0[0],phiOn) ### A0[r]
                                if max(phis)/min(phis) >= 100:
                                    plt.yscale('log')
                                ### Overlay instantaneous frequency
                                tsmooth = np.linspace(0,onD,10001)
                                axf = axS.twinx()
                                axf.set_yscale('log')
                                ft = self.f0*(self.fT/self.f0)**(tsmooth/onD)
                                axf.plot(tsmooth+delD,ft,'g')
                                axf.set_ylabel('$f\ \mathrm{[Hz]}$')
                            else:
                                axI = Ifig.add_subplot(111)
                                
                            plotLight(self.pulses,axI)
                        
                        elif protocol == 'sinusoid':
                            if nRuns >1: #len(phis) > 1: #nRuns???
                                gsSin = plt.GridSpec(2,3)
                                axIp = Ifig.add_subplot(gsSin[0,-1])
                                axIss = Ifig.add_subplot(gsSin[1,-1],sharex=axIp)
                                axI = Ifig.add_subplot(gsSin[:,:-1])
                            else:
                                axI = Ifig.add_subplot(111) # Combine with else condition below
                            plotLight(self.pulses, axI)
                            
                        elif protocol == 'inwardRect':
                            gsIR = plt.GridSpec(1,3)
                            axfV = Ifig.add_subplot(gsIR[0,-1])
                            axI = Ifig.add_subplot(gsIR[0,0:2],sharey=axfV)
                            plotLight(self.pulses,axI)
                            
                        elif protocol =='varyIPI': ############# Tidy up!!!
                            axI = Ifig.add_subplot(111)
                            for p in range(0, nPulses):
                                plt.axvspan(delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD,facecolor='y',alpha=0.2)
                                    
                        else: 
                            axI = Ifig.add_subplot(111)
                            plt.sca(axI)
                            plotLight(self.pulses,axI)
                            
                        plt.sca(axI)
                        plt.xlabel('$\mathrm{Time\ [ms]}$') #(r'\textbf{Time} [ms]')
                        plt.xlim((0,totT))
                        plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
                        if addTitles:
                            plt.title(figTitle) #'Photocurrent through time'
                    
                        
                    
                    else:
                        
                        ### Plot stimulus above photocurrent
                        if 'pulses' in self.__dict__: #verbose > 1 and 
                            if protocol == 'saturate': ### Generalise this!!!
                                self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,axS,None)
                            
                        
                        plt.figure(Ifig.number)
                        if protocol == "varyIPI" and (vInd == 0) and (phiInd == 0): ############# Tidy up!!!
                            #plotLight(self.pulses[1:,:]) # Plot all pulses but the first
                            
                            for p in range(1, nPulses): # Plot all pulses but the first
                                plt.axvspan(delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD,facecolor='y',alpha=0.2)
                        figTitle = "" # Could remove title from loops
                    
                        if protocol == 'ramp':
                            self.plotStimulus(self.phiFuncs[run][phiInd],self.pulses,self.totT,axS,light=None)
                            axS.set_xlabel('')
                    
                    
                    plt.sca(axI)
                    baseline, = axI.plot(t, I_RhO, color=col, linestyle=style, label=label)

                    if protocol == "varyPL":
                        axI.axvline(x=delD,linestyle=':',color='k')#plt.axvline(x=delD+(p*(onD+offD)),linestyle='--',color='k')
                        axI.axvline(x=delD+onD,linestyle=':',color=col)#baseline.get_color())   #plt.axvline(x=delD+(run*(onD+offD))+onD,linestyle='--',color='k')
        
                    # self.Ifig = Ifig
                    # self.ax = ax
            
                    ### PLOT
                    if self.plotPeakRecovery and peakInds and len(peakInds) > 0: ### Fit curves for recovery kinetics #(nPulses > 1): # and (nRuns == 1): #run+1 == nRuns:
                        plt.figure(Ifig.number)
                        self.fitPeaks(t[peakInds], I_RhO[peakInds], expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$',Ifig) #(1,100,1)
                        plt.legend(loc='best')
                    
                    
                    if self.plotKinetics:
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
                        popt, _, _=self.fitPeaks(t[peakInds[0]:onEndInd + 1], I_RhO[peakInds[0]:onEndInd + 1], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
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
                        popt, _, _=self.fitPeaks(t[onEndInd:], I_RhO[onEndInd:], biExpDecay, p0off, '$I_{{off}} = {:.3}e^{{-t/{:g}}} {:+.3}e^{{-t/{:g}}} {:+.3}$','')
                        ### Plot tau_off vs Irrad (for curves of V)
                        ### Plot tau_off vs V (for curves of Irrad)

                        # Draw boundary between ON and INACT phases
                        for p in peakInds:
                            plt.axvline(x=t[p],linestyle=':',color='m')
                            plt.axhline(y=I_RhO[peakInds[0]],linestyle=':',color='r')
                            plt.axhline(y=Iss,linestyle=':',color='b')

                        plt.legend(loc='best')

                    else:
                        if protocol == 'inwardRect':
                            #### Calculate Steady-state as the mean of the last 50ms of the On phase
                            onEndInd = np.searchsorted(t,onD+delD,side="left")
                            tFromOffInd = np.searchsorted(t,t[onEndInd]-tFromOff,side="left") # Generalise!!!
                            Iss = np.mean(I_RhO[tFromOffInd:onEndInd+1])
                            IssVals[run][phiInd][vInd] = Iss
                    
                    ### Plot additional elements for varying pulses or IPIs
                    col, style = self.getLineProps(run, vInd, phiInd)
                    if protocol == "varyPL":
        #                 plt.figure(Ifig.number)
                        # ymin, ymax = plt.ylim()
                        # pos = 0.02 * abs(ymax-ymin)
        
                        # axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)#'b')  #colours[c]
                        axI.plot(t[peakInds],I_RhO[peakInds],marker='*',color=col)
                        axLag.plot(pDs[run],(t[IpInds[run][phiInd][vInd]] - delD),marker='*',markersize=10,color=col)
                        axPeak.plot(pDs[run],I_RhO[peakInds],marker='*',markersize=10,color=col)
                        
                    #elif protocol == "varyIPI":
                        #plt.figure(Ifig.number)
                        ##plt.annotate('', (delD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col})
                        #plt.annotate('', (delD+onD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col})
                        ##arrow(delD, (run+1)*pos, onD+offD, 0, fc=colours[c], ec=colours[c], length_includes_head=True, head_width=0.01, head_length=1)
                    ##c=c+1
        #                 plt.axvline(x=delD,linestyle='--',color='k')plt.axvline(x=delD+(p*(onD+offD)),linestyle='--',color='k')
        #                 plt.axvline(x=delD+(run*(onD+offD))+onD,linestyle='--',color='k')
        #                 for p in range(0, nPulses):
                        #plt.hlines(y=(run+1)*1e-8,xmin=delD,xmax=delD+pulseCycles[run,0],linewidth=3,color='b')
                    
                    if (verbose > 2 or self.plotStateVars): # and Sim.simulator != 'NEURON': ##### Hack!!!!! #len(Vs) == 1: 
                        if protocol == 'varyPL':
                            pulses = np.array([[0,self.pDs[run]]])+self.delDs[run]
                        else:
                            pulses = self.pulses                        
                        soln = multisoln[run][phiInd][vInd]
                        RhO.plotStates(t,soln,pulses,RhO.labels,phiOn,IpInds[run][phiInd][vInd],'states{}s-{}-{}-{}'.format(RhO.nStates,run,phiInd,vInd))
                


                    #display(Ifig) # Show figure inline 
                    


                    
        # Freeze y-limits
        # ymin, ymax = ylim()
        # ylim(ymin, ymax)

        plt.sca(axI)
        plt.figure(Ifig.number)
        #plt.legend()
        if verbose > 1 and protocol == "varyIPI":
            print(self.tpIPI) # tPeaks
            print(self.IpIPI) # peaks
            #popt = fitPeaks(tPeaks, peaks, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$')
            #print("tau_R = {} ==> rate_R = {}".format(popt[1],1/popt[1]))

            
        if label:
            if protocol == 'custom' or protocol == 'step' or protocol == 'inwardRect':
                if len(Vs) == 1:
                    ncol=1
                else:
                    ncol=len(phis)
                lgd = plt.legend(loc='best', borderaxespad=0, ncol=ncol, fancybox=True) #, shadow=True , bbox_to_anchor=(1.02, 1)
            else:
                lgd = plt.legend(loc='best')

            
        if protocol == "saturate":
        
            plt.figure(Ifig.number)
            plt.axvline(x=t[peakInds[0]],linestyle=':',color='k')
            plt.axhline(y=I_RhO[peakInds[0]],linestyle=':',color='k')
            plt.text(1.05*t[peakInds[0]],1.05*I_RhO[peakInds[0]],'$I_{{peak}} = {:.3g}$'.format(I_RhO[peakInds[0]]),ha='left',va='bottom',fontsize=eqSize)

            
        if protocol == "inwardRect": # and RhO.useInwardRect: # Rewrite fitfV to handle other models
            plt.figure(Ifig.number) #IssVfig = plt.figure()
            ax = axfV #IssVfig.add_subplot(111)
            
            legLabels = [None for p in range(len(phis))]
            for phiInd, phiOn in enumerate(phis): 
                ### PLOT
                RhO.calcSteadyState(phiOn)
                #print(self.IssVals[run][phiInd][:])
                popt, pcov, eqString = self.fitfV(Vs,self.IssVals[run][phiInd][:],calcIssfromfV,p0fV,RhO,ax)#,eqString)
                
                # Add equations to legend
                if len(phis) > 1: 
                    legLabels[phiInd] = eqString + '$,\ \phi={:.3g}$'.format(phiOn)
                else:
                    legLabels[phiInd] = eqString
                
                ### Move this to fitting routines?
                # v0 = popt[0], v1 = popt[1], E = popt[2]
            
            #if len(phis) > 1:
            ax.legend(legLabels,loc='best')
            
            ax.spines['left'].set_position('zero')
            ax.spines['right'].set_color('none')
            ax.spines['bottom'].set_position('zero')
            ax.spines['top'].set_color('none')
            ax.spines['left'].set_smart_bounds(True)
            ax.spines['bottom'].set_smart_bounds(True)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            #ax.yaxis.set_major_formatter(mp.ticker.ScalarFormatter(useMathText=True))
            #ax.yaxis.set_minor_formatter(mp.ticker.ScalarFormatter(useMathText=True))
            
            
            ax.set_xlabel('$V_{clamp}$ $\mathrm{[mV]}$', position=(0.95,0.8)) #plt.xlabel
            #plt.xlim((min(Vs),max(Vs)))
            ax.set_ylabel('$I_{ss}$ $\mathrm{[nA]}$', position=(0.55,0.05))
            
            plt.tight_layout()
            
            #IssVfig.savefig(fDir+protocol+dataTag+"fV."+saveFigFormat, format=saveFigFormat)
        
        
        if protocol == 'sinusoid' and nRuns > 1: #len(phis) > 1:
            plt.figure(Ifig.number)
            axI.legend().set_visible(False)
            
            #if len(self.phis) > 1:
            fstars = np.zeros((len(phis),len(Vs)))
            Itemp = 0.0
            for phiInd, phiOn in enumerate(phis):
                for vInd, V in enumerate(Vs):
                    Ipeaks = np.zeros(nRuns) #[None for r in range(nRuns)]
                    for run in range(nRuns): 
                        Ipeaks[run] = max(abs(self.IpVals[run][phiInd][vInd])) # Maximum absolute value over all peaks from that trial
                        if Ipeaks[run] > Itemp:
                            fpMaxInd = np.argmax(abs(self.IpVals[run][phiInd][vInd]))
                            fpMaxSign = np.sign(self.IpVals[run][phiInd][vInd][fpMaxInd])
                            Itemp = Ipeaks[run]
                    col, style = self.getLineProps(run, vInd, phiInd)
                    axIp.plot(self.fs,Ipeaks,'x',color=col)
                    #intIp = UnivariateSpline(self.fs, Ipeaks)
                    intIp = InterpolatedUnivariateSpline(self.fs, Ipeaks)
                    #intIp = interp1d(self.fs, Ipeaks, kind='cubic')
                    fsmooth = np.logspace(np.log10(self.fs[0]), np.log10(self.fs[-1]), num=100)
                    axIp.plot(fsmooth,intIp(fsmooth))
                    fstar_p = self.fs[np.argmax(Ipeaks)]
                    fstars[phiInd,vInd] = fstar_p
                    Ap = max(Ipeaks)
                    #fpMaxInd = np.argmax(Ipeaks)
                    fpLabel = '$f^*_{{peak}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_p,3))
                    axIp.plot(fstar_p,Ap,'*',markersize=10)
                    #axIp.annotate(fpLabel, xy=(fstar_p,Ap), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})
            axIp.set_xscale('log')
            axIp.set_ylabel('$|A|_{peak}$ $\mathrm{[nA]}$')
            if addTitles:
                axIp.set_title('$\mathrm{|Amplitude|_{peak}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
            #axIp.set_aspect('auto')
                    
            # Calculate the time to allow for transition effects from the period of fstar_p
            buffer=2
            fstar_p=max(max(fstars))
            transD = buffer*np.ceil(1000/fstar_p) # [ms]
            transEndInd = round((delD+transD)/self.dt)
            if transEndInd >= (self.onDs[0])/self.dt: # If transition period is greater than the on period
                transEndInd = round((delD+self.onDs[0]/2)/self.dt) # Take the second half of the data
            tTransEnd = transEndInd*self.dt #ts[0][0][0]
            axI.axvline(x=tTransEnd,linestyle=':',color='k')
            #t_on = self.pulses[0,0]
            onBegInd = RhO.pulseInd[0,0]
            #t_off = self.pulses[0,1]
            onEndInd = RhO.pulseInd[0,1] # End of first pulse
            #print(self.totT,t[onEndInd])
            axI.annotate('', xy=(tTransEnd, Itemp*fpMaxSign), xytext=(t[onEndInd], Itemp*fpMaxSign), arrowprops={'arrowstyle':'<->','color':'black','shrinkA':0,'shrinkB':0}) ### Removed 'Search Zone' since text shifts the arrow slightly
        
            for phiInd, phiOn in enumerate(phis):
                for vInd, V in enumerate(Vs):
                    Iabs = np.zeros(nRuns) #[None for r in range(nRuns)]
                    for run in range(nRuns): 
                        t = ts[run][phiInd][vInd]
                        I_RhO = Is[run][phiInd][vInd]
                        #transEndInd = np.searchsorted(t,delD+transD,side="left") # Add one since upper bound is not included in slice
                        #onEndInd = np.searchsorted(t,self.PulseInds[run][phiInd][vInd][0,1],side="left") # End of first pulse
                        #onBegInd = RhO.pulseInd[0,0]
                        #onEndInd = RhO.pulseInd[0,1] # End of first pulse
                        #if transEndInd >= len(t): # If transition period is greater than the on period
                        #    transEndInd = round(len(t[onBegInd:onEndInd+1])/2) # Take the second half of the data
                        #print(fstar_p,'Hz --> ',transD,'ms;', transEndInd,':',onEndInd+1)
                        I_zone = I_RhO[transEndInd:onEndInd+1]
                        #print(I_zone)
                        try:
                            maxV=max(I_zone)
                        except ValueError:
                            maxV=0.0
                        try:
                            minV=min(I_zone)
                        except ValueError:
                            minV=0.0
                        Iabs[run] = abs(maxV-minV)
                    
                    #axI.axvline(x=t[transEndInd],linestyle=':',color='k')
                    #axI.annotate('Search zone', xy=(t[transEndInd], min(I_RhO)), xytext=(t[onEndInd], min(I_RhO)), arrowprops={'arrowstyle':'<->','color':'black'})
                    col, style = self.getLineProps(run, vInd, phiInd) ### Modify to match colours correctly
                    axIss.plot(self.fs,Iabs,'x',color=col)
                    #intIss = UnivariateSpline(self.fs, Iabs)
                    intIss = InterpolatedUnivariateSpline(self.fs, Iabs)
                    #intIss = interp1d(self.fs, Iabs, kind='cubic')
                    #fsmooth = np.logspace(self.fs[0], self.fs[-1], 100)
                    axIss.plot(fsmooth,intIss(fsmooth))
                    fstar_abs = self.fs[np.argmax(Iabs)]
                    fstars[phiInd,vInd] = fstar_abs
                    Aabs = max(Iabs)
                    fabsLabel = '$f^*_{{res}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_abs,3))
                    axIss.plot(fstar_abs,Aabs,'*',markersize=10,label=fabsLabel)
                    #axIss.legend(loc='best')
                    #axIss.annotate(fabsLabel, xy=(fstar_abs,Aabs), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})
            axIss.set_xscale('log')
            axIss.set_xlabel('$f$ $\mathrm{[Hz]}$')
            axIss.set_ylabel('$|A|_{ss}$ $\mathrm{[nA]}$')
            if addTitles:
                axIss.set_title('$\mathrm{|Amplitude|_{ss}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
            
            plt.tight_layout()
            
            self.fstars = fstars
            if len(phis)>1: # Multiple light amplitudes
                #for i, A0 in enumerate(self.A0):
                fstarAfig = plt.figure()
                for vInd, V in enumerate(self.Vs):
                    if self.A0[0] > 0: # A0[r]
                        plt.plot(np.array(phis)/self.A0[0],fstars[:,vInd])
                        plt.xlabel('$\mathrm{Modulating}\ \phi_1(t)/\phi_0(t)$')
                    else:
                        plt.plot(np.array(phis),fstars[:,vInd])
                        plt.xlabel('$\mathrm{Modulating}\ \phi_1(t)$')
                plt.xscale('log')
                #plt.xlabel('$\mathrm{Modulating}\ \phi_1(t)/\phi_0(t)\ \mathrm{[photons \cdot s^{-1} \cdot mm^{-2}]}$')
                #plt.xlabel('$\mathrm{Modulating}\ \phi_1(t)/\phi_0(t)$')
                plt.ylabel('$f^*\ \mathrm{[Hz]}$')
                if addTitles:
                    plt.title('$f^*\ vs.\ \phi_1(t).\ \mathrm{{Background\ illumination:}}\ \phi_0(t)={:.3g}$'.format(self.A0[0]))
            
        if protocol == "varyPL":
            # Freeze axis limits
            ymin, ymax = plt.ylim()
            plt.ylim(ymin, ymax)
            pos = 0.02 * abs(ymax-ymin)
            #plt.ylim(ax.get_ylim())
            for run in range(nRuns): 
                if nRuns > 1:
                    onD = pulseCycles[run,0]
                    offD = pulseCycles[run,1]
                    padD = padDs[run]
                # Loop over light intensity...
                for phiInd, phiOn in enumerate(phis): #for phiInd in range(0, len(phis)):
                    # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                    for vInd, V in enumerate(Vs): #range(0, len(Vs)):
                        col, style = self.getLineProps(run, vInd, phiInd)
                        peakInds = IpInds[run][phiInd][vInd]
                        axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)
        
        
                        #peakInds = IpInds[run][phiInd][vInd]
                        #axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)
            ### Plot figure to show time of Ipeak vs time of light off c.f. Nikolic et al. 2009 Fig 2b
        #     axLag = Ifig.add_subplot(gsPL[0,-1])
        #     tpeaks = [(t[IpInds[p][0][0]] - delD) for p in range(nRuns)]
        #     toffs = [(t[PulseInds[p][0][0][0][1]-1] - delD) for p in range(nRuns)] # pDs
        #     plot(toffs,tpeaks,marker='x',linestyle='')
            
        #     axLag.axis('equal')
            tmax = max(pDs)*1.25
            axLag.plot([0,tmax], [0,tmax], ls="--", c=".3")
            axLag.set_xlim(0,tmax)
            axLag.set_ylim(0,tmax)
            axLag.set_ylabel('$\mathrm{Time\ of\ peak\ [ms]}$')
        #     plt.tight_layout()
            axLag.set_aspect('auto')
        #     print(axLag.get_xlim())
        #     diag_line, = axLag.plot(axLag.get_xlim(), axLag.get_ylim(), ls="--", c=".3")
            
            
            ### Plot figure to show current peak vs time of light off c.f. Nikolic et al. 2009 Fig 2c
            axPeak.set_xlim(0,tmax)
            axPeak.set_xlabel('$\mathrm{Pulse\ duration\ [ms]}$')
            axPeak.set_ylabel('$\mathrm{Photocurrent\ peak\ [nA]}$')
            
            #plt.tight_layout()
        elif protocol == "varyIPI":
            # Freeze axis limits
            ymin, ymax = plt.ylim()
            plt.ylim(ymin, ymax)
            pos = 0.02 * abs(ymax-ymin)
            #plt.ylim(ax.get_ylim())
            for run in range(nRuns): 
                if nRuns > 1:
                    onD = pulseCycles[run,0]
                    offD = pulseCycles[run,1]
                    padD = padDs[run]
                # Loop over light intensity...
                for phiInd, phiOn in enumerate(phis): #for phiInd in range(0, len(phis)):
                    # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                    for vInd, V in enumerate(Vs): #range(0, len(Vs)):
                        col, style = self.getLineProps(run, vInd, phiInd)
                        #plt.figure(Ifig.number)
                        #plt.annotate('', (delD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col})
                        plt.annotate('', (delD+onD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col,'shrinkA':0,'shrinkB':0})
            
            ### PLOT
            popt, _, _ = self.fitPeaks(self.tpIPI, self.IpIPI, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$',Ifig) # tPeaks peaks
            if verbose > 0:
                print("tau_R = {} ==> rate_R = {}".format(popt[1],1/popt[1]))
            
        #plt.show()
        plt.tight_layout()
        
        if label:
            Ifig.savefig(fDir+protocol+dataTag+"."+saveFigFormat, bbox_extra_artists=(lgd,), bbox_inches='tight', format=saveFigFormat) # Use this to save figures when legend is beside the plot
        else:
            Ifig.savefig(fDir+protocol+dataTag+"."+saveFigFormat, format=saveFigFormat)
        
        return Ifig.number

            
class protCustom(Protocol):
    # Class attributes
    protocol = 'custom'
    squarePulse = False
    def __init__(self, params=protParams['custom'], saveData=True): #ProtParamsCustom #phis=[1e14,1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[10.,160.]], totT=200., dt=0.1): # , nRuns=1
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #self.phis = phis
        #self.Vs = Vs
        #if isinstance(pulses, (np.ndarray)): # , np.generic
        #    self.pulses = pulses
        #else:
        #self.pulses = np.array(pulses)
        #self.totT = totT
        #self.prepare() # Called in setParams() and run()
        
        #self.dt=dt

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        #self.pulseInds = np.array([[np.searchsorted(self.t, pulses[p,time]) for time in range(2)] for p in range(self.nPulses)])
        #pulses = np.array([[delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD] for p in range(nPulses)]) #np.array([[delD,onD]])
        self.nRuns = 1 #nRuns ### Reconsider this...
        self.phis.sort()
        self.Vs.sort(reverse=True)
        
        
class protStep(Protocol):
    # Heaviside Pulse
    protocol = 'step'
    squarePulse = True
    nRuns = 1
    def __init__(self, params=protParams['step'], saveData=True): #ProtParamsStep #phis=[1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[50.,200.]], totT=300., dt=0.1): # , nRuns=1
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #self.phis = phis
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.totT = totT
        #self.dt=dt
        # pass
        # delD = 25.0  # Delay before on phase [ms]
        # onD = 250.0  # Duration of on phase [ms]
        # offD = 0.0   # Duration of off phase [ms]
        # padD = 0.0   # Duration of padding after last off phase [ms]
        # onSt = delD
        # offSt = delD+onD
        # totT = delD+nPulses*(onD+offD)  # Total simulation time [ms]
        # nRuns = 1
        # dt = 0.1
        #self.prepare() # Called in setParams() and run()
        

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)        
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        self.nRuns = 1 #nRuns
        #phi_t = InterpolatedUnivariateSpline([start,delD-dt,delD,end], [0,0,A,A],k=1) # Heaviside
        self.phis.sort()
        self.Vs.sort(reverse=True)
   
   
   # def phi_t(t):
        # for row in self.nPulses:
            # if t > self.pulses[0] and t < self.pulses[1]:
                # return 1
            # else:
                # return 0
        
class protSinusoid(Protocol):
    protocol = 'sinusoid'
    squarePulse = False
    def __init__(self, params=protParams['sinusoid'], saveData=True): #ProtParamsSinusoid #phis=[1e14], A0=[1e12], Vs=[-70], fs=np.logspace(-1,3,num=9), pulses=[[50.,550.]], totT=600., dt=0.1):
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #self.phis = phis
        #self.A0 = A0 # Background illumination
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.prepare() # Called in setParams() and run()

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        
        #self.totT = totT
        #self.dt=dt
        self.fs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        self.ws = 2 * np.pi * self.fs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        #self.sr = min([(1000)/(10*max(self.fs)), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.sr = max([(10)*max(self.fs), 1000/self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.dt = 1000/self.sr
        self.nRuns = len(self.ws)
        self.pulseCycles=np.column_stack((self.onDs,self.offDs))
        #self.pulseCycles=np.tile(np.column_stack((self.onDs,self.offDs)),(self.nRuns,1))
        self.padDs = np.zeros(self.nRuns)
        
        #ws = 2 * np.pi * np.logspace(-4,10,num=7) # Frequencies [rads/s]
        #self.nRuns = len(freqs)
        if (1000)/min(self.fs) > min(self.onDs):
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
            #print('Warning: The period of the lowest frequency is longer than the total simulation time!')
        
        #figTitle = "Photocurrent through time "
        #self.phis.sort(reverse=True)
        self.phis.sort()
        self.Vs.sort(reverse=True)

        #self.fs.sort()
        #self.ws.sort()

class protDualTone(Protocol):
    protocol = 'dualTone'
    squarePulse = False
    # Change default parameter key to 'dualTone'!!!
    def __init__(self, params=protParams['custom'], saveData=True): #ProtParamsSinusoid #phis=[1e14], A0=[1e12], Vs=[-70], fs=np.logspace(-1,3,num=9), pulses=[[50.,550.]], totT=600., dt=0.1):
        
        self.saveData = saveData
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #self.phis = phis
        #self.A0 = A0 # Background illumination
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.prepare() # Called in setParams() and run()
        #fA
        #fB

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        
        #self.totT = totT
        #self.dt=dt
        self.fAs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        self.fBs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        self.wAs = 2 * np.pi * self.fAs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        self.wBs = 2 * np.pi * self.fBs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        #self.sr = min([(1000)/(10*max(self.fAs,self.fBs)), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.sr = max([(10)*max(self.fAs,self.fBs), 1000/self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.dt = 1000/self.sr
        for fA,fB in itertools.product(self.fAs,self.fBs):
            print(fA+fB)
        self.nRuns = len(self.ws) # Modify...
        self.pulseCycles=np.column_stack((self.onDs,self.offDs))
        #self.pulseCycles=np.tile(np.column_stack((self.onDs,self.offDs)),(self.nRuns,1))
        self.padDs = np.zeros(self.nRuns)
        
        if (1000)/min(self.fs) > min(self.onDs):
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
            #print('Warning: The period of the lowest frequency is longer than the total simulation time!')

        self.phis.sort()
        self.Vs.sort(reverse=True)


class protChirp(Protocol):
    protocol = 'chirp'
    squarePulse = False
    def __init__(self, params=protParams['chirp'], saveData=True): #ProtParamsSinusoid #phis=[1e14], A0=[1e12], Vs=[-70], fs=np.logspace(-1,3,num=9), pulses=[[50.,550.]], totT=600., dt=0.1):
        
        self.saveData = saveData
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #self.phis = phis
        #self.A0 = A0 # Background illumination
        #self.Vs = Vs
        #self.f0 = f0
        #self.fgrad = Hz/s
        #self.prepare() # Called in setParams() and run()
        
    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        
        #self.fs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        #self.ws = 2 * np.pi * self.fs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        #self.sr = min([1000/(10*max(self.fs)), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        #self.sr = min([(1000)/(100*self.fT), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.sr = max([(10)*max(self.f0,self.fT), 1000/self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.dt = 1000/self.sr
        #print(self.sr,self.dt)
        self.nRuns = 1 #len(self.ws)
        self.pulseCycles=np.column_stack((self.onDs,self.offDs))
        self.padDs = np.zeros(self.nRuns)
        
        #ws = 2 * np.pi * np.logspace(-4,10,num=7) # Frequencies [rads/s]
        #self.nRuns = len(freqs)
        if (1000)/self.f0 > min(self.onDs): #1/10**self.fs[0] > self.totT:
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
            #print('Warning: The period of the lowest frequency is longer than the total simulation time!')
        self.phis.sort()
        self.Vs.sort(reverse=True)
        


class protRamp(Protocol):
    protocol = 'ramp'
    squarePulse = False
    nRuns = 1
    def __init__(self, params=protParams['ramp'], saveData=True): # ProtParamsRamp #phis=[1e14,1e15,1e16,1e17,1e18], phi_ton = 0, Vs=[-70], pulses=[[25.,275.]], totT=300., dt=0.1): # , nRuns=1
        """Linearly increasing pulse"""
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #phi_t = InterpolatedUnivariateSpline([start,delD,end], [0,0,A],k=1)
        #self.phis = phis
        #self.phi_ton = phi_ton
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.totT = totT
        #self.dt = dt
        #self.prepare() # Called in setParams() and run()

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        self.nRuns = 1 #nRuns # Make len(phi_ton)?
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        self.pulseCycles=np.column_stack((self.onDs,self.offDs))
        #self.pulseCycles=np.tile(np.column_stack((self.onDs,self.offDs)),(self.nRuns,1))
        self.padDs = np.zeros(self.nRuns)
        self.phis.sort()
        self.Vs.sort(reverse=True)

class protSaturate(Protocol): 
    # One very short, saturation intensity pulse e.g. 10 ns @ 100 mW*mm^-2 for wild type ChR
    # Used to calculate gbar, assuming that O(1)-->1 as onD-->0 and phi-->inf
    protocol = 'saturate'
    squarePulse = True
    nRuns = 1
    def __init__(self, params=protParams['saturate'], saveData=True): # ProtParamsSaturate #phis=[irrad2flux(1000,470)], Vs=[-70], pulses=[[5.,5+1e-3]], totT=20., dt=1e-3): # delD=5, dt=1e-3, totT=20 # , nRuns=1
        if verbose > 0:
            print("Running saturation protocol to find the maximum conductance (bar{g})")
        
        #phis = [irrad2flux(1000,470)] # 100 mW*mm^-2
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = True #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        #self.plotStateVars = True
        self.setParams(params)
        #self.phis = phis
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.prepare() # Called in setParams() and run()

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        #self.totT = totT
        self.nRuns = 1 #nRuns
        #self.dt=dt
        if any(tp < self.dt for tp in self.onDs):
            warnings.warn('Warning: Time step is too large for the pulse width [pulse:{}]!'.format(p))
        #for p, t in enumerate(self.onDs):
        #    if self.onDs[p] < self.dt:
        #        warnings.warn('Warning: Time step is too large for the pulse width [pulse:{}; t={}]!'.format(p,t))
        self.phis.sort()
        self.Vs.sort(reverse=True)

class protInwardRect(Protocol):
    # Iss vs Vclamp
    protocol = 'inwardRect'
    squarePulse = True
    nRuns = 1
    def __init__(self, params=protParams['inwardRect'], saveData=True): # ProtParamsInwardRect #phis=[irrad2flux(1,470),irrad2flux(10,470)], Vs=[-100,-80,-60,-40,-20,0,20,40,60,80], pulses=[[50.,300.]], totT=400., dt=0.1): # , nRuns=1
        # Used to calculate v0 and v1
        
        if verbose > 0:
            print("Running inward rectification protocol to parameterise f(V)")
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False #plotStateVars
        self.plotKinetics = False #True
        
        self.setParams(params)
        #self.phis = phis
        #self.Vs = Vs
        #self.pulses = np.array(pulses)
        #self.totT = totT
        #self.dt=dt
        #self.prepare() # Called in setParams() and run()
    
    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pulses = np.asarray(self.pulses)
        self.nPulses = self.pulses.shape[0]
        self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        self.nRuns = 1 #nRuns
        self.phis.sort()
        self.Vs.sort(reverse=True)
        

class protVaryPL(Protocol):
    # Vary pulse length - See Nikolic+++++2009, Fig. 2 & 9
    #def __init__(self, phis=[1e12], Vs=[-70], delD=25, pulses=[[1,74],[2,73],[3,72],[5,70],[8,67],[10,65],[20,55]], nRuns=1, dt=0.1):
    protocol = 'varyPL'
    squarePulse = True
    nPulses = 1 #Fixed at 1
    def __init__(self, params=protParams['varyPL'], saveData=True): # ProtParamsVaryPL # phis=[1e12], Vs=[-70], delD=25, pDs=[1,2,3,5,8,10,20], totT=100, dt=0.1): #nPulses=1,
        if verbose > 0:
            print("Running variable pulse length protocol")
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False
        self.plotStateVars = True #plotStateVars
        self.plotKinetics = False #plotKinetics
        
        self.setParams(params)
        #delD = 25
        #IPI = 75 # Inter-Pulse-Interval
        #pDs = [1,2,3,5,8,10,20]
        #phis = [1e12]#[irrad2flux(0.65, 470)]#[1e50]#
        #Vs = [-70] # Could relax this?
        #delD = 25
        #nPulses = 1
        #totT = delD+IPI#delD+nPulses*(onD+offD)  # Total simulation time per run [ms]
        #self.totT = totT
        #self.phis = phis
        #self.Vs = Vs
        #self.pulses, _ = cycles2times(self.pulseCycles,self.delD)
        #print(self.pulses)
        #self.pulseCycles=np.column_stack((pDs,[IPI-pD for pD in pDs])) # [:,0] = on phase duration; [:,1] = off phase duration
        #self.nPulses = nPulses 
        #self.dt = dt
        #self.prepare() # Called in setParams() and run()

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.pDs = np.sort(np.array(self.pDs))
        self.nRuns = len(self.pDs)
        self.delDs = np.ones(self.nRuns)*self.delD
        self.onDs = self.pDs
        self.offDs = (np.ones(self.nRuns)*self.totT) - self.delDs - self.onDs
        self.pulseCycles = np.column_stack((self.onDs,self.offDs))
        self.padDs = np.zeros(self.nRuns)
        self.phis.sort()
        self.Vs.sort(reverse=True)
        
class protVaryIPI(Protocol):
    # Vary Inter-Pulse-Interval
    protocol = 'varyIPI'
    squarePulse = True
    nPulses = 2 # Fixed at 2 for this protocol
    def __init__(self, params=protParams['varyIPI'], saveData=True): # ProtParamsVaryIPI #phis=[1e14], Vs=[-70], delD=100, onD=200, IPIs=[500,1000,1500,2500,5000,7500,10000], dt=0.1): #nPulses=2, 
        if verbose > 0:
            print("Running S1-S2 protocol to solve for tau_R")
        
        self.saveData = saveData
        ###self.dataTag = str(RhO.nStates)+"s"
        #self.plotResults = plotResults
        self.plotPeakRecovery = False #plotPeakRecovery
        self.plotStateVars = False # plotStateVars
        self.plotKinetics = False #plotKinetics    
        
        self.setParams(params)
        # delD = 100
        # IPIs = [500,1000,1500,2500,5000,7500,10000]#[10,35,55,105,155] #IPIs = [10,20,30,50,80,100]#,200] #[55,105,155,305,1005]#
        # onD = 200                # ???
        # phis = [1e14]#[irrad2flux(0.65, 470)]#[irrad2flux(10, 470)]#
        # Vs = [-70]#,-50] # Could relax this?
        #self.phis = phis
        #self.Vs = Vs
        #self.nPulses = 2 # Fixed at 2 for this protocol
        #self.dt = dt
        # Make self. ?
        #self.plotPeakRecovery = True
        #plotStateVars = True
        #self.prepare() # Called in setParams() and run()

    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.IPIs = np.sort(np.asarray(self.IPIs)) #np.sort(np.array(IPIs))
        
        self.nRuns = len(self.IPIs)
        self.delDs = np.ones(self.nRuns)*self.delD
        self.onDs = np.ones(self.nRuns)*self.onD
        self.offDs = self.IPIs
        ##pulseCycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI-onD for IPI in IPIs])) # [:,0] = on phase duration; [:,1] = off phase duration
        #pulseCycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI for IPI in IPIs]))
        self.pulseCycles = np.column_stack((self.onDs,self.offDs))
        
        self.pulses, _ = cycles2times(self.pulseCycles,self.delD)
        
        #padDs = [totT-((onD+pOff)*nPulses)-delD for pOff in pulseCycles[:,1]] # Not necessary
        self.padDs = np.zeros(self.nRuns); self.padDs[-1] = -0.8*max(self.pulseCycles[:,1])
        
        self.totT = self.delD+self.nPulses*max(self.IPIs) -0.8*max(self.pulseCycles[:,1]) # Total simulation time per run [ms] ### Should IPI be one pulse cycle or time between end/start of two pulses?
        self.IpIPI = np.zeros(self.nRuns) #peaks
        self.tpIPI = np.zeros(self.nRuns) #tPeaks
        self.phis.sort()
        self.Vs.sort(reverse=True)



    



protocols = {'custom': protCustom, 'step': protStep, 'saturate': protSaturate, 'inwardRect': protInwardRect, 'varyPL': protVaryPL, 'varyIPI': protVaryIPI, 'sinusoid': protSinusoid, 'chirp': protChirp, 'ramp': protRamp}
# E.g. 
# protocols['varyPL']([1e12], [-70], 25, [1,2,3,5,8,10,20], 100, 0.1)

#squarePulses = [protocol for protocol in protocols if protocol.squarePulse]
#arbitraryPulses = [protocol for protocol in protocols if not protocol.squarePulse]
#squarePulses = {'custom': True, 'saturate': True, 'step': True, 'inwardRect': True, 'varyPL': True, 'varyIPI': True}
#arbitraryPulses = {'custom': True, 'sinusoid': True, 'chirp': True, 'ramp':True} # Move custom here
#smallSignalAnalysis = {'sinusoid': True, 'step': True, 'saturate': True} 
    
    
    
    
    
    
    
    
    
    
### Previous code ###    
    
#saveData = True
#dataTag = str(nStates)+"s"

# Solver parameters
#dt = 0.1 #0.01 #0.001
#StepsPerMS = 1
#extOrder = 5 #int(round(totT/(2*dt)))

# Plot settings
#plotPeakRecovery = False
#plotStateVars = False
#plotKinetics = False

# Initial guesses for kinetics curve fitting (tuples to be immutable)
#p0on = (-0.1,2,-1)
#p0inact = (-0.5,25,-0.5)
#p0off = (-0.1,7.5,-0.1,35,-0.1)
#p0fV = (40,4,RhO.E)
#p0IPI = (-1e-8,400,-1e-7)

