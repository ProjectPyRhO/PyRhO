import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl # for tick locators
from scipy.interpolate import InterpolatedUnivariateSpline #*
from scipy.optimize import curve_fit
from pyrho.parameters import *
from pyrho.utilities import * # times2cycles, cycles2times, plotLight, round_sig, expDecay, biExpDecay, findPeaks
from pyrho.loadData import * #import loadData
from pyrho.fitting import * #calcIssfromfV and fitPeaks...
from pyrho.models import *
from pyrho.simulators import * # For characterise()
from pyrho.config import * #verbose, saveFigFormat, eqSize, addTitles, addStimulus, colours, styles, dDir, fDir
import pickle
import warnings
import os


class Protocol(PyRhOobject): #object
    """Common base class for all protocols"""
    
    def __init__(self, params=None, saveData=True):
        if params is None:
            params = protParams[self.protocol]
        self.RhO = None
        self.dataTag = ""
        self.saveData = saveData
        self.plotPeakRecovery = False
        self.plotStateVars = False
        self.plotKinetics = False
        self.setParams(params)
        self.prepare()
        self.begT, self.endT = 0, self.totT
        self.phi_ts = None
        self.lam = 470 # Default wavelength [nm]
    
    def __str__(self):
        return self.protocol
    
    def __repr__(self):
        return "<PyRhO {} Protocol object (nRuns={}, nPhis={}, nVs={})>".format(self.protocol, self.nRuns, self.nPhis, self.nVs)
    
    def __iter__(self):
        """Iterator to return the pulse sequence for the next trial"""
        self.run = 0
        self.phiInd = 0
        self.vInd = 0
        return self
        
    def __next__(self):
        """Iterator to return the pulse sequence for the next trial"""
        self.run += 1
        if self.run > self.nRuns:
            raise StopIteration
        return self.getRunCycles(self.run - 1) #self.pulses#[self.run]
    
    '''
    def setCycles(self, cycles, delD=None):
        if delD is None:
            delD = self.delD
        else:
            self.delD = delD
        self.cycles = np.asarray(cycles, copy=True)
        #self.cycles = np.asarray([[self.onD, self.totT-self.delD]])
        self.nPulses = self.cycles.shape[0]
        self.pulses, self.totT = cycles2times(self.cycles, self.delD)
        self.delDs = np.array([row[0] for row in self.pulses], copy=True) # pulses[:,0]    # Delay Durations
        #self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        #self.offDs = np.append(self.pulses[1:,0], self.totT) - self.pulses[:,1]
        self.onDs = np.array(self.cycles[:,0])
        self.offDs = np.array(self.cycles[:,1])
    '''
    
    def prepare(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
        
        if np.isscalar(self.cycles): #self.cycles.shape[1] == 1: # Only on duration specified
            onD = self.cycles
            if hasattr(self, 'totT'):
                offD = self.totT - onD - self.delD
            else:
                offD = 0
            self.cycles = np.asarray([[onD, offD]]) #self.cycles[0]
        
        self.cycles = np.asarray(self.cycles)
        self.nPulses = self.cycles.shape[0]
        self.pulses, self.totT = cycles2times(self.cycles, self.delD)
        self.delDs = np.array([pulse[0] for pulse in self.pulses], copy=True) # pulses[:,0]    # Delay Durations #self.delDs = np.array([self.delD] * self.nRuns)
        self.onDs = np.array(self.cycles[:,0])  #self.onDs = np.array([cycle[0] for cycle in self.cycles])
        self.offDs = np.array(self.cycles[:,1]) #self.offDs = np.array([cycle[1] for cycle in self.cycles])
        
        if np.isscalar(self.phis):
            self.phis = [self.phis] #np.asarray([self.phis])
        self.phis.sort(reverse=True)
        self.nPhis = len(self.phis)
        
        if np.isscalar(self.Vs):
            self.Vs = [self.Vs] #np.asarray([self.Vs])
        self.Vs.sort(reverse=True)
        self.nVs = len(self.Vs)
        
        self.extraPrep()
        return
    
    def extraPrep(self):
        pass
    
    def genContainer(self):
        return [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
    
    def getShortestPeriod(self):
        return np.amin(self.cycles[self.cycles.nonzero()]) #min(self.delD, min(min(self.cycles)))
    
    def finish(self, PC, RhO):
        pass
    
    def getRunCycles(self, run):
        #return times2cycles(self.pulses, self.totT)
        return (self.cycles, self.delD) #self.delDs[run]
    
    '''    
    def run(self, Sim, RhO, verbose=verbose): 
        """Main routine to run the simulation protocol"""
        
        t0 = wallTime() #time.perf_counter()
        self.prepare()
        self.dt = Sim.prepare(self.getShortestPeriod())
        
        if verbose > 0:
            print("\n================================================================================")
            print("Running '{}' protocol with {} for the {}-state model... ".format(self.protocol, Sim.simulator, RhO.nStates))
            print("================================================================================\n")
            print("{{nRuns={}, nPhis={}, nVs={}}}".format(self.nRuns, self.nPhis, self.nVs))
        
        
        #Sim.prepare(min(self.onDs)/10)
        self.RhO = RhO
        self.Sim = Sim

        self.PD = ProtocolData(self.protocol, self.nRuns, self.phis, self.Vs)
        self.PD.peak_ = [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
        self.PD.ss_ = [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
        if hasattr(self, 'runLabels'):
            self.PD.runLabels = self.runLabels
        
        if verbose > 1: 
            self.printParams()
        
        # Loop over the number of runs...             ### Place within V & phi loops to test protocols at different V & phi?
        for run in range(self.nRuns):
            
            cycles, delD = self.getRunCycles(run)
            pulses, totT = cycles2times(cycles, delD)
            
            # Loop over light intensity...
            for phiInd, phiOn in enumerate(self.phis): 
                
                if verbose > 1 and (self.nPhis > 1 or (run == 0 and phiInd == 0)): # len(phis)>0
                    #print(RhO.phi); print(type(RhO.phi))
                    RhO.dispRates()
                
                # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                for vInd, V in enumerate(self.Vs): 
                    
                    if self.squarePulse: #protocol in squarePulses: ##### Change after changing 'custom'
                        #phi_t = InterpolatedUnivariateSpline([pStart,pEnd],[phiOn,phiOn], k=1, ext=1) 
                        #I_RhO,t,soln = Sim.runTrial(RhO, self.nPulses, V,phiOn,delD,onD,offD,padD,dt,verbose)
                        #cycles, delD = self.getRunCycles(run)
                        #cycles, delD = times2cycles(self.pulses,self.totT)
                        I_RhO, t, soln = Sim.runTrial(RhO, phiOn, V, delD, cycles, self.dt, verbose) #self.totT, 
                        
                    else: # Arbitrary functions of time: phi(t)
                        
                        #phi_t = self.genPulse(run, phiOn, delD, onD)
                        phi_ts = self.phi_ts[run][phiInd][:]
                        #I_RhO,t,soln = Sim.runTrialPhi_t(RhO, V, phi_t, delD, onD, self.totT, dt, verbose)
                        #runTrialPhi_t_new(self, RhO, phi_ts, V, delD, cycles, endT, dt, verbose=verbose):
                        I_RhO, t, soln = Sim.runTrialPhi_t(RhO, phi_ts, V, delD, cycles, self.totT, self.dt, verbose)
                        
                    # Save simulation results
                        # Change: Prot.pulses[run] := [[t_on1, t_off1],...]
                    #pulses = np.array([[delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD] for p in range(self.nPulses)]) #np.array([[delD,onD]])
                    
                    PC = PhotoCurrent(I_RhO, t, pulses, phiOn, V, self.protocol)
                    #PC.alignToTime()

                    PC.states = soln
                    #self.data[run][phiInd][vInd] = PC
                    #PD.trials.append(PC)
                    self.PD.trials[run][phiInd][vInd] = PC
                    self.PD.peak_[run][phiInd][vInd] = PC.peak_
                    self.PD.ss_[run][phiInd][vInd] = PC.ss_
                    
                    # self.labels[run][phiInd][vInd] = label
                    # label = ""
                    if verbose > 1:
                        #print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]; label=<{}>'.format(run,nRuns,phiInd,len(phis),vInd,len(Vs),PC.range_[0],PC.range_[1],label))
                        print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]'.format(run,self.nRuns,phiInd,len(self.phis),vInd,len(self.Vs),PC.range_[0],PC.range_[1]))
        
        self.finish(PC, RhO) # Sim
        
        if self.saveData:
            #from os import path
            self.dataTag = str(RhO.nStates)+"s"
            #pklFile = path.join(dDir, self.protocol+self.dataTag+".pkl")
            #fh = open(pklFile,"wb")
            #pickle.dump(self.PD, fh)
            #fh.close()
            saveData(self.PD, self.protocol+self.dataTag)
            #if verbose > 0:
            #    print("Protocol data saved to disk: {}".format(pklFile))
        
        self.runTime = wallTime() - t0 #time.perf_counter() - t0
        if verbose > 0:
            print("\nFinished '{}' protocol with {} for the {}-state model in {:.3g}s".format(self.protocol, Sim.simulator, RhO.nStates, self.runTime))
            print("--------------------------------------------------------------------------------\n")
            
        return self.PD
    '''
    
    def genPulseSet(self, genPulse=None):
        """Function to generate a set of spline functions to phi(t) simulations"""
        if genPulse is None: # Default to square pulse generator
            genPulse = self.genPulse
        phi_ts = [[[None for pulse in range(self.nPulses)] for phi in range(self.nPhis)] for run in range(self.nRuns)]
        for run in range(self.nRuns):
            cycles, delD = self.getRunCycles(run)
            pulses, totT = cycles2times(cycles, delD)
            for phiInd, phi in enumerate(self.phis):
                for pInd, pulse in enumerate(pulses):
                #for pulse, onD in enumerate(self.onDs):
                    #onD, offD = cycles[pulse]
                    phi_ts[run][phiInd][pInd] = genPulse(run, phi, pulse) #(run, phiOn, self.delDs[pulse], onD)
        self.phi_ts = phi_ts
        return phi_ts
    
    def genPulse(self, run, phi, pulse):
        """Default interpolation function for square pulses"""
        pStart, pEnd = pulse
        phi_t = InterpolatedUnivariateSpline([pStart,pEnd], [phi,phi], k=1, ext=1)
        return phi_t
    
    
    ### Move fitPeaks and fitfV to fitting.py ###
    def fitPeaks(self, t_peaks, I_peaks, curveFunc, p0, eqString, fig=None):
        #print(p0)
        shift = t_peaks[0] # ~ delD
    #     if protocol == 'recovery':
    #         plt.ylim(ax.get_ylim()) # Prevent automatic rescaling of y-axis
        popt, pcov = curve_fit(curveFunc, t_peaks-shift, I_peaks, p0=p0) #Needs ball-park guesses (0.3, 125, 0.5)
        peakEq = eqString.format(*[round_sig(p,3) for p in popt]) # *popt rounded to 3s.f.
        
        if fig:
            plt.figure(fig.number) # Select figure
    #     ext = 10 # Extend for ext ms either side
    #     xspan = t_peaks[-1] - t_peaks[0] + 2*ext 
    #     xfit=np.linspace(t_peaks[0]-ext-shift,t_peaks[-1]+ext-shift,xspan/dt)
            plt.plot(t_peaks, I_peaks, linestyle='', color='r', marker='*')
            #xfit=np.linspace(-shift,self.totT-shift,self.totT/self.dt) #totT
            nPoints = 10*int(round(self.totT-shift/self.dt))+1 # 1001
            xfit=np.linspace(-shift, self.totT-shift, nPoints) #totT
            yfit=curveFunc(xfit,*popt)
            
            plt.plot(xfit+shift,yfit,linestyle=':',color='#aaaaaa',linewidth=1.5*mp.rcParams['lines.linewidth'])#,label="$v={:+} \mathrm{{mV}}$, $\phi={:.3g}$".format(V,phiOn)) # color='#aaaaaa' 
            #ylower = copysign(1.0,I_peaks.min())*ceil(abs((I_peaks.min()*10**ceil(abs(log10(abs(I_peaks.min())))))))/10**ceil(abs(log10(abs(I_peaks.min()))))
            #yupper = copysign(1.0,I_peaks.max())*ceil(abs((I_peaks.max()*10**ceil(abs(log10(abs(I_peaks.max())))))))/10**ceil(abs(log10(abs(I_peaks.max()))))
        #     if (len(Vs) == 1) and (len(phis) == 1) and (nRuns == 1):
        #         x, y = 0.8, 0.9
        #     else:
            x = 0.8
            y = yfit[-1] #popt[2]
            
            plt.text(x*self.totT,y,peakEq,ha='center',va='bottom',fontsize=config.eqSize) #, transform=ax.transAxes)
        
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
    
    # def fitfV(self, Vs, Iss, curveFunc, p0, RhO, fig=None):#, eqString): =plt.gcf()
        # if fig==None:
            # fig=plt.gcf()
        # markerSize=40
        # eqString = r'$f(V) = \frac{{{v1:.3}}}{{V-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{V-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'
        # psi = RhO.calcPsi(RhO.steadyStates)
        ##sf = RhO.A * RhO.gbar * psi * 1e-6 # Six-state only
        # sf = RhO.g * psi * 1e-6 
        # fVs = np.asarray(Iss)/sf # np.asarray is not needed for the six-state model!!!
        # popt, pcov = curve_fit(curveFunc, Vs, fVs, p0=p0) # (curveFunc, Vs, Iss, p0=p0)
        # pFit = [round_sig(p,3) for p in popt]
        ##peakEq = eqString.format(pFit[0],pFit[2],pFit[2],pFit[1])
        # peakEq = eqString.format(v1=pFit[0],E=pFit[2],v0=pFit[1])
        
        # Vrange = max(Vs)-min(Vs)
        # xfit=np.linspace(min(Vs),max(Vs),Vrange/.1) #Prot.dt
        # yfit=curveFunc(xfit,*popt)*sf
        
        ##peakEq = eqString.format(*[round_sig(p,3) for p in popt])
        
        # fig.plot(xfit,yfit)#,label=peakEq)#,linestyle=':', color='#aaaaaa')
        ##col, = getLineProps(Prot, 0, 0, 0) #Prot, run, vInd, phiInd
        ##plt.plot(Vs,Iss,linestyle='',marker='x',color=col)
        # fig.scatter(Vs,Iss,marker='x',color=colours,s=markerSize)#,linestyle=''
        
        ##x = 1 #0.8*max(Vs)
        ##y = 1.2*yfit[-1]#max(IssVals[run][phiInd][:])
        ##plt.text(-0.8*min(Vs),y,peakEq,ha='right',va='bottom',fontsize=eqSize)#,transform=ax.transAxes)
        
        # if verbose > 1:
            # print(peakEq)
        # return popt, pcov, peakEq
    
    
    



    def plot(self, plotStateVars=False):
        
        Ifig = plt.figure() #plt.figure(figsize=(figWidth, figHeight))
        self.createLayout(Ifig)
        #self.genLabels()
        self.PD.plot(self.axI) #self.labels... #baseline, = axI.plot(t, I_RhO, color=col, linestyle=style, label=label)
        ### Add legend    
        # if label:
            # if protocol == 'custom' or protocol == 'step' or protocol == 'rectifier':
                # if len(self.Vs) == 1:
                    # ncol = 1
                # else:
                    # ncol = len(self.phis)
                # lgd = plt.legend(loc='best', borderaxespad=0, ncol=ncol, fancybox=True) #, shadow=True , bbox_to_anchor=(1.02, 1)
            # else:
                # lgd = plt.legend(loc='best')
        #protPulses = self.getProtPulses()
        
        # if self.squarePulse:
            # self.addStimulus = False
        # if self.addStimulus:
            # plotLight(protPulses, self.axI)
            # for run in range(self.nRuns):
                # for phiInd in range(self.nPhis):
                    # self.plotStimulus(self.phi_ts[run][phiInd], protPulses, self.totT, self.axS, None) #if protocol == 'delta':
        # else: # Create a separate stimulus figure
            # for run in range(self.nRuns):
                # for phiInd in range(self.nPhis):
                    # stimFig = self.plotStimulus(self.phi_ts[run][phiInd], protPulses, self.totT, ax=None, light='spectral')
            # if max(self.phis) / min(self.phis) >= 100:
                # plt.yscale('log')
        
        self.addAnnotations()
        self.plotExtras()
        
        self.plotStateVars = plotStateVars
        
        #animateStates = True # https://jakevdp.github.io/blog/2013/05/28/a-simple-animation-the-magic-triangle/
        if self.plotStateVars: 
            RhO = self.RhO
            for run in range(self.nRuns):
                #cycles, delD = self.getRunCycles(run)
                #pulses, totT = cycles2times(cycles, delD)
                for phiInd, phi in enumerate(self.phis):
                    for vInd in range(self.nVs):
                        # write getPulseSeries() to return the correct set for shortPulses and IPIs?
                        # if self.protocol == 'shortPulse':
                            # pulses = np.array([[0,self.pDs[run]]])+self.delDs[run]
                        # else:
                            # pulses = self.pulses
                        pc = self.PD.trials[run][phiInd][vInd]
                        # soln = pc.states #multisoln[run][phiInd][vInd]
                        # t = pc.t
                        #RhO.plotStates(t,soln,pulses,RhO.labels,phiOn,IpInds[run][phiInd][vInd],'states{}s-{}-{}-{}'.format(RhO.nStates,run,phiInd,vInd))
                        fileName = '{}States{}s-{}-{}-{}'.format(self.protocol,RhO.nStates,run,phiInd,vInd)#; print(fileName)
                        RhO.plotStates(pc.t, pc.states, pc.pulses, RhO.stateLabels, phi, pc.peakInds_, fileName)
        
        plt.figure(Ifig.number)
        plt.sca(self.axI)
        self.axI.set_xlim(self.PD.begT, self.PD.endT)
        # if addTitles:
            # figTitle = self.genTitle()
            # plt.title(figTitle) #'Photocurrent through time'
            
        #plt.show()
        plt.tight_layout()
        
        externalLegend = False
        #from os import path
        figName = os.path.join(fDir, self.protocol+self.dataTag+"."+config.saveFigFormat)
        #plt.figure(Ifig.number)
        if externalLegend:
            Ifig.savefig(figName, bbox_extra_artists=(lgd,), bbox_inches='tight', format=config.saveFigFormat) # Use this to save figures when legend is beside the plot
        else:
            Ifig.savefig(figName, format=config.saveFigFormat)
        
        return #Ifig.number
    
    
    def genPlottingStimuli(self, genPulse=None, vInd=0):
        """Redraw stimulus functions in case data has been realigned"""
        if genPulse is None:
            genPulse = self.genPulse
        
        #self.addStimulus = config.addStimulus # Necessary?
        #if self.addStimulus: # Redraw stimulus functions in case data has been realigned
            #vInd = 0
            
            # # for delD in len(self.delDs):
                # # self.delDs -= self.PD.trials[run][phiInd][vInd]
        phi_ts = [[[None for pulse in range(self.nPulses)] for phi in range(self.nPhis)] for run in range(self.nRuns)]
        for run in range(self.nRuns):
            #cycles, delD = self.getRunCycles(run)
            #pulses, totT = cycles2times(cycles, delD)
            for phiInd, phi in enumerate(self.phis):
                pc = self.PD.trials[run][phiInd][vInd]
                # if pc.pulseAligned:
                for p, pulse in enumerate(pc.pulses):
                    phi_ts[run][phiInd][p] = genPulse(run, pc.phi, pulse)
        #self.phi_ts = self.genPulseSet()
        return phi_ts
        #else:
        #    return None
    
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()
        
        # Default layout
        self.axI = Ifig.add_subplot(111)
        plt.sca(self.axI)
        #plotLight(self.pulses, self.axI)
    
    '''
    def genLabels(self):
        
        self.labels = [[[None for v in range(self.nVs)] for p in range(self.nPhis)] for r in range(self.nRuns)]
        figTitle = "Photocurrent through time "
        if self.protocol == "shortPulse":
            figTitle += "for varying pulse length "
        elif self.protocol == "recovery":
            figTitle += "for varying inter-pulse-interval "
        else:
            figTitle += "\n "
        
        if len(self.phis) == 1:
            figTitle += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
        if len(self.Vs) == 1:
            figTitle += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V) # v=Vclamp
            
        for run in range(self.nRuns): 
            
            cycles, delD = self.getRunCycles(run)
            onD, offD = cycles[0]
            pulses, totT = cycles2times(cycles, delD)
            
            for phiInd, phiOn in enumerate(self.phis):

                for vInd, V in enumerate(self.Vs):
                    if verbose > 1:
                        print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}'.format(run,self.nRuns,phiInd,len(self.phis),vInd,len(self.Vs)))
                    
                    #label = self.labels[run][phiInd][vInd]
                    #label = ""
                    
                    if self.protocol == "shortPulse":
                        label = "$\mathrm{{Pulse}}={}\mathrm{{ms}}$ ".format(cycles[0][0]) #onD
                        #figTitle += "for varying pulse length "
                    elif self.protocol == "recovery":
                        label = "$\mathrm{{IPI}}={}\mathrm{{ms}}$ ".format(self.IPIs[run])
                        #figTitle += "for varying inter-pulse-interval "
                    elif self.protocol == 'sinusoid':
                        label = "$f={}\mathrm{{Hz}}$ ".format(round_sig(self.fs[run],3))
                    else:
                        label = ""                        
                        #figTitle += "\n "
                        
                    if len(self.phis) > 1:
                        label += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    #else:                    
                        #figTitle += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(phiOn)
                    
                    if len(self.Vs) > 1:
                        label += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V) # v=Vclamp
                    #else:
                        #figTitle += "$\mathrm{{V}} = {:+}\ \mathrm{{mV}}$ ".format(V)
                    
                    self.labels[run][phiInd][vInd] = label
                    
        return figTitle    
    
    
    def genTitle(self):
        figTitle = "Photocurrent through time "
        if self.protocol == "shortPulse":
            figTitle += "for varying pulse length \n"
        elif self.protocol == "recovery":
            figTitle += "for varying inter-pulse-interval "
        else:
            figTitle += "\n "
        
        if len(self.phis) == 1:
            figTitle += "$\phi = {:.3g}\ \mathrm{{photons \cdot s^{{-1}} \cdot mm^{{-2}}}}$ ".format(self.phis[0])
        if len(self.Vs) == 1:
            figTitle += "$\mathrm{{v}} = {:+}\ \mathrm{{mV}}$ ".format(self.Vs[0])
        
        return figTitle
    '''

    def getLineProps(self, run, vInd, phiInd):
        
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
    
        
    def plotExtras(self):
        pass
    
    def addAnnotations(self):
        pass
    
    
    def plotStimulus(self, phi_ts, begT, pulses, endT, ax=None, light='shade', col=None, style=None):
        
        nPulses = pulses.shape[0]
        assert(nPulses == len(phi_ts))
        #t = np.linspace(0,totT,10*int(round(totT/self.dt))+1) #10001) # 
        #t = np.linspace(begT,endT,10*int(round(endT-begT/self.dt))+1)
        nPoints = 10*int(round(endT-begT/self.dt))+1
        t = np.linspace(begT, endT, nPoints) #100001
        
        if ax == None:
            #fig = plt.figure()    
            ax = plt.gca()
        else:
            #plt.figure(fig.number)
            plt.sca(ax)
            
        if col is None:
            for p in range(nPulses):
                plt.plot(t, phi_ts[p](t))
        else:
            if style is None:
                style = '-'
            for p in range(nPulses):
                plt.plot(t, phi_ts[p](t), color=col, linestyle=style)

                
        if light == 'spectral':
            plotLight(pulses, ax=ax, light='spectral', lam=self.lam)
        else:
            plotLight(pulses, ax=ax, light=light)
        
        plt.xlabel('$\mathrm{Time\ [ms]}$') #(r'\textbf{Time} [ms]')
        #plt.xlim((0,totT))
        plt.xlim((begT,endT))
        #plt.ylabel('$\mathrm{\phi\ [photons \cdot mm^{-2} \cdot s^{-1}]}$')
        plt.ylabel('$\mathrm{\phi\ [ph./mm^{2}/s]}$')
        
        return ax
    
    
    def plotKinetics(self):
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
    



    

            
class protCustom(Protocol):
    
    # Class attributes
    protocol = 'custom'
    squarePulse = False
    #custPulseGenerator = None
    phi_ft = None
    
    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.nRuns = 1 #nRuns ### Reconsider this...
        
        #self.custPulseGenerator = self.phi_ft
        if not hasattr(self, 'phi_ts') or self.phi_ts == None:
            #self.phi_ts = self.genPulseSet()
            #self.genPulseSet(self.custPulseGenerator)
            self.genPulseSet(self.phi_ft)
        
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        
        if self.addStimulus: 
            phi_ts = self.genPlottingStimuli(self.phi_ft) #self.genPlottingStimuli(self.custPulseGenerator)
            
            gsStim = plt.GridSpec(4,1)
            self.axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
            self.axI = Ifig.add_subplot(gsStim[1:,:],sharex=self.axS) # Photocurrent axes
            pc = self.PD.trials[0][0][0]
            plotLight(pc.pulses, ax=self.axS, light='spectral', lam=470, alpha=0.2)
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    pc = self.PD.trials[run][phiInd][vInd]
                    col, style = self.getLineProps(run, vInd, phiInd)
                    self.plotStimulus(phi_ts[run][phiInd], pc.begT, self.pulses, pc.endT, self.axS, light=None, col=col, style=style) #light='spectral'
            plt.setp(self.axS.get_xticklabels(), visible=False)
            self.axS.set_xlabel('')
        else:
            self.axI = Ifig.add_subplot(111)
    
        
    
    
    def plotExtras(self):
        pass
        
class protStep(Protocol):
    # Heaviside Pulse
    protocol = 'step'
    squarePulse = True
    nRuns = 1
    

    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        
        self.nRuns = 1 #nRuns
        
        self.phi_ts = self.genPulseSet()
        #self.genPulseSet()
   
   # def phi_t(t):
        # for row in self.nPulses:
            # if t > self.pulses[0] and t < self.pulses[1]:
                # return 1
            # else:
                # return 0
                
    
    def addAnnotations(self):
        self.axI.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        self.axI.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        self.axI.grid(b=True, which='minor', axis='both', linewidth=.2)
        self.axI.grid(b=True, which='major', axis='both', linewidth=1)


class protSinusoid(Protocol):
    
    protocol = 'sinusoid'
    squarePulse = False
    
    def extraPrep(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
        
        self.fs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        self.ws = 2 * np.pi * self.fs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        self.sr = max(10000, int(round(10*max(self.fs)))) # Nyquist frequency - sampling rate (10*f) >= 2*f >= 10/ms
        #self.dt = 1000/self.sr # dt is set by simulator but used for plotting
        self.nRuns = len(self.ws)
        
        if (1000)/min(self.fs) > min(self.onDs):
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
        
        if isinstance(self.phi0, (int, float, complex)):
            self.phi0 = np.ones(self.nRuns) * self.phi0
        elif isinstance(self.phi0, (list, tuple, np.ndarray)):
            if len(self.phi0) != self.nRuns:
                self.phi0 = np.ones(self.nRuns) * self.phi0[0]
        else:
            warning.warn('Unexpected data type for phi0: ', type(self.phi0))
        
        assert(len(self.phi0) == self.nRuns)
        
        self.begT, self.endT = 0, self.totT
        self.phi_ts = self.genPulseSet()
        self.runLabels = ["$f={}\mathrm{{Hz}}$ ".format(round_sig(f,3)) for f in self.fs]
        
    def getShortestPeriod(self):
        return 1000/self.sr # dt [ms]
        
    def genPulse(self, run, phi, pulse):
        pStart, pEnd = pulse
        onD = pEnd - pStart
        t = np.linspace(0.0, onD, int(round((onD*self.sr/1000))+1), endpoint=True) # Create smooth series of time points to interpolate between
        if self.startOn: # Generalise to phase offset
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[run] + 0.5*phi*(1+np.cos(self.ws[run]*t)), ext=1, k=5) # phi0[r]
        else:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[run] + 0.5*phi*(1-np.cos(self.ws[run]*t)), ext=1, k=5) # phi0[r]
        
        return phi_t
    
    '''
    def getPhiFunc_orig(self, run, phiOn, delD, onD):
        pStart, pEnd = delD, (delD+onD)
        t = np.linspace(0.0, onD, (onD*self.sr/1000)+1, endpoint=True) # Create smooth series of time points to interpolate between
        if self.startOn:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[0] + 0.5*phiOn*(1+np.cos(self.ws[run]*t)), ext=1) # phi0[r]
        else:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[0] + 0.5*phiOn*(1-np.cos(self.ws[run]*t)), ext=1) # phi0[r]
        
        return phi_t
    '''    
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()
        
        if self.nRuns > 1: #len(phis) > 1: #nRuns???
            gsSin = plt.GridSpec(2,3)
            self.axIp = Ifig.add_subplot(gsSin[0,-1])
            self.axIss = Ifig.add_subplot(gsSin[1,-1], sharex=self.axIp)
            self.axI = Ifig.add_subplot(gsSin[:,:-1])
        else:
            self.addStimulus = config.addStimulus

            if self.addStimulus: 
                phi_ts = self.genPlottingStimuli()
                
                gsStim = plt.GridSpec(4,1)
                self.axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
                self.axI = Ifig.add_subplot(gsStim[1:,:], sharex=self.axS) # Photocurrent axes
                for run in range(self.nRuns):
                    for phiInd in range(self.nPhis):
                        pc = self.PD.trials[run][phiInd][vInd]
                        col, style = self.getLineProps(run, vInd, phiInd)
                        self.plotStimulus(phi_ts[run][phiInd], pc.begT, pc.pulses, pc.endT, self.axS, light='spectral', col=col, style=style)
                plt.setp(self.axS.get_xticklabels(), visible=False)
                self.axS.set_xlabel('') #plt.xlabel('')
                self.axS.set_ylim(self.phi0[0], max(self.phis)) ### phi0[r]
                if max(self.phis) / min(self.phis) >= 100:
                    self.axS.set_yscale('log') #plt.yscale('log')
            else:
                self.axI = Ifig.add_subplot(111)
            #self.axI = Ifig.add_subplot(111) # Combine with else condition below
        #plotLight(self.pulses, axI)

    
    def plotExtras(self):
        
        splineOrder = 2     #[1,5]
        
        trim = 0.1
        transEndInd = int(self.delDs[0] + round(self.onDs[0]*trim/self.dt))
        
        if self.nRuns > 1:
            #plt.figure(Ifig.number)
            #axI.legend().set_visible(False)
            
            #if len(self.phis) > 1:
            fstars = np.zeros((self.nPhis, self.nVs))
            for phiInd, phiOn in enumerate(self.phis): ### These loops need reconsidering...!!!
                for vInd, V in enumerate(self.Vs):
                    Ipeaks = np.zeros(self.nRuns)
                    for run in range(self.nRuns): 
                        PC = self.PD.trials[run][phiInd][vInd]
                        Ipeaks[run] = abs(PC.peak_) # Maximum absolute value over all peaks from that trial
                        Ip = self.PD.trials[np.argmax(Ipeaks)][phiInd][vInd].peak_
                    col, style = self.getLineProps(run, vInd, phiInd)
                    self.axIp.plot(self.fs, Ipeaks, 'x', color=col)
                    try:
                        intIp = InterpolatedUnivariateSpline(self.fs, Ipeaks, k=splineOrder)
                        #nPoints = 10*int(round(abs(np.log10(self.fs[-1])-np.log10(self.fs[0]))+1))
                        fsmooth = np.logspace(np.log10(self.fs[0]), np.log10(self.fs[-1]), num=1001)
                        self.axIp.plot(fsmooth, intIp(fsmooth))
                    except:
                        if verbose > 0:
                            print('Unable to plot spline for current peaks!')
                    fstar_p = self.fs[np.argmax(Ipeaks)]
                    fstars[phiInd,vInd] = fstar_p
                    Ap = max(Ipeaks)
                    fpLabel = '$f^*_{{peak}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_p,3))
                    self.axIp.plot(fstar_p, Ap, '*', markersize=10)
                    #axIp.annotate(fpLabel, xy=(fstar_p,Ap), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})
            
            self.axIp.set_xscale('log')
            self.axIp.set_ylabel('$|A|_{peak}$ $\mathrm{[nA]}$')
            if config.addTitles:
                #self.axIp.set_title('$\mathrm{|Amplitude|_{peak}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
                self.axIp.set_title('$f^*:=arg\,max_f(|A|_{peak})$')
            #axIp.set_aspect('auto')
            
            # Calculate the time to allow for transition effects from the period of fstar_p
            # buffer = 3
            # fstar_p = max(max(fstars))
            # transD = buffer * np.ceil(1000/fstar_p) # [ms]
            # transEndInd = round((self.delDs[0]+transD)/self.dt)
            # if transEndInd >= (self.onDs[0])/self.dt: # If transition period is greater than the on period
                # transEndInd = round((self.delDs[0]+self.onDs[0]/2)/self.dt) # Take the second half of the data
            
            tTransEnd = transEndInd*self.dt #ts[0][0][0]
            self.axI.axvline(x=tTransEnd, linestyle=':', color='k')
            for phiInd, phiOn in enumerate(self.phis): ### These loops need reconsidering...!!!
                for vInd, V in enumerate(self.Vs):
                    PC = self.PD.trials[np.argmax(Ipeaks)][phiInd][vInd]
                    onBegInd, onEndInd = PC.pulseInds[0]
                    t = PC.t
                    self.axI.annotate('', xy=(tTransEnd, Ip), xytext=(t[onEndInd], Ip), 
                                        arrowprops={'arrowstyle':'<->','color':'black','shrinkA':0,'shrinkB':0})
        
            for phiInd, phiOn in enumerate(self.phis):
                for vInd, V in enumerate(self.Vs):
                    Iabs = np.zeros(self.nRuns) #[None for r in range(nRuns)]
                    for run in range(self.nRuns): 
                        PC = self.PD.trials[run][phiInd][vInd]
                        onBegInd, onEndInd = PC.pulseInds[0]
                        t = PC.t #t = ts[run][phiInd][vInd]
                        I_RhO = PC.I #I_RhO = Is[run][phiInd][vInd]
                        #transEndInd = np.searchsorted(t,delD+transD,side="left") # Add one since upper bound is not included in slice
                        #if transEndInd >= len(t): # If transition period is greater than the on period
                        #    transEndInd = round(len(t[onBegInd:onEndInd+1])/2) # Take the second half of the data
                        #print(fstar_p,'Hz --> ',transD,'ms;', transEndInd,':',onEndInd+1)
                        I_zone = I_RhO[transEndInd:onEndInd+1]
                        try:
                            maxV = max(I_zone)
                        except ValueError:
                            maxV = 0.0
                        try:
                            minV = min(I_zone)
                        except ValueError:
                            minV = 0.0
                        Iabs[run] = abs(maxV-minV)
                    
                    #axI.axvline(x=t[transEndInd],linestyle=':',color='k')
                    #axI.annotate('Search zone', xy=(t[transEndInd], min(I_RhO)), xytext=(t[onEndInd], min(I_RhO)), arrowprops={'arrowstyle':'<->','color':'black'})
                    col, style = self.getLineProps(run, vInd, phiInd) ### Modify to match colours correctly
                    self.axIss.plot(self.fs, Iabs, 'x', color=col)
                    try:
                        intIss = InterpolatedUnivariateSpline(self.fs, Iabs, k=splineOrder)
                        #fsmooth = np.logspace(self.fs[0], self.fs[-1], 100)
                        self.axIss.plot(fsmooth, intIss(fsmooth))
                    except:
                        if verbose > 0:
                            print('Unable to plot spline for current steady-states!')
                    fstar_abs = self.fs[np.argmax(Iabs)]
                    fstars[phiInd,vInd] = fstar_abs
                    Aabs = max(Iabs)
                    fabsLabel = '$f^*_{{res}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_abs,3))
                    self.axIss.plot(fstar_abs, Aabs, '*', markersize=10, label=fabsLabel)
                    self.axIss.legend(loc='best')
                    #axIss.annotate(fabsLabel, xy=(fstar_abs,Aabs), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})
                    if verbose > 0:
                        print('Resonant frequency (phi={}; V={}) = {} Hz'.format(phiOn, V, fstar_abs))
            self.axIss.set_xscale('log')
            self.axIss.set_xlabel('$f$ $\mathrm{[Hz]}$')
            self.axIss.set_ylabel('$|A|_{ss}$ $\mathrm{[nA]}$')
            if config.addTitles:
                #axIss.set_title('$\mathrm{|Amplitude|_{ss}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
                self.axIss.set_title('$f^*:=arg\,max_f(|A|_{ss})$')
            
            plt.tight_layout()
            
            self.fstars = fstars
            if len(self.phis) > 1: # Multiple light amplitudes
                #for i, phi0 in enumerate(self.phi0):
                fstarAfig = plt.figure()
                for vInd, V in enumerate(self.Vs):
                    if self.phi0[0] > 0: # phi0[r]
                        plt.plot(np.array(self.phis)/self.phi0[0], fstars[:,vInd])
                        plt.xlabel('$\mathrm{Modulating}\ \phi(t)/\phi_0$') #\phi_1(t)/\phi_0(t)
                    else:
                        plt.plot(np.array(self.phis), fstars[:,vInd])
                        plt.xlabel('$\mathrm{Modulating}\ \phi(t)$') #\phi_1(t)
                plt.xscale('log')
                plt.ylabel('$f^*\ \mathrm{[Hz]}$')
                if config.addTitles:
                    plt.title('$f^*\ vs.\ \phi_1(t).\ \mathrm{{Background\ illumination:}}\ \phi_0(t)={:.3g}$'.format(self.phi0[0]))
            
        
class protDualTone(Protocol):
    # http://uk.mathworks.com/products/demos/signaltlbx/dtmf/dtmfdemo.html
    # http://dspguru.com/sites/dspguru/files/Sum_of_Two_Sinusoids.pdf
    protocol = 'dualTone'
    squarePulse = False
    
    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        # self.pulses = np.asarray(self.pulses)
        # self.nPulses = self.pulses.shape[0]
        # self.delDs = [row[0] for row in self.pulses] # pulses[:,0]    # Delay Durations
        # self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        # self.offDs = np.append(self.pulses[1:,0],self.totT) - self.pulses[:,1]
        
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
        self.cycles=np.column_stack((self.onDs,self.offDs))
        #self.cycles=np.tile(np.column_stack((self.onDs,self.offDs)),(self.nRuns,1))
        #self.padDs = np.zeros(self.nRuns)
        
        if (1000)/min(self.fs) > min(self.onDs):
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
            #print('Warning: The period of the lowest frequency is longer than the total simulation time!')

        # self.phis.sort(reverse=True)
        # self.Vs.sort(reverse=True)
        # self.nPhis = len(self.phis)
        # self.nVs = len(self.Vs)
        self.phi_ts = self.genPulseSet()

        self.runLabels = ["$\omega={}\mathrm{{rads/ms}}$ ".format(round_sig(w,3)) for w in self.ws]
        


class protChirp(Protocol):
    # http://en.wikipedia.org/wiki/Chirp
    protocol = 'chirp'
    squarePulse = False
    
    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
    
        #self.fs = np.sort(np.array(self.fs)) # Frequencies [Hz] 
        #self.ws = 2 * np.pi * self.fs / (1000) # Frequencies [rads/ms] (scaled from /s to /ms
        #self.sr = min([1000/(10*max(self.fs)), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        #self.sr = min([(1000)/(100*self.fT), self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        #self.dt = Sim.prepare(self.getShortestPeriod())
        
        #self.sr = max([(10)*max(self.f0,self.fT), 1000/self.dt]) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.sr = max(10000, int(round(10 * max(self.f0, self.fT)))) # Nyquist frequency - sampling rate (10*f) >= 2*f
        #self.dt = 1000/self.sr # Set by simulator but used for plotting
        #print(self.sr,self.dt)
        self.nRuns = 1 #len(self.ws)
        self.cycles = np.column_stack((self.onDs,self.offDs))
        #self.padDs = np.zeros(self.nRuns)
        
        #ws = 2 * np.pi * np.logspace(-4,10,num=7) # Frequencies [rads/s]
        #self.nRuns = len(freqs)
        if (1000)/self.f0 > min(self.onDs): #1/10**self.fs[0] > self.totT:
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')
            #print('Warning: The period of the lowest frequency is longer than the total simulation time!')
        
        if isinstance(self.phi0, (int, float, complex)):
            self.phi0 = np.ones(self.nRuns) * self.phi0
        elif isinstance(self.phi0, (list, tuple, np.ndarray)):
            if len(self.phi0) != self.nRuns:
                self.phi0 = np.ones(self.nRuns) * self.phi0[0]
        else:
            warning.warn('Unexpected data type for phi0: ', type(self.phi0))
        
        assert(len(self.phi0) == self.nRuns)
        
        self.phi_ts = self.genPulseSet()

    def getShortestPeriod(self):
        return 1000/self.sr
    
    def genPulse(self, run, phi, pulse):
        pStart, pEnd = pulse
        onD = pEnd - pStart
        t = np.linspace(0.0, onD, (onD*self.sr/1000)+1, endpoint=True) # Create smooth series of time points to interpolate between
        if self.linear: # Linear sweep
            ft = self.f0 + (self.fT-self.f0)*(t/onD)
        else:           # Exponential sweep
            ft = self.f0 * (self.fT/self.f0)**(t/onD)
        ft /= 1000 # Convert to frequency in ms
        if self.startOn:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[run] + 0.5*phi*(1+np.cos(ft*t)), ext=1, k=5)
        else:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[run] + 0.5*phi*(1-np.cos(ft*t)), ext=1, k=5)
        return phi_t
    
    '''    
    def getPhiFunc_orig(self, run, phiOn, delD, onD):
        pStart, pEnd = delD, (delD+onD)
        t = np.linspace(0.0, onD, (onD*self.sr/1000)+1, endpoint=True) # Create smooth series of time points to interpolate between
        if self.linear:
            ft = self.f0 + (self.fT-self.f0)*(t/pEnd)
        else: # Exponential sweep
            ft = self.f0 * (self.fT/self.f0)**(t/pEnd)
        ft /= 1000 # Convert to frequency in ms
        if self.startOn:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[0] + 0.5*phiOn*(1+np.cos(ft*t)), ext=1)
        else:
            phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi0[0] + 0.5*phiOn*(1-np.cos(ft*t)), ext=1)
        
        return phi_t        
    '''
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus

        if self.addStimulus: 
            phi_ts = self.genPlottingStimuli()
            
            gsStim = plt.GridSpec(4,1)
            self.axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
            self.axI = Ifig.add_subplot(gsStim[1:,:],sharex=self.axS) # Photocurrent axes
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    pc = self.PD.trials[run][phiInd][vInd]
                    col, style = self.getLineProps(run, vInd, phiInd)
                    self.plotStimulus(phi_ts[run][phiInd], pc.begT, pc.pulses, pc.endT, self.axS, light='spectral', col=col, style=style)
            plt.setp(self.axS.get_xticklabels(), visible=False)
            self.axS.set_xlabel('') #plt.xlabel('')
            
            self.axS.set_ylim(self.phi0[0], max(self.phis)) ### phi0[r]

            if max(self.phis) / min(self.phis) >= 100:
                self.axS.set_yscale('log') #plt.yscale('log')

            ### Overlay instantaneous frequency
            self.axF = self.axS.twinx()
            if not self.linear:
                self.axF.set_yscale('log')
            pc = self.PD.trials[0][0][0]
            for p in range(self.nPulses):
                pStart, pEnd = self.PD.trials[0][0][0].pulses[p]
                #pEnd = self.onDs[p]
                #delD = self.delDs[p]
                onD = pEnd - pStart
                nPoints = 10*int(round(onD/self.dt))+1 # 10001
                tsmooth = np.linspace(0, onD, nPoints)

                if self.linear:
                    ft = self.f0 + (self.fT-self.f0)*(tsmooth/onD)
                else: # Exponential
                    ft = self.f0 * (self.fT/self.f0)**(tsmooth/onD)
                #print(ft)
                self.axF.plot(tsmooth+pStart, ft, 'g')
            self.axF.set_ylabel('$f\ \mathrm{[Hz]}$')
        else:
            self.axI = Ifig.add_subplot(111)
            
        #plotLight(self.pulses, self.axI)
        

class protRamp(Protocol):
    """Linearly increasing pulse"""
    
    protocol = 'ramp'
    squarePulse = False
    nRuns = 1

    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.nRuns = 1 #nRuns # Make len(phi_ton)?
        self.cycles = np.column_stack((self.onDs,self.offDs))
        #self.cycles=np.tile(np.column_stack((self.onDs,self.offDs)),(self.nRuns,1))
        #self.padDs = np.zeros(self.nRuns)
        
        self.phi_ts = self.genPulseSet()# [[[None for pulse in range(self.nPulses)] for phi in range(self.nPhis)] for run in range(self.nRuns)]
        
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        
        if self.addStimulus: 
            phi_ts = self.genPlottingStimuli()
            
            gsStim = plt.GridSpec(4,1)
            self.axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
            self.axI = Ifig.add_subplot(gsStim[1:,:],sharex=self.axS) # Photocurrent axes
            pc = self.PD.trials[0][0][0]
            plotLight(pc.pulses, ax=self.axS, light='spectral', lam=470, alpha=0.2)
            #for p in range(self.nPulses):
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    pc = self.PD.trials[run][phiInd][vInd]
                    col, style = self.getLineProps(run, vInd, phiInd)
                    self.plotStimulus(phi_ts[run][phiInd], pc.begT, self.pulses, pc.endT, self.axS, light=None, col=col, style=style) #light='spectral'
            plt.setp(self.axS.get_xticklabels(), visible=False)
            #plt.xlabel('')
            self.axS.set_xlabel('')
            #if phis[-1]/phis[0] >= 100:
            #    plt.yscale('log')
        else:
            self.axI = Ifig.add_subplot(111)
                
    
    
    def genPulse(self, run, phi, pulse):
        """Generate spline for a particular pulse. phi0 is the offset so decreasing ramps can be created with negative phi values. """
        pStart, pEnd = pulse
        phi_t = InterpolatedUnivariateSpline([pStart, pEnd], [self.phi0, self.phi0+phi], k=1, ext=1) #self.phi_ton
        return phi_t
    
    """
    def getPhiFunc_orig(self, run, phiOn, delD, onD):
        pStart, pEnd = delD, (delD+onD)
        phi_t = InterpolatedUnivariateSpline([pStart,pEnd], [self.phi_ton,phiOn], k=1, ext=1) #[start,delD,end,totT], [0,self.phi_ton,phiOn,0] 
        #phi_t = InterpolatedUnivariateSpline(pStart + t, self.phi_ton + phiOn*(t/onD), k=1, ext=1) #[start,delD,end,totT], [0,self.phi_ton,phiOn,0] 
        return phi_t
    """
        
class protDelta(Protocol): 
    # One very short, saturation intensity pulse e.g. 10 ns @ 100 mW*mm^-2 for wild type ChR
    # Used to calculate gbar, assuming that O(1)-->1 as onD-->0 and phi-->inf
    protocol = 'delta'
    squarePulse = True
    nRuns = 1
    
    def prepare(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
        self.cycles = np.asarray([[self.onD, self.totT-self.delD]])
        self.nPulses = self.cycles.shape[0]
        self.pulses, self.totT = cycles2times(self.cycles, self.delD)
        self.delDs = np.array([row[0] for row in self.pulses], copy=True) # pulses[:,0]    # Delay Durations
        self.onDs = [row[1]-row[0] for row in self.pulses] # pulses[:,1] - pulses[:,0]   # Pulse Durations
        self.offDs = np.append(self.pulses[1:,0], self.totT) - self.pulses[:,1]
        
        if np.isscalar(self.phis):
            self.phis = np.asarray([self.phis])
        self.phis.sort(reverse=True)
        self.nPhis = len(self.phis)
        
        if np.isscalar(self.Vs):
            self.Vs = np.asarray([self.Vs])
        self.Vs.sort(reverse=True)
        self.nVs = len(self.Vs)
        
        self.extraPrep()
        return

    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        
        self.nRuns = 1 #nRuns
        
        self.phi_ts = self.genPulseSet()
        
        
    def finish(self, PC, RhO):
        # Take the max over all runs, phis and Vs?
        # Ipmax = minmax(self.IpVals[run][phiInd][vInd][:])# = I_RhO[peakInds]
        if PC.V is None:
            return
            
        try: #if V != RhO.E:
            Gmax = PC.peak_ / (PC.V - RhO.E) #Ipmax / (V - RhO.E) # Assuming [O_p] = 1 ##### Should fV also be used?
        except ZeroDivisionError: #else:
            print("The clamp voltage must be different to the reversal potential!")
        
        gbar_est = Gmax * 1e6
        
        if verbose > 0:
            print("Estimated maximum conductance (g) = {} uS".format(round_sig(gbar_est,3)))

    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
            
        if self.addStimulus: 
            phi_ts = self.genPlottingStimuli()
            
            gsStim = plt.GridSpec(4,1)
            self.axS = Ifig.add_subplot(gsStim[0,:]) # Stimulus axes
            self.axI = Ifig.add_subplot(gsStim[1:,:],sharex=self.axS) # Photocurrent axes
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis):
                    pc = self.PD.trials[run][phiInd][vInd]
                    col, style = self.getLineProps(run, vInd, phiInd)
                    self.plotStimulus(phi_ts[run][phiInd], pc.begT, pc.pulses, pc.endT, self.axS, light='spectral', col=col, style=style)
            plt.setp(self.axS.get_xticklabels(), visible=False)
            self.axS.set_xlabel('') #plt.xlabel('')
            if max(self.phis) / min(self.phis) >= 100:
                self.axS.set_yscale('log')
        else:
            self.axI = Ifig.add_subplot(111)
        
        #plotLight(self.pulses, self.axI)
        
    
    def addAnnotations(self):
        #plt.figure(Ifig.number)
        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
                    pc = self.PD.trials[run][phiInd][vInd]
                    # Maximum only...
                    #Ip = pc.peak_
                    #tp = pc.tpeak_
                    for p in range(self.nPulses):
                        Ip = pc.peaks_[p]
                        tp = pc.tpeaks_[p]
                        tlag = pc.lags_[p]
                        self.axI.axvline(x=tp, linestyle=':', color='k')
                        #plt.axhline(y=I_RhO[peakInds[0]], linestyle=':', color='k')
                        label = '$I_{{peak}} = {:.3g}\mathrm{{nA;}}\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(Ip, tlag)
                        plt.text(1.05*tp, 1.05*Ip, label, ha='left', va='bottom', fontsize=config.eqSize)

            
            
class protRectifier(Protocol):
    """Protocol to determine the rectification parameters of rhodopsins. Typically they are inward rectifiers where current is more easily passed into the cell than out. """
    # Iss vs Vclamp
    # http://en.wikipedia.org/wiki/Inward-rectifier_potassium_ion_channel
    protocol = 'rectifier'
    squarePulse = True
    nRuns = 1
        
    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.nRuns = 1 #nRuns
        
        self.phi_ts = self.genPulseSet()
        
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()
                                                        
        self.gsIR = plt.GridSpec(2,3)
        self.axI = Ifig.add_subplot(self.gsIR[:,0:2])
        self.axVI = Ifig.add_subplot(self.gsIR[0,-1])#, sharey=self.axI)
        self.axfV = Ifig.add_subplot(self.gsIR[-1,-1], sharex=self.axVI)

        
    def plotExtras(self):
        # TODO: Refactor!!!
        #plt.figure(Ifig.number) #IssVfig = plt.figure()
        ax = self.axVI #IssVfig.add_subplot(111)
        
        legLabels = [None for p in range(self.nPhis)]
        #eqString = r'$f(v) = \frac{{{v1:.3}}}{{v-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{v-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'
        Vs = self.Vs
        for run in range(self.nRuns):
            for phiInd, phiOn in enumerate(self.phis): 
                #RhO.calcSteadyState(phiOn) ##################################### Is this necessary? Only for adjusting the gain (g)
                #print(self.IssVals[run][phiInd][:])
                #popt, pcov, eqString = self.fitfV(Vs,self.IssVals[run][phiInd][:],calcIssfromfV,p0fV,RhO,ax)#,eqString)
                #popt, pcov, eqString = self.fitfV(self.Vs, self.PD.ss_[run][phiInd][:], calcIssfromfV, p0fV, RhO, ax)#,eqString)
                
                Iss = self.PD.ss_[run][phiInd][:]
                
                ### Original routines
                ##popt, pcov, eqString = fitFV(self.Vs, Iss, p0FV, ax=ax)
                #p0FV = (35, 15, 0)
                E_i = 0
                v0_i = 35
                g0 = 25000 ### Rethink...
                pseudoV1 = calcV1(E_i, v0_i) * (g0 * 1e-6 * 0.5) # g0*f(phi)*v1 (assuming I in nA and f(phi)=0.5)
                p0FV = (v0_i, pseudoV1, E_i)
                
                poptI, poptg = fitFV(Vs, Iss, p0FV) #, ax=ax)
                '''
                ### From fitfV() and fitFV() --> poptI
                def calcRect(V, v0, v1, E): #, gpsi):
                    if type(V) != np.ndarray:
                        V = np.array(V)
                    fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
                    fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
                    return fV * (V - E) # * gpsi
                poptI, pcov = curve_fit(calcRect, Vs, Iss, p0FV)
                '''
                
                ### New routines
                pfV = Parameters()
                pfV.add_many(
                ('E',   0,           True, -100 , 100,  None),
                ('v0',  50,          True, -1e12, 1e12, None),
                ('v1',  calcV1(0,50),True, -1e9 , 1e9,  None))
                
                pfV = fitfV(Vs, Iss, pfV)
                #print(pfV)
                
                Vrange = max(Vs) - min(Vs)
                Vsmooth = np.linspace(min(Vs), max(Vs), 1+Vrange/.1) #Prot.dt
                
                E = poptI[2]
            #    E = pfV['E'].value
                #v0 = pfV['v0'].value
                #v1 = pfV['v1'].value
                
                ### Top plot
                # self.RhO.v0 = poptI[0]
                # self.RhO.v1 = poptI[1]
                # self.RhO.E = E #poptI[2]
                # fVsmooth = self.RhO.calcfV(Vsmooth)
                
                pfV['E'].value = E
                pfV['v0'].value = poptI[0] #v0
                pfV['v1'].expr = ''
                pfV['v1'].value = poptI[1] #v1 # poptg[1]?
                
                FVsmooth = errFV(pfV, Vsmooth)
                
                '''
                pFV = Parameters()
                copyParam(['E', 'v0', 'v1'], pfV, pFV)
                pFV['v1'].expr = ''
                pFV['v0'].value, pFV['v1'].value, pFV['E'].value = poptI
                FVsmooth = errFV(pFV, Vsmooth)
                '''
                
                '''
                gNorm = getNormGs(Vs, Iss, E, v=-70)
                
                fVmin = minimize(errfV, pfV, args=(Vs, gNorm), method=method)#, tol=1e-12)
                pfVfinal = fVmin.params
                v0 = pfVfinal['v0'].value
                v1 = pfVfinal['v1'].value
                
                zeroErrs = np.isclose(Vs, np.ones_like(Vs)*E)
                gNorm[zeroErrs] = v1/v0
                
                ### From fitFV() --> poptg
                def calcScale(V, v0, v1): #, gpsi):
                    if type(V) != np.ndarray:
                        V = np.array(V)
                    fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
                    fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
                    return fV
                
                poptg, pcov = curve_fit(calcScale, Vs, gNorm, p0=(v0, v1))
                '''
                

                
                ax.plot(Vsmooth, FVsmooth) #*(Vsmooth-E)) #,label=peakEq)#,linestyle=':', color='#aaaaaa')
                #FVsmooth = errfV(pfV, Vsmooth)
                #ax.plot(Vsmooth, FVsmooth*(Vsmooth-E))
                
                #col, = getLineProps(Prot, 0, 0, 0) #Prot, run, vInd, phiInd
                #plt.plot(Vs,Iss,linestyle='',marker='x',color=col)
                markerSize=40
                ax.scatter(Vs, Iss, marker='x', color=colours, s=markerSize)#,linestyle=''
                
                
                ### Bottom plot
                # self.RhO.v0 = poptg[0]
                # self.RhO.v1 = poptg[1]
                # self.RhO.E = E #poptI[2]
                # fVsmooth = self.RhO.calcfV(Vsmooth)
                
                #pfV['E'].value = E
                pfV['v0'].value = poptg[0] #v0
                pfV['v1'].value = poptg[1] #v1
                fVsmooth = errfV(pfV, Vsmooth)
                
                self.axfV.plot(Vsmooth, fVsmooth)
            #    fVstring = eqString.format(v0=poptg[0], E=poptI[2], v1=poptg[1])
                v0 = pfV['v0'].value
                v1 = pfV['v1'].value
                if np.isclose(E, 0, atol=0.005):
                    eqString = r'$f(v) = \frac{{{v1:.3}}}{{v-{E:.0f}}} \cdot \left[1-\exp\left({{-\frac{{v-{E:.0f}}}{{{v0:.3}}}}}\right)\right]$'
                    fVstring = eqString.format(E=np.abs(E), v0=v0, v1=v1)
                else:
                    eqString = r'$f(v) = \frac{{{v1:.3}}}{{v-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{v-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'
                
                    fVstring = eqString.format(E=E, v0=v0, v1=v1)
                
                #v0 = poptg[0]
                #v1 = poptg[1]
                
                #E = poptI[2]
                #vInd = np.searchsorted(self.Vs, (-70 - E))
                #sf = Iss[vInd]
                
                #sf = Iss[Vs.index(-70)]
                #g0 = Iss / (Vs - E)
                #gNorm = g0 / (sf / (-70 - E))
                
                gs = Iss / (np.asarray(Vs) - E) # 1e6 * 
                gm70 = Iss[Vs.index(-70)] / (-70 - E)# * -70 # 1e6 * 
                if verbose > 0:
                    print('g(v=-70) = ', gm70)
                #g0[(Vs - E)==0] = None #(v1/v0)
                gNorm = gs / gm70 # Normalised conductance relative to V=-70
                
                self.axfV.scatter(Vs, gNorm, marker='x', color=colours, s=markerSize)#,linestyle=''
                if verbose > 1:
                    print(gm70)
                    if verbose > 2:
                        print(np.c_[Vs, np.asarray(Vs)-E, Iss, gs, gNorm])
                
                
                # Add equations to legend
                if self.nPhis > 1: 
                    legLabels[phiInd] = fVstring + '$,\ \phi={:.3g}$'.format(phiOn)
                else:
                    legLabels[phiInd] = fVstring
                
                ### Move this to fitting routines?
                # v0 = popt[0], v1 = popt[1], E = popt[2]
        
        #if len(phis) > 1:
        #ax.legend(legLabels, loc='best')
        
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
        ax.set_xlim(min(Vs), max(Vs))
        
        ax.set_ylabel('$I_{ss}$ $\mathrm{[nA]}$')#, position=(0.95,0.8)) #plt.xlabel
        
        ax = self.axfV
        ax.set_ylabel('$f(v)$ $\mathrm{[1]}$')#, position=(0.95,0.8)) #plt.xlabel
        
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlim(min(Vs), max(Vs))
        
        # yticks = ax.get_yticklabels()
        # ax.set_ylim(0, float(yticks[-1].get_text()))
        
        useLegend = True
        if useLegend:
            #ax.legend(legLabels, bbox_to_anchor=(0., 1.01, 1., .101), loc=3, mode="expand", borderaxespad=0., prop={'size':mp.rcParams['font.size']})
            ax.legend(legLabels, loc='best')
        else:
            ymin, ymax = ax.get_ylim()
            #ax.set_ylim(ymin, ymax)
            
            ax.text(min(Vs), 0.98*ymax, legLabels[phiInd], ha='left', va='top')#, fontsize=eqSize) #, transform=ax.transAxes)
        
        # ax.axvline(x=-70, linestyle=':', color='k')
        # yind = np.searchsorted(Vsmooth, -70)
        # ax.axhline(y=fVsmooth[yind], linestyle=':', color='k')
        
        ax.vlines(x=-70, ymin=0, ymax=1, linestyle=':', color='k')
        ax.hlines(y=1, xmin=-70, xmax=0, linestyle=':', color='k')
        
        #ax.set_xlabel('$V_{clamp}$ $\mathrm{[mV]}$', position=(0.8,0.8)) #plt.xlabel
        ax.set_xlabel('$V_{clamp}$ $\mathrm{[mV]}$', position=(xLabelPos,0), ha='right')
        #plt.xlim((min(Vs),max(Vs)))
        #ax.set_ylabel('$I_{ss}$ $\mathrm{[nA]}$', position=(0.55,0.05)) # Shared axis
        
        
        
        self.axI.grid(b=True, which='minor', axis='both', linewidth=.2)
        self.axI.grid(b=True, which='major', axis='both', linewidth=1)
        
        plt.tight_layout()
        


        
class protShortPulse(Protocol):
    # Vary pulse length - See Nikolic+++++2009, Fig. 2 & 9
    #def __init__(self, phis=[1e12], Vs=[-70], delD=25, pulses=[[1,74],[2,73],[3,72],[5,70],[8,67],[10,65],[20,55]], nRuns=1, dt=0.1):
    protocol = 'shortPulse'
    squarePulse = True
    nPulses = 1 #Fixed at 1
    
    # def __next__(self):
        # if self.run >= self.nRuns:
            # raise StopIteration
        # #return cycles2times(self.cycles[run], self.delDs[run]) #np.asarray[self.pulses[self.run]]
        # #return self.cycles[run], self.totT
        # return self.getRunCycles(self, self.run)
    
    def prepare(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
        self.pDs = np.sort(np.array(self.pDs))
        self.nRuns = len(self.pDs)
        self.delDs = np.ones(self.nRuns)*self.delD
        self.onDs = self.pDs
        self.offDs = (np.ones(self.nRuns)*self.totT) - self.delDs - self.onDs
        self.cycles = np.column_stack((self.onDs,self.offDs))
        #self.padDs = np.zeros(self.nRuns)
        self.phis.sort(reverse=True)
        self.Vs.sort(reverse=True)
        self.nPhis = len(self.phis)
        self.nVs = len(self.Vs)
        
        self.phi_ts = self.genPulseSet()
        
        self.runLabels = ["$\mathrm{{Pulse}}={}\mathrm{{ms}}$ ".format(pD) for pD in self.pDs]
            

    def getRunCycles(self, run):
        return np.asarray([[self.onDs[run],self.offDs[run]]]), self.delDs[run]
    
    
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()
        
        gsPL = plt.GridSpec(2,3)
        self.axLag = Ifig.add_subplot(gsPL[0,-1])
        self.axPeak = Ifig.add_subplot(gsPL[1,-1], sharex=self.axLag)
        self.axI = Ifig.add_subplot(gsPL[:,:-1])
        
    
    def addAnnotations(self):
        # Freeze axis limits
        ymin, ymax = plt.ylim()
        #plt.ylim(ymin, ymax)
        pos = 0.02 * abs(ymax-ymin)
        #plt.ylim(ax.get_ylim())
        plt.ylim(ymin, pos*(self.nRuns+1)) # Allow extra space for thick lines
        for run in range(self.nRuns): 
            
            # Loop over light intensity...
            for phiInd, phiOn in enumerate(self.phis): #for phiInd in range(0, len(phis)):
                # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                for vInd, V in enumerate(self.Vs): #range(0, len(Vs)):
                    col, style = self.getLineProps(run, vInd, phiInd)
                    #peakInds = IpInds[run][phiInd][vInd]
                    
                    PC = self.PD.trials[run][phiInd][vInd]
                    t_on, t_off = PC.pulses[0,:]
                    
                    # self.axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)
                    # self.axI.axvline(x=delD,linestyle=':',color='k')
                    # self.axI.axvline(x=delD+onD,linestyle=':',color=col)
                    
                    self.axI.hlines(y=(run+1)*pos, xmin=t_on, xmax=t_off, linewidth=4, color=col)
                    self.axI.axvline(x=t_on, linestyle=':', color='k')
                    self.axI.axvline(x=t_off, linestyle=':', color=col)
                    
                    
                    #plt.figure(Ifig.number)
                    # ymin, ymax = plt.ylim()
                    # pos = 0.02 * abs(ymax-ymin)
    
                    # axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)#'b')  #colours[c]
                    #axI.plot(t[peakInds],I_RhO[peakInds],marker='*',color=col)
                    #axLag.plot(pDs[run],(t[IpInds[run][phiInd][vInd]] - delD),marker='*',markersize=10,color=col)
                    #axPeak.plot(pDs[run],I_RhO[peakInds],marker='*',markersize=10,color=col)
                    self.axI.plot(PC.tpeaks_, PC.peaks_, marker='*', color=col)
                    #self.axLag.plot(self.pDs[run], (PC.tpeaks_ - delD), marker='*', markersize=10, color=col)
                    self.axLag.plot(self.pDs[run], PC.lags_[0], marker='*', markersize=10, color=col)
                    self.axPeak.plot(self.pDs[run], PC.peaks_, marker='*', markersize=10, color=col)
    
                    #peakInds = IpInds[run][phiInd][vInd]
                    #axI.hlines(y=(run+1)*pos,xmin=delD,xmax=delD+onD,linewidth=4,color=col)
        ### Plot figure to show time of Ipeak vs time of light off c.f. Nikolic et al. 2009 Fig 2b
    #     axLag = Ifig.add_subplot(gsPL[0,-1])
    #     tpeaks = [(t[IpInds[p][0][0]] - delD) for p in range(nRuns)]
    #     toffs = [(t[PulseInds[p][0][0][0][1]-1] - delD) for p in range(nRuns)] # pDs
    #     plot(toffs,tpeaks,marker='x',linestyle='')
        
    #     axLag.axis('equal')
        tmax = max(self.pDs)*1.25
        self.axLag.plot([0,tmax], [0,tmax], ls="--", c=".3")
        self.axLag.set_xlim(0,tmax)
        self.axLag.set_ylim(0,tmax)
        self.axLag.set_ylabel('$\mathrm{Time\ of\ peak\ [ms]}$')
    #     plt.tight_layout()
        self.axLag.set_aspect('auto')
    #     print(axLag.get_xlim())
    #     diag_line, = axLag.plot(axLag.get_xlim(), axLag.get_ylim(), ls="--", c=".3")
        
        
        ### Plot figure to show current peak vs time of light off c.f. Nikolic et al. 2009 Fig 2c
        self.axPeak.set_xlim(0,tmax)
        self.axPeak.set_xlabel('$\mathrm{Pulse\ duration\ [ms]}$')
        self.axPeak.set_ylabel('$\mathrm{Photocurrent\ peak\ [nA]}$')
        
class protRecovery(Protocol):
    # Vary Inter-Pulse-Interval
    protocol = 'recovery'
    squarePulse = True
    nPulses = 2 # Fixed at 2 for this protocol
    
    # def __next__(self):
        # if self.run >= self.nRuns:
            # raise StopIteration
        # return np.asarray[self.pulses[self.run]]
        
    def prepare(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.IPIs = np.sort(np.asarray(self.IPIs)) #np.sort(np.array(IPIs))
        
        self.nRuns = len(self.IPIs)
        self.delDs = np.ones(self.nRuns)*self.delD
        self.onDs = np.ones(self.nRuns)*self.onD
        self.offDs = self.IPIs
        ##cycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI-onD for IPI in IPIs])) # [:,0] = on phase duration; [:,1] = off phase duration
        #cycles=np.column_stack((onD*np.ones(len(IPIs)),[IPI for IPI in IPIs]))
        self.cycles = np.column_stack((self.onDs, self.offDs))
        
        self.pulses, _ = cycles2times(self.cycles, self.delD)
        self.runCycles = np.zeros((self.nPulses, 2, self.nRuns))
        for run in range(self.nRuns):
            self.runCycles[:,:,run] = np.asarray([[self.onDs[run],self.offDs[run]],[self.onDs[run],self.offDs[run]]])
        #padDs = [totT-((onD+pOff)*nPulses)-delD for pOff in cycles[:,1]] # Not necessary
        
        #self.padDs = np.zeros(self.nRuns); self.padDs[-1] = -0.8*max(self.cycles[:,1])
        #self.totT = self.delD+self.nPulses*max(self.IPIs) -0.8*max(self.cycles[:,1]) # Total simulation time per run [ms] ### Should IPI be one pulse cycle or time between end/start of two pulses?
        self.begT = 0
        self.endT = self.totT
        IPIminD = max(self.delDs) + (2*max(self.onDs)) + max(self.IPIs)
        if self.endT < IPIminD:
            warnings.warn("Insufficient run time for all stimulation periods!")
        else:
            #self.cycles[:,1] = self.totT - IPIminD
            self.runCycles[-1,1,:] = self.totT - IPIminD
            
        
        self.IpIPI = np.zeros(self.nRuns) #peaks
        self.tpIPI = np.zeros(self.nRuns) #tPeaks
        
        if np.isscalar(self.phis):
            self.phis = np.asarray([self.phis])
        self.phis.sort(reverse=True)
        self.nPhis = len(self.phis)
        
        if np.isscalar(self.Vs):
            self.Vs = np.asarray([self.Vs])
        self.Vs.sort(reverse=True)
        self.nVs = len(self.Vs)
        
        #self.protPulses = self.getProtPulses()
        
        self.phi_ts = self.genPulseSet()

        self.runLabels = ["$\mathrm{{IPI}}={}\mathrm{{ms}}$ ".format(IPI) for IPI in self.IPIs]
        

    def getRunCycles(self, run):
        #return np.asarray([[self.onDs[run],self.offDs[run]],[self.onDs[run],self.offDs[run]-self.padDs[run]]]), self.delDs[run]
        #return np.asarray([[self.onDs[run],self.offDs[run]],[self.onDs[run],self.offDs[run]]]), self.delDs[run]
        return self.runCycles[:,:,run], self.delDs[run]
    
    def finish(self, PC, RhO):
        ### Build array of second peaks
        self.PD.IPIpeaks_ = np.zeros((self.nRuns, self.nPhis, self.nVs)) #[[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.PD.tIPIpeaks_ = np.zeros((self.nRuns, self.nPhis, self.nVs))
        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
            #for phiInd, phiOn in enumerate(self.phis): 
                #for vInd, V in enumerate(self.Vs): 
                    PC = self.PD.trials[run][phiInd][vInd]
                    PC.alignToPulse(pulse=0, alignPoint=2) # End of the first pulse
                    self.PD.IPIpeaks_[run][phiInd][vInd] = PC.peaks_[1]
                    self.PD.tIPIpeaks_[run][phiInd][vInd] = PC.tpeaks_[1]
                    ### Search only within the on phase of the second pulse
                    # startInd = self.PulseInds[run][phiInd][vInd][1,0]
                    # endInd = self.PulseInds[run][phiInd][vInd][1,1]
                    # extOrder = int(1+endInd-startInd) #100#int(round(len(I_RhO)/5))
                    # #peakInds = findPeaks(I_RhO[:endInd+extOrder+1],minmax,startInd,extOrder)
                    # peakInds = findPeaks(I_RhO[:endInd+extOrder+1],startInd,extOrder)
                    # if len(peakInds) > 0: # Collect data at the (second) peak
                        # self.IpIPI[run] = I_RhO[peakInds[0]] #-1 peaks
                        # self.tpIPI[run] = t[peakInds[0]] # tPeaks

        if verbose > 1:
            print(self.PD.tIPIpeaks_) # tPeaks
            print(self.PD.IPIpeaks_) # peaks
            #popt = fitPeaks(tPeaks, peaks, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$')
            #print("tau_R = {} ==> rate_R = {}".format(popt[1],1/popt[1]))

    '''
    def createLayout(self, Ifig=None, vInd=0):
    
        if Ifig == None:
            Ifig = plt.figure()
        
        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()
        
        # Default layout
        self.axI = Ifig.add_subplot(111)
        #plt.sca(self.axI)
        #plotLight(self.pulses, self.axI)
    '''
        
    def addAnnotations(self):
    

        # Freeze axis limits
        ymin, ymax = plt.ylim()
        pos = 0.02 * abs(ymax-ymin)
        #plt.ylim(ymin, ymax)
        plt.ylim(ymin, pos*self.nRuns)
        #plt.ylim(ax.get_ylim())
        
        xmin, xmax = plt.xlim()
        plt.xlim(xmin, xmax)
        
        for run in range(self.nRuns): 
            if self.nRuns > 1:
                delD = self.delDs[run]
                onD = self.cycles[run,0]
                offD = self.cycles[run,1]
                #padD = self.padDs[run]
            else:
                delD = self.delD
                #...
            # Loop over light intensity...
            for phiInd, phiOn in enumerate(self.phis): #for phiInd in range(0, len(phis)):
                # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                for vInd, V in enumerate(self.Vs): #range(0, len(Vs)):
                    col, style = self.getLineProps(run, vInd, phiInd)
                    #plt.figure(Ifig.number)
                    #plt.annotate('', (delD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col})
                    #plt.annotate('', (delD+onD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col,'shrinkA':0,'shrinkB':0})
                    
                    #plt.annotate('', (delD+onD, (run+1)*pos), (delD+onD+offD, (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col,'shrinkA':0,'shrinkB':0})
                    
                    pulses = self.PD.trials[run][phiInd][vInd].pulses
                    plt.annotate('', (pulses[0,1], (run+1)*pos), (pulses[1,0], (run+1)*pos), arrowprops={'arrowstyle':'<->','color':col,'shrinkA':0,'shrinkB':0})
                    if run == 0:
                        
                        ### Fitting
                        #popt, _, _ = fitRecovery(self.PD.tIPIpeaks_[:,phiInd,vInd], self.PD.IPIpeaks_[:,phiInd,vInd], self.endT, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$',self.axI) # tPeaks peaks
                        
                        # Prepend t_off0 and Iss0
                        # tss0 = self.PD.trials[run][phiInd][vInd].pulses[0,1]
                        # Iss0 = self.PD.trials[run][phiInd][vInd].sss_[0]
                        # Ipeak0 = self.PD.trials[run][phiInd][vInd].peaks_[0]
                        # t_peaks = np.r_[tss0, self.PD.tIPIpeaks_[:,phiInd,vInd]]
                        # I_peaks = np.r_[Iss0, self.PD.IPIpeaks_[:,phiInd,vInd]]
                        t_peaks, I_peaks, Ipeak0, Iss0 = getRecoveryPeaks(self.PD, phiInd, vInd, usePeakTime=True)
                        params = Parameters()
                        params.add('Gr0', value=0.002, min=0.0001, max=0.1)
                        #self.RhO.exportParams(params)
                        params = fitRecovery(t_peaks, I_peaks, params, Ipeak0, Iss0, self.axI)
                        if verbose > 0:
                            #print("tau_R = {} ==> rate_R = {}".format(popt[1],1/popt[1]))
                            print("tau_r0 = {} ==> G_r0 = {}".format(1/params['Gr0'].value, params['Gr0'].value))
        
        #plt.legend(loc='upper right')
            
from collections import OrderedDict
protocols = OrderedDict([('step', protStep), ('delta', protDelta), ('sinusoid', protSinusoid), ('chirp', protChirp), ('ramp', protRamp), ('recovery', protRecovery), ('rectifier', protRectifier), ('shortPulse', protShortPulse), ('custom', protCustom)])

# E.g. 
# protocols['shortPulse']([1e12], [-70], 25, [1,2,3,5,8,10,20], 100, 0.1)

#squarePulses = [protocol for protocol in protocols if protocol.squarePulse]
#arbitraryPulses = [protocol for protocol in protocols if not protocol.squarePulse]
#squarePulses = {'custom': True, 'delta': True, 'step': True, 'rectifier': True, 'shortPulse': True, 'recovery': True}
#arbitraryPulses = {'custom': True, 'sinusoid': True, 'chirp': True, 'ramp':True} # Move custom here
#smallSignalAnalysis = {'sinusoid': True, 'step': True, 'delta': True} 



def selectProtocol(protocol, params=None, saveData=True):
    """Protocol selection function"""
    if protocol in protList:
        if params:
            return protocols[protocol](params, saveData=saveData)
        else:
            return protocols[protocol](params=protParams[protocol], saveData=saveData)
    else:
        raise NotImplementedError(protocol)
        #print("Error in selecting protocol - please choose from 'custom', 'delta', 'rectifier', 'shortPulse' or 'recovery'")
        

### Protocols to be included in the next version:
### - Temperature (Q10)
### - pH (intracellular and extracellular)
### - Wavelength (lambda)


def characterise(RhO):
    """Run small signal analysis on Rhodopsin"""
    for protocol in smallSignalAnalysis:
        RhO.setLight(0.0)
        Prot = protocols[protocol]()
        Sim = simulators['Python'](Prot, RhO)
        Sim.run()
        Sim.plot()
    return
