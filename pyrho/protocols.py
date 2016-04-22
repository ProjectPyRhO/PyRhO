"""
Stimulation protocols to run on the opsin models
    * Neuro-engineering stimuli: ``step``, ``sinusoid``, ``chirp``, ``ramp``, ``delta``
    * Opsin-specific protocols: ``rectifier``, ``shortPulse``, ``recovery``.
    * The ``custom`` protocol can be used with arbitrary interpolation fuctions
"""

from __future__ import print_function, division
import warnings
import logging
import os
import abc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl # for tick locators
from scipy.interpolate import InterpolatedUnivariateSpline as spline
#from scipy.optimize import curve_fit
from lmfit import Parameters

from pyrho.parameters import *
from pyrho.parameters import PyRhOobject
from pyrho.utilities import * # times2cycles, cycles2times, plotLight, round_sig, expDecay, biExpDecay, findPeaks
from pyrho.expdata import * #import loadData
from pyrho.fitting import fitFV, errFV, fitfV, errfV, getRecoveryPeaks, fitRecovery
from pyrho.models import *
from pyrho.simulators import * # For characterise()
from pyrho.config import *
from pyrho import config

__all__ = ['protocols', 'selectProtocol', 'characterise']


class Protocol(PyRhOobject): #, metaclass=ABCMeta
    """Common base class for all protocols"""

    __metaclass__ = abc.ABCMeta

    protocol = None
    nRuns = None
    delD = None
    cycles = None
    totT = None
    dt = None
    phis = None
    Vs = None

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
        self.PD = None
        self.Ifig = None

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
        return self.getRunCycles(self.run - 1)

    def prepare(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""

        if np.isscalar(self.cycles): # Only 'on' duration specified
            onD = self.cycles
            if hasattr(self, 'totT'):
                offD = self.totT - onD - self.delD
            else:
                offD = 0
            self.cycles = np.asarray([[onD, offD]])

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
        return (self.cycles, self.delD)

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
                    phi_ts[run][phiInd][pInd] = genPulse(run, phi, pulse)
        self.phi_ts = phi_ts
        return phi_ts

    def genPulse(self, run, phi, pulse):
        """Default interpolation function for square pulses"""
        pStart, pEnd = pulse
        phi_t = spline([pStart,pEnd], [phi,phi], k=1, ext=1)
        return phi_t

    def genPlottingStimuli(self, genPulse=None, vInd=0):
        """Redraw stimulus functions in case data has been realigned"""
        if genPulse is None:
            genPulse = self.genPulse

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


    def getStimArray(self, run, phiInd, dt) : #phi_ts, delD, cycles, dt):
        """Return a stimulus array (not spline) with the same sampling rate as the photocurrent"""

        cycles, delD = self.getRunCycles(run)
        phi_ts = self.phi_ts[run][phiInd][:]

        nPulses = cycles.shape[0]
        assert(len(phi_ts) == nPulses)

        #start, end = RhO.t[0], RhO.t[0]+delD #start, end = 0.00, delD
        start, end = 0, delD
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)
        phi_tV = np.zeros_like(t)
        #pulseInds = np.empty([0,2],dtype=int) # Light on and off indexes for each pulse

        for p in range(nPulses):
            start = end
            onD, offD = cycles[p,0], cycles[p,1]
            end = start + onD + offD
            nSteps = int(round(((end-start)/dt)+1))
            tPulse = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]
            phiPulse = phi_t(tPulse) # -tPulse[0] # Align time vector to 0 for phi_t to work properly

            #onInd = len(t) - 1 # Start of on-phase
            #offInd = onInd + int(round(onD/dt))
            #pulseInds = np.vstack((pulseInds, [onInd,offInd]))

            #t = np.r_[t, tPulse[1:]]

            phi_tV = np.r_[phi_tV, phiPulse[1:]]

        phi_tV[np.ma.where(phi_tV < 0)] = 0 # Safeguard for negative phi values
        return phi_tV #, t, pulseInds


    def plot(self, plotStateVars=False):
        """Plot protocol"""
        self.Ifig = plt.figure()
        self.createLayout(self.Ifig)
        self.PD.plot(self.axI)
        self.addAnnotations()
        self.plotExtras()

        self.plotStateVars = plotStateVars

        #animateStates = True # https://jakevdp.github.io/blog/2013/05/28/a-simple-animation-the-magic-triangle/
        if self.plotStateVars:
            #RhO = self.RhO
            for run in range(self.nRuns):
                for phiInd in range(self.nPhis): #, phi in enumerate(self.phis):
                    for vInd in range(self.nVs):
                        pc = self.PD.trials[run][phiInd][vInd]
                        fileName = '{}States{}s-{}-{}-{}'.format(self.protocol, pc.nStates, run, phiInd, vInd) # RhO.nStates
                        #RhO.plotStates(pc.t, pc.states, pc.pulses, RhO.stateLabels, phi, pc.peakInds_, fileName)
                        pc.plotStates(name=fileName)
                        logging.info('Plotting states to: {}'.format(fileName))

        plt.figure(self.Ifig.number)
        plt.sca(self.axI)
        self.axI.set_xlim(self.PD.begT, self.PD.endT)
        # if addTitles:
            # figTitle = self.genTitle()
            # plt.title(figTitle) #'Photocurrent through time'

        #plt.show()
        plt.tight_layout()
        plt.show()

        figName = os.path.join(config.fDir, self.protocol+self.dataTag+"."+config.saveFigFormat)
        #externalLegend = False
        #if externalLegend:
        #    self.Ifig.savefig(figName, bbox_extra_artists=(lgd,), bbox_inches='tight', format=config.saveFigFormat) # Use this to save figures when legend is beside the plot
        #else:
        self.Ifig.savefig(figName, format=config.saveFigFormat)

        return


    def createLayout(self, Ifig=None, vInd=0):
        """Create axes for protocols with multiple subplots"""
        if Ifig == None:
            Ifig = plt.figure()

        self.addStimulus = config.addStimulus
        #phi_ts = self.genPlottingStimuli()

        # Default layout
        self.axI = Ifig.add_subplot(111)
        plt.sca(self.axI)
        #plotLight(self.pulses, self.axI)

    # TODO: Refactor multiple getLineProps
    def getLineProps(self, run, vInd, phiInd):

        colours = config.colours
        styles = config.styles
        
        if config.verbose > 1 and (self.nRuns>len(colours) or len(self.phis)>len(colours) or len(self.Vs)>len(colours)):
            warnings.warn("Warning: only {} line colours are available!".format(len(colours)))
        if config.verbose > 0 and self.nRuns>1 and len(self.phis)>1 and len(self.Vs)>1:
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


    def plotExtras(self):
        pass

    def addAnnotations(self):
        pass


    def plotStimulus(self, phi_ts, begT, pulses, endT, ax=None, light='shade', col=None, style=None):

        nPulses = pulses.shape[0]
        assert(nPulses == len(phi_ts))
        nPoints = 10*int(round(endT-begT/self.dt))+1
        t = np.linspace(begT, endT, nPoints)

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

        plt.xlabel(r'$\mathrm{Time\ [ms]}$')
        plt.xlim((begT, endT))
        plt.ylabel(r'$\mathrm{\phi\ [ph./mm^{2}/s]}$')

        return ax







class protCustom(Protocol):

    # Class attributes
    protocol = 'custom'
    squarePulse = False
    #custPulseGenerator = None
    phi_ft = None

    def extraPrep(self):
        '''Function to set-up additional variables and make parameters consistent after any changes'''
        self.nRuns = 1 #nRuns ### Reconsider this...

        #self.custPulseGenerator = self.phi_ft
        if not hasattr(self, 'phi_ts') or self.phi_ts is None:
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
        """Function to set-up additional variables and make parameters consistent after any changes"""
        self.nRuns = 1
        self.phi_ts = self.genPulseSet()

    def addAnnotations(self):
        self.axI.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        self.axI.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        self.axI.grid(b=True, which='minor', axis='both', linewidth=.2)
        self.axI.grid(b=True, which='major', axis='both', linewidth=1)


class protSinusoid(Protocol):

    protocol = 'sinusoid'
    squarePulse = False

    startOn = False
    phi0 = 0

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
            warnings.warn('Unexpected data type for phi0: ', type(self.phi0))

        assert(len(self.phi0) == self.nRuns)

        self.begT, self.endT = 0, self.totT
        self.phi_ts = self.genPulseSet()
        self.runLabels = [r'$f={}\mathrm{{Hz}}$ '.format(round_sig(f,3)) for f in self.fs]

    def getShortestPeriod(self):
        return 1000/self.sr # dt [ms]

    def genPulse(self, run, phi, pulse):
        pStart, pEnd = pulse
        onD = pEnd - pStart
        t = np.linspace(0.0, onD, int(round((onD*self.sr/1000))+1), endpoint=True) # Create smooth series of time points to interpolate between
        if self.startOn: # Generalise to phase offset
            phi_t = spline(pStart + t, self.phi0[run] + 0.5*phi*(1+np.cos(self.ws[run]*t)), ext=1, k=5)
        else:
            phi_t = spline(pStart + t, self.phi0[run] + 0.5*phi*(1-np.cos(self.ws[run]*t)), ext=1, k=5)

        return phi_t

    def createLayout(self, Ifig=None, vInd=0):

        if Ifig == None:
            Ifig = plt.figure()

        self.addStimulus = config.addStimulus

        if self.nRuns > 1: #len(phis) > 1:
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
                        intIp = spline(self.fs, Ipeaks, k=splineOrder)
                        #nPoints = 10*int(round(abs(np.log10(self.fs[-1])-np.log10(self.fs[0]))+1))
                        fsmooth = np.logspace(np.log10(self.fs[0]), np.log10(self.fs[-1]), num=1001)
                        self.axIp.plot(fsmooth, intIp(fsmooth))
                    except:
                        if config.verbose > 0:
                            print('Unable to plot spline for current peaks!')
                    fstar_p = self.fs[np.argmax(Ipeaks)]
                    fstars[phiInd,vInd] = fstar_p
                    Ap = max(Ipeaks)
                    #fpLabel = r'$f^*_{{peak}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_p,3))
                    self.axIp.plot(fstar_p, Ap, '*', markersize=10)
                    #axIp.annotate(fpLabel, xy=(fstar_p,Ap), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})

            self.axIp.set_xscale('log')
            self.axIp.set_ylabel(r'$|A|_{peak}$ $\mathrm{[nA]}$')
            if config.addTitles:
                #self.axIp.set_title('$\mathrm{|Amplitude|_{peak}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
                self.axIp.set_title(r'$f^*:=arg\,max_f(|A|_{peak})$')
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
                        intIss = spline(self.fs, Iabs, k=splineOrder)
                        #fsmooth = np.logspace(self.fs[0], self.fs[-1], 100)
                        self.axIss.plot(fsmooth, intIss(fsmooth))
                    except:
                        if config.verbose > 0:
                            print('Unable to plot spline for current steady-states!')
                    fstar_abs = self.fs[np.argmax(Iabs)]
                    fstars[phiInd,vInd] = fstar_abs
                    Aabs = max(Iabs)
                    fabsLabel = r'$f^*_{{res}}={}$ $\mathrm{{[Hz]}}$'.format(round_sig(fstar_abs,3))
                    self.axIss.plot(fstar_abs, Aabs, '*', markersize=10, label=fabsLabel)
                    self.axIss.legend(loc='best')
                    #axIss.annotate(fabsLabel, xy=(fstar_abs,Aabs), xytext=(0.7, 0.9), textcoords='axes fraction', arrowprops={'arrowstyle':'->','color':'black'})
                    if config.verbose > 0:
                        print('Resonant frequency (phi={}; V={}) = {} Hz'.format(phiOn, V, fstar_abs))
            self.axIss.set_xscale('log')
            self.axIss.set_xlabel(r'$f$ $\mathrm{[Hz]}$')
            self.axIss.set_ylabel(r'$|A|_{ss}$ $\mathrm{[nA]}$')
            if config.addTitles:
                #axIss.set_title('$\mathrm{|Amplitude|_{ss}\ vs.\ frequency}.\ f^*:=arg\,max_f(|A|)$')
                self.axIss.set_title(r'$f^*:=arg\,max_f(|A|_{ss})$')

            plt.tight_layout()

            self.fstars = fstars
            if len(self.phis) > 1: # Multiple light amplitudes
                #for i, phi0 in enumerate(self.phi0):
                fstarAfig = plt.figure()
                for vInd, V in enumerate(self.Vs):
                    if self.phi0[0] > 0: # phi0[r]
                        plt.plot(np.array(self.phis)/self.phi0[0], fstars[:,vInd])
                        plt.xlabel(r'$\mathrm{Modulating}\ \phi(t)/\phi_0$') #\phi_1(t)/\phi_0(t)
                    else:
                        plt.plot(np.array(self.phis), fstars[:,vInd])
                        plt.xlabel(r'$\mathrm{Modulating}\ \phi(t)$') #\phi_1(t)
                plt.xscale('log')
                plt.ylabel(r'$f^*\ \mathrm{[Hz]}$')
                if config.addTitles:
                    plt.title(r'$f^*\ vs.\ \phi_1(t).\ \mathrm{{Background\ illumination:}}\ \phi_0(t)={:.3g}$'.format(self.phi0[0]))

#TODO: Finish Dual Tone protocol
'''
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
'''


class protChirp(Protocol):
    """Sweep through a range of frequencies from f0 to fT either linearly or exponentially"""
    protocol = 'chirp'
    squarePulse = False

    f0 = 0
    fT = 0
    linear = True
    startOn = False
    phi0 = 0

    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'

        self.sr = max(10000, int(round(10 * max(self.f0, self.fT)))) # Nyquist frequency - sampling rate (10*f) >= 2*f
        self.nRuns = 1 #len(self.ws)
        #self.cycles = np.column_stack((self.onDs,self.offDs))
        #ws = 2 * np.pi * np.logspace(-4,10,num=7) # Frequencies [rads/s]

        if (1000)/self.f0 > min(self.onDs): #1/10**self.fs[0] > self.totT:
            warnings.warn('Warning: The period of the lowest frequency is longer than the stimulation time!')

        if isinstance(self.phi0, (int, float, complex)):
            self.phi0 = np.ones(self.nRuns) * self.phi0
        elif isinstance(self.phi0, (list, tuple, np.ndarray)):
            if len(self.phi0) != self.nRuns:
                self.phi0 = np.ones(self.nRuns) * self.phi0[0]
        else:
            warnings.warn('Unexpected data type for phi0: ', type(self.phi0))

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
            phi_t = spline(pStart + t, self.phi0[run] + 0.5*phi*(1+np.cos(ft*t)), ext=1, k=5)
        else:
            phi_t = spline(pStart + t, self.phi0[run] + 0.5*phi*(1-np.cos(ft*t)), ext=1, k=5)
        return phi_t

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
                onD = pEnd - pStart
                nPoints = 10*int(round(onD/self.dt))+1 # 10001
                tsmooth = np.linspace(0, onD, nPoints)

                if self.linear:
                    ft = self.f0 + (self.fT-self.f0)*(tsmooth/onD)
                else: # Exponential
                    ft = self.f0 * (self.fT/self.f0)**(tsmooth/onD)
                self.axF.plot(tsmooth+pStart, ft, 'g')
            self.axF.set_ylabel(r'$f\ \mathrm{[Hz]}$')
        else:
            self.axI = Ifig.add_subplot(111)

        #plotLight(self.pulses, self.axI)


class protRamp(Protocol):
    """Linearly increasing pulse"""

    protocol = 'ramp'
    squarePulse = False
    nRuns = 1

    phi0 = 0

    def extraPrep(self):
        'Function to set-up additional variables and make parameters consistent after any changes'
        self.nRuns = 1 #nRuns # Make len(phi_ton)?
        self.cycles = np.column_stack((self.onDs,self.offDs))
        self.phi_ts = self.genPulseSet()

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
        phi_t = spline([pStart, pEnd], [self.phi0, self.phi0+phi], k=1, ext=1)
        return phi_t


class protDelta(Protocol):
    # One very short, saturation intensity pulse e.g. 10 ns @ 100 mW*mm^-2 for wild type ChR
    # Used to calculate gbar, assuming that O(1)-->1 as onD-->0 and phi-->inf
    protocol = 'delta'
    squarePulse = True
    nRuns = 1

    onD = 0

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

        if config.verbose > 0:
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
                        label = r'$I_{{peak}} = {:.3g}\mathrm{{nA;}}\ t_{{lag}} = {:.3g}\mathrm{{ms}}$'.format(Ip, tlag)
                        plt.text(1.05*tp, 1.05*Ip, label, ha='left', va='bottom', fontsize=config.eqSize)



class protRectifier(Protocol):
    """Protocol to determine the rectification parameters of rhodopsins. Typically they are inward rectifiers where current is more easily passed into the cell than out. """
    # Iss vs Vclamp
    # http://en.wikipedia.org/wiki/Inward-rectifier_potassium_ion_channel
    protocol = 'rectifier'
    squarePulse = True
    nRuns = 1

    def extraPrep(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
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
        colours = config.colours
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
                # TODO: Set from config
                markerSize = 40
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
                if config.verbose > 0:
                    print('g(v=-70) = ', gm70)
                #g0[(Vs - E)==0] = None #(v1/v0)
                gNorm = gs / gm70 # Normalised conductance relative to V=-70

                self.axfV.scatter(Vs, gNorm, marker='x', color=colours, s=markerSize)#,linestyle=''
                if config.verbose > 1:
                    print(gm70)
                    if config.verbose > 2:
                        print(np.c_[Vs, np.asarray(Vs)-E, Iss, gs, gNorm])


                # Add equations to legend
                if self.nPhis > 1:
                    legLabels[phiInd] = fVstring + r'$,\ \phi={:.3g}$'.format(phiOn)
                else:
                    legLabels[phiInd] = fVstring

                ### Move this to fitting routines?
                # v0 = popt[0], v1 = popt[1], E = popt[2]

        #if len(phis) > 1:
        #ax.legend(legLabels, loc='best')

        #ax.spines['left'].set_position('zero')
        setCrossAxes(ax, zeroY=True)
        ax.set_xlim(min(Vs), max(Vs))

        ax.set_ylabel(r'$I_{ss}$ $\mathrm{[nA]}$')#, position=(0.95,0.8)) #plt.xlabel

        ax = self.axfV
        ax.set_ylabel(r'$f(v)$ $\mathrm{[1]}$')#, position=(0.95,0.8)) #plt.xlabel

        #ax.spines['left'].set_position('zero')
        setCrossAxes(ax, zeroY=True)
        ax.set_xlim(min(Vs), max(Vs))

        # yticks = ax.get_yticklabels()
        # ax.set_ylim(0, float(yticks[-1].get_text()))

        useLegend = True
        if useLegend:
            #ax.legend(legLabels, bbox_to_anchor=(0., 1.01, 1., .101), loc=3, mode="expand", borderaxespad=0., prop={'size':mpl.rcParams['font.size']})
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

        ax.set_xlabel(r'$V_{clamp}$ $\mathrm{[mV]}$', position=(config.xLabelPos,0), ha='right')
        #plt.xlim((min(Vs),max(Vs)))

        self.axI.grid(b=True, which='minor', axis='both', linewidth=.2)
        self.axI.grid(b=True, which='major', axis='both', linewidth=1)

        plt.tight_layout()




class protShortPulse(Protocol):
    # Vary pulse length - See Nikolic+++++2009, Fig. 2 & 9
    protocol = 'shortPulse'
    squarePulse = True
    nPulses = 1 # Fixed at 1

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
        self.phis.sort(reverse=True)
        self.Vs.sort(reverse=True)
        self.nPhis = len(self.phis)
        self.nVs = len(self.Vs)
        self.phi_ts = self.genPulseSet()
        self.runLabels = [r'$\mathrm{{Pulse}}={}\mathrm{{ms}}$ '.format(pD) for pD in self.pDs]


    def getRunCycles(self, run):
        return np.asarray([[self.onDs[run],self.offDs[run]]]), self.delDs[run]


    def createLayout(self, Ifig=None, vInd=0):

        if Ifig == None:
            Ifig = plt.figure()
        self.addStimulus = config.addStimulus
        gsPL = plt.GridSpec(2,3)
        self.axLag = Ifig.add_subplot(gsPL[0,-1])
        self.axPeak = Ifig.add_subplot(gsPL[1,-1], sharex=self.axLag)
        self.axI = Ifig.add_subplot(gsPL[:,:-1])


    def addAnnotations(self):
        # Freeze axis limits
        ymin, ymax = plt.ylim()
        pos = 0.02 * abs(ymax-ymin)
        plt.ylim(ymin, pos*(self.nRuns+1)) # Allow extra space for thick lines
        lightBarWidth = 2 * mpl.rcParams['lines.linewidth']
        peakMarkerSize = 1.5 * mpl.rcParams['lines.markersize']

        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
                    colour, style = self.getLineProps(run, vInd, phiInd)
                    PC = self.PD.trials[run][phiInd][vInd]
                    t_on, t_off = PC.pulses[0,:]

                    self.axI.hlines(y=(run+1)*pos, xmin=t_on, xmax=t_off, lw=lightBarWidth, color=colour)
                    self.axI.axvline(x=t_on, linestyle=':', c='k')
                    self.axI.axvline(x=t_off, linestyle=':', c=colour)

                    self.axI.plot(PC.tpeaks_, PC.peaks_, marker='*', ms=peakMarkerSize, c=colour)

                    ### Plot figure to show time of Ipeak vs time of light off c.f. Nikolic et al. 2009 Fig 2b
                    self.axLag.plot(self.pDs[run], PC.lags_[0], marker='*', ms=peakMarkerSize, c=colour)

                    ### Plot figure to show current peak vs time of light off c.f. Nikolic et al. 2009 Fig 2c
                    self.axPeak.plot(self.pDs[run], PC.peaks_, marker='*', ms=peakMarkerSize, c=colour)

        # axLag.axis('equal')
        tmax = max(self.pDs)*1.25
        self.axLag.plot([0,tmax], [0,tmax], ls="--", c=".3")
        self.axLag.set_xlim(0,tmax)
        self.axLag.set_ylim(0,tmax)
        self.axLag.set_ylabel(r'$\mathrm{Time\ of\ peak\ [ms]}$')
        self.axLag.set_aspect('auto')

        self.axPeak.set_xlim(0,tmax)
        self.axPeak.set_xlabel(r'$\mathrm{Pulse\ duration\ [ms]}$')
        self.axPeak.set_ylabel(r'$\mathrm{Photocurrent\ peak\ [nA]}$')


class protRecovery(Protocol):
    '''Two pulse stimulation protocol with varying inter-pulse interval to determine the dark recovery rate'''
    # Vary Inter-Pulse-Interval
    protocol = 'recovery'
    squarePulse = True
    nPulses = 2 # Fixed at 2 for this protocol

    onD = 0
    # def __next__(self):
        # if self.run >= self.nRuns:
            # raise StopIteration
        # return np.asarray[self.pulses[self.run]]

    def prepare(self):
        """Function to set-up additional variables and make parameters consistent after any changes"""
        self.IPIs = np.sort(np.asarray(self.IPIs))

        self.nRuns = len(self.IPIs)
        self.delDs = np.ones(self.nRuns)*self.delD
        self.onDs = np.ones(self.nRuns)*self.onD
        self.offDs = self.IPIs
        self.cycles = np.column_stack((self.onDs, self.offDs)) # [:,0] = on phase duration; [:,1] = off phase duration

        self.pulses, _ = cycles2times(self.cycles, self.delD)
        self.runCycles = np.zeros((self.nPulses, 2, self.nRuns))
        for run in range(self.nRuns):
            self.runCycles[:,:,run] = np.asarray([[self.onDs[run],self.offDs[run]],[self.onDs[run],self.offDs[run]]])

        self.begT = 0
        self.endT = self.totT
        IPIminD = max(self.delDs) + (2*max(self.onDs)) + max(self.IPIs)
        if self.endT < IPIminD:
            warnings.warn("Insufficient run time for all stimulation periods!")
        else:
            self.runCycles[-1,1,:] = self.totT - IPIminD

        self.IpIPI = np.zeros(self.nRuns)
        self.tpIPI = np.zeros(self.nRuns)

        if np.isscalar(self.phis):
            self.phis = np.asarray([self.phis])
        self.phis.sort(reverse=True)
        self.nPhis = len(self.phis)

        if np.isscalar(self.Vs):
            self.Vs = np.asarray([self.Vs])
        self.Vs.sort(reverse=True)
        self.nVs = len(self.Vs)

        self.phi_ts = self.genPulseSet()
        self.runLabels = [r'$\mathrm{{IPI}}={}\mathrm{{ms}}$ '.format(IPI) for IPI in self.IPIs]


    def getRunCycles(self, run):
        return self.runCycles[:,:,run], self.delDs[run]


    def finish(self, PC, RhO):
        ### Build array of second peaks
        self.PD.IPIpeaks_ = np.zeros((self.nRuns, self.nPhis, self.nVs)) #[[[None for v in range(len(Vs))] for p in range(len(phis))] for r in range(nRuns)]
        self.PD.tIPIpeaks_ = np.zeros((self.nRuns, self.nPhis, self.nVs))
        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
                    PC = self.PD.trials[run][phiInd][vInd]
                    PC.alignToPulse(pulse=0, alignPoint=2) # End of the first pulse
                    self.PD.IPIpeaks_[run][phiInd][vInd] = PC.peaks_[1]
                    self.PD.tIPIpeaks_[run][phiInd][vInd] = PC.tpeaks_[1]

        if config.verbose > 1:
            print(self.PD.tIPIpeaks_)
            print(self.PD.IPIpeaks_)


    def addAnnotations(self):

        # Freeze axis limits
        ymin, ymax = plt.ylim()
        pos = 0.02 * abs(ymax-ymin)
        plt.ylim(ymin, pos*self.nRuns)
        xmin, xmax = plt.xlim()
        plt.xlim(xmin, xmax)

        for run in range(self.nRuns):
            for phiInd in range(self.nPhis):
                for vInd in range(self.nVs):
                    col, style = self.getLineProps(run, vInd, phiInd)
                    pulses = self.PD.trials[run][phiInd][vInd].pulses
                    plt.annotate('', (pulses[0,1], (run+1)*pos), (pulses[1,0], (run+1)*pos), 
                                 arrowprops={'arrowstyle':'<->', 'color':col, 'shrinkA':0, 'shrinkB':0})

                    if run == 0:    ### Fit peak recovery
                        t_peaks, I_peaks, Ipeak0, Iss0 = getRecoveryPeaks(self.PD, phiInd, vInd, usePeakTime=True)
                        params = Parameters()
                        params.add('Gr0', value=0.002, min=0.0001, max=0.1)
                        params = fitRecovery(t_peaks, I_peaks, params, Ipeak0, Iss0, self.axI)
                        if config.verbose > 0:
                            print("tau_r0 = {} ==> G_r0 = {}".format(1/params['Gr0'].value, params['Gr0'].value))
        return



from collections import OrderedDict
protocols = OrderedDict([('step', protStep), ('delta', protDelta), ('sinusoid', protSinusoid), 
                         ('chirp', protChirp), ('ramp', protRamp), ('recovery', protRecovery), 
                         ('rectifier', protRectifier), ('shortPulse', protShortPulse), ('custom', protCustom)])

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
