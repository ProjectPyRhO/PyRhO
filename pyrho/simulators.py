# simulators.py

from pyrho.parameters import *
from pyrho.utilities import * # cycles2times, plotLight
from pyrho.loadData import * 
from pyrho.models import *
from pyrho.config import * #verbose
from pyrho import config
import numpy as np
import warnings
from os import path
import copy

#from brian2 import *

class Simulator(PyRhOobject): #object
    """Common base class for all simulators"""
    
    def __init__(self, Prot, RhO, simulator='Python'):
        self.simulator = simulator
        # Simulator is now initialised according to a particular protocol
    
    def __str__(self):
        return self.simulator
    
    def __repr__(self):
        return "<PyRhO {} Simulator object>".format(self.simulator)
    
    def checkDt(self, dt):
        """Function to compare simulator's timestep to the timestep required by the protocol"""
        if dt < self.dt:
            self.dt = dt #min(self.h.dt, dt)
            if config.verbose > 0:
                print('Time step reduced to {}ms by protocol!'.format(self.dt))
        return # self.dt
        
    def prepare(self, Prot):
        """Function to prepare the simulator according to the protocol and rhodopsin"""
        Prot.prepare()
        dt = Prot.getShortestPeriod()
        self.checkDt(dt)
        return # self.dt
        
    def run(self, verbose=verbose): 
        """Main routine to run the simulation protocol"""
        
        t0 = wallTime() #time.perf_counter()
        
        RhO = self.RhO
        #self.Sim = Sim
        Prot = self.Prot
        
        #Prot.prepare()
        self.prepare(Prot)
        
        if verbose > 0:
            print("\n================================================================================")
            print("Running '{}' protocol with {} for the {} model... ".format(Prot, self, RhO))
            print("================================================================================\n")
            print("{{nRuns={}, nPhis={}, nVs={}}}".format(Prot.nRuns, Prot.nPhis, Prot.nVs))
        
        # Change to self.PD ?
        Prot.PD = ProtocolData(Prot.protocol, Prot.nRuns, Prot.phis, Prot.Vs)
        Prot.PD.peak_ = [[[None for v in range(Prot.nVs)] for p in range(Prot.nPhis)] for r in range(Prot.nRuns)]
        Prot.PD.ss_ = [[[None for v in range(Prot.nVs)] for p in range(Prot.nPhis)] for r in range(Prot.nRuns)]
        if hasattr(Prot, 'runLabels'):
            Prot.PD.runLabels = Prot.runLabels
        
        if verbose > 1: 
            Prot.printParams()
        
        for run in range(Prot.nRuns):                   # Loop over the number of runs...       ### Place within V & phi loops to test protocols at different V & phi?
            
            cycles, delD = Prot.getRunCycles(run)
            pulses, totT = cycles2times(cycles, delD)
            
            for phiInd, phiOn in enumerate(Prot.phis):  # Loop over light intensity...
                
                if verbose > 1 and (Prot.nPhis > 1 or (run == 0 and phiInd == 0)): # len(phis)>0
                    RhO.dispRates()
                
                for vInd, V in enumerate(Prot.Vs):      # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                    
                    if Prot.squarePulse: #protocol in squarePulses: ##### Change after changing 'custom'
                        I_RhO, t, soln = self.runTrial(RhO, phiOn, V, delD, cycles, self.dt, verbose) #self.totT, 

                    else: # Arbitrary functions of time: phi(t)
                        phi_ts = Prot.phi_ts[run][phiInd][:]
                        I_RhO, t, soln = self.runTrialPhi_t(RhO, phi_ts, V, delD, cycles, Prot.totT, self.dt, verbose)
                        
                    # Save simulation results
                        # Change: Prot.pulses[run] := [[t_on1, t_off1],...]
                    #pulses = np.array([[delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD] for p in range(self.nPulses)]) #np.array([[delD,onD]])
                    
                    PC = PhotoCurrent(I_RhO, t, pulses, phiOn, V, Prot.protocol)
                    #PC.alignToTime()

                    PC.states = soln
                    #self.data[run][phiInd][vInd] = PC
                    #PD.trials.append(PC)
                    Prot.PD.trials[run][phiInd][vInd] = PC
                    Prot.PD.peak_[run][phiInd][vInd] = PC.peak_
                    Prot.PD.ss_[run][phiInd][vInd] = PC.ss_
                    
                    self.saveExtras(run, phiInd, vInd)
                    
                    # self.labels[run][phiInd][vInd] = label
                    # label = ""
                    if verbose > 1:
                        #print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]; label=<{}>'.format(run,nRuns,phiInd,len(phis),vInd,len(Vs),PC.range_[0],PC.range_[1],label))
                        print('Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]'.format(run,self.nRuns, phiInd,len(self.phis), vInd,len(self.Vs), PC.range_[0],PC.range_[1]))
        
        Prot.finish(PC, RhO) # Sim
        
        if Prot.saveData:
            Prot.dataTag = str(RhO.nStates)+"s"
            saveData(Prot.PD, Prot.protocol+Prot.dataTag)
        
        self.runTime = wallTime() - t0 #time.perf_counter() - t0
        if verbose > 0:
            print("\nFinished '{}' protocol with {} for the {} model in {:.3g}s".format(Prot, self, RhO, self.runTime))
            print("--------------------------------------------------------------------------------\n")
            
        return Prot.PD
    
    def saveExtras(self, run, phiInd, vInd):
        pass
    
    def plot(self):
        self.Prot.plot()
        self.plotExtras()
    
    """
    def plot(self, plotStateVars=False):
        
        Prot = self.Prot
        Ifig = plt.figure() #plt.figure(figsize=(figWidth, figHeight))
        Prot.createLayout(Ifig)
        #self.genLabels()
        Prot.PD.plot(Prot.axI) #self.labels... #baseline, = axI.plot(t, I_RhO, color=col, linestyle=style, label=label)

        protPulses = Prot.getProtPulses()
        
        Prot.addAnnotations()
        Prot.plotExtras()
        
        Prot.plotStateVars = plotStateVars
        
        #animateStates = True # https://jakevdp.github.io/blog/2013/05/28/a-simple-animation-the-magic-triangle/
        if Prot.plotStateVars: 
            RhO = self.RhO
            for run in range(Prot.nRuns):
                #cycles, delD = self.getRunCycles(run)
                #pulses, totT = cycles2times(cycles, delD)
                for phiInd, phi in enumerate(Prot.phis):
                    for vInd in range(Prot.nVs):
                        pc = self.PD.trials[run][phiInd][vInd]
                        fileName = '{}States{}s-{}-{}-{}'.format(self.protocol,RhO.nStates,run,phiInd,vInd)#; print(fileName)
                        RhO.plotStates(pc.t, pc.states, pc.pulses, RhO.stateLabels, phi, pc.peakInds_, fileName)
        
        plt.figure(Ifig.number)
        plt.sca(Prot.axI)
        Prot.axI.set_xlim(Prot.PD.begT, Prot.PD.endT)

        #plt.show()
        plt.tight_layout()
        
        externalLegend = False
        figName = os.path.join(fDir, self.protocol+self.dataTag+"."+config.saveFigFormat)
        #plt.figure(Ifig.number)
        if externalLegend:
            Ifig.savefig(figName, bbox_extra_artists=(lgd,), bbox_inches='tight', format=config.saveFigFormat) # Use this to save figures when legend is beside the plot
        else:
            Ifig.savefig(figName, format=config.saveFigFormat)
        
        return #Ifig.number
    """
        
    def plotExtras(self):
        pass
        
        
class simPython(Simulator):
    simulator = 'Python'
    """Class for channel level simulations with Python"""
    
    def __init__(self, Prot, RhO, params=simParams['Python']):
        self.dt = params['dt'].value
        
        self.Prot = Prot
        self.RhO = RhO

    '''
    # Add this into runTrial and fitting routines...
    def run(RhO, t):
        # P = RhO.P; Gd = RhO.Gd; Gr = RhO.Gr
        # if not RhO.useAnalyticSoln or 2*(P*Gd + P*Gr + Gd*Gr) > (P**2 + Gd**2 + Gr**2):
            # soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
        # else:
            # soln = RhO.calcSoln(t, RhO.states[-1,:])
            
        try:
            soln = RhO.calcSoln(t, RhO.states[-1,:])
        except: # Any exception e.g. NotImplementedError or ValueError
            soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
        
        # if RhO.useAnalyticSoln:
            # soln = RhO.calcSoln(t, RhO.states[-1,:])
            # if np.any(np.isnan(soln)):
                # soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
        # else:
            # soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
        return soln
    '''
    # runTrial(self, RhO, nPulses, V,phiOn,delD,onD,offD,padD,dt,verbose=verbose): #dt; a1,a3,b2,b4,I_RhO
    def runTrial(self, RhO, phiOn, V, delD, cycles, dt, verbose=verbose): 
        """
        Main routine for simulating a square pulse train
        
        
        Returns
            I_RhO   := [I_t0, I_t1, ..., I_tn]      (1 x) nSamples row vector
            t       := [t0, t1, ..., tn]            (1 x) nSamples row vector
            states  := [[S1_t0, S2_t0, ... Sk_t0],  nSamples x nStates array
                        [S1_t1, S2_t1, ... Sk_t1],
                                ...
                        [S1_tn, S2_tn, ... Sk_tn]]
        """
        
        nPulses = cycles.shape[0]
        
        if verbose > 1:
            if V is not None:
                Vstr = 'V = {:+}mV, '.format(V)
            else:
                Vstr = ''
            info = "Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: [delD={:.4g}ms".format(phiOn,Vstr,delD)
            for p in range(nPulses):
                info += "; [onD={:.4g}ms; offD={:.4g}ms]".format(cycles[p,0],cycles[p,1])
            info += "]"
            print(info)
        
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+delD
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True) #t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
            if verbose > 2:
                print("Simulating t_del = [{},{}]".format(start,end))
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(t, RhO.s0)
        else:
            soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian) #t_del #delay
        #t = t_del
        RhO.storeStates(soln[1:],t[1:])
        
        for p in range(0, nPulses):
            
            ### Light on phase
            RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            #end = start + onD
            end = start + cycles[p,0] # onD
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True)
            onInd = len(RhO.t) - 1  # Start of on-phase
            offInd = onInd + len(t) - 1 # Start of off-phase
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            # Turn on light and set transition rates
            phi = phiOn  # Light flux
            RhO.setLight(phi)
            if verbose > 1:
                print("On-phase initial conditions:{}".format(RhO.s_on))
                if verbose > 2:
                    print("Simulating t_on = [{},{}]".format(start,end))
            if RhO.useAnalyticSoln:
                soln = RhO.calcSoln(t, RhO.s_on)
            else:
                soln = odeint(RhO.solveStates, RhO.s_on, t, args=(None,), Dfun=RhO.jacobian)
            RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            ### Light off phase
            RhO.s_off = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            #end = start + offD
            end = start + cycles[p,1] # offD
            #if (p+1) == nPulses: # Add (or subtract) extra time after (during) the off phase
            #    end += padD
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True) #t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # endpoint=True
            # Turn off light and set transition rates
            phi = 0  # Light flux
            RhO.setLight(phi)
            if verbose > 1:
                print("Off-phase initial conditions:{}".format(RhO.s_off))
                if verbose > 2:
                    print("Simulating t_off = [{},{}]".format(start,end))
            if RhO.useAnalyticSoln:
                soln = RhO.calcSoln(t, RhO.s_off)
            else:
                soln = odeint(RhO.solveStates, RhO.s_off, t, args=(None,), Dfun=RhO.jacobian)
            RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            if verbose > 1:
                print('t_pulse{} = [{}, {}]'.format(p,RhO.t[onInd],RhO.t[offInd]))
            
        ### Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states,t = RhO.getStates()
        
        return I_RhO, t, states
        
    
    def runTrialPhi_t(self, RhO, phi_ts, V, delD, cycles, endT, dt, verbose=verbose):
        """Main routine for simulating a pulse train"""
        # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
        #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
        
        ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
        ### Combine this with the original runTrial
        
        nPulses = cycles.shape[0]
        assert(len(phi_ts) == nPulses)
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+delD #start, end = 0.00, delD
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
        soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:],t[1:])

        ### Stimulation phases
        for p in range(0, nPulses):
            RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            onD, offD = cycles[p,0], cycles[p,1]
            if p < nPulses - 1:
                end = start + onD + offD #totT #start + stimD
            else:
                end = endT
            
            onInd = len(RhO.t) - 1  # Start of on-phase
            #offInd = onInd + len(t) - 1 # Start of off-phase
            #offInd = onInd + int(round(stimD/dt)) # Consider rounding issues...
            offInd = onInd + int(round(onD/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]
            
            if verbose > 1:
                print("Pulse initial conditions:{}".format(RhO.s_on))
            
            if verbose > 2:
                soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian, full_output=True)
                print(out)
            else:
                soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian)
            
            RhO.storeStates(soln[1:], t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            if verbose > 1:
                print('t_pulse{} = [{}, {}]'.format(p,RhO.t[onInd],RhO.t[offInd]))
        
        ### Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states, t = RhO.getStates()
        
        return I_RhO, t, states

    
        
    def runTrialPhi_t_old(self, RhO,V,phi_t,delD,stimD,totT,dt,verbose=verbose): #dt; a1,a3,b2,b4,I_RhO
        """Main routine for simulating a pulse train"""
        # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
        #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
        
        ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
        ### Combine this with the original runTrial
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+delD #start, end = 0.00, delD
        t = np.linspace(start, end, int(round(((end-start)/dt)+1)), endpoint=True) #t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
        soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian) #t_del#delay # odeint(RhO.solveStates, RhO.s_0, t_del, Dfun=RhO.jacobian)
        #t = t_del
        RhO.storeStates(soln[1:],t[1:])

        ### Stimulation phase
        RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
        start = end
        end = totT #start + stimD
        t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) #t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
        onInd = len(RhO.t) - 1  # Start of on-phase
        #offInd = onInd + len(t) - 1 # Start of off-phase
        offInd = onInd + int(round(stimD/dt)) # Consider rounding issues...
        RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
        # Turn on light and set transition rates
        #phi = phiOn  # Light flux
        #RhO.setLight(phi)
        if verbose > 1:
            print("On-phase initial conditions:{}".format(RhO.s_on))
            #print('Pulse = [{}, {}]'.format(onInd,offInd))
        #print(RhO.pulseInd)
        
        if verbose > 1:
            soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian, full_output=True)
            print(out)
        else:
            soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
        if verbose > 1:
            print('Pulse = [{}, {}]'.format(RhO.t[onInd],RhO.t[offInd]))
        
        ### Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states,t = RhO.getStates()
        
        # if verbose:
            # print("Run={}, I_min={}, I_max={}, label={}".format(run,I_RhO.min(),I_RhO.max(),label))
        
        return I_RhO, t, states

        
        
        
        
class simNEURON(Simulator):
    
    """Class for cellular level simulations with NEURON"""
    simulator = 'NEURON'
    mechanisms = {3:'RhO3', 4:'RhO4', 6:'RhO6'}
    #mods = {3:'RhO3.mod', 4:'RhO4.mod', 6:'RhO6.mod'}
    #paramExceptions = ['useIR']
    
    def __init__(self, Prot, RhO, params=simParams['NEURON'], recInd=0): #v_init=-70, integrator='fixed'):
        
        ### Model Specification
        # Topology
        
        # Geometry
        
        # Biophysics
        
        ### Instrumentation
        # Synaptic Input
        
        # Graphing - try using NEURON's facilities
        
        ### Simulation control
        
        #RhO = selectModel(nStates=6)
        #p = Parameters()
        #p.add_many(('phis',[irrad2flux(20)]),('Vs',[-70]),('pulses',[[100.,200.]]),('totT',300),('dt',0.1))
        #Prot = protocols['step'](p)
        
        
        ### To load mod files:
        # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
        
        self.Prot = Prot
        self.RhO = RhO
        
        #import neuron as nrn
        #self.nrn = nrn
        #self.nrn.load_mechanisms(os.environ['NRN_NMODL_PATH'])
        
        from neuron import h
        #self.neuron = neuron
        self.h = h
        
        self.h.load_file('stdrun.hoc')
        
        #self.reset()
        
        ### Simulation control
        self.CVode = params['CVode'].value
        
        if self.CVode == True: #self.integrator == 'variable': #... if dt <= 0:
            self.h('objref cvode')
            self.h.cvode = self.h.CVode()
            self.h.cvode.active(1)
        else: #elif 'fixed':
            #self.h.cvode.active(0)
            self.h.dt = params['dt'].value
            self.dt = self.h.dt
        #else:
        #    raise ValueError('Unknown integrator!')
        if verbose > 0:
            print('Integrator tolerances: absolute=',self.h.cvode.atol(),' relative=',self.h.cvode.rtol()) ### Set as parameters
        #self.h.cvode.atol(0.000001)
        #self.h.cvode.rtol(0.000001)
        
        self.h.v_init = params['v_init'].value #-60
        
        #expdist = 600 ### Redundant
        
        
        ### Move into prepare() in order to insert the correct type of mechanism (continuous or discrete) according to the protocol
        
        self.mod = self.mechanisms[self.RhO.nStates] ### Use this to select the appropriate mod file for insertion ##### How to set this for continuous mods?
        if not Prot.squarePulse:
            self.mod += 'c'
        
        self.buildCell(params['cell'].value)
        self.h.topology() # Print topology
        #mod = self.mechanisms[self.RhO.nStates] ### Use this to select the appropriate mod file for insertion ##### How to set this for continuous mods?
        self.transduce(self.RhO, expProb=params['expProb'].value)
        self.setRhodopsinParams(self.rhoList, self.RhO, modelParams[str(self.RhO.nStates)])
        
        #print(recInd)
        #for sec in self.cell:
        #    print(sec)
        #print(self.cell.count())
        #print(self.rhoList)
        #print(self.rhoList.count())
        
        self.rhoRec = self.rhoList[recInd] #self.compList.o(recInd) # Choose a rhodopsin to record
        self.setRecords(self.rhoRec, params['Vcomp'].value, self.RhO)
        self.Vclamp = params['Vclamp'].value
        #if self.Vclamp == True:
        #    self.addVclamp() # params['Vhold'].value
        
        ### else skip variations in Vs for protocol...
        
        
        #if self.cell has params['Vcomp'].value # Check for existence of section
        #comp='soma'
        #memRec = self.h.Section(name=comp)
        
        
        # Run then plot
        
        
    def reset(self):
        ### Clear any previous models
        # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=2&t=2367
        # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2655137/
        # http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html
        
        pc = self.h.ParallelContext()
        pc.gid_clear()
        
        #unref all the cell objects
        self.h("forall delete_section()")
        for sec in self.h.allsec():
            #self.h("%s{delete_section()}"%sec.name())
            self.h("{}{{delete_section()}}".format(sec.name()))
        #self.h('objectref compList')
        #self.h('objectref rhoList')
        
        #unref all the NetCon
        
        #Other top level hoc objects (Vectors, Graphs etc.)
        #pc.done()
        return

    # self.dt should point to self.h.dt
    
    def prepare(self, Prot):
        """Function to prepare the simulator according to the protocol and rhodopsin"""
        Prot.prepare()
        
        if self.CVode is False:
            dt = Prot.getShortestPeriod()
            self.checkDt(dt)
            self.h.dt = self.dt
        
        # Include code to set Vclamp = False if Prot.Vs=[None]
        if self.Vclamp is False:
            Prot.Vs = [None]
            Prot.nVs = 1
        else:
            self.addVclamp()
        
        self.Vms = Prot.genContainer()
        return #self.h.dt
    
    
    def buildCell(self, scriptList=[]):
        """Pass a list of hoc or python files in the order in which they are to be processed"""
        import os
        import subprocess
        ### Check for script in cwd and then fallback to looking in NRN_NMODL_PATH
        cwd = os.getcwd()
        for script in scriptList:
            #name, ext = os.path.splitext(script)
            ext = script.split(".")[-1]
            if not os.path.isfile(os.path.join(cwd, script)):
                script = os.path.join(os.environ['NRN_NMODL_PATH'], script)
            #script = 'PyRhO/NEURON/'+script #.join(script)
            
            #if verbose > 0:
            print('Loading ',script)
            if ext == 'hoc' or ext == 'HOC':
                self.h.load_file(script)
            elif ext == 'py':
                #exec(script) # Is this safe?!
                #import script
                #os.system("""python {}""".format(script))
                subprocess.Popen(["python", script])
                #self.cell = subprocess.Popen(["python", script])
            else:
                raise ValueError('Unknown file type!')
        
        #cell = 'L5PC'
        self.cell = self.h.SectionList()
        self.nSecs = 0
        for sec in self.h.allsec():
            self.cell.append(sec)
            self.nSecs += 1
        if verbose > 0:
            print('Total sections: ',self.nSecs) #self.cell.count())
        ### Model Specification (Topology, Geometry, Biophysics)
        #h.load_file("import3d.hoc")
        #h.load_file("models/L5PCbiophys3.hoc") # Biophysics
        #h.load_file("models/L5PCtemplate.hoc")
        
        ########################################################################## Place in another script!!!!!
        #self.h('objref L5PC') #h('objref '+cell)
        #self.h.L5PC = self.h.L5PCtemplate("morphologies/cell1.asc")
        #self.cell = self.h.L5PC
        
        
        #self.cell = self.h.L5PCtemplate("morphologies/cell1.asc")
        


    def transduce(self, RhO, expProb=1):
        import random
        #print(expProb)
        ### Insert Rhodopsins
        #filename = 'ChR-in-apical.hoc'
        #self.h.load_file(filename)
        
        # Separate from protocol variables?
        #ChRexp = RhO.g/15000000
        # The following should be set during the protocol
        #Irrad = 0
        #delay = 0
        #onD = 0
        #offD = 0
        #expdist = 600
        
        self.compList = self.h.List() #self.h('new List()')
        self.rhoList = self.h.List()
        self.h('objectvar rho')
        mech = self.mechanisms[RhO.nStates]
        for sec in self.cell: #self.h.allsec(): # Loop over every section in the cell
            #print(sec)
            self.compList.append(sec)
            #print(random.random())
            #print(expProb)
            if random.random() <= expProb: # Insert a rhodopsin and append it to a rhoList
                #self.compList.append(sec)
                self.h('rho = new {}(0.5)'.format(mech))
                self.rhoList.append(self.h.rho)
        
        #print(self.compList.count())
        #print(self.rhoList.count())
        
        #self.h.ChR_in_apical(ChRexp,Irrad,delay,onD,offD,expdist) # Sets Er, del, ton, toff, num, gbar # 0,108
        # ChR_in_basal 0,83
        # ChR_in_axon 0,1
        # ChR_in_soma 0,0
        #self.compList = self.h.ChR_apic
        
        # proc ChR_in_apical() {
            # ChRexp = $1
            # pd = $2         // 10 mW/mm^2
            # del_light = $3  // 650 ms  delay time before light ON
            # t_on = $4       // 300 ms  
            # t_off = $5       // 300 ms  
            # expdist = $6
            # ChR_apic = new List()
            # ChR2apic=ChRexp  
            # i=0
            # print "ChR distribution [index, section, segment, area, gbar, distance]:"  
            # for k=0,108 {
                # L5PC.apic[k] {
                    # print x
                    # for(x,0){
                    # //if(distance(x)>600){
                        # ChR_apic.append(new ChR(x))
                        # ChR_apic.o(i).gbar= ChR2apic*area(x)
                        # ChR_apic.o(i).Er=pd
                        # ChR_apic.o(i).del=del_light
                        # ChR_apic.o(i).ton=t_on
                        # ChR_apic.o(i).toff=t_off
                        # ChR_apic.o(i).num=n_pulses
                        # i=i+1      
                        # //}
                    # }
                # }
            # }
            # print "Number of segments populated (ChR, apical): ",i
        # }

    ### Set Rhodopsin parameters
    def setRhodopsinParams(self, rhoList, RhO, pSet): #compList
        # for i in range(int(compList.count())): #377
            # #for p in compList.o(i).__dict__: ### Make the mod files have exactly the same parameters as the models
            # compList.o(i).e1 = RhO.a10
            # compList.o(i).b10 = RhO.b1
            # compList.o(i).a2dark = RhO.a30
            # compList.o(i).a2light = RhO.a31
            # compList.o(i).b2dark = RhO.b20
            # compList.o(i).b2light = RhO.b21
            # compList.o(i).a30 = RhO.a3
            # compList.o(i).e3 = RhO.b40
            # compList.o(i).a40 = RhO.a4
            # compList.o(i).gama = RhO.gam
            
        for rho in rhoList: #range(int(rhoList.count())): # not self.rhoList so that subsets can be passed
            for p in pSet: #modelParams[str(RhO.nStates)]:
                #if p not in self.paramExceptions:
                setattr(rho, p, pSet[p].value)
            
            ### Should this be set in transduce?
            #compList.o(i).gbar = RhO.g/150000000 # ChRexp
            
            ### These need to be made accessible in the mod
            #h.ChR_apic.o(i).Er0 = flux2irrad(RhO.phi0)
            #h.ChR_apic.o(i).v0 = RhO.v0
            #h.ChR_apic.o(i).v1 = RhO.v1
            # Add e_rev
    
    def getRhodopsinParams():
        # Use nrnpython to set Python variables from within hoc?
        pass
    
    
    def setRecords(self,rhoRec,Vcomp,RhO): #memRec
    
        #comp='soma'
        #memRec = self.h.Section(name=comp)
        
        self.h('objref Vm')
        self.h.Vm = self.h.Vector()
        #self.h.Vm.record(memRec(0.5)._ref_v) #cvode.record(&v(0.5),vsoma,tvec)
        self.h('Vm.record(&{}.v(0.5))'.format(Vcomp))
        
        #self.h('objref Vsoma') #objref vsoma, i_ChR, tvec
        self.h('objref Iphi')
        self.h('objref tvec')
        #self.h.Vsoma = self.h.Vector() #vsoma = new Vector()
        self.h.Iphi = self.h.Vector() #i_ChR = new Vector()
        self.h.tvec = self.h.Vector() #tvec = new Vector()
        
        ### Record time vector (since the solver uses variable time steps
        #h.Vsoma.record(soma(0.5)._ref_v) #cvode.record(&v(0.5),vsoma,tvec)
        #h('cvode.record(&rho.iChR,Iphi,tvec)')
        self.h.tvec.record(self.h._ref_t) # Record time points
        
        #h('Iphi.record(h.ChR_apic.o(77).iChR)')
        #h.Iphi.record(h.ref(rhoRec.iChR)) # Only records initial value
        #h.Iphi.record(h.ChR_apic.o(7)._ref_iChR) # Works!
        self.h.Iphi.record(rhoRec._ref_i)
        
        # h.setpointer(_ref_hocvar, 'POINTER_name', point_proces_object)
        # h.setpointer(_ref_hocvar, 'POINTER_name', nrn.Mechanism_object)
        # For example if a mechanism with suffix foo has a POINTER bar and you want it to point to t use
        # h.setpointer(_ref_t, 'bar', sec(x).foo)
        
        
        ### Save state variables according to RhO model
        #states = getStateVars(RhO)
        for s in RhO.stateVars:
            self.h('objref {}Vec'.format(s))
            #self.h('objref tmpVec')
            self.h.rho = rhoRec                             ########### Check this works with multiple sections/rhodopsins
            self.h('{}Vec = new Vector()'.format(s))
            #self.h('{}Vec.record(rhoRec._ref_{})'.format(s,s))
            #self.h('{}Vec.record(&rhoRec.{})'.format(s,s))
            #self.h.{}Vec.record(rhoRec._ref_{})
            #self.h('{}Vec.record(rho, &{})'.format(s,s))
            self.h('{}Vec.record(&rho.{})'.format(s,s))
            #print(s,':',self.h('{}Vec.printf()'.format(s)))
        
        # vec['v_pre'].record(pre(0.5)._ref_v)
        # vec['v_post'].record(post(0.5)._ref_v)
        # vec['i_syn'].record(syn._ref_i)
        # vec['t'].record(h._ref_t)
        
        ### Record cell membrane potential from the soma
        #self.h('access L5PC.soma')
        #self.h('cvode.record(&v(0.5),Vsoma,tvec)') # Works
        
        #h.L5PC.soma[0].push()
        #h.Vsoma.record(h.L5PC.soma[0](0.5)._ref_v)
        #h.pop_section()
    
    
    def addVclamp(self):
        #Vhold=-70
        self.h('objref Vcl')
        self.h.Vcl = self.h.SEClamp(0.5)#('new SEClamp(.5)') # voltage-clamp at Vcl.amp1
        #self.h.Vcl.dur1 = self.h.tstop #-10
        #self.h.Vcl.amp1 = Vhold
    
    def setVclamp(self, Vhold=-70):
        self.h.Vcl.dur1 = self.h.tstop
        self.h.Vcl.amp1 = Vhold
    
    def setPulses(self, phiOn, delD, onD, offD, nPulses):
        # for i in range(int(self.rhoList.count())): #377
            # #for p in compList.o(i).__dict__: ### Make the mod files have exactly the same parameters as the models
            # self.rhoList[i].phiOn = phiOn #self.compList.o(i).Er = flux2irrad(phiOn)
            # self.rhoList[i].delD = delD #setattr(self.compList.o(i), 'del', delD)
            # self.rhoList[i].onD = onD #self.compList.o(i).ton = onD
            # self.rhoList[i].offD = offD #self.compList.o(i).toff = offD
            # self.rhoList[i].nPulses = nPulses #self.compList.o(i).num = nPulses
            # # Add V
        for rho in self.rhoList:
            rho.phiOn = phiOn
            rho.delD = delD
            rho.onD = onD
            rho.offD = offD
            rho.nPulses = nPulses
    
    # runTrial(self, RhO, nPulses, V, phiOn, delD, onD, offD, padD, dt, verbose=verbose):
    def runTrial(self, RhO, phiOn, V, delD, cycles, dt, verbose=verbose): 
    
        # Notes on the integrator
        # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=8&t=1330
        # http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html
        
        # Model events: http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html#modeldescriptionissues-events
        
                
        nPulses = cycles.shape[0]
        times, totT = cycles2times(cycles,delD)
        self.times = times
        
        # phi_ts = self.Prot.phi_ts[0][0][:]#[run][phiInd][:]
        # self.runTrialPhi_t(self, RhO, self.Prot.genPhiFuncs(), V, delD, cycles, totT, dt)
        # return
        
        
        RhO.initStates(phi=0) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial states for plotting
        
        # Set simulation run time
        self.h.tstop = totT #delD + np.sum(cycles) #nPulses*(onD+offD) + padD
        
        if self.Vclamp == True:
            self.setVclamp(V)
        
        fixedPulses = False # True
        
        if fixedPulses:
            for p in range(1, nPulses):
                if (cycles[p,0] != cycles[p-1,0]) or cycles[p,1] != cycles[p-1,1]:
                    warnings.warn('Warning: Pulse cycles must be identical for NEURON simulations - replicating pulse!')
                    #self.Prot.squarePulse = False
                    #self.runTrialPhi_t(self, RhO, self.Prot.genPhiFuncs(), V, delD, cycles, self.Prot.totT, dt)
                    #return
            onD, offD, padD = cycles[0,0], cycles[0,1], 0   # HACK!!! Use first pulse timing only...
            
            if verbose > 0:
                Vstr = '' if V is None else 'V = {:+}mV, '.format(V)
                info = "Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: [delD={:.4g}ms; onD={:.4g}ms; offD={:.4g}ms]".format(phiOn,Vstr,delD,onD,offD+padD)
            
            self.setPulses(phiOn, delD, onD, offD, nPulses)
            self.h.init() #self.neuron.init()
            
            #self.h.finitialize()
            
            self.h.run() #self.neuron.run(self.h.tstop)
        
        else: # Work in progress to simulate variable pulses
            
            # progress = delD
            # padD = 0
            # for p in range(nPulses):
                # delay = delD if p==0 else 0
                # onD, offD = cycles[p]
                # self.setPulses(phiOn, delay, onD, offD, 1)
                # progress += (onD + offD)
                # #tstop = delay + np.sum(cycles[p])
                # self.h.tstop = progress # tstop #delay + np.sum(cycles[p])
                # #while self.h.t<tstop:
                # #    self.h.fadvance()
                # self.h.run() 
        

            padD = 0
            progress = delD
            self.h.init()
            
            if verbose > 0:
                Vstr = '' if V is None else 'V = {:+}mV, '.format(V)
                print("Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: delD={:.4g}ms;".format(phiOn,Vstr,delD))
            
            #t = np.asarray([])
            #I_RhO = np.asarray([])
            #self.Vm = np.asarray([])
            firstRun = True
            p = 0
            while p < nPulses:
                delay = delD if p==0 else 0
                identPulses = 1
                onD, offD = cycles[p]
                #prevOffD = offD
                for pCheck in range(p+1, nPulses):
                    if (cycles[pCheck,0] == onD) and (cycles[pCheck,1] == offD):
                        identPulses += 1
                        p += 1
                    else:
                        break
                #elapsed = progress
                progress += ((onD + offD) * identPulses)
                if verbose > 0:
                    print(" [onD={:.4g}ms; offD={:.4g}ms] x {}".format(onD, offD, identPulses))
                
                
                self.setPulses(phiOn, delay, onD, offD, nPulses)#identPulses+1)
                
                #self.h.run()
                #while self.h.t < delay + (onD + offD) * identPulses:
                #    self.h.fadvance()
                if firstRun:
                    #self.setPulses(phiOn, delay, onD, offD, identPulses+1)
                    self.h.tstop = progress #- offD
                    self.h.run()
                    firstRun = False
                else:
                    #self.setPulses(phiOn, prevOffD, onD, offD, identPulses+1)
                    self.h.continuerun(progress) # - offD)
                p+=1
                # for rho in self.rhoList:
                    # print(rho.phiOn)
                    # print(rho.delD)
                    # print(rho.onD)
                    # print(rho.offD)
                    # print(rho.nPulses)
            
            #I_RhO = np.r_[I_RhO, np.array(self.h.Iphi.to_python(), copy=True)]
            #t = np.r_[t, np.array(self.h.tvec.to_python(), copy=True)]
            #self.Vm = np.r_[self.Vm, np.array(self.h.Vm.to_python(), copy=True)]
            #self.t = t #np.r_[self.t, t]
            #print(t)
        
        I_RhO = np.array(self.h.Iphi.to_python(), copy=True)
        t = np.array(self.h.tvec.to_python(), copy=True)
        self.Vm = np.array(self.h.Vm.to_python(), copy=True)
        self.t = t
        
        #print(t)
        #print(self.Vm)
        
        for p in range(nPulses):
            ##onInd = np.searchsorted(t,delD+p*(onD+offD),side="left")
            #onInd = np.searchsorted(t,times[p,0],side="left")
            ##offInd = np.searchsorted(t,delD+p*(onD+offD)+onD,side="left")
            #offInd = np.searchsorted(t,times[p,1],side="left")
            #RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            #[onInd,offInd] = np.searchsorted(t,times[p,:],side="left")
            pInds = np.searchsorted(t,times[p,:],side="left") # Combine searches!
            RhO.pulseInd = np.vstack((RhO.pulseInd,pInds))
        #print(RhO.pulseInd)
        
        ### Get solution variables
        soln = np.zeros((len(t),RhO.nStates)) ################ Hack!!!!!
        for sInd, s in enumerate(RhO.stateVars):
            #print(s,':',self.h('{}Vec.printf()'.format(s)))
            self.h('objref tmpVec')
            self.h('tmpVec = {}Vec'.format(s))
            #soln[:,sInd] = np.array(self.h('{}Vec.to_python()'.format(s)), copy=True)
            soln[:,sInd] = np.array(self.h.tmpVec.to_python(), copy=True)
            self.h('{}Vec.clear()'.format(s))
            self.h('{}Vec.resize(0)'.format(s))
        
        # if self.Vclamp is False:
            # self.plotVm(times)
        
        # Reset vectors
        self.h.tvec.clear() # Necessary?
        self.h.tvec.resize(0)
        self.h.Iphi.clear()
        self.h.Iphi.resize(0)
        self.h.Vm.clear()
        self.h.Vm.resize(0)
        
        
        return I_RhO, t, soln


    def runTrialPhi_t(self, RhO, phi_ts, V, delD, cycles, endT, dt, verbose=verbose): 
        # Make RhO.mod a continuous function of light - 9.12: Discontinuities 9.13: Time-dependent parameters p263
            #- Did cubic spline interpolation ever get implemented for vec.play(&rangevar, tvec, 1)?
            #- http://www.neuron.yale.edu/neuron/faq#vecplay
            #- http://www.neuron.yale.edu/neuron/static/new_doc/programming/math/vector.html#Vector.play
            #- Ramp example: http://www.neuron.yale.edu/phpbb/viewtopic.php?f=15&t=2602
        #raise NotImplementedError('Error: Light as a continuous function of time has not been implemented for NEURON yet!')
        
        # vsrc.play(&var, Dt)
        # vsrc.play(&var, tvec)
        # vsrc.play("stmt involving $1", optional Dt or tvec arg)
        # vsrc.play(index)
        # vsrc.play(&var or stmt, tvec, continuous)
        # vsrc.play(&var or stmt, tvec, indices_of_discontinuities_vector)
        # vsrc.play(point_process_object, &var, ...)
        
        # The vsrc vector values are assigned to the "var" variable during a simulation.
        # The same vector can be played into different variables.

        # If the "stmt involving $1" form is used, that statement is executed with the appropriate value of the $1 arg. This is not as efficient as the pointer form but is useful for playing a value into a set of variables as in "forall g_pas = $1"
        # The index form immediately sets the var (or executes the stmt) with the value of vsrc.x[index]

        # v.x[i] -> var(t) where t(i) is Dt*i or tvec.x[i]
        
        # The discrete event delivery system is used to determine the precise time at which values are copied from vsrc to var. 
        # Note that for variable step methods, unless continuity is specifically requested, the function is a step function. 
        # Also, for the local variable dt method, var MUST be associated with the cell that contains the currently accessed section (but see the paragraph below about the use of a point_process_object inserted as the first arg).
        # For the fixed step method transfers take place on entry to finitialize() and on entry to fadvance(). 
        # At the beginning of finitialize(), var = v.x[0]. On fadvance() a transfer will take place if t will be (after the fadvance increment) equal or greater than the associated time of the next index. For the variable step methods, transfers take place exactly at the times specified by the Dt or tvec arguments.
        
        # If the end of the vector is reached, no further transfers are made (var becomes constant)
        # c.f. UnivariateSpline ext=3
        
        # Note well: for the fixed step method, if fadvance exits with time equal to t (ie enters at time t-dt), then on entry to fadvance, var is set equal to the value of the vector at the index appropriate to time t. Execute tests/nrniv/vrecord.hoc to see what this implies during a simulation. ie the value of var from t-dt to t played into by a vector is equal to the value of the vector at index(t). If the vector was meant to serve as a continuous stimulus function, this results in a first order correct simulation with respect to dt. If a second order correct simulation is desired, it is necessary (though perhaps not sufficient since all other equations in the system must also be solved using methods at least second order correct) to fill the vector with function values at f((i-.5)*dt).
        
        # When continuous is 1 then linear interpolation is used to define the values between time points. # [c.f. UnivariateSpline k=1]
        # However, events at each Dt or tvec are still used and that has beneficial performance implications for variable step methods since vsrc is equivalent to a piecewise linear function and variable step methods can excessively reduce dt as one approaches a discontinuity in the first derivative. 
        # Note that if there are discontinuities in the function itself, then tvec should have adjacent elements with the same time value. As of version 6.2, when a value is greater than the range of the t vector, linear extrapolation of the last two points is used instead of a constant last value. 
        # c.f. UnivariateSpline ext=0
        # If a constant outside the range is desired, make sure the last two points have the same y value and have different t values (if the last two values are at the same time, the constant average will be returned). (note: the 6.2 change allows greater variable time step efficiency as one approaches discontinuities.)
        
        # The indices_of_discontinuities_vector argument is used to specifying the indices in tvec of the times at which discrete events should be used to notify that a discontinuity in the function, or any derivative of the function, occurs. Presently, linear interpolation is used to determine var(t) in the interval between these discontinuities (instead of cubic spline) so the length of steps used by variable step methods near the breakpoints depends on the details of how the parameter being played into affects the states.

        # For the local variable timestep method, CVode.use_local_dt() and/or multiple threads, ParallelContext.nthread() , it is often helpful to provide specific information about which cell the var pointer is associated with by inserting as the first arg some POINT_PROCESS object which is located on the cell. This is necessary if the pointer is not a RANGE variable and is much more efficient if it is. The fixed step and global variable time step method do not need or use this information for the local step method but will use it for multiple threads. It is therefore a good idea to supply it if possible.
        
        # vec.x[i] to access values
        
        """Main routine for simulating a pulse train"""
        # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
        #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
        
        ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
        ### Combine this with the original runTrial
        
        
        
        if self.Vclamp == True:
            self.setVclamp(V)
        
        nPulses = cycles.shape[0]
        assert(len(phi_ts) == nPulses)
        
        times, totT = cycles2times(cycles,delD)
        self.times = times
        
        # Set simulation run time
        self.h.tstop = totT #delD + np.sum(cycles) #nPulses*(onD+offD) + padD
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        
        start, end = RhO.t[0], RhO.t[0]+delD #start, end = 0.00, delD
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)
        #t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) # Time vector
        phi_tV = np.zeros_like(t)
        
        discontinuities = np.asarray([len(t) - 1]) # -1?
        for p in range(nPulses):
            start = end
            onD, offD = cycles[p,0], cycles[p,1]
            if p < nPulses - 1:
                end = start + onD + offD
            else:
                end = endT
            nSteps = int(round(((end-start)/dt)+1))
            tPulse = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]
            phiPulse = phi_t(tPulse) # -tPulse[0] # Align time vector to 0 for phi_t to work properly
            discontinuities = np.r_[discontinuities, len(tPulse) - 1]
            
            #onInd = len(RhO.t) - 1 # Start of on-phase
            onInd = len(t) - 1 # Start of on-phase
            offInd = onInd + int(round(onD/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            
            t = np.r_[t, tPulse[1:]]
            phi_tV = np.r_[phi_tV, phiPulse[1:]]
            
            #RhO.storeStates(soln[1:], t[1:])
        
        tVec = self.h.Vector(t) # tVec.from_python(t)
        tVec.label('Time [ms]')
        # phi_t = InterpolatedUnivariateSpline([pStart, pEnd], [self.phi_ton, phi], k=1, ext=1) # Ramp
        phiVec = self.h.Vector(phi_tV) #phi_t(t)
        phiVec.label('phi [ph./mm^2/s]')
        
        
        #phiVec.play(&phi, tVec, continuous) # indices_of_discontinuities_vector
        phiVec.play(self.rhoRec._ref_phi, tVec, 1, discontinuities) # continuous, indices_of_discontinuities_vector
        
        #phiVec.play_remove()
        
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
        #soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian)
        #RhO.storeStates(soln[1:],t[1:])
        
        #self.setPulses(phiOn, delD, onD, offD, nPulses)
        self.h.init() #self.neuron.init()
        #self.h.finitialize()
        
        self.h.run() #self.neuron.run(self.h.tstop)
        
        I_RhO = np.array(self.h.Iphi.to_python(), copy=True)
        t = np.array(self.h.tvec.to_python(), copy=True)
        self.Vm = np.array(self.h.Vm.to_python(), copy=True)
        self.t = t
        
        ### Get solution variables
        soln = np.zeros((len(t),RhO.nStates)) ################ Hack!!!!!
        for sInd, s in enumerate(RhO.stateVars):
            #print(s,':',self.h('{}Vec.printf()'.format(s)))
            self.h('objref tmpVec')
            self.h('tmpVec = {}Vec'.format(s))
            #soln[:,sInd] = np.array(self.h('{}Vec.to_python()'.format(s)), copy=True)
            soln[:,sInd] = np.array(self.h.tmpVec.to_python(), copy=True)
            self.h('{}Vec.clear()'.format(s))
            self.h('{}Vec.resize(0)'.format(s))
        
        # if self.Vclamp is False:
            # self.plotVm(times)
        
        # Reset vectors
        self.h.tvec.clear() # Necessary?
        self.h.tvec.resize(0)
        #self.h.tVec.clear() # Necessary?
        #self.h.tVec.resize(0)
        self.h.Iphi.clear()
        self.h.Iphi.resize(0)
        self.h.Vm.clear()
        self.h.Vm.resize(0)
        # phiVec ?
        
        RhO.storeStates(soln[1:],t[1:])
        
        '''
        ### Stimulation phases
        for p in range(0, nPulses):
            RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            onD, offD = cycles[p,0], cycles[p,1]
            if p < nPulses - 1:
                end = start + onD + offD #totT #start + stimD
            else:
                end = endT
            
            onInd = len(RhO.t) - 1  # Start of on-phase
            #offInd = onInd + len(t) - 1 # Start of off-phase
            #offInd = onInd + int(round(stimD/dt)) # Consider rounding issues...
            offInd = onInd + int(round(onD/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            
            t = np.linspace(start, end, int(round(((end-start)/dt)+1)), endpoint=True)
            phi_t = phi_ts[p]
            
            tVec = self.h.Vector(t) # tVec.from_python(t)
            tVec.label('Time [ms]')
            # phi_t = InterpolatedUnivariateSpline([pStart, pEnd], [self.phi_ton, phi], k=1, ext=1) # Ramp
            phiVec = self.h.Vector(phi_t(t))
            phiVec.label('phi [ph./mm^2/s]')
            
            #discontinuities = np.asarray([])
            #phiVec.play(&phi, tVec, continuous) # indices_of_discontinuities_vector
            phiVec.play(rhoRec._ref_phi, tVec, 1) # continuous, indices_of_discontinuities_vector
        
            if verbose > 1:
                print("Pulse initial conditions:{}".format(RhO.s_on))
                
            
            if verbose > 2:
                soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian, full_output=True)
                print(out)
            else:
                soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian)
            
            RhO.storeStates(soln[1:], t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            if verbose > 1:
                print('t_pulse{} = [{}, {}]'.format(p,RhO.t[onInd],RhO.t[offInd]))
        '''
        
        ### Calculate photocurrent
        #I_RhO = RhO.calcI(V, RhO.states) # Recorded directly
        states, t = RhO.getStates()
        
        return I_RhO, t, soln
    
    def saveExtras(self, run, phiInd, vInd):
        ### HACK!!!
        self.Vms[run][phiInd][vInd] = copy.copy(self.Vm)
        return
    
    def plotExtras(self): ### REVISE!!!
        Prot = self.Prot
        RhO = self.RhO
        if not self.Vclamp:
            Vfig = plt.figure()
            axV = Vfig.add_subplot(111)
            for run in range(Prot.nRuns):                   # Loop over the number of runs...   ### Place within V & phi loops to test protocols at different V & phi?
                cycles, delD = Prot.getRunCycles(run)
                pulses, totT = cycles2times(cycles, delD)
                for phiInd, phiOn in enumerate(Prot.phis):  # Loop over light intensity...
                    for vInd, V in enumerate(Prot.Vs):      # Loop over clamp voltage ### N.B. solution variables are not currently dependent on V
                        #Vfig = self.plotVm(pulses)
                        col, style = Prot.getLineProps(run, vInd, phiInd)
                        Vm = self.Vms[run][phiInd][vInd]
                        t = self.t - delD ### HACK Change me!!!
                        plt.plot(t, Vm, color=col, ls=style) #plt.plot(t, Vm, color='g')
                        #if times is not None:
                        #    plotLight(times)
                        if run == 0 and phiInd == 0 and vInd == 0:
                            plotLight(pulses-delD) ### HACK
                            
            plt.ylabel('$\mathrm{Membrane\ Potential\ [mV]}$') #axV.set_ylabel('Voltage [mV]')
            plt.xlim((-delD, self.h.tstop - delD)) ### HACK
            #plt.xlabel('$\mathrm{Time\ [ms]}$')
            plt.xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
            #axI.set_xlim(0,tstop)
            
            axV.spines['right'].set_color('none')
            axV.spines['bottom'].set_position('zero') # x-axis
            axV.spines['top'].set_color('none')
            axV.spines['left'].set_smart_bounds(True)
            axV.spines['bottom'].set_smart_bounds(True)
            axV.xaxis.set_ticks_position('bottom')
            axV.yaxis.set_ticks_position('left')
            
            axV.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
            axV.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
            axV.grid(b=True, which='minor', axis='both', linewidth=.2)
            axV.grid(b=True, which='major', axis='both', linewidth=1)
            
            plt.legend(Prot.PD.legLabels)
            #print(axV.get_xlim())
            #print(axV.get_ylim())
            #axV.set_ylim(axV.get_ylim())
            ymin, ymax = axV.get_ylim()
            if ymax < 0: ### HACK to fix matplotlib bug
                axV.set_ylim((ymin,0))
            plt.tight_layout() ### Bugged when min and max are negative?
            #figName = '{}Vm{}s-{}-{}-{}'.format(Prot.protocol,RhO.nStates,run,phiInd,vInd)
            figName = '{}Vm{}s'.format(Prot.protocol,RhO.nStates)
            fileName = os.path.join(fDir, figName+"."+config.saveFigFormat)
            Vfig.savefig(fileName, format=config.saveFigFormat)
    
    def plotVm(self, times=None):

        # Count number of spikes
        # Find other properties e.g. spike delay, duration
        
        Vfig = plt.figure()
        axV = Vfig.add_subplot(111)
        #axI = fig.add_subplot(111)
        #axV = axI.twinx()
        #axI.plot(h.tvec,h.Iphi,color='b')
        #axI.set_ylabel('Current [nA]')
        #print(self.h.tvec.to_python())
        #print(self.h.Vm.to_python())
        #plt.plot(self.h.tvec.to_python(), self.h.Vm.to_python(), color='g') #axV.plot(h.tvec,h.Vm,color='g')
        plt.plot(self.t, self.Vm, color='g')
        plt.ylabel('$\mathrm{Membrane\ Potential\ [mV]}$') #axV.set_ylabel('Voltage [mV]')
        if times is not None:
            plotLight(times)
        plt.xlim((0, self.h.tstop))
        #plt.xlabel('$\mathrm{Time\ [ms]}$')
        plt.xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
        #axI.set_xlim(0,tstop)
        
        axV.spines['right'].set_color('none')
        axV.spines['bottom'].set_position('zero') # x-axis
        axV.spines['top'].set_color('none')
        axV.spines['left'].set_smart_bounds(True)
        axV.spines['bottom'].set_smart_bounds(True)
        axV.xaxis.set_ticks_position('bottom')
        axV.yaxis.set_ticks_position('left')
        
        axV.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.grid(b=True, which='minor', axis='both', linewidth=.2)
        axV.grid(b=True, which='major', axis='both', linewidth=1)
        
        plt.tight_layout()
        
        return Vfig
    
    
    
    ### Unused functions
        
    '''
    def plot(self):
        fig = plt.figure()
        axI = fig.add_subplot(111)
        axV = axI.twinx()
        axI.plot(self.h.tvec,self.h.Iphi,color='b')
        axI.set_ylabel('Current [nA]')
        axV.plot(self.h.tvec,self.h.Vsoma,color='g')
        axV.set_ylabel('Voltage [mV]')
        #plotLight()#Prot.pulses
        axI.set_xlim(0,self.h.tstop)
        axI.set_xlabel('Time [ms]')
    '''
    
    
    def go(self,tstop):
        h.finitialize()

        #g.begin()
        while h.t<tstop:
            h.fadvance()
            # update the graph while simulating
            #g.plot(h.t)

        #g.flush()
    



        
class simBrian(Simulator):
    simulator = 'Brian'
    """Class for network level simulations with Brian"""
    
    def __init__(self, Prot, RhO, params=simParams['Brian'], network=None, netParams=None, monitors=None): # G_RhO='G0'):
        
        #from brian2 import *
        # from brian2.only import * # Do not import pylab etc
        # Brian2 has its own version of numpy... http://brian2.readthedocs.org/en/latest/user/import.html
        #import brian2.numpy_ as np
        #import brian2.only as br2
        
        self.Prot = Prot
        self.RhO = RhO
        
        import brian2.only as br
        self.br = br
        
        self.varIrho = 'I'
        
        self.dt = params['dt'].value
        #self.method_choice = params['method_choice'].value
        # Threshold
        #self.threshold = params['threshold'].value # threshold='v > -50*mV'
        # Reset
        #self.reset = params['reset'].value # reset='v = v_r'
        
        # refractory='(1 + 2*rand())*ms'
        # reset='''refractory += 1*ms''' # Make refractory dynamic to model adaptation
        #G = NeuronGroup(N, '''dv/dt = -(v + w)/ tau_v : 1 (unless refractory)
        #              dw/dt = -w / tau_w : 1''',
        #        threshold='v > 1', reset='v=0; w+=0.1', refractory=2*ms)
        
        # Prepend S_ to state variables to avoid nameclash with MSVC under Windows
        self.stateVars = RhO.brianStateVars # ['S_'+s for s in RhO.stateVars]
        
        self.namespace = self.setParams(RhO)
        self.namespace.update(netParams)
        self.netParams = netParams
        
        #if network is None:
        #    network = self.buildNetwork(self.namespace)
        
        self.net = network
        
        #setRecords(self, GroupRec=None, IndRec=None)
        self.monitors = monitors
        
        self.G_RhO = 'Inputs' #'Retina'
        self.net[self.G_RhO].__dict__[self.stateVars[0]] = 1. # Set rhodopsins to dark-adapted state
    
    def setParams(self, RhO):
        params = {} #dict(RhO.paramsList) #Parameters()        
        #params = RhO.exportParams(params)
        #namespace = params.valuesdict()
        #for k,v in namespace:
        #    namespace[k] = v * modelUnits[k]
        for p in RhO.paramsList:
            params[p] = RhO.__dict__[p] * modelUnits[p]
        
        return params
    
    
    def prepare(self, Prot):
        """Function to compare simulator's timestep to the timestep required by the protocol"""
        Prot.prepare()
        dt = Prot.getShortestPeriod()
        self.checkDt(dt)
        self.br.defaultclock.dt = self.dt*ms
        self.rasters = Prot.genContainer()
        return # self.dt
    
    
    def buildNetwork(self, network, namespace, G_RhO='G0', varIrho='I_RhO'):
        

        self.G_RhO = G_RhO
        #G_RhO = br.NeuronGroup(N, eqs+RhO.eqs, threshold='V>V_th', reset='V=V_r', refractory=self.t_r*ms)
        
        #self.syns = []
        #S1 = br.Synapses(G_RhO, G1, pre='V_post += 1.75*mV')
        #S1.connect(True, p=0.2)
        #S1.delay = 0.5*ms
        self.namespace = namespace
        self.varIrho = varIrho
        
        # photoNeurons = neuronModel + Equations(RhO.brian, I=varIrho, g=vargrho, V=varV)
        self.net = network
        self.setRecords(self.G_RhO, 0)        
        
    
    def setRecords(self, GroupRec=None, IndRec=None): #StateRec=('I'),
        #if GroupRec == None:
        #    pass
            #Record all groups...
        #else:
        #    M = br.SpikeMonitor(GroupRec)
        #    R = br.PopulationRateMonitor(GroupRec) # Redundant with SpikeMonitor?
            
        
        #States = br.StateMonitor(GroupRec, StateRec, record=IndRec)
        #States = StateMonitor(GroupRec, ('I', 'V'), record=True)
        #States = StateMonitor(GroupRec, ('I', 'V'), record=[0, 10, 20])
        
        
        
        ### http://brian2.readthedocs.org/en/latest/user/running.html
        # If placing monitors in a list...
        self.monitors = [self.br.StateMonitor(self.net[GroupRec], self.stateVars, record=IndRec),   # States
                         self.br.StateMonitor(self.net[GroupRec], self.varIrho, record=IndRec),     # I
                         self.br.StateMonitor(self.net[GroupRec], 'V', record=IndRec),              # V
                         self.br.SpikeMonitor(self.net[GroupRec])]                                  # Spikes
        # Make this a dictionary... 'states', 'spikes', 'rates'
        self.net = self.br.Network(collect())  # automatically include G and S
        self.net.add(self.monitors)  # manually add the monitors
        #self.stateMon
        #self.currentMon
        #self.spikeMon
    
    def runTrial(self, RhO, phiOn, V, delD, cycles, dt, verbose=verbose): 
        """Main routine for simulating a square pulse train"""
        
        
        # phim = 3.54*10**17*mm**-2*second**-1
        # p = 0.985
        # q = 1.58
        # E = 0*mvolt
        # v = -70*mvolt
        # v0 = 43*mvolt
        # v1 = 4.1*mvolt
        # g = 1.1e5*psiemens
        # gam = 0.0161
        # k1 = 13.4*ms**-1
        # k2 = 2.71*ms**-1
        # kf = 0.103*ms**-1
        # kb = 0.139*ms**-1
        # Gd1 = 0.112*ms**-1
        # Gd2 = 0.0185*ms**-1
        # Gf0 = 0.0389*ms**-1
        # Gb0 = 0.0198*ms**-1
        # Gr0 = 0.00163*ms**-1
        # Go1 = 2*ms**-1
        # Go2 = 0.0567*ms**-1
        
        
        
        # N = 100
        # tau = 10*ms
        # vr = -70*mV
        # vt0 = -50*mV
        # delta_vt0 = 5*mV
        # tau_t = 75*ms
        # sigma = 0.5*(vt0-vr)
        
        # Rm = 70*Mohm
        
        #'sigma':'0.5*(vt0-vr)'
        
        #brianNamespace = {'N':100, 'tau':10*ms, 'vr':-70*mV, 'vt0':-50*mV, 'delta_vt0':5*mV, 'tau_t':75*ms, 'sigma':0.5*(-50--70)*mV, 'duration':100*ms, 'Rm':70.*Mohm,
        #'phim':3.54*10**17*mm**-2*second**-1, 'p':0.985, 'q':1.58, 'E':0*mvolt, 'v':-70*mvolt, 'v0':43*mvolt, 'v1':4.1*mvolt, 'g':1.1e5*psiemens, 'gam':0.0161, 
        #'k1':13.4*ms**-1, 'k2':2.71*ms**-1, 'kf':0.103*ms**-1, 'kb':0.139*ms**-1, 
        #'Gd1':0.112*ms**-1, 'Gd2':0.0185*ms**-1, 'Gf0':0.0389*ms**-1, 'Gb0':0.0198*ms**-1, 'Gr0':0.00163*ms**-1, 'Go1':2*ms**-1, 'Go2':0.0567*ms**-1}
        
        #'N':100, 
        #'duration':100*ms,
        #'v':-70*mvolt,
        #brianNamespace = {'tau':10*ms, 'V_r':-70*mV, 'V_th':-50*mV, 'delta_V':5*mV, 'tau_t':75*ms, 'R_m':70.*Mohm,
        #'phim':3.54*10**17*mm**-2*second**-1, 'p':0.985, 'q':1.58, 'E':0*mvolt, 'v0':43*mvolt, 'v1':4.1*mvolt, 'g':1.1e5*psiemens, 'gam':0.0161, 
        #'k1':13.4*ms**-1, 'k2':2.71*ms**-1, 'kf':0.103*ms**-1, 'kb':0.139*ms**-1, 
        #'Gd1':0.112*ms**-1, 'Gd2':0.0185*ms**-1, 'Gf0':0.0389*ms**-1, 'Gb0':0.0198*ms**-1, 'Gr0':0.00163*ms**-1, 'Go1':2*ms**-1, 'Go2':0.0567*ms**-1}
        
        nPulses = cycles.shape[0]
        #delD *= self.br.ms
        #cycles *= self.br.ms
        
        if verbose > 1:
            if V is not None:
                Vstr = 'V = {:+}mV, '.format(V)
            else:
                Vstr = ''
            info = "Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: [delD={:.4g}ms".format(phiOn,Vstr,delD)
            #info = "Simulating experiment at phi = {}, V = {:+}mV, pulse cycles: [delD={:.4g}ms".format(phiOn,V,delD) # {:.3g}photons/s/mm^2
            for p in range(nPulses):
                info += "; [onD={:.4g}ms; offD={:.4g}ms]".format(cycles[p,0],cycles[p,1])
            info += "]"
            print(info)
            
            report = 'text' #'stdout', 'stderr', function
        else:
            report = None
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        
        #self.namespace.update({'phi':phi*modelUnits['phim'], 'stimulus':bool(phi>0)})
        self.net[self.G_RhO].phi = phi*modelUnits['phim']
        self.net[self.G_RhO].stimulus = False
        # Skip last state value since this is defined as 1 - sum(states) not as an ODE
        #self.net[G_RhO].set_states({s: o for s, o in zip(self.stateVars[:-1], RhO.s_0[:-1])})
        
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        #start, end = RhO.t[0], RhO.t[0]+delD
        #t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) #t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
            #if verbose > 2:
            #    print("Simulating t_del = [{},{}]".format(start,end))
        #if RhO.useAnalyticSoln:
        #    soln = RhO.calcSoln(t, RhO.s0)
        #else:
        #    soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian) #t_del #delay
        
        #self.br.Network.run(delD*ms, report)
        self.net.run(duration=delD*ms, namespace=self.namespace, report=report) #brianNamespace
        #self.net.run(duration=delD*ms, report=report) #brianNamespace
        #self.br.run(duration=delD*ms, report=report)
        
        #t = t_del
        #RhO.storeStates(soln[1:],t[1:])
        #RhO.storeStates(self.monitors[0].S1[0][1:], self.monitors[0].t[1:]/ms)
        
        elapsed = delD
        
        for p in range(0, nPulses):
            
            ### Light on phase
            #RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            #start = end
            #end = start + cycles[p,0] # onD
            #t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) #t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
            #onInd = len(RhO.t) - 1  # Start of on-phase
            #offInd = onInd + len(t) - 1 # Start of off-phase
            #RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            
            elapsed += cycles[p,0]
            onInd = int(round(elapsed/dt))
            
            # Turn on light and set transition rates
            phi = phiOn  # Light flux
            #RhO.setLight(phi)
            # Update G_RhO parameters?
            
            #self.namespace.update({'phi':phi*modelUnits['phim'], 'stimulus':bool(phi>0)})
            self.net[self.G_RhO].phi = phi*modelUnits['phim']
            self.net[self.G_RhO].stimulus = True
            
            if verbose > 1:
                print(self.namespace)
            
            #if verbose > 1:
                #print("On-phase initial conditions:{}".format(RhO.s_on))
                #if verbose > 2:
                #    print("Simulating t_on = [{},{}]".format(start,end))
            #if RhO.useAnalyticSoln:
            #    soln = RhO.calcSoln(t, RhO.s_on)
            #else:
            #    soln = odeint(RhO.solveStates, RhO.s_on, t, args=(None,), Dfun=RhO.jacobian)
            #RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            
            self.net.run(duration=cycles[p,0]*ms, namespace=self.namespace, report=report)
            #self.net.run(duration=cycles[p,0]*ms, report=report)
            
            ### Light off phase
            #RhO.s_off = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            #start = end
            #end = start + cycles[p,1] # offD
            #if (p+1) == nPulses: # Add (or subtract) extra time after (during) the off phase
            #    end += padD
            #t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) #t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # endpoint=True
            # Turn off light and set transition rates
            phi = 0  # Light flux
            #RhO.setLight(phi)
            # Update G_RhO parameters?
            
            self.net[self.G_RhO].phi = phi*modelUnits['phim']
            self.net[self.G_RhO].stimulus = False
            #self.namespace.update({'phi':phi*modelUnits['phim'], 'stimulus':bool(phi>0)})
            
            #if verbose > 1:
                #print("Off-phase initial conditions:{}".format(RhO.s_off))
                #if verbose > 2:
                    #print("Simulating t_off = [{},{}]".format(start,end))
            #if RhO.useAnalyticSoln:
            #    soln = RhO.calcSoln(t, RhO.s_off)
            #else:
            #    soln = odeint(RhO.solveStates, RhO.s_off, t, args=(None,), Dfun=RhO.jacobian)
            #RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            
            self.net.run(duration=cycles[p,1]*ms, namespace=self.namespace, report=report)
            #self.net.run(duration=cycles[p,1]*ms, report=report)
            
            elapsed += cycles[p,1]
            offInd = int(round(elapsed/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            
        ### Calculate photocurrent
        ### Assumes that only one neuron is recorded from...
        
        assert(self.monitors['states'].n_indices == 1)
        self.monitors['states'].record_single_timestep()
        soln = np.hstack([self.monitors['states'].variables[s].get_value() for s in self.stateVars])
        
        #self.monitors[0].record_single_timestep()
        #soln = np.hstack([self.monitors[0].variables[s].get_value() for s in self.stateVars])
        
        #soln = np.hstack([self.monitors['states'].variables['_recorded_'+s].get_value() for s in self.stateVars])
        #assert(self.monitors[0].n_indices == 1)
        #soln = np.hstack([self.monitors[0].variables['_recorded_'+s].get_value() for s in self.stateVars])
        
        #RhO.storeStates(self.monitors[0].S1[0][1:], self.monitors[0].t[1:]/ms)
        

        
        # Print last value is not recorded!
        # https://github.com/brian-team/brian2/issues/452
        
        #t = self.monitors[0].t/ms
        t = self.monitors['states'].t/ms
        #t = self.monitors[0].t/ms
        #t = np.append(t, t[-1]+dt)          ### HACK!!!
        
        
        #RhO.storeStates(soln, t[1:]) # 1000 ### HACK!!! (not soln[1:])
        RhO.storeStates(soln[1:], t[1:]) # 1000 ### HACK!!! (not soln[1:])
        
        states, t = RhO.getStates()
        
        if V != None: # and 'I' in self.monitors:
            I_RhO = RhO.calcI(V, RhO.states)
        else:
            self.monitors['I'].record_single_timestep()
            I_RhO = self.monitors['I'].variables[self.varIrho].get_value() * 1e9# / nA
            #self.monitors[1].record_single_timestep()
            #I_RhO = self.monitors[1].variables[self.varIrho].get_value() * 1e9# / nA
            
            #I_RhO = self.monitors['I'].variables['_recorded_'+self.varIrho].get_value() * 1e9# / nA
            #I_RhO = self.monitors[1].variables['_recorded_'+self.varIrho].get_value() * 1e9# / nA
            
            #I_RhO = np.asarray(I_RhO.T.tolist()[0])
        
        
        
        
        #I_RhO = np.append(0., I_RhO)        ### HACK!!!
        #print(self.monitors[1].t[0])
        #print(len(self.monitors[1].t))
        
        # print(len(I_RhO))
        # print(I_RhO)
        # print(len(t))
        # print(t)
        # print(states.shape)
        # print(states)
        
        times, totT = cycles2times(cycles, delD)
        self.plotRasters(times, totT)
        self.raster = self.monitors['spikes']
        
        return I_RhO, t, states #I_RhO, t[1:], states
        
    

    
    def runTrialPhi_t(self, RhO, phi_ts, V, delD, cycles, endT, dt, verbose=verbose):
        """Main routine for simulating a pulse train"""
        # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
        #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
        
        ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
        ### Combine this with the original runTrial
        
        raise NotImplementedError('Error: Light as a continuous function of time has not been implemented for Brian yet!')
        
        nPulses = cycles.shape[0]
        assert(len(phi_ts) == nPulses)
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1,:]   # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+delD #start, end = 0.00, delD
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True) # int(round(((end-start)/dt)+1)), endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
        soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:],t[1:])

        ### Stimulation phases
        for p in range(0, nPulses):
            RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            onD, offD = cycles[p,0], cycles[p,1]
            if p < nPulses - 1:
                end = start + onD + offD #totT #start + stimD
            else:
                end = endT
            
            onInd = len(RhO.t) - 1  # Start of on-phase
            #offInd = onInd + len(t) - 1 # Start of off-phase
            #offInd = onInd + int(round(stimD/dt)) # Consider rounding issues...
            offInd = onInd + int(round(onD/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True) #int(round(((end-start)/dt)+1))
            phi_t = phi_ts[p]
            
            if verbose > 1:
                print("Pulse initial conditions:{}".format(RhO.s_on))
                if verbose > 2:
                    print('Pulse = [{}, {}]'.format(RhO.t[onInd],RhO.t[offInd]))
            
            if verbose > 2:
                soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian, full_output=True)
                print(out)
            else:
                soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian)
            
            RhO.storeStates(soln[1:], t[1:]) # Skip first values to prevent duplicating initial conditions and times
        
        
        ### Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states, t = RhO.getStates()
        
        return I_RhO, t, states
        
    
    def saveExtras(self, run, phiInd, vInd):
        ### HACK!!!
        self.rasters[run][phiInd][vInd] = copy.copy(self.monitors['spikes'])
    
    def plotRasters(self, times=None, totT=None):
        # Plot rasters
        
        nLayers = len(self.monitors['spikes'])
        
        Rfig = plt.figure() #figsize=(10, 6)
        gs = plt.GridSpec(nLayers,1)
        axes = [None for lay in range(nLayers)]
        
        for lay in range(nLayers):
            
            axes[lay] = Rfig.add_subplot(gs[lay])
            gInd = nLayers-lay-1
            plt.plot(self.monitors['spikes'][gInd].t/ms, self.monitors['spikes'][gInd].i, '.')
            plt.ylabel('$\mathrm{' + self.monitors['spikes'][gInd].name + '}$')
            if gInd > 0:
                plt.setp(axes[lay].get_xticklabels(), visible=False)
            else:
                plt.xlabel('$\mathrm{Time\ [ms]}$')
            if totT is not None:
                plt.xlim((0,totT))
            plt.ylim((0, len(self.monitors['spikes'][gInd].spike_trains())))
            if times is not None:
                #plotLight(times)
                for pulse in times:
                    axes[lay].axvspan(pulse[0], pulse[1], facecolor='y', alpha=0.2)
            
        # axV1 = fig.add_subplot(gs[0])
        # plt.plot(self.monitors['spikes'][2].t/ms, self.monitors['spikes'][2].i, '.')
        # #plt.plot(self.monitors[3][2].t/ms, self.monitors[3][2].i, '.')
        # plt.setp(axV1.get_xticklabels(), visible=False)
        # plt.ylabel('V1 Neuron')

        # axLGN = fig.add_subplot(gs[1], sharex=axV1)
        # plt.plot(self.monitors['spikes'][1].t/ms, self.monitors['spikes'][1].i, '.')
        # #plt.plot(self.monitors[3][1].t/ms, self.monitors[3][1].i, '.')
        # plt.setp(axLGN.get_xticklabels(), visible=False)
        # plt.ylabel('LGN Neuron')
        
        # axRet = fig.add_subplot(gs[2], sharex=axV1)
        # plt.plot(self.monitors['spikes'][0].t/ms, self.monitors['spikes'][0].i, '.')
        # #plt.plot(self.monitors[3][0].t/ms, self.monitors[3][0].i, '.')
        # plt.ylabel('Retinal Neuron')
        # plt.xlabel('Time/ms')
        

        # Add stimulus shading
        #times, totT = cycles2times(cycles, delD)
        # for pulse in times:
            # axRet.axvspan(pulse[0], pulse[1], facecolor='y', alpha=0.2)
            # axLGN.axvspan(pulse[0], pulse[1], facecolor='y', alpha=0.2)
            # axV1.axvspan(pulse[0], pulse[1], facecolor='y', alpha=0.2)

        plt.tight_layout()
        
        #from os import path
        figName = path.join(fDir, "raster."+config.saveFigFormat)
        Rfig.savefig(figName, format=config.saveFigFormat)
    
    def plotRaster(self, group=None, addHist=True):
        if group==None:
            group = all
        
        # Plot all groups on one figure or raster + hist for each group in separate figures?
        if addHist:
            pass
            # Add to subplot or create separate figure?
            
        Rfig = plt.figure() # Raster figure
        
        plt.plot(M.t/ms, M.i, '.b')
        xlabel('Time [ms]')
        ylabel('Neuron Index')
        
        return
    
    
    
    
from collections import OrderedDict
simulators = OrderedDict([('Python', simPython), ('NEURON', simNEURON), ('Brian', simBrian)])#{'Python': simPython, 'NEURON': simNEURON, 'Brian': simBrian}
