"""
Classes to wrap simulation engines for a uniform interface
    * ``Python`` for simulations at the `channel` level
    * ``NEURON`` for simulations at the `cellular` level
    * ``Brian2`` for simulations at the `network` level
"""

import warnings
import logging
import os
import copy
import abc
from collections import OrderedDict

import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from matplotlib import pyplot as plt

from pyrho.parameters import *
from pyrho.parameters import PyRhOobject, modelUnits, ms
from pyrho.utilities import *  # cycles2times, plotLight
from pyrho.expdata import *
from pyrho.models import *
from pyrho.config import *  # verbose
from pyrho.config import wall_time, _DASH_LINE, _DOUB_DASH_LINE
from pyrho import config

__all__ = ['simulators']

logger = logging.getLogger(__name__)


class Simulator(PyRhOobject):  # object
    """Common base class for all simulators."""

    __metaclass__ = abc.ABCMeta

    Prot = None
    RhO = None
    simulator = None
    clampProts = ['rectifier']

    @abc.abstractmethod
    def __init__(self, Prot, RhO, params=None):
        pass

    @abc.abstractmethod
    def runTrial(self, RhO, phiOn, V, Dt_delay, cycles, dt, verbose=config.verbose):
        pass

    @abc.abstractmethod
    def runTrialPhi_t(self, RhO, phi_ts, V, Dt_delay, cycles, dt, verbose=config.verbose):
        pass

    def __str__(self):
        return self.simulator

    def __repr__(self):
        return "<PyRhO {} Simulator object>".format(self.simulator)

    def checkDt(self, dt):
        """Compare simulator's timestep to the timestep required by the protocol."""
        if dt < self.dt:
            self.dt_prev, self.dt = self.dt, dt  # min(self.h.dt, dt)
            if config.verbose > 0:
                info = 'Time step reduced to {}ms by protocol!'.format(self.dt)
                print(info)
                logger.info(info)
        return self.dt

    def prepare(self, Prot):
        """Prepare the simulator according to the protocol and opsin."""
        Prot.prepare()
        dt = Prot.getShortestPeriod()
        Prot.dt = self.checkDt(dt)
        assert self.dt > 0
        assert Prot.dt > 0
        return  # self.dt

    def initialise(self):
        pass

    def run(self, verbose=config.verbose):
        """Main routine to run the simulation protocol."""
        # TODO: Use joblib to cache results https://pythonhosted.org/joblib/

        t0 = wall_time()
        RhO = self.RhO
        Prot = self.Prot

        self.prepare(Prot)

        if verbose > 0:
            #print("\n================================================================================")
            print("\n"+_DOUB_DASH_LINE)
            prot_info = "Running '{}' protocol with {} " \
                        "for the {} model... ".format(Prot, self, RhO)
            print(prot_info)
            logger.info(prot_info)
            #print("================================================================================\n")
            print(_DOUB_DASH_LINE+"\n")
            prot_details = "{{nRuns={}, nPhis={}, nVs={}}}".format(Prot.nRuns, Prot.nPhis, Prot.nVs)
            print(prot_details)
            logger.info(prot_details)

        Prot.PD = ProtocolData(Prot.protocol, Prot.nRuns, Prot.phis, Prot.Vs)
        Prot.PD.I_peak_ = [[[None for v in range(Prot.nVs)] for p in range(Prot.nPhis)] for r in range(Prot.nRuns)]
        Prot.PD.I_ss_ = [[[None for v in range(Prot.nVs)] for p in range(Prot.nPhis)] for r in range(Prot.nRuns)]
        if hasattr(Prot, 'runLabels'):
            Prot.PD.runLabels = Prot.runLabels

        if verbose > 1:
            Prot.printParams()
            Prot.logParams()

        for run in range(Prot.nRuns):

            cycles, Dt_delay = Prot.getRunCycles(run)
            pulses, Dt_total = cycles2times(cycles, Dt_delay)

            for phiInd, phiOn in enumerate(Prot.phis):

                if verbose > 1 and (Prot.nPhis > 1 or (run == 0 and phiInd == 0)):
                    RhO.dispRates()

                for vInd, V in enumerate(Prot.Vs):
                    # N.B. solution variables are not currently dependent on V

                    # Reset simulation environment...
                    self.initialise()

                    # TODO: Deprecate special square pulse fucntions
                    if Prot.squarePulse and self.simulator == 'Python':
                        I_RhO, t, soln = self.runTrial(RhO, phiOn, V, Dt_delay, cycles, self.dt, verbose)
                    else:  # Arbitrary functions of time: phi(t)
                        phi_ts = Prot.phi_ts[run][phiInd][:]
                        I_RhO, t, soln = self.runTrialPhi_t(RhO, phi_ts, V, Dt_delay, cycles, self.dt, verbose)

                    stim = Prot.getStimArray(run, phiInd, self.dt)
                    PC = PhotoCurrent(I_RhO, t, pulses, phiOn, V, stimuli=stim,
                                      states=soln, stateLabels=RhO.stateLabels,
                                      label=Prot.protocol)
                    #PC.alignToTime()

                    PC.ssInf = np.array(RhO.ssInf)
                    Prot.PD.trials[run][phiInd][vInd] = PC
                    Prot.PD.I_peak_[run][phiInd][vInd] = PC.I_peak_
                    Prot.PD.I_ss_[run][phiInd][vInd] = PC.I_ss_

                    self.saveExtras(run, phiInd, vInd)

                    if verbose > 1:
                        prot_details = 'Run=#{}/{}; phiInd=#{}/{}; vInd=#{}/{}; Irange=[{:.3g},{:.3g}]'.format(run, Prot.nRuns, phiInd, Prot.nPhis, vInd, Prot.nVs, PC.I_range_[0], PC.I_range_[1])
                        print(prot_details)
                        logger.debug(prot_details)

        Prot.finish(PC, RhO)
        # self.finish() # Reset dt and Vclamp

        if Prot.saveData:
            Prot.dataTag = str(RhO.nStates)+"s"
            saveData(Prot.PD, Prot.protocol+Prot.dataTag)

        # Reset any variables overridden by the protocol
        if hasattr(self, 'dt_prev') and self.dt_prev is not None:
            self.dt, self.dt_prev = self.dt_prev, None

        # TODO: Move Vclamp reset to NEURON simulator in self.finish()
        if hasattr(self, 'Vclamp_prev') and self.Vclamp_prev is not None:
            self.Vclamp, self.Vclamp_prev = self.Vclamp_prev, None
            self.h.Vcl.rs = 1e9  # TODO; Find a way to fully remove SEClamp

        self.runTime = wall_time() - t0
        if verbose > 0:
            prot_info = "\nFinished '{}' protocol with {} for the {} model in"\
                        " {:.3g}s".format(Prot, self, RhO, self.runTime)
            print(prot_info)
            #print("--------------------------------------------------------------------------------\n")
            print(_DASH_LINE+"\n")
            logger.info(prot_info)

        return Prot.PD

    def saveExtras(self, run, phiInd, vInd):
        pass

    def plot(self):
        self.Prot.plot()
        self.plotExtras()

    def plotExtras(self):
        pass


class simPython(Simulator):
    """Class for channel level simulations with Python."""

    simulator = 'Python'

    def __init__(self, Prot, RhO, params=simParams['Python']):
        self.dt = params['dt'].value
        self.Prot = Prot
        self.RhO = RhO

    def runSoln(self, RhO, t):

        if RhO.useAnalyticSoln:  # Leave the ability to manually override
            logger.debug('Calculating solution analytically.')
            try:
                soln = RhO.calcSoln(t, RhO.states[-1, :])
            except:  # Any exception e.g. NotImplementedError or ValueError
                logger.debug('Falling back to numerical integration.')
                soln = odeint(RhO.solveStates, RhO.states[-1, :], t,
                              args=(None,), Dfun=RhO.jacobian)
        else:
            logger.debug('Calculating solution numerically.')
            soln = odeint(RhO.solveStates, RhO.states[-1, :], t,
                          args=(None,), Dfun=RhO.jacobian)

        if np.isnan(np.sum(soln)):  # np.any(np.isnan(soln)):
            warnings.warn('The state solution is undefined!')

        return soln

    def runTrial(self, RhO, phiOn, V, Dt_delay, cycles, dt, verbose=config.verbose):
        """
        Main routine for simulating a square pulse train.


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
            info = "Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: [Dt_delay={:.4g}ms".format(phiOn, Vstr, Dt_delay)
            for p in range(nPulses):
                info += "; [Dt_on={:.4g}ms; Dt_off={:.4g}ms]".format(cycles[p, 0], cycles[p, 1])
            info += "]"
            print(info)

        # Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi)         # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1, :]  # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+Dt_delay
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)

        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
            if verbose > 2:
                print("Simulating t_del = [{},{}]".format(start, end))

        soln = self.runSoln(RhO, t)
        RhO.storeStates(soln[1:], t[1:])

        for p in range(0, nPulses):

            # Light on phase
            RhO.s_on = soln[-1, :]
            start = end
            end = start + cycles[p, 0]
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True)
            onInd = len(RhO.t) - 1       # Start of on-phase
            offInd = onInd + len(t) - 1  # Start of off-phase
            RhO.pulseInd = np.vstack((RhO.pulseInd, [onInd, offInd]))
            # Turn on light and set transition rates
            phi = phiOn  # Light flux
            RhO.setLight(phi)

            if verbose > 1:
                print("On-phase initial conditions:{}".format(RhO.s_on))
                if verbose > 2:
                    print("Simulating t_on = [{},{}]".format(start, end))

            soln = self.runSoln(RhO, t)
            RhO.storeStates(soln[1:], t[1:])  # Skip first values to prevent duplicating initial conditions and times
            RhO.ssInf.append(RhO.calcSteadyState(phi))


            # Light off phase
            RhO.s_off = soln[-1, :]
            start = end
            end = start + cycles[p, 1]
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True)
            # Turn off light and set transition rates
            phi = 0  # Light flux
            RhO.setLight(phi)

            if verbose > 1:
                print("Off-phase initial conditions:{}".format(RhO.s_off))
                if verbose > 2:
                    print("Simulating t_off = [{},{}]".format(start, end))

            soln = self.runSoln(RhO, t)
            RhO.storeStates(soln[1:], t[1:])  # Skip first values to prevent duplicating initial conditions and times

            if verbose > 1:
                print('t_pulse{} = [{}, {}]'.format(p, RhO.t[onInd], RhO.t[offInd]))

        # Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states, t = RhO.getStates()
        return I_RhO, t, states

    def runTrialPhi_t(self, RhO, phi_ts, V, Dt_delay, cycles, dt, verbose=config.verbose):
        """Main routine for simulating a pulse train."""

        nPulses = cycles.shape[0]
        assert len(phi_ts) == nPulses

        # Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi)                     # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1, :]               # Store initial state used
        start, end = RhO.t[0], RhO.t[0]+Dt_delay
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)  # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
        soln = odeint(RhO.solveStates, RhO.s0, t, args=(None,), Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], t[1:])

        # Stimulation phases
        for p in range(0, nPulses):
            RhO.s_on = soln[-1, :]
            start = end
            Dt_on, Dt_off = cycles[p, 0], cycles[p, 1]
            end = start + Dt_on + Dt_off

            onInd = len(RhO.t) - 1                 # Start of on-phase
            offInd = onInd + int(round(Dt_on/dt))  # Start of off-phase
            RhO.pulseInd = np.vstack((RhO.pulseInd, [onInd, offInd]))
            nSteps = int(round(((end-start)/dt)+1))
            t = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]

            if verbose > 1:
                print("Pulse initial conditions:{}".format(RhO.s_on))

            if verbose > 2:
                soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,),
                                   Dfun=RhO.jacobian, full_output=True)
                print(out)
            else:
                soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,),
                              Dfun=RhO.jacobian)

            RhO.storeStates(soln[1:], t[1:])  # Skip first values to prevent duplicating initial conditions and times
            RhO.ssInf.append(RhO.calcSteadyState(phi_t(end-Dt_off)))

            if verbose > 1:
                print('t_pulse{} = [{}, {}]'.format(p, RhO.t[onInd], RhO.t[offInd]))

        # Calculate photocurrent
        I_RhO = RhO.calcI(V, RhO.states)
        states, t = RhO.getStates()

        return I_RhO, t, states






class simNEURON(Simulator):
    """Class for cellular level simulations with NEURON."""

    simulator = 'NEURON'
    mechanisms = {3: 'RhO3c', 4: 'RhO4c', 6: 'RhO6c'}

    def __init__(self, Prot, RhO, params=simParams['NEURON'], recInd=0):  #v_init=-70, integrator='fixed'):

        ### Model Specification
        # Topology
        # Geometry
        # Biophysics

        ### Instrumentation
        # Synaptic Input
        # Graphing - try using NEURON's facilities

        self.Prot = Prot
        self.RhO = RhO

        import copy
        from neuron import h
        #self.neuron = neuron
        self.h = h  # Enables self.h('any hoc statement')
                    # c.f. nrnpython('any python statement')

        self.h.load_file('stdrun.hoc')
        ### To load mod files:
        # Add os.environ['NRN_NMODL_PATH'] to environment variables. See $NEURONPATH/nrn/lib/python/neuron/__init__.py
        #import neuron as nrn
        #self.nrn = nrn
        #self.nrn.load_mechanisms(os.environ['NRN_NMODL_PATH'])

        #self.reset()

        ### Simulation control
        self.CVode = params['CVode'].value

        if self.CVode is True:  # self.integrator == 'variable': #... if dt <= 0:
            self.h('objref cvode')
            self.h.cvode = self.h.CVode()
            self.h.cvode.active(1)
        else:  # 'fixed' N.B. NEURON overrides this, so it is reset at the end of runTrial
            # self.h.cvode.active(0)
            self.h.dt = params['dt'].value
            self.dt = self.h.dt

        # TODO: Set as parameters
        # self.h.cvode.atol(0.000001)
        # self.h.cvode.rtol(0.000001)
        if config.verbose > 0:
            self.cvode = h.CVode()
            print('Integrator tolerances: absolute=', self.cvode.atol(),
                                        ' relative=', self.cvode.rtol())

        # self.h.v_init = params['v_init'].value
        self.h.finitialize(params['v_init'].value)

        ### TODO: Move into prepare() in order to insert the correct type of mechanism (continuous or discrete) according to the protocol
        #for states, mod in self.mechanisms.items():
        #    self.mechanisms[states] += 'c'
        self.mod = self.mechanisms[self.RhO.nStates]  # Use this to select the appropriate mod file for insertion
        #if not Prot.squarePulse:
        #self.mod += 'c'

        # self.buildCell(params['cell'].value)
        self.build_cell()
        if config.verbose > 0:
            self.h.topology()  # Print topology
        self.transduce(self.RhO, expProb=params['expProb'].value)
        self.rhoParams = copy.deepcopy(modelParams[str(self.RhO.nStates)])
        self.RhO.exportParams(self.rhoParams)
        self.setOpsinParams(self.rhoList, self.rhoParams)  # self.RhO, modelParams[str(self.RhO.nStates)])
        print(recInd)
        print(self.rhoList)
        self.rhoRec = self.rhoList[recInd]  # Choose a rhodopsin to record
        self.setRecords(self.rhoRec, params['Vcomp'].value, self.RhO)
        self.Vclamp = params['Vclamp'].value

    def reset(self):
        # TODO
        """Clear any previous models"""
        # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=2&t=2367
        # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2655137/
        # http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html

        pc = self.h.ParallelContext()
        pc.gid_clear()

        # unref all the cell objects
        self.h("forall delete_section()")
        for sec in self.h.allsec():
            self.h("{}{{delete_section()}}".format(sec.name()))
        #self.h('objectref compList')
        #self.h('objectref rhoList')

        # unref all the NetCon

        # Other top level hoc objects (Vectors, Graphs etc.)
        #pc.done()
        return

    def prepare(self, Prot):
        """Function to prepare the simulator according to the protocol and rhodopsin"""
        Prot.prepare()

        #if self.CVode is False: # ~ Always set Prot.dt
        dt = Prot.getShortestPeriod()
        Prot.dt = self.checkDt(dt)
        self.h.dt = self.dt

        # Include code to set Vclamp = False if Prot.Vs=[None]
        if Prot.protocol in self.clampProts:
            self.Vclamp_prev, self.Vclamp = self.Vclamp, True
            warnings.warn('Voltage clamp set for protocol: {}'.format(Prot.protocol))
        if self.Vclamp is False:
            Prot.Vs = [None]
            Prot.nVs = 1
        else:
            self.addVclamp()

        self.Vms = Prot.genContainer()
        return  # self.h.dt

    def build_cell(self):
        self.cell = [self.h.Section(name='soma')]
        self.nSecs = len(self.cell)
        if config.verbose > 0:
            print(self.h.topology())
            print(f'Total sections: {self.nSecs}')

    def buildCell(self, scriptList=[]):
        """Pass a list of hoc or python files in the order in which they are to
        be processed."""

        #import subprocess
        # Check for script in path cwd and then fallback to NRN_NMODL_PATH
        cwd = os.getcwd()
        for script in scriptList:
            #name, ext = os.path.splitext(script)
            ext = script.split('.')[-1]
            if not os.path.isfile(os.path.join(cwd, script)):
                # script = os.path.join(os.environ['NRN_NMODL_PATH'], script)
                script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'NEURON', script)

            if config.verbose > 0:
                print('Loading ', script)

            if ext.lower() == 'hoc':
                self.h.load_file(script)
                '''
                elif ext == 'py':
                    with open(script, "r") as fh:
                        exec(fh.read(), globals(), locals()) # Is this safe?!
                    #import script
                    #os.system("""python {}""".format(script))
                    subprocess.Popen(["python", script])
                    #self.cell = subprocess.Popen(["python", script])
                '''
            else:
                raise ValueError('Unknown file type!')

        self.cell = self.h.SectionList()
        self.nSecs = 0
        for sec in self.h.allsec():
            self.cell.append(sec)
            self.nSecs += 1

        if config.verbose > 0:
            print('Total sections: ', self.nSecs)  # self.cell.count())

    def transduce(self, RhO, expProb=1):
        import random

        self.compList = self.h.List()
        self.rhoList = self.h.List()
        self.h('objectvar rho')
        mech = self.mechanisms[RhO.nStates]
        for sec in self.cell:  # self.h.allsec(): # Loop over every section in the cell
            self.compList.append(sec)
            if random.random() <= expProb:  # Insert a rhodopsin and append it to a rhoList
                self.h('rho = new {}(0.5)'.format(mech))
                self.rhoList.append(self.h.rho)

    def setOpsinParams(self, rhoList, pSet):  # compList
        for rho in rhoList:  # not self.rhoList so that subsets can be passed
            for p in pSet:
                setattr(rho, p, pSet[p].value)

    def getOpsinParams(self):
        # Use nrnpython to set Python variables from within hoc?
        pass

    def setRecords(self, rhoRec, Vcomp, RhO):

        self.h('objref Vm')
        self.h.Vm = self.h.Vector()
        self.h('Vm.record(&{}.v(0.5))'.format(Vcomp))   # self.h.Vm.record(memRec(0.5)._ref_v) #cvode.record(&v(0.5),vsoma,tvec)
        self.h('objref Iphi')
        self.h('objref tvec')
        self.h.Iphi = self.h.Vector()
        self.h.tvec = self.h.Vector()

        # Record time vector (since the solver may use variable time steps
        self.h.tvec.record(self.h._ref_t)  # Record time points

        # Record photocurrent
        self.h.Iphi.record(rhoRec._ref_i)

        # h.setpointer(_ref_hocvar, 'POINTER_name', point_proces_object)
        # h.setpointer(_ref_hocvar, 'POINTER_name', nrn.Mechanism_object)
        # For example if a mechanism with suffix foo has a POINTER bar and you want it to point to t use
        # h.setpointer(_ref_t, 'bar', sec(x).foo)

        # Save state variables according to RhO model
        for s in RhO.stateVars:
            self.h('objref {}Vec'.format(s))
            self.h.rho = rhoRec                             # TODO: Check this works with multiple sections/rhodopsins
            self.h('{}Vec = new Vector()'.format(s))
            self.h('{}Vec.record(&rho.{})'.format(s, s))

    def addVclamp(self):
        self.h('objref Vcl')
        self.h.Vcl = self.h.SEClamp(0.5)
        self.h.Vcl.rs = 1  # Set to default in case it has been increased to disable it

    def setVclamp(self, Vhold=-70):
        self.h.Vcl.dur1 = self.h.tstop
        self.h.Vcl.amp1 = Vhold

    def setPulses(self, phiOn, Dt_delay, Dt_on, Dt_off, nPulses):
        for rho in self.rhoList:
            rho.phiOn = phiOn
            rho.Dt_delay = Dt_delay
            rho.Dt_on = Dt_on
            rho.Dt_off = Dt_off
            rho.nPulses = nPulses

    def runTrial(self, RhO, phiOn, V, Dt_delay, cycles, dt, verbose=config.verbose):

        # Notes on the integrator
        # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=8&t=1330
        # http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html

        # Model events: http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html#moDt_delayescriptionissues-events

        nPulses = cycles.shape[0]
        times, Dt_total = cycles2times(cycles,Dt_delay)

        # phi_ts = self.Prot.phi_ts[0][0][:]#[run][phiInd][:]
        # self.runTrialPhi_t(self, RhO, self.Prot.genPhiFuncs(), V, Dt_delay, cycles, dt)
        # return

        RhO.initStates(phi=0) # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1, :]   # Store initial states for plotting

        # Set simulation run time
        self.h.tstop = Dt_total

        if self.Vclamp is True:
            self.setVclamp(V)

        fixedPulses = False

        if fixedPulses:
            for p in range(1, nPulses):
                if (cycles[p, 0] != cycles[p-1, 0]) or cycles[p, 1] != cycles[p-1, 1]:
                    warnings.warn('Warning: Pulse cycles must be identical for NEURON simulations - replicating pulse!')
                    #self.Prot.squarePulse = False
                    #self.runTrialPhi_t(self, RhO, self.Prot.genPhiFuncs(), V, Dt_delay, cycles, dt)
                    #return
            Dt_on, Dt_off, padD = cycles[0, 0], cycles[0, 1], 0   # HACK!!! Use first pulse timing only...

            if verbose > 0:
                Vstr = '' if V is None else 'V = {:+}mV, '.format(V)
                info = ("Simulating experiment at phi = {:.3g}photons/mm^2/s, "
                        "{}pulse cycles: [Dt_delay={:.4g}ms; Dt_on={:.4g}ms; "
                        "Dt_off={:.4g}ms]".format(phiOn, Vstr, Dt_delay, Dt_on, Dt_off+padD))

            self.setPulses(phiOn, Dt_delay, Dt_on, Dt_off, nPulses)
            self.h.init()  # self.neuron.init()
            # self.h.finitialize()

            self.h.run()  # self.neuron.run(self.h.tstop)

        else:  # Work in progress to simulate variable pulses

            # progress = Dt_delay
            # padD = 0
            # for p in range(nPulses):
                # delay = Dt_delay if p==0 else 0
                # Dt_on, Dt_off = cycles[p]
                # self.setPulses(phiOn, delay, Dt_on, Dt_off, 1)
                # progress += (Dt_on + Dt_off)
                # #tstop = delay + np.sum(cycles[p])
                # self.h.tstop = progress # tstop #delay + np.sum(cycles[p])
                # #while self.h.t<tstop:
                # #    self.h.fadvance()
                # self.h.run()

            if verbose > 0:
                Vstr = '' if V is None else 'V = {:+}mV, '.format(V)
                print("Simulating experiment at phi = {:.3g}photons/mm^2/s, {}pulse cycles: Dt_delay={:.4g}ms;".format(phiOn, Vstr, Dt_delay))

            progress = Dt_delay
            self.h.init()
            firstRun = True
            p = 0
            while p < nPulses:
                delay = Dt_delay if p == 0 else 0
                identPulses = 1
                Dt_on, Dt_off = cycles[p]
                #prevOffD = Dt_off
                for pCheck in range(p+1, nPulses):
                    #if (cycles[pCheck, 0] == Dt_on) and (cycles[pCheck, 1] == Dt_off):
                    if np.isclose(cycles[pCheck, 0], Dt_on) and np.isclose(cycles[pCheck, 1], Dt_off):
                        identPulses += 1
                        p += 1
                    else:
                        break
                progress += ((Dt_on + Dt_off) * identPulses)
                if verbose > 0:
                    print(" [Dt_on={:.4g}ms; Dt_off={:.4g}ms] x {}".format(Dt_on, Dt_off, identPulses))

                self.setPulses(phiOn, delay, Dt_on, Dt_off, nPulses)

                if firstRun:
                    # self.setPulses(phiOn, delay, Dt_on, Dt_off, identPulses+1)
                    self.h.tstop = progress
                    self.h.run()
                    firstRun = False
                else:
                    # self.setPulses(phiOn, prevOffD, Dt_on, Dt_off, identPulses+1)
                    self.h.continuerun(progress)
                p += 1

        I_RhO = np.array(self.h.Iphi.to_python(), copy=True)
        t = np.array(self.h.tvec.to_python(), copy=True)
        self.Vm = np.array(self.h.Vm.to_python(), copy=True)
        self.t = t

        for p in range(nPulses):
            pInds = np.searchsorted(t, times[p, :], side="left")
            RhO.pulseInd = np.vstack((RhO.pulseInd, pInds))
            RhO.ssInf.append(RhO.calcSteadyState(phiOn))

        # Get solution variables
        soln = np.zeros((len(t), RhO.nStates))
        for sInd, s in enumerate(RhO.stateVars):
            self.h('objref tmpVec')
            self.h('tmpVec = {}Vec'.format(s))
            soln[:, sInd] = np.array(self.h.tmpVec.to_python(), copy=True)
            self.h('{}Vec.clear()'.format(s))
            self.h('{}Vec.resize(0)'.format(s))

        # Reset vectors
        self.h.tvec.clear()  # Necessary?
        self.h.tvec.resize(0)
        self.h.Iphi.clear()
        self.h.Iphi.resize(0)
        self.h.Vm.clear()
        self.h.Vm.resize(0)

        # NEURON changes the timestep! Set the actual timestep for plotting stimuli
        self.dt = self.h.dt
        self.Prot.dt = self.h.dt

        return I_RhO, t, soln

    def runTrialPhi_t(self, RhO, phi_ts, V, Dt_delay, cycles, dt, verbose=config.verbose):
        """Main routine for simulating a pulse train."""
        # TODO: Rename RhOnc.mod to RhOn.mod and RhOn.mod to RhOnD.mod (discrete) and similar for runTrialPhi_t
        # Make RhO.mod a continuous function of light - 9.12: Discontinuities 9.13: Time-dependent parameters p263
            #- Did cubic spline interpolation ever get implemented for vec.play(&rangevar, tvec, 1)?
            #- http://www.neuron.yale.edu/neuron/faq#vecplay
            #- http://www.neuron.yale.edu/neuron/static/new_doc/programming/math/vector.html#Vector.play
            #- Ramp example: http://www.neuron.yale.edu/phpbb/viewtopic.php?f=15&t=2602

        ### Some notes from the NEURON documentation
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

        if self.Vclamp is True:
            self.setVclamp(V)

        nPulses = cycles.shape[0]
        assert len(phi_ts) == nPulses

        times, Dt_total = cycles2times(cycles, Dt_delay)

        if verbose > 1:
            if V is not None:
                Vstr = 'V = {:+}mV, '.format(V)
            else:
                Vstr = ''
            info = "Simulating experiment {}pulse cycles: [Dt_delay={:.4g}ms".format(Vstr, Dt_delay)
            for p in range(nPulses):
                info += "; [Dt_on={:.4g}ms; Dt_off={:.4g}ms]".format(cycles[p, 0], cycles[p, 1])
            info += "]"
            print(info)

        # Set simulation run time
        self.h.tstop = Dt_total  # Dt_delay + np.sum(cycles) #nPulses*(Dt_on+Dt_off) + padD

        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi)  # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1, :]   # Store initial state used

        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))

        start, end = RhO.t[0], RhO.t[0]+Dt_delay  # start, end = 0.00, Dt_delay
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)
        phi_tV = np.zeros_like(t)

        discontinuities = np.asarray([len(t) - 1])
        for p in range(nPulses):
            start = end
            Dt_on, Dt_off = cycles[p, 0], cycles[p, 1]
            end = start + Dt_on + Dt_off
            nSteps = int(round(((end-start)/dt)+1))
            tPulse = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]
            phiPulse = phi_t(tPulse) # -tPulse[0] # Align time vector to 0 for phi_t to work properly
            discontinuities = np.r_[discontinuities, len(tPulse) - 1]

            onInd = len(t) - 1  # Start of on-phase
            offInd = onInd + int(round(Dt_on/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd, [onInd, offInd]))

            t = np.r_[t, tPulse[1:]]
            phi_tV = np.r_[phi_tV, phiPulse[1:]]

            RhO.ssInf.append(RhO.calcSteadyState(phi_t(end-Dt_off)))

        tvec = self.h.Vector(t)
        tvec.label('Time [ms]')
        phi_tV[np.ma.where(phi_tV < 0)] = 0  # Safeguard for negative phi values # Necessary? Handled in mods
        phiVec = self.h.Vector(phi_tV)
        phiVec.label('phi [ph./mm^2/s]')

        phiVec.play(self.rhoRec._ref_phi, tvec, 1, discontinuities)
        #phiVec.play_remove()

        self.h.init()
        self.h.run()

        # Collect data
        I_RhO = np.array(self.h.Iphi.to_python(), copy=True)
        t = np.array(self.h.tvec.to_python(), copy=True)
        self.Vm = np.array(self.h.Vm.to_python(), copy=True)
        self.t = t

        # Get solution variables N.B. NEURON changes the sampling rate
        soln = np.zeros((len(t), RhO.nStates))
        for sInd, s in enumerate(RhO.stateVars):
            self.h('objref tmpVec')
            self.h('tmpVec = {}Vec'.format(s))
            soln[:, sInd] = np.array(self.h.tmpVec.to_python(), copy=True)
            self.h('{}Vec.clear()'.format(s))
            self.h('{}Vec.resize(0)'.format(s))

        # Reset vectors
        self.h.tvec.clear()  # Necessary?
        self.h.tvec.resize(0)
        self.h.Iphi.clear()
        self.h.Iphi.resize(0)
        self.h.Vm.clear()
        self.h.Vm.resize(0)
        # phiVec ?

        RhO.storeStates(soln[1:], t[1:])

        # NEURON changes the timestep! Set the actual timestep for plotting stimuli
        self.dt = self.h.dt
        self.Prot.dt = self.h.dt

        ### Calculate photocurrent
        #I_RhO = RhO.calcI(V, RhO.states) # Recorded directly
        states, t = RhO.getStates()

        return I_RhO, t, soln

    def saveExtras(self, run, phiInd, vInd):
        # TODO: Clean this up
        self.Vms[run][phiInd][vInd] = copy.copy(self.Vm)
        return

    def plotExtras(self):  # TODO: REVISE!!!

        if self.Vclamp:
            return

        Prot = self.Prot
        RhO = self.RhO

        Vfig = plt.figure()
        axV = Vfig.add_subplot(111)
        for run in range(Prot.nRuns):
            for phiInd, phiOn in enumerate(Prot.phis):
                for vInd, V in enumerate(Prot.Vs):
                    colour, style = Prot.getLineProps(run, vInd, phiInd)
                    Vm = self.Vms[run][phiInd][vInd]
                    pc = Prot.PD.trials[run][phiInd][vInd]
                    t = pc.t
                    pulses = pc.pulses
                    axV.plot(t, Vm, color=colour, ls=style)

                    if Prot.protocol == 'shortPulse':
                        if run == 0 and phiInd == 0 and vInd == 0:
                            ymin, ymax = axV.get_ylim()
                            pos = 0.02 * abs(ymax-ymin)
                            axV.set_ylim(ymin, ymax + pos*(Prot.nRuns+1), auto=True)  # Allow extra space for thick lines
                        plot_light_bar(axV, pulses, y=ymax+(run+1)*pos, colour=colour)

                    else:
                        if run == 0 and phiInd == 0 and vInd == 0:
                            plotLight(pulses)
        # TODO: Plot extra pulses and refactor
        #self.t_start, self.t_end = min(t_starts), max(t_ends)
        # Add stimuli
        #for p in range(self.nPulses):
        #    sameStart, sameEnd = False, False
        #    if np.allclose(pulseSet[p, 0, :], np.tile(pulseSet[p, 0, 0], (1, 1, self.nRuns))): #pth t_on are the same
        #        sameStart = True
        #    if np.allclose(pulseSet[p, 1, :], np.tile(pulseSet[p, 1, 0], (1, 1, self.nRuns))): #pth t_off are the same
        #        sameEnd = True

        #    if sameStart and sameEnd: #np.allclose(pulseSet[p,:,run], np.tile(pulseSet[p,:,0], (1,1,self.nRuns))): #pth pulses are the same
        #        plotLight(np.asarray([pulseSet[p, :, 0]]), ax=ax, light=light, lam=470, alpha=0.2)

        #    elif not sameStart and not sameEnd: # No overlap
        #        for run in range(self.nRuns):
        #            plotLight(np.asarray([pulseSet[p, :, run]]), ax=ax, light=light, lam=470, alpha=0.2)

        #    else: #not (sameStart and sameEnd): # One or the other - xor
        #        pass # This applies to shortPulse only at present - do not shade!

        plt.ylabel('$\mathrm{Membrane\ Potential\ [mV]}$')  # axV.set_ylabel('Voltage [mV]')
        axV.set_xlim((t[0], t[-1]))
        plt.xlabel('$\mathrm{Time\ [ms]}$', position=(config.xLabelPos, 0), ha='right')

        #setCrossAxes(axV, zeroX=False)  # Zeroing x-axis caused the plot to become too big!
        # #axV.spines['bottom'].set_position('zero') # x-axis # Caused the plot to be too big!
        setCrossAxes(axV, zeroX=False, zeroY=True)
        axV.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.grid(visible=True, which='minor', axis='both', linewidth=.2)
        axV.grid(visible=True, which='major', axis='both', linewidth=1)

        if len(Prot.PD.legLabels) > 0 and Prot.PD.legLabels[0] != '':
            plt.legend(Prot.PD.legLabels)

        #axV.set_ylim(axV.get_ylim())
        ymin, ymax = axV.get_ylim()
        plt.tight_layout()
        plt.show()

        #figName = '{}Vm{}s-{}-{}-{}'.format(Prot.protocol, RhO.nStates, run, phiInd, vInd)
        figName = '{}Vm{}s'.format(Prot.protocol, RhO.nStates)
        fileName = os.path.join(config.fDir, figName+"."+config.saveFigFormat)
        logger.info('Saving Membrane potential figure to: {}'.format(fileName))
        Vfig.savefig(fileName, format=config.saveFigFormat)
        self.axV = axV


class simBrian(Simulator):
    """Class for network level simulations with Brian."""

    simulator = 'Brian'

    def __init__(self, Prot, RhO, params=simParams['Brian'],
                 network=None, netParams=None, monitors=None,
                 G_RhO='Inputs', varIrho='I', varV='v'):

        #from brian2 import *
        # from brian2.only import * # Do not import pylab etc
        # Brian2 has its own version of numpy... http://brian2.readthedocs.org/en/latest/user/import.html
        #import brian2.numpy_ as np
        #import brian2.only as br2

        self.Prot = Prot
        self.RhO = RhO

        import brian2.only as br
        self.br = br

        self.dt = params['dt'].value
        # http://brian2.readthedocs.org/en/2.0b4/advanced/state_update.html
        # http://brian2.readthedocs.org/en/2.0b4/user/models.html
        #self.method_choice = params['method_choice'].value
        # {'linear', 'euler', 'milstein', 'independent', 'exponential_euler', 'rk2', 'rk4'}
        # Threshold
        #self.threshold = params['threshold'].value # threshold='v > -50*mV'
        # Reset
        #self.reset = params['reset'].value # reset='v = v_r'
        # Refractory
        # self.refractory = params['refractory'].value #'(1 + 2*rand())*ms'

        # Avoid state variables nameclash with MSVC under Windows
        self.stateVars = RhO.brianStateVars  # ['S_'+s for s in RhO.stateVars]

        self.namespace = self.setParams(RhO)
        self.namespace.update(netParams)
        self.netParams = netParams

        #if network is None:
        #    network = self.buildNetwork(self.namespace)

        self.net = network

        #setRecords(self, GroupRec=None, IndRec=None)
        self.monitors = monitors

        self.G_RhO = G_RhO      # 'Inputs' ### Change to match RhO type?
        self.varIrho = varIrho  # 'I'
        self.varV = varV        # 'v'
        self.net[self.G_RhO].__dict__[self.stateVars[0]] = 1.  # Set rhodopsins to dark-adapted state


    def setParams(self, RhO):
        params = {}  # dict(RhO.paramsList)
        for p in RhO.paramsList:
            params[p] = RhO.__dict__[p] * modelUnits[p]
        return params


    def prepare(self, Prot):
        """Function to prepare everything for a simulation accounting for
        changes in protocol parameters.
        """
        Prot.prepare()
        dt = Prot.getShortestPeriod()
        Prot.dt = self.checkDt(dt)
        self.br.defaultclock.dt = self.dt * ms
        self.rasters = Prot.genContainer()
        self.Vms = Prot.genContainer()
        # Skip last state value since this is defined as 1 - sum(states) not as an ODE
        self.net[self.G_RhO].set_states({s: o for s, o in zip(self.stateVars[:-1], self.RhO.s_0[:-1])})  # Necessary?
        self.net.store()  # http://brian2.readthedocs.org/en/latest/user/running.html
        return  # self.dt

    def initialise(self):
        self.net.restore()

    def runTrial(self, RhO, phiOn, V, Dt_delay, cycles, dt, verbose=config.verbose):
        """Main routine for simulating a square pulse train."""

        nPulses = cycles.shape[0]

        if verbose > 1:
            if V is not None:
                Vstr = 'V = {:+}mV, '.format(V)
            else:
                Vstr = ''
            info = ("Simulating experiment at phi = {:.3g}photons/mm^2/s, "
                   "{}pulse cycles: [Dt_delay={:.4g}ms".format(phiOn, Vstr, Dt_delay))
            for p in range(nPulses):
                info += "; [Dt_on={:.4g}ms; Dt_off={:.4g}ms]".format(cycles[p, 0], cycles[p, 1])
            info += "]"
            print(info)

            report = 'text'  # 'stdout', 'stderr', function
        else:
            report = None

        # Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi)  # Reset state and time arrays from previous runs

        #self.namespace.update({'phi':phi*modelUnits['phi_m'], 'stimulus':bool(phi>0)})
        self.net[self.G_RhO].phi = phi * modelUnits['phi_m']
        #self.net[self.G_RhO].stimulus = False

        RhO.s0 = RhO.states[-1, :]   # Store initial state used
        #start, end = RhO.t[0], RhO.t[0]+Dt_delay
        #t = np.linspace(start,end,int(round(((end-start)/dt)+1)), endpoint=True) #t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))
            #if verbose > 2:
            #    print("Simulating t_del = [{},{}]".format(start,end))

        if verbose > 1:
            print(self.namespace)

        #self.br.Network.run(Dt_delay*ms, report)
        self.net.run(duration=Dt_delay*ms, namespace=self.namespace, report=report) #brianNamespace
        #self.br.run(duration=Dt_delay*ms, report=report)

        elapsed = Dt_delay

        for p in range(nPulses):

            ### Light on phase
            elapsed += cycles[p, 0]
            onInd = int(round(elapsed/dt))
            # Turn on light and set transition rates
            phi = phiOn if (phiOn > 0) else 0  # Light flux
            #self.namespace.update({'phi':phi*modelUnits['phi_m'], 'stimulus':bool(phi>0)})
            self.net[self.G_RhO].phi = phi * modelUnits['phi_m']
            #self.net[self.G_RhO].stimulus = True if (phiOn>0) else False #True
            self.net.run(duration=cycles[p, 0]*ms, namespace=self.namespace, report=report)
            RhO.ssInf.append(RhO.calcSteadyState(phi))

            ### Light off phase
            # Turn off light and set transition rates
            phi = 0  # Light flux
            self.net[self.G_RhO].phi = phi * modelUnits['phi_m']
            #self.net[self.G_RhO].stimulus = False
            #self.namespace.update({'phi':phi*modelUnits['phi_m'], 'stimulus':bool(phi>0)})
            self.net.run(duration=cycles[p, 1]*ms, namespace=self.namespace, report=report)
            elapsed += cycles[p, 1]
            offInd = int(round(elapsed/dt))

            RhO.pulseInd = np.vstack((RhO.pulseInd, [onInd, offInd]))

        # Calculate photocurrent
        '''
        There are two subtly different ways to get the values for specific neurons:
        you can either index the 2D array stored in the attribute with the variable name (as in the example above)
        or you can index the monitor itself. The former will use an index relative to the recorded neurons
        (e.g. M.v[1] will return the values for the second recorded neuron which is the neuron with the index 10
        whereas M.v[10] would raise an error because only three neurons have been recorded),
        whereas the latter will use an absolute index corresponding to the recorded group
        (e.g. M[1].v will raise an error because the neuron with the index 1 has not been recorded
        and M[10].v will return the values for the neuron with the index 10).
        If all neurons have been recorded (e.g. with record=True) then both forms give the same result.
        '''

        assert(self.monitors['states'].n_indices == 1)  # Assumes that only one neuron is recorded from...
        # Last value is not automatically recorded!
        # https://github.com/brian-team/brian2/issues/452
        self.monitors['states'].record_single_timestep()
        soln = np.hstack([self.monitors['states'].variables[s].get_value()
                          for s in self.stateVars])
        t = self.monitors['states'].t/ms
        RhO.storeStates(soln[1:], t[1:])
        states, t = RhO.getStates()

        if V is not None:  # and 'I' in self.monitors:
            I_RhO = RhO.calcI(V, RhO.states)
        else:
            self.monitors['I'].record_single_timestep()
            I_RhO = self.monitors['I'].variables[self.varIrho].get_value() * 1e9  # / nA
            I_RhO = np.squeeze(I_RhO)  # Replace with record indexing?
        self.monitors['V'].record_single_timestep()
        return I_RhO, t, states

    def runTrialPhi_t(self, RhO, phi_ts, V, Dt_delay, cycles, dt, verbose=config.verbose):
        """Main routine for simulating a pulse train."""
        # TimedArray([x1, x2, ...], dt=my_dt), the value x1 will be returned for all 0<=t<my_dt, x2 for my_dt<=t<2*my_dt etc.

        nPulses = cycles.shape[0]
        assert(len(phi_ts) == nPulses)

        times, Dt_total = cycles2times(cycles, Dt_delay)

        if verbose > 1:
            if V is not None:
                Vstr = 'V = {:+}mV, '.format(V)
            else:
                Vstr = ''
            info = "Simulating experiment {}pulse cycles: [Dt_delay={:.4g}ms".format(Vstr, Dt_delay)
            for p in range(nPulses):
                info += "; [Dt_on={:.4g}ms; Dt_off={:.4g}ms]".format(cycles[p, 0], cycles[p, 1])
            info += "]"
            print(info)

            report = 'text'  # 'stdout', 'stderr', function
        else:
            report = None

        # Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi)  # Reset state and time arrays from previous runs
        RhO.s0 = RhO.states[-1, :]   # Store initial state used

        if verbose > 1:
            print("Trial initial conditions:{}".format(RhO.s0))

        if verbose > 1:
            print(self.namespace)

        duration = Dt_total  # self.Dt_total instead?
        start, end = RhO.t[0], RhO.t[0]+Dt_delay
        nSteps = int(round(((end-start)/dt)+1))
        t = np.linspace(start, end, nSteps, endpoint=True)
        phi_tV = np.zeros_like(t)

        for p in range(nPulses):
            start = end
            Dt_on, Dt_off = cycles[p, 0], cycles[p, 1]
            end = start + Dt_on + Dt_off

            nSteps = int(round(((end-start)/dt)+1))
            tPulse = np.linspace(start, end, nSteps, endpoint=True)
            phi_t = phi_ts[p]
            phiPulse = phi_t(tPulse)  # -tPulse[0] # Align time vector to 0 for phi_t to work properly

            onInd = len(t) - 1  # Start of on-phase
            offInd = onInd + int(round(Dt_on/dt))
            RhO.pulseInd = np.vstack((RhO.pulseInd, [onInd, offInd]))

            t = np.r_[t, tPulse[1:]]
            phi_tV = np.r_[phi_tV, phiPulse[1:]]

            RhO.ssInf.append(RhO.calcSteadyState(phi_t(end-Dt_off)))

        phi_tV[np.ma.where(phi_tV < 0)] = 0  # Safeguard for negative phi values # np.clip(phi_tV, 0, phiMax)?
        phi = self.br.TimedArray(phi_tV * modelUnits['phi_m'], dt=dt*ms, name='phi')  # 'phi_t'
        self.namespace.update({'phi': phi})
        self.net.run(duration=duration*ms, namespace=self.namespace, report=report)

        # Calculate photocurrent
        assert(self.monitors['states'].n_indices == 1)  # Assumes that only one neuron is recorded
        # Last value is not automatically recorded!
        # https://github.com/brian-team/brian2/issues/452
        self.monitors['states'].record_single_timestep()
        soln = np.hstack([self.monitors['states'].variables[s].get_value()
                          for s in self.stateVars])
        t = self.monitors['states'].t/ms
        RhO.storeStates(soln[1:], t[1:])
        states, t = RhO.getStates()

        if V is not None:  # and 'I' in self.monitors:
            I_RhO = RhO.calcI(V, RhO.states)
        else:
            self.monitors['I'].record_single_timestep()
            I_RhO = self.monitors['I'].variables[self.varIrho].get_value() * 1e9  # / nA
            I_RhO = np.squeeze(I_RhO)  # Replace with record indexing?
        self.monitors['V'].record_single_timestep()

        #for mon in self.monitors:
        #    mon.record_single_timestep()
        return I_RhO, t, states

    def saveExtras(self, run, phiInd, vInd):
        # TODO: Rethink this HACK!!!
        #self.rasters[run][phiInd][vInd] = copy.copy(self.monitors['spikes'])
        self.rasters[run][phiInd][vInd] = [{'name': self.monitors['spikes'][lay].name,
                                            't': self.monitors['spikes'][lay].t/ms,
                                            'i': self.monitors['spikes'][lay].i,
                                            'n': len(self.monitors['spikes'][lay].spike_trains())} for lay in range(len(self.monitors['spikes']))]
        self.Vms[run][phiInd][vInd] = copy.copy(self.monitors['V'])

    def plotExtras(self):
        Prot = self.Prot
        RhO = self.RhO
        for run in range(Prot.nRuns):
            cycles, Dt_delay = Prot.getRunCycles(run)
            pulses, Dt_total = cycles2times(cycles, Dt_delay)
            for phiInd, phiOn in enumerate(Prot.phis):
                for vInd, V in enumerate(Prot.Vs):
                    col, style = Prot.getLineProps(run, vInd, phiInd)
                    Vm = self.Vms[run][phiInd][vInd]
                    spikes = self.rasters[run][phiInd][vInd]  # TODO: Change name

                    figName = '{}Vm{}s-{}-{}-{}'.format(Prot.protocol, RhO.nStates, run, phiInd, vInd)
                    Vmonitor = self.Vms[run][phiInd][vInd]
                    self.plotVm(Vmonitor=Vmonitor, times=pulses,
                                Dt_total=Dt_total, offset=Dt_delay,
                                figName=figName)

                    figName = '{}Spikes{}s-{}-{}-{}'.format(Prot.protocol, RhO.nStates, run, phiInd, vInd)
                    spikeMonitors = self.rasters[run][phiInd][vInd]
                    self.plotRasters(spikeSets=spikeMonitors, times=pulses,
                                     Dt_total=Dt_total, offset=Dt_delay,
                                     figName=figName)
        return

    def plotVm(self, Vmonitor, times=None, Dt_total=None, offset=0, figName=None):  # TODO: REVISE (esp. offset)!!!
        Prot = self.Prot
        RhO = self.RhO
        Vfig = plt.figure()
        axV = Vfig.add_subplot(111)
        Vm = Vmonitor.variables[self.varV].get_value() * 1000  # Convert to mV # HACK #self.monitors['V']
        t = Vmonitor.t/ms - offset
        # times -= offset # Do not edit in place or it causes errors in subsequent functions!
        plt.plot(t, Vm, 'g')
        plotLight(times - offset, axV)
        axV.set_xlim((-offset, Dt_total - offset))

        plt.ylabel(r'$\mathrm{Membrane\ Potential\ [mV]}$')
        if Dt_total is None:
            Dt_total = Prot.Dt_total

        plt.xlabel(r'$\mathrm{Time\ [ms]}$', position=(config.xLabelPos, 0), ha='right')

        #setCrossAxes(axV, zeroX=False)  # Zeroing x-axis caused the plot to become too big!
        # #axV.spines['bottom'].set_position('zero') # x-axis # Caused the plot to be too big!
        setCrossAxes(axV, zeroX=False, zeroY=True)
        axV.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axV.grid(visible=True, which='minor', axis='both', linewidth=.2)
        axV.grid(visible=True, which='major', axis='both', linewidth=1)

        if len(Prot.PD.legLabels) > 0 and Prot.PD.legLabels[0] != '':
            plt.legend(Prot.PD.legLabels)

        # axV.set_ylim(axV.get_ylim())
        # ymin, ymax = axV.get_ylim()

        plt.tight_layout()
        plt.show()

        if figName is None:
            figName = '{}Vm{}s'.format(Prot.protocol, RhO.nStates)
        fileName = os.path.join(config.fDir, figName+'.'+config.saveFigFormat)
        logger.info('Saving Membrane potential figure to: {}'.format(fileName))
        Vfig.savefig(fileName, format=config.saveFigFormat)

    def plotRasters(self, spikeSets, times=None, Dt_total=None, offset=0, figName=None):
        """Plot spike rasters for each layer."""

        nLayers = len(spikeSets)  # self.monitors['spikes']

        if Dt_total is None:
            Dt_total = self.Prot.Dt_total

        Rfig = plt.figure()
        gs = plt.GridSpec(nLayers, 1)
        axes = [None for lay in range(nLayers)]

        for lay in range(nLayers):

            axes[lay] = Rfig.add_subplot(gs[lay])
            gInd = nLayers-lay-1
            plt.plot(spikeSets[gInd]['t'] - offset, spikeSets[gInd]['i'], '.')
            plt.ylabel(r'$\mathrm{' + spikeSets[gInd]['name'] + '}$')
            if gInd > 0:
                plt.setp(axes[lay].get_xticklabels(), visible=False)
            else:
                plt.xlabel(r'$\mathrm{Time\ [ms]}$')

            if times is not None:
                plotLight(times - offset, axes[lay])

            axes[lay].set_ylim((0, spikeSets[gInd]['n']))
            axes[lay].set_xlim((-offset, Dt_total - offset))

        plt.tight_layout()
        plt.show()

        if figName is None:
            figName = '{}Spikes{}s'.format(self.Prot.protocol, self.RhO.nStates)
        fileName = os.path.join(config.fDir, figName+'.'+config.saveFigFormat)
        logger.info('Saving Raster plot figure to: {}'.format(fileName))
        Rfig.savefig(fileName, format=config.saveFigFormat)


simulators = OrderedDict([('Python', simPython),
                          ('NEURON', simNEURON),
                          ('Brian', simBrian)])

# TODO: Use this for loading components of the GUI
# TODO: Check on initialising sim objects and raise an error if unavailable
simAvailable = OrderedDict([(sim, simAvailable(sim)) for sim in simulators])
