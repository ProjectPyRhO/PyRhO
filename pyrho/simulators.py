# simulators.py

from .parameters import *
from .models import *
from .config import verbose
import numpy as np

class Simulator(object):
    """Common base class for all simulators"""
    
    def __init__(self, simulator='Python'):
        self.simulator = simulator

    
    def __str__(self):
        return "Simulator type: "+self.simulator
        

class simPython(Simulator):
    simulator = 'Python'
    """Class for channel level simulations with Python"""
    
    def __init__(self, RhO, params=simParams['Python']):
        self.dt = params['dt'].value
    
    def prepare(self, dt):
        """Function to compare simulator's timestep to the timestep required by the protocol"""
        if dt < self.dt:
            self.dt = dt #min(self.h.dt, dt)
            print('Time step set to {}ms by protocol!'.format(self.dt))
    
    def runTrial(self, RhO, nPulses, V,phiOn,delD,onD,offD,padD,dt): #dt; a1,a3,b2,b4,I_RhO
        """Main routine for simulating a square pulse train"""
        
        if verbose >= 0:
            print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse: [delD={:.4g}ms; onD={:.4g}ms; offD={:.4g}ms]".format(V,phiOn,delD,onD,offD+padD))
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        start, end = 0.00, delD
        t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        soln = odeint(RhO.solveStates, RhO.s_0, t_del, args=(None,), Dfun=RhO.jacobian) #delay
        t = t_del
        RhO.storeStates(soln,t)
        
        for p in range(0, nPulses):
            
            ### Light on phase
            RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            end = start + onD
            t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
            onInd = len(RhO.t) - 1  # Start of on-phase
            offInd = onInd + len(t) - 1 # Start of off-phase
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
            # Turn on light and set transition rates
            phi = phiOn  # Light flux
            RhO.setLight(phi)
            if verbose > 1:
                print("On-phase initial conditions:{}".format(RhO.s_on))
            soln = odeint(RhO.solveStates, RhO.s_on, t, args=(None,), Dfun=RhO.jacobian)
            RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
            ### Light off phase
            RhO.s_off = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
            start = end
            end = start + offD
            if (p+1) == nPulses: # Add (or subtract) extra time after (during) the off phase
                end += padD
            t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # endpoint=True
            # Turn off light and set transition rates
            phi = 0  # Light flux
            RhO.setLight(phi)
            if verbose > 1:
                print("Off-phase initial conditions:{}".format(RhO.s_off))
            soln = odeint(RhO.solveStates, RhO.s_off, t, args=(None,), Dfun=RhO.jacobian)
            RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
            
        ### Calculate photocurrent
        I_RhO = RhO.calcPhotocurrent(V, RhO.states)
        states,t = RhO.getStates()
        
        return I_RhO, t, states
        
        
        
    def runTrialPhi_t(self, RhO,V,phi_t,delD,stimD,totT,dt): #dt; a1,a3,b2,b4,I_RhO
        """Main routine for simulating a pulse train"""
        # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
        #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
        
        ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
        ### Combine this with the original runTrial
        
        ### Delay phase (to allow the system to settle)
        phi = 0
        RhO.initStates(phi) # Reset state and time arrays from previous runs
        start, end = 0.00, delD
        t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
        soln = odeint(RhO.solveStates, RhO.s_0, t_del, args=(None,), Dfun=RhO.jacobian) #delay # odeint(RhO.solveStates, RhO.s_0, t_del, Dfun=RhO.jacobian)
        t = t_del
        RhO.storeStates(soln,t)

        ### Stimulation phase
        RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
        start = end
        end = totT #start + stimD
        t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
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
        I_RhO = RhO.calcPhotocurrent(V, RhO.states)
        states,t = RhO.getStates()
        
        # if verbose:
            # print("Run={}, I_min={}, I_max={}, label={}".format(run,I_RhO.min(),I_RhO.max(),label))
        
        return I_RhO, t, states

        
        
        
        
class simNEURON(Simulator):
    
    """Class for cellular level simulations with NEURON"""
    simulator = 'NEURON'
    mechanisms = {3:'RhO3', 4:'RhO4', 6:'RhO6'}
    #mods = {3:'RhO3.mod', 4:'RhO4.mod', 6:'RhO6.mod'}
    
    def __init__(self, RhO, params=simParams['NEURON'], recInd=0): #v_init=-70, integrator='fixed'):
        
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
        
        
        
        #import neuron
        from neuron import h
        #self.neuron = neuron
        self.h = h
        
        self.h.load_file('stdrun.hoc')
        
        ### Simulation control
        self.CVODE = params['CVODE'].value
        
        if self.CVODE == True: #self.integrator == 'variable': #... if dt <= 0:
            self.h('objref cvode')
            self.h.cvode = self.h.CVode()
            self.h.cvode.active(1)
        else: #elif 'fixed':
            #self.h.cvode.active(0)
            self.h.dt = params['dt'].value
        #else:
        #    raise ValueError('Unknown integrator!')
        print('Integrator tolerances: absolute=',self.h.cvode.atol(),' relative=',self.h.cvode.rtol()) ### Set as parameters
        #self.h.cvode.atol(0.000001)
        #self.h.cvode.rtol(0.000001)
        
        self.h.v_init = params['v_init'].value #-60
        
        #expdist = 600 ### Redundant
        
        self.buildCell(params['cell'].value)
        mod = self.mechanisms[RhO.nStates] ### Use this to select the appropriate mod file for insertion
        self.transduce(RhO, expProb=params['expProb'].value)
        self.setRhodopsinParams(self.rhoList, RhO, modelParams[str(RhO.nStates)])
        #print(recInd)
        #for sec in self.cell:
        #    print(sec)
        #print(self.cell.count())
        #print(self.rhoList.count())
        self.rhoRec = self.rhoList[recInd] #self.compList.o(recInd) # Choose a rhodopsin to record
        self.setRecords(self.rhoRec, RhO)
        self.Vclamp = params['Vclamp'].value
        if self.Vclamp == True:
            self.addVclamp() # params['Vhold'].value
        # Run then plot
        
        
    def prepare(self, dt):
        """Function to compare simulator's timestep to the timestep required by the protocol"""
        if dt < self.h.dt:
            self.h.dt = dt #min(self.h.dt, dt)
            print('Time step set to {}ms by protocol!'.format(self.h.dt))
        
        
    def buildCell(self, scriptList=[]):
        """Pass a list of hoc or python files in the order in which they are to be processed"""
        #import os
        import subprocess
        for script in scriptList:
            #name, ext = os.path.splitext(script)
            ext = script.split(".")[-1]
            script = 'PyRhO/pyrho/NEURON/'+script #.join(script)
            if ext == 'hoc':
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
        for sec in self.h.allsec():
            self.cell.append(sec)
        #print('Total sections:',self.cell.count())
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
            if random.random() <= expProb: # Insert a rhodopsin and append it to a rhoList
                self.compList.append(sec)
                self.h('rho = new {}(0.5)'.format(mech))
                self.rhoList.append(self.h.rho)
        
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
    
    
    def setRecords(self,rhoRec,RhO):
        self.h('objref Vsoma') #objref vsoma, i_ChR, tvec
        self.h('objref Iphi')
        self.h('objref tvec')
        self.h.Vsoma = self.h.Vector() #vsoma = new Vector()
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
    
    def setVclamp(self,Vhold):
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
    
    def runTrial(self, RhO, nPulses, V, phiOn, delD, onD, offD, padD, dt):
        
        # Notes on the integrator
        # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=8&t=1330
        # http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html
        
        # Model events: http://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html#modeldescriptionissues-events
        
        RhO.initStates(phi=0) # Reset state and time arrays from previous runs
        
        # Set simulation run time
        self.h.tstop = delD + nPulses*(onD+offD) + padD
        
        if self.Vclamp == True:
            self.setVclamp(V)
        
        self.setPulses(phiOn, delD, onD, offD, nPulses)

        if verbose >= 0:
            print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse: [delD={:.4g}ms; onD={:.4g}ms; offD={:.4g}ms]".format(V,phiOn,delD,onD,offD+padD))        
                
        self.h.init() #self.neuron.init()
        
        self.h.run() #self.neuron.run(self.h.tstop)
            
        I_RhO = np.array(self.h.Iphi.to_python(), copy=True)
        t = np.array(self.h.tvec.to_python(), copy=True)
        
        for p in range(nPulses):
            onInd = np.searchsorted(t,delD+p*(onD+offD),side="left")
            offInd = np.searchsorted(t,delD+p*(onD+offD)+onD,side="left")
            RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
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
        
        # Reset vectors
        self.h.tvec.clear() # Necessary?
        self.h.tvec.resize(0)
        self.h.Iphi.clear()
        self.h.Iphi.resize(0)
        
        
        
        
        return I_RhO, t, soln


    def runTrialPhi_t(self, RhO,V,phi_t,delD,stimD,totT,dt): 
        # Make RhO.mod a continuous function of light - 9.12: Discontinuities 9.13: Time-dependent parameters p263
            #- Did cubic spline interpolation ever get implemented for vec.play(&rangevar, tvec, 1)?
            #- http://www.neuron.yale.edu/neuron/faq#vecplay
            #- http://www.neuron.yale.edu/neuron/static/new_doc/programming/math/vector.html#Vector.play
            #- Ramp example: http://www.neuron.yale.edu/phpbb/viewtopic.php?f=15&t=2602
        raise NotImplementedError('Error: Light as a continuous function of time has not been implemented for NEURON yet!')
        return I_RhO,t,soln
    
    
    
    
    
    
    ### Unused functions
    
    def go():
        h.finitialize()

        g.begin()
        while h.t<tstop:
            h.fadvance()
            # update the graph while simulating
            g.plot(h.t)

        g.flush()
    
        
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




        
class simBrian(Simulator):
    simulator = 'Brian'
    """Class for network level simulations with Brian"""
    
    def __init__(self, RhO, params=simParams['Brian']):
        pass
        

        
simulators = {'Python': simPython, 'NEURON': simNEURON, 'Brian': simBrian}