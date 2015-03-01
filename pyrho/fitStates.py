# Call this script in IPython
# %run fitStates.py
# import (as a module)
# %import
# %load

# import lmfit
from lmfit import minimize, Parameters, fit_report
import numpy as np
import matplotlib.pyplot as plt
from .parameters import *
from .loadData import * #import loadData
from .protocols import * # for SSA and numerical solution
#import models # for calculating analytic solution
#global phi0
from .models import * # for fitPeaks
from .config import verbose, saveFigFormat, addTitles

#import protocols

methods=('leastsq','nelder','lbfgsb','powell','cg','newton','cobyla','tnc','trust-ncg','dogleg','slsqp')
meth = methods[3] #3
# http://scipy-lectures.github.io/advanced/mathematical_optimization/#choosing-a-method

### Alternative optimisation toolkits
# Stimfit: http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00016/full
# NeuroTools: http://neuralensemble.org/NeuroTools/
# http://mdolab.engin.umich.edu/sites/default/files/pyOpt.pdf
# CVXopt
# OpenOpt
# Optimizer
# scikit-learn




# ### Temporary code...
# if ChR.nStates != 6:
    # gmax = ChR.g
# else:
    # gmax = ChR.gbar * A
# #######################


# gmax = dataset['saturate'].gmax
# Gr0 = 1/dataset['varyIPI'].tau_r # Gr,dark
# I
# RhO = fit3states(I,t,onInd,offInd,dataset,gmax,Ipeak)

def calc3on(p,t):
    """Fit a biexponential curve to the on-phase to find lambdas"""
    a0 = p['a0'].value
    a1 = p['a1'].value
    a2 = p['a2'].value
    tau_act = p['tau_act'].value
    tau_deact = p['tau_deact'].value
    return -(a0 + a1*(1-np.exp(-t/tau_act)) + a2*np.exp(-t/tau_deact))
    
def resid3on(p,I,t):
    return I - calc3on(p,t)
    
def calc3off(p,t):
    A = p['A'].value
    Gd = p['Gd'].value
    return -(A*np.exp(-Gd*t))

def resid3off(p,I,t):
    return I - calc3off(p,t)

def calc3off2exp(p,t):
    A = p['A'].value
    B = p['B'].value
    C = p['C'].value
    Gd1 = p['Gd1'].value
    Gd2 = p['Gd2'].value
    return -(A*np.exp(-Gd1*t)+B*np.exp(-Gd2*t)+C)

def resid3off2exp(p,I,t):
    return I - calc3off2exp(p,t)
    
# def calc3offLin(p,t):
    # B = p['B'].value
    # Gd = p['Gd'].value
    # return B -Gd*t
    
# def resid3offLin(p,I,t):
    # return I - calc3offLin(p,t)
    
def calcG(photocurrent, E):
    """Function to calculate a lower bound on the cell's maximum conductance from its peak current
    Ipmax :=  Peak current [nA]
    V     :=  Clamp Voltage [mV]
    return gmax [pS]"""
    ### This explicitly depends upon the reversal potential
    ### Also depends on fV? Or just has to be done at -70mV?
    ### gmax = Ipeak/([O_ss]*(V-E)) # assuming [O_ss] = 1 ### Should this be O_p?
    ### This assumption is an underestimate for 4 & 6 state models: [O_ss] =~ 0.71 (depending on rates)
    if hasattr(photocurrent, 'gmax'):
        gmax = photocurrent.gmax
    else: # findPeaks
        if (photocurrent.V < E): # Negative potential - Find Minima ### if abs(min(photocurrent.I)) > abs(max(photocurrent.I)): 
            Ipmax = min(photocurrent.I)
        else:       # Positive potential - Find Maxima
            Ipmax = max(photocurrent.I)
        gmax = Ipmax/(photocurrent.V-E)
        photocurrent.gmax = gmax
    return gmax * (1e6) # nA / mV = pS * (1e6)

def calcGr0(dataset):
    if hasattr(dataset['varyIPI'], 'tau_r'):
        Gr0 = 1/dataset['varyIPI'].tau_r # Gr,dark
    else:
        print("Extract the peaks and fit an exponential...")
        
    return Gr0
    
def fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,Iss=None): ### Modify to pass in p3s
    """
    I       := Photocurrent to fit [nA]
    t       := Timepoints corresponding to I samples [ms]
    onInd   := Index for I and t arrays corresponding to the start of the light pulse
    offInd  := Index for I and t arrays corresponding to the end of the light pulse
    phi     := Flux intensity at which the photocurrent was recorded [photons s^-1 mm^-2]
    V       := Clamp voltage at which the photocurrent was recorded [mV]
    Gr0     := Dark recovery transition rate [ms^-1]
    gmax    := Maximum conductance for rhodopsins in the neuron [nS]
    Ipmax   := Maximum Peak photocurrent [nA]
    Iss     := 
    
    The limitations of the 3-state model are such that it only holds for a particular value of phi"""
    
    p3s = Parameters()
    
    # Three-state model parameters
    # self.Gr_dark = 1/5000      # [ms^-1] tau_r,dark = 5-10s p405 Nikolic et al. 2009
    # self.Gr_light = 1/165      # [ms^-1] tau_r,light
    # self.Gd = 1/11             # [ms^-1] @ 1mW mm^-2
    # self.eF = self.set_eF(phi) # [ms^-1] Quantum efficiency * number of photons absorbed by a ChR2 molecule per unit time
    # self.g = 100e-3 * 1.67e5   # [pS] g1 (100fS): Conductance of an individual ion channel * N (100000)
    ## self.c = 1/165             # [ms^-1] (Gr,light @ 1mW mm^-2)
    # self.e = 0.5*1.2e-8*1e-6/1.1 # Coefficient for the activation rate

    #  set_eF(self, phi):
    #   return self.e * phi #0.5 * phi * (1.2e-8 * 1e-6) / 1.1  # eps * phi * sigma_ret (mm^-2) / wloss

    #  set_Gr(self, phi):
    #   return self.Gr_dark + self.Gr_light
    

    
    ### 1. Load data for 'saturate' protocol to find Ipmax in order to calculate gmax
    ### gmax = Ipmax/([O_p]*(V-E)) = Ipmax/(V-E) # with the assumption [O_p] = 1
    ### This assumption is an underestimate for 4 & 6 state models: [O_p] =~ 0.71 (depending on rates)
    # if hasattr(dataset['saturate'], 'gmax'):
        # gmax = dataset['saturate'].gmax
    # else: # findPeaks
        # if (dataset['saturate'].V < E): # Find Minima
            # Ipmax = min(dataset['saturate'].I_phi)
        # else:       # Find Maxima
            # Ipmax = max(dataset['saturate'].I_phi)
        # gmax = Ipmax/(dataset['saturate'].V-E)
        # dataset['saturate'].gmax = gmax
    
    p3s.add('g', value=gmax, vary=False) # derived['gmax'] = True
    # Change the model to be consistent so that g = gbar * A
    
    
    ### 2. Fit exponential to peak recovery plots
    # if hasattr(dataset['varyIPI'], 'tau_r'):
        # Gr0 = 1/dataset['varyIPI'].tau_r # Gr,dark
    # else:
        # print("Extract the peaks and fit an exponential...")
    
    p3s.add('Gr0', value=Gr0, vary=False)
    
    ### 3. Use Optimisation algorithms to fit the remaining parameters
    
    # Divide curve into on and off phases
    # Replace indices with time-value matching for real data
    # t[onInd] <= t_on < t[offInd+1]
    # t[offInd] <= t_off 
    Ion, Ioff = I[onInd:offInd+1], I[offInd:]
    ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    ### 3a. Fit exponential to off curve to find Gd
    ### Fit off curve
    pOff = Parameters() # Create parameter dictionary
    #print(-Ioff[0])
    
    ### Original single exponential fit
    #pOff.add('A', value=-Ioff[0])#vary=False) #-1 # Constrain fitting to start at the beginning of the experimental off curve
    #pOff.add('Gd', value=0.1, min=0)
    #offP = minimize(resid3off,pOff,args=(Ioff,toff),method=meth)
    
    pOff.add('A', value=-1) # Islow
    pOff.add('Gd1', value=0)#, min=0))
    pOff.add('B', value=-1) # Ifast
    pOff.add('Gd2', value=100)#, min=0))
    pOff.add('C', value=0)
    #offP = minimize(resid3off,pOff,args=(Ioff[0:100],toff[0:100]),method=meth)
    offP = minimize(resid3off2exp,pOff,args=(Ioff,toff),method=meth)
    #pOff.add('Gd', value=max(pOff['Gd1'].value,pOff['Gd2'].value), min=0)
    nf = pOff['A'].value + pOff['B'].value
    pOff.add('Gd', value=pOff['A'].value*pOff['Gd1'].value/nf+pOff['B'].value*pOff['Gd2'].value/nf, min=0)
    
    print(">>> Off-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOff)) # offP?
    print(offP.message)
    print("Chi^2 (Off-phase): ",offP.chisqr)
    p3s.add('Gd', value=pOff['Gd'].value, vary=False)
    
    #print(pOff['Gd'].value)
    #print(p3s['Gd'].value)
    
    ### 3b. Optimise to find eF and Gr at a particular phi from the on curve
    ### Fit on curve
    pOn = Parameters()
    pOn.add('a0', value=-1, expr='-a2')
    pOn.add('a1', value=1)
    pOn.add('a2', value=1)
    pOn.add('tau_act', value=2, min=0)
    pOn.add('tau_deact', value=10, min=0)
    onP = minimize(resid3on,pOn,args=(Ion,ton),method=meth)
    
    # Calculate model parameters from fitting parameters
    tau_act = pOn['tau_act'].value
    tau_deact = pOn['tau_deact'].value
    lam1 = 1/tau_deact
    lam2 = 1/tau_act
    
    Gd = pOff['Gd'].value
    
    # findPeaks in experimental photocurrent
    # if (V < E): # Find Minima
        # Ip = min(I)
    # else:       # Find Maxima
        # Ip = max(I)
    # Oss = Iplat/(gmax*(Vsaturate-E)) # f(V) here?
    Iplat = findPlateauCurrent(Ion,p=0.1)
    ###Ipmax = gmax*(V-E) # This is Ipmax (from saturate)
    Oss = Iplat/Ipmax
    alpha = (lam1+lam2)/Gd - 1
    beta = alpha/((1/Oss)-1)
    
    if (alpha**2 - 4*beta) < 0:
        print('\n\nEarly Warning! No real solution exists!\n\n')
        ### Rerun off-curve fitting with single exponential
    
    x1 = (alpha + np.sqrt(alpha**2 - 4*beta))/2
    Gr1 = x1 * Gd
    P1 = lam1+lam2 - (Gr1+Gd)
    print("Solution 1: eF = {}; Gd = {}; Gr = {}".format(P1,Gd,Gr1))
    
    x2 = (alpha - np.sqrt(alpha**2 - 4*beta))/2
    Gr2 = x2 * Gd
    P2 = lam1+lam2 - (Gr2+Gd)
    print("Solution 2: eF = {}; Gd = {}; Gr = {}".format(P2,Gd,Gr2))
    
    if P1 > Gr1:
        print("==> Solution 1 (assuming that eF > Gr): {{eF={:.3}, Gd={:.3}, Gr={:.3}}}".format(P1,Gd,Gr1))
        P = P1
        Gr = Gr1
    else:
        print("==> Solution 2 (assuming that eF > Gr): {{eF={:.3}, Gd={:.3}, Gr={:.3}}}".format(P2,Gd,Gr2))
        P = P2
        Gr = Gr2
    #print("==> Original values: {{eF={:.3}, Gd={:.3}, Gr={:.3}}}".format(  ,Gd,ChR.Gr))
    
    if (Gr + Gd - 2*np.sqrt(Gr*Gd)) < P < (Gr + Gd + 2*np.sqrt(Gr*Gd)):
        print('\n\nWarning! No real solution exists!\n\n')
    
    print("\n\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOn)) # onP?
    print(onP.message)
    print("Chi^2 (On-phase): ",onP.chisqr)
    
    k = P/phi # To find the coefficient which scales phi to calculate the activation rate
    p3s.add('k', value=k, vary=False)
    # N.B. The three-state model does not generalise well to other values of phi
    
    GrPhi = Gr - Gr0 # To find Grlight for a particular phi
    p3s.add('GrPhi', value=GrPhi, vary=False)

#     ci = lmfit.conf_interval(pOff)
#     lmfit.printfuncs.report_ci(ci)
    

    ### 4. Optionally fit f(V) parameters with inwardRect data
    # if 'inwardRect' in dataset:
        # if hasattr(dataset['inwardRect'], 'Iss'): # Use extracted values
            # print("Finish me!")
        # else: # Extract steady state values
            # print("Finish me!")
        # # Fit Curve of V vs Iss
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?
    
    if verbose > 1:
        for key, value in p3s.items() :
            print (key, value)
    
    RhO = selectModel(3)
    #RhO.setParams(d3sp)
    RhO.g = p3s['g'].value              # 16700     [pS]
    RhO.k = p3s['k'].value              # 0.545e-14 [ms^-1 * photons^-1 * s * mm^2]
    RhO.Gr0 = p3s['Gr0'].value          # 1/5000    [ms^-1] #Gr_dark
    RhO.Gr1 = p3s['GrPhi'].value        # 1/165     [ms^-1] #Gr_light
    RhO.Gd = p3s['Gd'].value            # 1/11      [ms^-1]
    RhO.phiFit = phi                    # Flux intensity at which the parameters were fit. 
    RhO.setLight(0)                     # Re-initialise model to dark state
    RhO.reportParams()
    
    ### Plot curves
    totT = max(t)
    Ifig = plt.figure()
    gsPL = plt.GridSpec(4,1)
    axFit = Ifig.add_subplot(gsPL[:-1,:])
    
    #ax = Ifig.add_subplot(111)
    # for p in range(0, nPulses):
        # plt.axvspan(delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD,facecolor='y',alpha=0.2)
    plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    plt.xlim((0,totT))
    #plt.xlabel('$\mathrm{Time\ [ms]}$')
    plt.setp(axFit.get_xticklabels(), visible=False)
    plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    plt.plot(t,I,color='b',label='$\mathrm{Experimental\ Data}$')
    RhO.setLight(phi)
    if addTitles:
        plt.title('Three-state model fit to data (phi={:.3g}): [eF={:.3g}; Gd={:.3g}; Gr={:.3g}]'.format(phi,RhO.eF,RhO.Gd,RhO.Gr))
    RhO.setLight(0)
    
    
    if 1: # Use pseudo parameters
        IfitOn = calc3on(pOn, ton) # onP
        pOff['A'].value = -IfitOn[-1] # Constrain off curve to start at the end of the fitted on curve
        IfitOff = calc3off(pOff, toff) # offP
        ###pOff.add('A',value=np.exp(pOff['B'].value),vary=False)
    elif 0: # Use parameters
        I_RhO, t, states = runTrial(RhO, 1, V,phi,t[onInd],t[offInd]-t[onInd],totT-t[offInd],0,0.01) # Problem matching time arrays
        #I_RhO, t, states = runTrial(RhO, 1, V,phi,round(t[onInd]),round(t[offInd]-t[onInd]),round(totT-t[offInd]),0,0.01)
        IfitOn = I_RhO[onInd:offInd+1]
        IfitOff = I_RhO[offInd:]
        #Prot = protocols['custom']([phi], [-70], [[t[onInd],t[offInd]]], totT, 1, 0.1)
        #Prot.runProtocol(RhO)
        #Prot.plotProtocol()
        #IfitOn = Prot.Is[0][0][0][onInd:offInd]
        #IfitOff = Prot.Is[0][0][0][offInd:]
    else: # Use analytic solution
        #print(ton)
        RhO.setLight(phi)
        RhO.reportState()
        onStates = RhO.calcSoln(ton)
        IfitOn = RhO.calcPhotocurrent(V, onStates)
        print(onStates)
        RhO.setLight(0)
        #print(toff)
        IfitOff = RhO.calcPhotocurrent(V, RhO.calcSoln(toff, onStates[-1,:]))
    
#     tfit = np.append(ton,toff)
#     plot(tfit,Ifit,color='r')
    plt.plot(ton+t[onInd],IfitOn,color='g',label='$\mathrm{Three-state\ model\ fit}$')
    plt.plot(toff+t[offInd],IfitOff,color='g')
    plt.legend(loc='best')
    

    axRes = Ifig.add_subplot(gsPL[-1,:],sharex=axFit)
    
    # axLag.set_aspect('auto')
                    
    # Rfig = plt.figure() # Replace this with a subplot...
    #print(IfitOn)
    #print(Ion)
    #print((IfitOn-Ion)/Ion)
    plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    #plt.plot(toff+t[offInd],Ioff-IfitOff)
    plt.ylabel('$\mathrm{Residuals}$') # % relative error')
    plt.xlabel('$\mathrm{Time\ [ms]}$')
    #plt.setp(axRes.get_xticklabels(), visible=False)
    plt.xlim((0,totT))
    plt.tight_layout()
    

    print("Parameters have been fit for the three-state model at a flux of {:.3g} [photons * s^-1 * mm^-2]".format(phi))
    
    return RhO


    
def calc4off(p,t):
    """Fit a biexponential curve to the off-phase to find lambdas"""
    a0 = p['a0'].value
    a1 = p['a1'].value
    a2 = p['a2'].value
    lam1 = p['lam1'].value
    lam2 = p['lam2'].value #####*pOff
    return -(a0 + a1*(1-np.exp(-lam1*t)) + a2*np.exp(-lam2*t))
    #return -(a0 + a1*np.exp(-lam1*t) + a2*np.exp(-lam2*t)) #
    
def resid4off(p,I,t):
    return I - calc4off(p,t)
    
def calc4offPP(p,t,RhO,V):
    
    s_off = RhO.states[-1,:]
    RhO.initStates(0)
    # Override light-sensitive transition rates
    RhO.Ga1 = 0
    RhO.Ga2 = 0
    RhO.e12 = 0.01
    RhO.e21 = 0.01
    
    RhO.Gd1 = p['Gd1'].value
    RhO.Gd2 = p['Gd2'].value
    soln = odeint(RhO.solveStates, s_off, t[:offInd+1], Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:],t[1:])
    
    I_RhO = RhO.calcPhotocurrent(V, RhO.states)
    
    return I_RhO
    
    
def calc4PP(p,t,offInd,RhO,V,phi):
    print('.', end="") # sys.stdout.write('.')
    
    #print(t.shape, offInd)
    RhO.initStates(phi)
    #RhO.setLight(phi)
    # Override light-sensitive transition rates
    RhO.Ga1 = p['Ga1'].value
    RhO.Ga2 = p['Ga2'].value
    RhO.e12 = p['e12'].value
    RhO.e21 = p['e21'].value
    
    RhO.Gd1 = p['Gd1'].value
    RhO.Gd2 = p['Gd2'].value
    soln = odeint(RhO.solveStates, RhO.s_0, t[:offInd+1], Dfun=RhO.jacobian)
    RhO.storeStates(soln,t)
    #Ion = RhO.calcPhotocurrent(V, soln)
    #print(soln.shape)
    
    
    RhO.setLight(0)
    # Override light-sensitive transition rates
    RhO.Ga1 = 0
    RhO.Ga2 = 0
    RhO.e12 = 0.01
    RhO.e21 = 0.01
    
    RhO.s_off = soln[-1,:]
    soln = odeint(RhO.solveStates, RhO.s_off, t[offInd:], Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:],t[1:])
    #Ioff = RhO.calcPhotocurrent(V, soln)
    #print(soln.shape)
    
    
    #I_RhO = np.concatenate((Ion, Ioff[1:]))
    I_RhO = RhO.calcPhotocurrent(V, RhO.states)
    
    return I_RhO
    
def resid4PP(p,I,t,offInd,RhO,V,phi):
    #print(I.shape,t.shape,offInd)
    #print(p)
    #print(RhO.reportParams())
    #print(V,phi)
    #print(calc4PP(p,t,offInd,RhO,V,phi).shape)
    return I - calc4PP(p,t,offInd,RhO,V,phi)
    
    
def calc4onPP(p,t,RhO,V,phi):
    if verbose > 1:
        print('.', end="") # sys.stdout.write('.')
    #RhO.setLight(phi)
    RhO.setParams(p)
    # Override light-sensitive transition rates
    RhO.Ga1 = p['Ga1'].value
    RhO.Ga2 = p['Ga2'].value
    RhO.e12 = p['e12'].value
    RhO.e21 = p['e21'].value
    
    RhO.Gd1 = p['Gd1'].value
    RhO.Gd2 = p['Gd2'].value
    soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
    #RhO.storeStates(soln, t) # Assumes no previous storage (i.e. not [-1,:]
    I_RhO = RhO.calcPhotocurrent(V, soln)
    return I_RhO
    
def resid4onPP(p,I,t,RhO,V,phi):
    return I - calc4onPP(p,t,RhO,V,phi)
    
    

def calc4on(p,t,RhO,V,phi):
    print('.', end="") # sys.stdout.write('.')
    #RhO.setParams(p)
    RhO.k1 = p['k1'].value # Ga1 = k1*phi
    RhO.k2 = p['k2'].value # Ga2 = k2*phi
    RhO.Gd1 = p['Gd1'].value
    RhO.Gd2 = p['Gd2'].value
    RhO.e12d = p['e12d'].value
    RhO.e21d = p['e21d'].value
    RhO.c1 = p['c1'].value
    RhO.c2 = p['c2'].value
    # nPulses = 1
    # delD = 0
    # onD = max(t) # ton
    # offD = 0
    # padD = 0
    # I_RhO, _, _ = runTrial(RhO,nPulses,V,phi,delD,onD,offD,padD)
    soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
    I_RhO = RhO.calcPhotocurrent(V, soln)
    return I_RhO
    
def resid4on(p,I,t,RhO,V,phi):
    return I - calc4on(p,t,RhO,V,phi)
    
    
def fit4states(I,t,onInd,offInd,phi,V,Gr,gmax):#,Iss): # ,Ipeak
    # Four-state model parameters
    # self.Gr = 400 # [ms^-1] ==> tau_r = 2.5s
    # self.Gd1 = 0.11 # [ms^-1]
    # self.Gd2 = 0.025 # [ms^-1]
    # self.c1 = 0.03
    # self.c2 = 0.0115
    # self.e12d = 0.01 # [ms^-1]
    # self.e21d = 0.015 # [ms^-1]
    # self.g = 100e-3 * 1.67e5   # [pS] g1 (100fS): Conductance of an individual ion channel * N (~150000)
    
    #  set_Ga1(self, phi):
    # N.B. making Ga a function of time (as in Appendix 1) results in the Six-state model
    # Gai = ei * F * f(t,tChR) See App 1
    # Ga = 1/tauChR
    #    e = 0.5
    #    sigma_ret = 1.2e-8 * 1e-6 # Convert from m^2 to mm^2
    #    w_loss = 1.1
    #    return e*phi*sigma_ret / w_loss
    
    #  set_Ga2(self, phi):
    #    e = 0.15
    #    sigma_ret = 1.2e-8  * 1e-6 # Convert from m^2 to mm^2
    #    w_loss = 1.1
    #    return e*phi*sigma_ret / w_loss
    
    #  set_e12(self, phi):
    #    return self.e12d + self.c1*np.log(1+(phi/phi0))
    
    #  set_e21(self, phi):
    #    return self.e21d + self.c2*np.log(1+(phi/phi0))
    
    
    
    
    p4s = Parameters()
    
    ### 0. Optionally fit f(V) parameters with inwardRect data
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?
    
    
    ### phi0, gamma and E...
    
    ### 1. Load data for 'saturate' protocol to find Ipeak in order to calculate gmax : g
    p4s.add('g', value=gmax, vary=False) # derived['gmax'] = True
    # Change the model to be consistent so that g = gbar * A
    
    
    ### 2. Fit exponential to peak recovery plots : Gr
    p4s.add('Gr', value=Gr, vary=False)
    # Is this a valid assumption for the 4-state model?
    
    
    ### 3. Use Optimisation algorithms to fit the remaining parameters
    
    # Divide curve into on and off phases
    # Replace indices with time-value matching for real data
    # t[onInd] <= t_on < t[offInd+1]
    # t[offInd] <= t_off 
    Ion, Ioff = I[onInd:offInd+1], I[offInd:]
    ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    
    fitPseudoParams = True
    fitWholeCurve = False
    frac = 1    
    
    
    if not fitWholeCurve:
        ### 3a. Fit biexponential to off curve to find lambdas
        ### Fit off curve - if initial conditions can be calculated, it might be better to use an analytic solution relating directly to model parameters c.f. analysis_4state_off_new.m

        pOff = Parameters() # Create parameter dictionary
        I0 = abs(I[offInd])
        pOff.add('a0', value=I0, expr='{}-a2'.format(I0))#vary=False) #'{}-a1-a2'
        pOff.add('a1', value=-1) # Islow
        pOff.add('lam1', value=0.1, min=0)
        pOff.add('a2', value=-1) # Ifast
        pOff.add('lam2', value=10, min=0)
        offP = minimize(resid4off,pOff,args=(Ioff,toff),method=meth)
        print(">>> Off-phase fitting summary: <<<")
        if verbose > 0:
            print(fit_report(pOff))
        print('lambda1 = {:.5g}; lambda2 = {:.5g}'.format(pOff['lam1'].value, pOff['lam2'].value))
        print(offP.message)
        print("Chi^2: ",offP.chisqr)

        #p3s.add('Gd', value=pOff['Gd'].value, vary=False)
        
        ### Replace the off-curve fitting with the analytic expression and pass remaining parameters to interactive tools...
        
        ### 3b. Feed the constraints from Equations 29, 30 and 35 into fitting the on curve
        ### ==> Find Gd1 Gd2 e12d e21d c1 c2
        # runTrial(ChR, nPulses, V,phiOn,delD,onD,offD,padD)
        
        sumLam = pOff['lam1'].value + pOff['lam2'].value
        prodLam = pOff['lam1'].value * pOff['lam2'].value
        
        ### Try fitting the actual model parameters to the off curve before fitting the on curve
    
    RhO = selectModel(4) ### phi0, gam and E initialised to defaults
    RhO.setLight(0)
        
    pOn = Parameters()
    pOn.add('g', value=gmax, vary=False) # derived['gmax'] = True
    pOn.add('Gr', value=Gr, vary=False)
    
    if fitPseudoParams:   # Fit Pseudo parameters
        pOn.add('Ga1',value=1,min=0.01,max=20) # Find good starting values
        pOn.add('Ga2',value=1,min=0.01,max=10) # Find good starting values
        pOn.add('Gd1',value=0.12,min=0.01,vary=True)#,expr='sL-Gd2-e12d-e21d')
        pOn.add('Gd2',value=0.02, min=0.01,vary=True)
        pOn.add('e12',value=0.2, min=0.001)
        pOn.add('e21',value=0.1, min=0.001)
        
        RhO.setLight(phi) # Calculate transition rates for phi then override within resid4onPP
        #RhO.setParams(pOn)
        print('Optimising',end='')
        if not fitWholeCurve:
            onP = minimize(resid4onPP,pOn,args=(Ion[0:int((offInd-onInd)/frac)],ton[0:int((offInd-onInd)/frac)],RhO,V,phi),method=meth)
        else:
            #print(len(I[onInd:]),len(t[onInd:]),offInd-len(I[:onInd]))
            onP = minimize(resid4PP,pOn,args=(I[onInd:],t[onInd:],offInd-len(I[:onInd]),RhO,V,phi),method=meth)
    else: # Fit parameters


        
        pOn.add('sL',value=sumLam,vary=False)
        pOn.add('pL',value=prodLam,vary=False)
        #pOn.add('Ga1',value=0.111,min=0.05*phi) # Find good starting values
        #pOn.add('Ga2',value=0.111,min=0.05*phi) # Find good starting values
        pOn.add('k1',value=0.5,min=0) # Ga1 = k1 * phi
        pOn.add('k2',value=0.2,min=0) # Ga2 = k2 * phi
        pOn.add('Gd1',value=0.1,min=0.01)#,expr='sL-Gd2-e12d-e21d')
        pOn.add('Gd2',value=0.02, min=0.01)
        pOn.add('e12d',value=0.01, min=0)#, expr='(pL-(Gd1*Gd2)-(Gd1*e21d))/Gd2')
        pOn.add('e21d',value=0.01, min=0)
        pOn.add('c1',value=0.05, min=0)
        pOn.add('c2',value=0.01, min=0)
        # Place RhO,V,phi in onP?
        
        ### Trim down ton? Take 10% of data or one point every ms?
        ### Instead try looping over a coarse grid of parameter values and saving RMSE for each combination c.f. analysis_4state_on_new.m
        RhO.setLight(phi) # Calculate transition rates for phi
        print('Optimising',end='')
        onP = minimize(resid4on,pOn,args=(Ion[0::5],ton[0::5],RhO,V,phi),method=meth)
    
    print("\n\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOn))
    print(onP.message)
    print("Chi^2: ",onP.chisqr)
    
#     ci = lmfit.conf_interval(pOff)
#     lmfit.printfuncs.report_ci(ci)
    
    # (I,t,onInd,offInd,V,gmax,Gr,phi,Ipeak,Iss)
    p4s.add('Gd1',value=pOn['Gd1'].value,vary=False)
    p4s.add('Gd2',value=pOn['Gd2'].value,vary=False)
    if fitPseudoParams:
        #p4s.add('k1',value=pOn['Ga1'].value/phi,vary=False)
        #p4s.add('k2',value=pOn['Ga2'].value/phi,vary=False)
        p4s.add('Ga1',value=pOn['Ga1'].value,vary=False)
        p4s.add('Ga2',value=pOn['Ga2'].value,vary=False)
        p4s.add('e12',value=pOn['e12'].value,vary=False)
        p4s.add('e21',value=pOn['e21'].value,vary=False)
        
    else:
        p4s.add('k1',value=pOn['k1'].value,vary=False)
        p4s.add('k2',value=pOn['k2'].value,vary=False)
        p4s.add('e12d',value=pOn['e12d'].value,vary=False)
        p4s.add('e21d',value=pOn['e21d'].value,vary=False)
        p4s.add('c1',value=pOn['c1'].value,vary=False)
        p4s.add('c2',value=pOn['c2'].value,vary=False)
    

    
        # Generate Rhodopsin model
        RhO.g = p4s['g'].value              # 16700     [pS]
        RhO.Gr = p4s['Gr'].value            # 1/5000    [ms^-1]
        RhO.k1 = p4s['k1'].value            #           [ms^-1 * photons^-1 * s * mm^2]
        RhO.k2 = p4s['k2'].value            #           [ms^-1 * photons^-1 * s * mm^2]
        RhO.Gd1 = p4s['Gd1'].value          #           [ms^-1]
        RhO.Gd2 = p4s['Gd2'].value          #           [ms^-1]
        RhO.e12d = p4s['e12d'].value        #           [ms^-1]
        RhO.e21d = p4s['e21d'].value        #           [ms^-1]
        RhO.c1 = p4s['c1'].value            #           [ms^-1 * photons^-1 * s * mm^2]
        RhO.c2 = p4s['c2'].value            #           [ms^-1 * photons^-1 * s * mm^2]
        RhO.phiFit = phi                    # Flux intensity at which the parameters were fit. 
        RhO.setLight(0)                     # Re-initialise model to dark state

    
    ### Plot curves
    totT = max(t)
    Ifig = plt.figure()
    gsPL = plt.GridSpec(4,1)
    axFit = Ifig.add_subplot(gsPL[:-1,:])
    
    #ax = Ifig.add_subplot(111)
    plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    plt.xlim((0,totT))
    #plt.xlabel('$\mathrm{Time\ [ms]}$')
    plt.setp(axFit.get_xticklabels(), visible=False)
    plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    plt.plot(t,I,color='b',label='$\mathrm{Experimental\ Data}$')
    if addTitles:
        plt.title('Four-state model fit to data (phi={:.3g}): \n[Ga1={:.3g}; Ga2={:.3g}; e12={:.3g}; e21={:.3g}; Gd1={:.3g}; Gd2={:.3g}]'.format(phi,RhO.Ga1,RhO.Ga2,RhO.e12,RhO.e21,RhO.Gd1,RhO.Gd2))
    
    # IfitOn = calc4on(pOn, ton) # onP
    # IfitOff = calc4off(pOff, toff) # offP
#     tfit = np.append(ton,toff)
#     plot(tfit,Ifit,color='r')
    
    if fitPseudoParams:
        RhO.setLight(phi) # Calculate transition rates for phi then override within resid4onPP
        if not fitWholeCurve:
            IfitOn = calc4onPP(p4s,ton,RhO,V,phi) ##### p4s
            plt.plot(ton+t[onInd],IfitOn,color='g')
            #IfitOff = calc4offPP(p4s,toff,RhO,V)
            #plt.plot(toff+t[offInd],IfitOff,color='g',linestyle=':')
            IfitOff = calc4off(pOff,toff)
            plt.plot(toff+t[offInd],IfitOff,color='g',linestyle=':',label='$\mathrm{Four-state\ model\ fit}$')#,linewidth=3)
        else:
            Ifit = calc4PP(p4s,t[onInd:]-t[onInd],offInd-len(I[:onInd]),RhO,V,phi)
            plt.plot(t[onInd:],Ifit,color='g',label='$\mathrm{Four-state\ model\ fit}$')
    else:
        nPulses = 1
        dt=t[1]-t[0]
        delD = t[onInd]
        onD = t[offInd]-t[onInd] # ton
        offD = t[-1]-t[offInd]
        padD = 0
        Ifit, tfit, _ = runTrial(RhO,nPulses,V,phi,delD,onD,offD,padD,dt)
        plt.plot(tfit,Ifit,color='g',label='$\mathrm{Four-state\ model\ fit}$')
        # plt.plot(ton+t[onInd],IfitOn,color='g')
        # plt.plot(toff+t[offInd],IfitOff,color='g')
    
    plt.legend(loc='best')
    
    axRes = Ifig.add_subplot(gsPL[-1,:],sharex=axFit)
    
    # axLag.set_aspect('auto')
                    
    # Rfig = plt.figure() # Replace this with a subplot...
    #print(IfitOn)
    #print(Ion)
    #print((IfitOn-Ion)/Ion)
    plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    #plt.plot(toff+t[offInd],Ioff-IfitOff)
    plt.ylabel('$\mathrm{Residuals}$')# % relative error')
    plt.xlabel('$\mathrm{Time\ [ms]}$')
    #plt.setp(axRes.get_xticklabels(), visible=False)
    plt.xlim((0,totT))
    plt.tight_layout()
    
    print("Parameters have been fit for the four-state model")# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
    
    return RhO


    
def fit6states(I,t,onInd,offInd,gbar,phi): # ,Ipeak
    pass
    
    
#def runSSA(RhO):
     #for protocol in []: # List of protocols for characterisation
        # Run protocols...
    #from .protocols import *
    #smallSignal = { 'saturate': protSaturate, 'step': protStep, 'sinusoid': protSinusoid }
#    for key in smallSignalAnalysis: ### for key, value in smallSignal.iteritems():
#        P = smallSignalAnalysis[key]()
#        P.runProtocol(RhO)
#        P.plotProtocol()

        
def fitModels(dataSet, highestState=3):
    """Routine to fit as many models as possible and select between them according to some parsimony criterion"""
    
    ### Trim data slightly to remove artefacts from light on/off transition ramps?
    
    ### Could use precalculated lookup tables to find the values of steady state O1 & O2 occupancies?

    if isinstance(dataSet['custom'], PhotoCurrent): # Single photocurrent
        #targetPC = dataSet['custom']
        nRuns = 1
        nPhis = 1
        nVs = 1
        #setPC = [[[dataSet['custom']]]]
        setPC = ProtocolData('custom',nRuns,[dataSet['custom'].phi],[dataSet['custom'].V])
        setPC.trials[0][0][0] = dataSet['custom']
    elif isinstance(dataSet['custom'], ProtocolData): # Set of photocurrents
        setPC = dataSet['custom']
        nRuns = setPC.nRuns
        nPhis = setPC.nPhis
        nVs = setPC.nVs
        #if (setPC.nPhis * setPC.nVs) > 1:
    else:
        print(type(dataSet['custom']))
        print(dataSet['custom'])
        raise TypeError("dataSet['custom']")
    
    
    ### Extract the parameters relevant to all models - move inside loop for varyIPI protocol?
    if 'params' in dataSet:
        params = dataSet['params']
    else: 
        params = None
        
    # Now check for the following: E,[phi0,gam,A]; Gr0,gmax,Ipeak,Iss
    
    ### E
    if 'E' in params:
        E = params['E'].value
    else:
        #global E # Set from global parameters
        E = dataSet['E']
    print('E = {}'.format(E))
    
    ### Gr0
    if 'Gr' in params:
        Gr0 = params['Gr'].value
    elif 'Gr0' in params:
        Gr0 = params['Gr0'].value
    elif 'Gr_dark' in params:
        Gr0 = params['Gr_dark'].value
    elif 'a6' in params:
        Gr0 = params['a6'].value
    else: ### 2. Fit exponential to peak recovery plots
        if hasattr(dataSet['varyIPI'], 'tau_r'):
            Gr0 = 1/dataSet['varyIPI'].tau_r # Gr,dark
        else:
            print("Extract the peaks and fit an exponential...")
            if not (hasattr(dataSet['varyIPI'], 'tpIPI') and hasattr(dataSet['varyIPI'], 'IpIPI')):
                # Extract peaks
                dataSet['varyIPI'].IpIPI = np.zeros(dataSet['varyIPI'].nRuns)
                dataSet['varyIPI'].tpIPI = np.zeros(dataSet['varyIPI'].nRuns)
                for r in range(dataSet['varyIPI'].nRuns):
                    ### Search only within the on phase of the second pulse
                    I_RhO = dataSet['varyIPI'].Is[run][0][0] # phiInd=0 and vInd=0 Run for each phi and V?
                    startInd = dataSet['varyIPI'].PulseInds[run][0][0][1,0]
                    endInd = dataSet['varyIPI'].PulseInds[run][0][0][1,1]
                    extOrder = int(1+endInd-startInd) #100#int(round(len(I_RhO)/5))
                    #peakInds = findPeaks(I_RhO[:endInd+extOrder+1],minmax,startInd,extOrder)
                    peakInds = findPeaks(I_RhO[:endInd+extOrder+1],startInd,extOrder)
                    if len(peakInds) > 0: # Collect data at the (second) peak
                        dataSet['varyIPI'].IpIPI[run] = I_RhO[peakInds[0]] #-1 peaks
                        dataSet['varyIPI'].tpIPI[run] = t[peakInds[0]] # tPeaks
            # Fit exponential
            popt, _, _ = fitPeaks(dataSet['varyIPI'].tpIPI, dataSet['varyIPI'].IpIPI, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$')
            Gr0 = 1/popt[1]
            ### calcGr0()
        params.add('Gr0', value=Gr0, vary=False)
    print('Gr0 = {}'.format(Gr0))
    
    ### Ipmax        
    if 'saturate' in dataSet:
        if hasattr(dataSet['saturate'], 'Ipmax'): 
            Ipmax = dataSet['saturate'].Ipmax
        else: # Find maximum peak for saturate protocol
            # peakInd = findPeaks(I_phi,startInd=0,extOrder=5) 
            if (dataSet['saturate'].V < E): # Find Minima
                Ipmax = min(dataSet['saturate'].I)
            else:       # Find Maxima
                Ipmax = max(dataSet['saturate'].I)
            dataSet['saturate'].Ipmax = Ipmax
    else: #hasattr(dataSet['custom'], 'Ipeak'): # Use peak of sample photocurrent as an estimate
        Ipmax, inds = setPC.getIpmax()
        Vsat = setPC[inds[0]][inds[1]][inds[2]].V
        #Ipmax = dataSet['custom'].Ipeak
        #Vsat = dataSet['custom'].V
    print('Ipmax = {}'.format(Ipmax))
    
    ### g        
    if 'g' in params: ###Ipeak
        gmax = params['g'].value        
    elif 'saturate' in dataSet:
        if hasattr(dataSet['saturate'], 'gbar_est'):
            gmax = dataSet['saturate'].gbar_est
        Vsat = dataSet['saturate'].V
    else: ### 1. Load data for 'saturate' protocol to find Ipmax in order to calculate gmax
        ### gmax = Ipmax/([O_p]*(V-E)) = Ipmax/(V-E) # with the assumption [O_p] = 1
        ### This assumption is an underestimate for 4 & 6 state models: [O_p] =~ 0.71 (depending on rates)
        assert(Vsat != E) #if dataSet['saturate'].V != E:
        gmax = Ipmax/(Vsat-E) # Assuming [O_p] = 1
        dataSet['saturate'].gbar_est = gmax
        ### calcG()
    print('g = {}'.format(gmax))
        
    # Change the model to be consistent so that g = gbar * A
    

    ##### FINISH THIS!!! #####
    # if hasattr(dataSet['custom'], 'Iss'):
        # Iss = dataSet['custom'].Iss
    # else: 

    ### Optionally fit f(V) parameters with inwardRect data - MUST MEASURE E AND FIT AFTER OTHER PARAMETERS
    if 'v0' in params and 'v1' in params:
        v0 = params['v0'].value
        v1 = params['v1'].value
    else:
        if 'inwardRect' in dataSet:
            if hasattr(dataSet['inwardRect'], 'Iss'): # Use extracted values
                Iss = dataSet['inwardRect'].Iss
                Vs = dataSet['inwardRect'].Vs
            else: # Extract steady state values
                print("Finish f(V) fitting!")
                Iss = None
        elif setPC.nVs > 1:
            IssSet, VsSet = setPC.getIRdata()
            for phiInd, phiOn in enumerate(phis): 
                ### PLOT
                RhO.calcSteadyState(phiOn)
                popt, pcov, eqString = fitfV(Vs,self.IssVals[run][phiInd][:],calcIssfromfV,p0fV)#,eqString)
                
                # Add equations to legend
                if len(phis) > 1: 
                    legLabels[phiInd] = eqString + '$,\ \phi={:.3g}$'.format(phiOn)
                else:
                    legLabels[phiInd] = eqString
                
                ### Move this to fitting routines?
                # v0 = popt[0], v1 = popt[1], E = popt[2]
            # Fit Curve of V vs Iss
        ###else: # Assume f(V) = (V-E)
    
        
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?
    
    
    #Models = {'3':[[None for v in len(Vs)] for p in len(phis)]}
    
    ### Loop over phi and extrapolate for parameters - skip nRuns
    # for phiInd in range(nPhis):
        # for vInd in range(nVs):
            # targetPC = setPC.trials[0][phiInd][vInd] # Take the first run only
            #< Start fitting here...
    targetPC = setPC.trials[0][0][0]
    I = targetPC.I
    t = targetPC.t
    onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
    offInd = targetPC.pulseInds[0,1]
    V = targetPC.V
    phi = targetPC.phi
    Iss = targetPC.Iss # Iplat ############################# Change name to avoid clash with inwardRect!    
    ###Iplat is only required for the 3-state fitting procedure
            
    if highestState == 3:
        RhO3 = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,Iss)
        ###RhO3.E
        #RhO3.k
        #RhO3.Gd
        ###RhO3.Gr_dark
        ###RhO3.Gr_light
        #Models['3'][phiInd][vInd] = RhO3
        RhO = RhO3
        
        
    elif highestState == 4:
        #Models['4'] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,Iss)
        RhO4 = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax)#,Iss)
        RhO = RhO4
    elif highestState == 6:
        Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret)#,Iss)#...
    else:
        raise Exception('Invalid choice for highestState: {}!'.format(highestState))
            
            

    
    # Compare fit vs computational complexity...
    # RhO = selectModel(nStates)
    # Calculate chisqr for each model and select between them. 
    ###RhO = RhO3
    
    # Run small signal analysis
    #runSSA(RhO)
    #characterise(RhO)
    
    return RhO #Models # # [RhO3,RhO4,RhO6]