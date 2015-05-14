# Call this script in IPython
# %run fitStates.py
# import (as a module)
# %import
# %load

# import lmfit
from lmfit import minimize, Parameters, fit_report
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from .parameters import *
from .loadData import * #import loadData
from .protocols import * # for SSA and numerical solution
#import models # for calculating analytic solution
#global phi0
from .models import * # for fitPeaks
from .config import verbose, saveFigFormat, addTitles, fDir, dDir
#from .config import verbose, saveFigFormat, eqSize, addTitles, addStimulus, colours, styles, dDir, fDir
import os
import pickle
import time

#import protocols

methods=('leastsq','nelder','lbfgsb','powell','cg','newton','cobyla','tnc','trust-ncg','dogleg','slsqp')
#method = methods[3] #3
# Notes
# =====
# Empirically, 'powell', 'lbfgsb' and 'nelder' typically provide the best results
# 'newton', 'trust-ncg' and 'dogleg' require a Jacobian 
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
# Gr0 = 1/dataset['recovery'].tau_r # Gr,dark
# I
# RhO = fit3states(I,t,onInd,offInd,dataset,gmax,Ipeak)

def copyParam(name, source, target):
    if name not in target:
        target.add(name, value=source[name].value, vary=source[name].vary, min=source[name].min, max=source[name].max, expr=source[name].expr)
    else:
        target[name].set(value=source[name].value, vary=source[name].vary, min=source[name].min, max=source[name].max, expr=source[name].expr)
    return #target


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
    if hasattr(dataset['recovery'], 'tau_r'):
        Gr0 = 1/dataset['recovery'].tau_r # Gr,dark
    else:
        print("Extract the peaks and fit an exponential...")
        
    return Gr0
    
def fit3statesIndiv(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=None,method=methods[3]): #Iss=None, ### Modify to pass in p3s
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
    # self.P = self.set_P(phi) # [ms^-1] Quantum efficiency * number of photons absorbed by a ChR2 molecule per unit time
    # self.g = 100e-3 * 1.67e5   # [pS] g1 (100fS): Conductance of an individual ion channel * N (100000)
    ## self.c = 1/165             # [ms^-1] (Gr,light @ 1mW mm^-2)
    # self.e = 0.5*1.2e-8*1e-6/1.1 # Coefficient for the activation rate

    #  set_P(self, phi):
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
    # if hasattr(dataset['recovery'], 'tau_r'):
        # Gr0 = 1/dataset['recovery'].tau_r # Gr,dark
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
    #offPmin = minimize(resid3off,pOff,args=(Ioff,toff),method=meth)
    
    pOff.add('A', value=-1) # Islow
    pOff.add('Gd1', value=0)#, min=0))
    pOff.add('B', value=-1) # Ifast
    pOff.add('Gd2', value=100)#, min=0))
    pOff.add('C', value=0)
    #offPmin = minimize(resid3off,pOff,args=(Ioff[0:100],toff[0:100]),method=meth)
    offPmin = minimize(resid3off2exp,pOff,args=(Ioff,toff),method=method)
    #pOff.add('Gd', value=max(pOff['Gd1'].value,pOff['Gd2'].value), min=0)
    nf = pOff['A'].value + pOff['B'].value
    pOff.add('Gd', value=pOff['A'].value*pOff['Gd1'].value/nf+pOff['B'].value*pOff['Gd2'].value/nf, min=0)
    
    print(">>> Off-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOff)) # offPmin?
    print(offPmin.message)
    print("Chi^2 (Off-phase): ",offPmin.chisqr)
    p3s.add('Gd', value=pOff['Gd'].value, vary=False)
    
    #print(pOff['Gd'].value)
    #print(p3s['Gd'].value)
    
    ### 3b. Optimise to find P and Gr at a particular phi from the on curve
    ### Fit on curve
    pOn = Parameters()
    pOn.add('a0', value=-1, expr='-a2')
    pOn.add('a1', value=1)
    pOn.add('a2', value=1)
    pOn.add('tau_act', value=2, min=0)
    pOn.add('tau_deact', value=10, min=0)
    onPmin = minimize(resid3on,pOn,args=(Ion,ton),method=method)
    
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
    print("Solution 1: P = {}; Gd = {}; Gr = {}".format(P1,Gd,Gr1))
    
    x2 = (alpha - np.sqrt(alpha**2 - 4*beta))/2
    Gr2 = x2 * Gd
    P2 = lam1+lam2 - (Gr2+Gd)
    print("Solution 2: P = {}; Gd = {}; Gr = {}".format(P2,Gd,Gr2))
    
    if P1 > Gr1:
        print("==> Solution 1 (assuming that P > Gr): {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(P1,Gd,Gr1))
        P = P1
        Gr = Gr1
    else:
        print("==> Solution 2 (assuming that P > Gr): {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(P2,Gd,Gr2))
        P = P2
        Gr = Gr2
    #print("==> Original values: {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(  ,Gd,ChR.Gr))
    
    if (Gr + Gd - 2*np.sqrt(Gr*Gd)) < P < (Gr + Gd + 2*np.sqrt(Gr*Gd)):
        print('\n\nWarning! No real solution exists!\n\n')
    
    print("\n\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOn)) # onPmin?
    print(onPmin.message)
    print("Chi^2 (On-phase): ",onPmin.chisqr)
    
    k = P/phi # To find the coefficient which scales phi to calculate the activation rate
    p3s.add('k', value=k, vary=False)
    # N.B. The three-state model does not generalise well to other values of phi
    
    GrPhi = Gr - Gr0 # To find Grlight for a particular phi
    p3s.add('GrPhi', value=GrPhi, vary=False)

#     ci = lmfit.conf_interval(pOff)
#     lmfit.printfuncs.report_ci(ci)
    

    ### 4. Optionally fit f(V) parameters with rectifier data
    # if 'rectifier' in dataset:
        # if hasattr(dataset['rectifier'], 'Iss'): # Use extracted values
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
        plt.title('Three-state model fit to data (phi={:.3g}): [P={:.3g}; Gd={:.3g}; Gr={:.3g}]'.format(phi,RhO.P,RhO.Gd,RhO.Gr))
    RhO.setLight(0)
    
    
    if 1: # Use pseudo parameters
        IfitOn = calc3on(pOn, ton) # onPmin
        pOff['A'].value = -IfitOn[-1] # Constrain off curve to start at the end of the fitted on curve
        IfitOff = calc3off(pOff, toff) # offPmin
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
        IfitOn = RhO.calcI(V, onStates)
        print(onStates)
        RhO.setLight(0)
        #print(toff)
        IfitOff = RhO.calcI(V, RhO.calcSoln(toff, onStates[-1,:]))
    
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
    
    return p3s # RhO


    
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
    
    I_RhO = RhO.calcI(V, RhO.states)
    
    return I_RhO
    
    
def calc4PP(p,t,offInd,RhO,V,phi):
    if verbose > 1:
        print('.', end="") # sys.stdout.write('.')
    
    #print(t.shape, offInd)
    RhO.initStates(phi)
    #RhO.setLight(phi)
    # Override light-sensitive transition rates
    v = p.valuesdict()
    RhO.Ga1 = v['Ga1']#p['Ga1'].value
    RhO.Ga2 = v['Ga2']#p['Ga2'].value
    RhO.e12 = v['e12']#p['e12'].value
    RhO.e21 = v['e21']#p['e21'].value
    
    RhO.Gd1 = v['Gd1']#p['Gd1'].value
    RhO.Gd2 = v['Gd2']#p['Gd2'].value
    soln = odeint(RhO.solveStates, RhO.s_0, t[:offInd+1], Dfun=RhO.jacobian)
    RhO.storeStates(soln,t)
    #Ion = RhO.calcI(V, soln)
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
    #Ioff = RhO.calcI(V, soln)
    #print(soln.shape)
    
    
    #I_RhO = np.concatenate((Ion, Ioff[1:]))
    I_RhO = RhO.calcI(V, RhO.states)
    
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
    RhO.updateParams(p) #RhO.setParams(p)
    # Override light-sensitive transition rates
    RhO.Ga1 = p['Ga1'].value
    RhO.Ga2 = p['Ga2'].value
    RhO.e12 = p['e12'].value
    RhO.e21 = p['e21'].value
    
    RhO.Gd1 = p['Gd1'].value
    RhO.Gd2 = p['Gd2'].value
    soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
    #RhO.storeStates(soln, t) # Assumes no previous storage (i.e. not [-1,:]
    I_RhO = RhO.calcI(V, soln)
    return I_RhO
    
def resid4onPP(p,I,t,RhO,V,phi):
    return I - calc4onPP(p,t,RhO,V,phi)
    
    

def calc4on(p,t,RhO,V,phi):
    """Simulate the on-phase from base parameters for the 4-state model"""
    if verbose > 1:
        print('.', end="") # sys.stdout.write('.')
    RhO.initStates(0)
    #RhO.setLight(phi) # Calculate transition rates for phi
    
    #RhO.setParams(p)
    # RhO.k1 = p['k1'].value # Ga1 = k1*phi
    # RhO.k2 = p['k2'].value # Ga2 = k2*phi
    # RhO.Gd1 = p['Gd1'].value
    # RhO.Gd2 = p['Gd2'].value
    # RhO.e12d = p['e12d'].value
    # RhO.e21d = p['e21d'].value
    # RhO.c1 = p['c1'].value
    # RhO.c2 = p['c2'].value
    
    RhO.updateParams(p)
    
    # if 1:
        # # New activation function
        # RhO.p = p['p'].value
        # RhO.phim = p['phim'].value
    
    # nPulses = 1
    # delD = 0
    # onD = max(t) # ton
    # offD = 0
    # padD = 0
    # I_RhO, _, _ = runTrial(RhO,nPulses,V,phi,delD,onD,offD,padD)
    
    RhO.setLight(phi) # Calculate transition rates for phi
    
    soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
    # soln,out = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian, full_output=True)
    # if out['message'] != 'Integration successful.':
        # #print(out)
        # print(RhO.reportParams())
    I_RhO = RhO.calcI(V, soln)
    return I_RhO
    
def resid4on(p,I,t,RhO,V,phi):
    return I - calc4on(p,t,RhO,V,phi)
    
    
def fit4statesIndiv(I,t,onInd,offInd,phi,V,Gr,gmax,params=None,method=methods[3]):#,Iss): # ,Ipeak
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
    #    return e*phi*sigma_ret / w_loss
    
    #  set_Ga2(self, phi):
    #    return e*phi*sigma_ret / w_loss
    
    #  set_e12(self, phi):
    #    return self.e12d + self.c1*np.log(1+(phi/phi0))
    
    #  set_e21(self, phi):
    #    return self.e21d + self.c2*np.log(1+(phi/phi0))
    
    p4s = Parameters()    
    ### 0. Optionally fit f(V) parameters with rectifier data
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?

    ### phi0, gamma and E...
    
    ### 1. Load data for 'saturate' protocol to find Ipeak in order to calculate gmax : g # Make model consistent so that g = gbar * A
    p4s.add('g', value=gmax, vary=False) # derived['gmax'] = True
            
    ### 2. Fit exponential to peak recovery plots : Gr # Is this a valid assumption for the 4-state model?
    p4s.add('Gr', value=Gr, vary=False)

    ### 3. Use Optimisation algorithms to fit the remaining parameters
    
    # Divide curve into on and off phases
    # Replace indices with time-value matching for real data
    # t[onInd] <= t_on < t[offInd+1]
    # t[offInd] <= t_off 
    Ion, Ioff = I[onInd:offInd+1], I[offInd:]
    ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    
    fitPseudoParams = True  # Fit by light dependent parameters
    fitByPhase = True   # Fit on and off curves separately
    frac = 1    
    
    
    if fitByPhase:
        ### 3a. Fit biexponential to off curve to find lambdas
        ### Fit off curve - if initial conditions can be calculated, it might be better to use an analytic solution relating directly to model parameters c.f. analysis_4state_off_new.m

        pOff = Parameters() # Create parameter dictionary
        I0 = abs(I[offInd])
        pOff.add('a0', value=I0, expr='{}-a2'.format(I0))#vary=False) #'{}-a1-a2'
        pOff.add('a1', value=-1) # Islow
        pOff.add('lam1', value=0.1, min=0)
        pOff.add('a2', value=-1) # Ifast
        pOff.add('lam2', value=10, min=0)
        offPmin = minimize(resid4off,pOff,args=(Ioff,toff),method=method)
        print(">>> Off-phase fitting summary: <<<")
        if verbose > 1:
            print(fit_report(pOff))
            print("Error bars: ",onPmin.errorbars)
        print('lambda1 = {:.5g}; lambda2 = {:.5g}'.format(pOff['lam1'].value, pOff['lam2'].value))
        print(offPmin.message)
        print("Chi^2: ",offPmin.chisqr)

        #p3s.add('Gd', value=pOff['Gd'].value, vary=False)
        
        ### Replace the off-curve fitting with the analytic expression and pass remaining parameters to interactive tools...
        
        ### 3b. Feed the constraints from Equations 29, 30 and 35 into fitting the on curve
        ### ==> Find Gd1 Gd2 e12d e21d c1 c2
        
        sumLam = pOff['lam1'].value + pOff['lam2'].value
        prodLam = pOff['lam1'].value * pOff['lam2'].value

        ### Try fitting the actual model parameters to the off curve before fitting the on curve
    
    RhO = selectModel(4) ### phi0, gam and E initialised to defaults
    RhO.setLight(0)
    
    if params is not None:
        pOn = params #p4s
    else:
        pOn = Parameters()
    if 'g' not in pOn:
        pOn.add('g', value=gmax, vary=False) # derived['gmax'] = True
    if 'Gr' not in pOn:
        pOn.add('Gr', value=Gr, vary=False)
    
    if fitPseudoParams:   # Fit Pseudo parameters
        if 'Ga1' not in pOn: #params is None:
            pOn.add('Ga1',value=1,min=0.01,max=5) # Find good starting values   #20
        if 'Ga2' not in pOn: 
            pOn.add('Ga2',value=1,min=0.01,max=5) # Find good starting values   #10
        if 'delta' not in pOn:
            pOn.add('delta',value=sumLam/2,min=0.1,max=sumLam,vary=True)
        if 'Gd1' not in pOn:     
            pOn.add('Gd1',value=0.12,min=0.01,max=1,vary=True,expr='delta-Gd2')#,expr='sL-Gd2-e12d-e21d')
        if 'Gd2' not in pOn:     
            pOn.add('Gd2',value=0.02, min=0.01,max=1,vary=True)
        if 'e12' not in pOn: 
            pOn.add('e12',value=0.2, min=0.001, max=0.1)
        if 'e21' not in pOn: 
            pOn.add('e21',value=0.1, min=0.001, max=0.1)
        
        
        RhO.setLight(phi) # Calculate transition rates for phi then override within resid4onPP
        #RhO.setParams(pOn)
        if verbose > 1:
            print('Optimising',end='')
        if fitByPhase:
            onPmin = minimize(resid4onPP,pOn,args=(Ion[0:int((offInd-onInd)/frac)],ton[0:int((offInd-onInd)/frac)],RhO,V,phi),method=method)
        else:
            #print(len(I[onInd:]),len(t[onInd:]),offInd-len(I[:onInd]))
            onPmin = minimize(resid4PP,pOn,args=(I[onInd:],t[onInd:],offInd-len(I[:onInd]),RhO,V,phi),method=method)
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
        # Place RhO,V,phi in onPmin?
        
        ### Trim down ton? Take 10% of data or one point every ms?
        ### Instead try looping over a coarse grid of parameter values and saving RMSE for each combination c.f. analysis_4state_on_new.m
        RhO.setLight(phi) # Calculate transition rates for phi
        if verbose > 1:
            print('Optimising',end='')
        onPmin = minimize(resid4on,pOn,args=(Ion[0::5],ton[0::5],RhO,V,phi),method=method)
    
    print("\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOn))
        print("Error bars: ",onPmin.errorbars)
    if not onPmin.success:
        print(onPmin.success)
        print(onPmin.message)
        print(onPmin.lmdif_message)
    print("Chi^2: ",onPmin.chisqr)
    
    
#     ci = lmfit.conf_interval(pOff)
#     lmfit.printfuncs.report_ci(ci)
    
    # (I,t,onInd,offInd,V,gmax,Gr,phi,Ipeak,Iss)
    #if params is None:
    if 'Gd1' not in p4s:
        copyParam('Gd1',pOn,p4s) #p4s.add('Gd1',value=pOn['Gd1'].value,vary=False)
    if 'Gd2' not in p4s:
        copyParam('Gd2',pOn,p4s) #p4s.add('Gd2',value=pOn['Gd2'].value,vary=False)
        
        
    
    
    if fitPseudoParams:
        #p4s.add('k1',value=pOn['Ga1'].value/phi,vary=False)
        #p4s.add('k2',value=pOn['Ga2'].value/phi,vary=False)
        #if params is None:
        copyParam('Ga1',pOn,p4s) #p4s.add('Ga1',value=pOn['Ga1'].value,vary=False)
        copyParam('Ga2',pOn,p4s) #p4s.add('Ga2',value=pOn['Ga2'].value,vary=False)
        copyParam('e12',pOn,p4s) #p4s.add('e12',value=pOn['e12'].value,vary=False)
        copyParam('e21',pOn,p4s) #p4s.add('e21',value=pOn['e21'].value,vary=False)
        
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
    
    # IfitOn = calc4on(pOn, ton) # onPmin
    # IfitOff = calc4off(pOff, toff) # offPmin
#     tfit = np.append(ton,toff)
#     plot(tfit,Ifit,color='r')
    
    if fitPseudoParams:
        RhO.setLight(phi) # Calculate transition rates for phi then override within resid4onPP
        if fitByPhase:
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
    plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff))
    #plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    #plt.plot(toff+t[offInd],Ioff-IfitOff)
    plt.ylabel('$\mathrm{Residuals}$')# % relative error')
    plt.xlabel('$\mathrm{Time\ [ms]}$')
    plt.axhline(y=0,linestyle=':',color='k')
    #plt.setp(axRes.get_xticklabels(), visible=False)
    plt.xlim((0,totT))
    plt.tight_layout()
    
    print("Parameters have been fit for the four-state model")# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
    
    return p4s #RhO


    
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


        
def fitCurve(dataSet, nStates=3, params=None):
    """Fit a single photocurrent to flux-dependent transition rates"""
    # E.g.  P, Gr [Gd]
    #       Ga1, Ga2, e12, e21, [Gd1, Gd2, Gr]
    #       a1, a3, b2, b4, [a2, a4, a6, b1, b3]
    # Metaparameters are then fit to the trends of these in fitModels()
    # E.g.  k
    
    
    ### Check contents of dataSet and produce report on model features which may be fit. 
    # e.g. if not 'rectifier': f(V)=1
    
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
    
    
    ### Extract the parameters relevant to all models - move inside loop for recovery protocol?
    if params is not None: # Allow you to pass paramaters outside of dataSet
        params = params
    elif 'params' in dataSet:
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
        if hasattr(dataSet['recovery'], 'tau_r'):
            Gr0 = 1/dataSet['recovery'].tau_r # Gr,dark
        else:
            print("Extract the peaks and fit an exponential...")
            if not (hasattr(dataSet['recovery'], 'tpIPI') and hasattr(dataSet['recovery'], 'IpIPI')):
                # Extract peaks
                dataSet['recovery'].IpIPI = np.zeros(dataSet['recovery'].nRuns)
                dataSet['recovery'].tpIPI = np.zeros(dataSet['recovery'].nRuns)
                for r in range(dataSet['recovery'].nRuns):
                    ### Search only within the on phase of the second pulse
                    I_RhO = dataSet['recovery'].Is[run][0][0] # phiInd=0 and vInd=0 Run for each phi and V?
                    startInd = dataSet['recovery'].PulseInds[run][0][0][1,0]
                    endInd = dataSet['recovery'].PulseInds[run][0][0][1,1]
                    extOrder = int(1+endInd-startInd) #100#int(round(len(I_RhO)/5))
                    #peakInds = findPeaks(I_RhO[:endInd+extOrder+1],minmax,startInd,extOrder)
                    peakInds = findPeaks(I_RhO[:endInd+extOrder+1],startInd,extOrder)
                    if len(peakInds) > 0: # Collect data at the (second) peak
                        dataSet['recovery'].IpIPI[run] = I_RhO[peakInds[0]] #-1 peaks
                        dataSet['recovery'].tpIPI[run] = t[peakInds[0]] # tPeaks
            # Fit exponential
            popt, _, _ = fitPeaks(dataSet['recovery'].tpIPI, dataSet['recovery'].IpIPI, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$')
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

    ### Optionally fit f(V) parameters with rectifier data - MUST MEASURE E AND FIT AFTER OTHER PARAMETERS
    if 'v0' in params and 'v1' in params:
        v0 = params['v0'].value
        v1 = params['v1'].value
    else:
        if 'rectifier' in dataSet:
            if hasattr(dataSet['rectifier'], 'Iss'): # Use extracted values
                Iss = dataSet['rectifier'].Iss
                Vs = dataSet['rectifier'].Vs
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
    Iss = targetPC.Iss # Iplat ############################# Change name to avoid clash with rectifier!    
    ###Iplat is only required for the 3-state fitting procedure
            
    if nStates == 3:
        fitParams = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,Iss)
        ###RhO3.E
        #RhO3.k
        #RhO3.Gd
        ###RhO3.Gr_dark
        ###RhO3.Gr_light
        #Models['3'][phiInd][vInd] = RhO3
        #RhO = RhO3
        
        
    elif nStates == 4:
        #Models['4'] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,Iss)
        #RhO4 = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax)#,Iss)
        fitParams = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax)
        #RhO = RhO4

    elif nStates == 6:
        Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret)#,Iss)#...
    else:
        raise Exception('Invalid choice for nStates: {}!'.format(nStates))
            
            
    RhO = models[nStates]()
    RhO.setParams(fitParams)
    
    # Compare fit vs computational complexity...
    # RhO = selectModel(nStates)
    # Calculate chisqr for each model and select between them. 
    ###RhO = RhO3
    
    # Run small signal analysis
    #runSSA(RhO)
    #characterise(RhO)
    
    
    return fitParams #RhO #Models # # [RhO3,RhO4,RhO6]

    
def aggregateFits(phis,phiFits,nStates=3):
    nTrials = len(phis)
    
    if nStates==3:
        if verbose > 0: # Report optimised parameters from each run
            print('phi     |P      |Gd     |Gr      \n=================================')
        for trial in range(nTrials):
            v = phiFits[trial].valuesdict()
            #Models[trial].setLight(phis[trial])
            if verbose > 0: # Report optimised parameters from each run
                print("{:.2g}\t|{:.3g}\t|{:.3g}\t|{:.3g}".format(phis[trial],v['P'],v['Gd'],v['Gr']))
            Ps[trial] = v['P']
            Gds[trial] = v['Gd']
            Grs[trial] = v['Gr']
            
        print('Mean (Median): \tGd = {:.3g} ({:.3g})'.format(np.mean(Gds),np.median(Gds)))
        aggregates = {'P':Ps, 'Gd':Gds, 'Gr':Grs}
    
    elif nStates==4:
        Ga1s = [None for trial in range(nTrials)]
        Ga2s = [None for trial in range(nTrials)]
        e12s = [None for trial in range(nTrials)]
        e21s = [None for trial in range(nTrials)]
        Gd1s = [None for trial in range(nTrials)]
        Gd2s = [None for trial in range(nTrials)]
        
        if verbose > 0: # Report optimised parameters from each run
            print('phi     |Ga1    |Ga2    |e12    |e21    |Gd1    |Gd2     ')
            print('=========================================================')
        for trial in range(nTrials):
            #Models[trial].setLight(phis[trial])
            v = phiFits[trial].valuesdict()
            if verbose > 0: # Report optimised parameters from each run
                print("{:.2g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}".format(phis[trial],v['Ga1'],v['Ga2'],v['e12'],v['e21'],v['Gd1'],v['Gd2']))
            Ga1s[trial] = v['Ga1']
            Ga2s[trial] = v['Ga2']
            e12s[trial] = v['e12']
            e21s[trial] = v['e21']
            Gd1s[trial] = v['Gd1']
            Gd2s[trial] = v['Gd2']
        
        print('Means (Medians): \tGd1 = {:.3g} ({:.3g}); Gd2 = {:.3g} ({:.3g})'.format(np.mean(Gd1s),np.median(Gd1s), np.mean(Gd2s),np.median(Gd2s)))
        #print('Medians: \tGd1 = {:.3g}; Gd2 = {:.3g}'.format(np.median(Gd1s),np.median(Gd2s)))
        aggregates = {'Ga1':Ga1s, 'Ga2':Ga2s, 'e12':e12s, 'e21':e21s, 'Gd1':Gd1s, 'Gd2':Gd2s}
        
    elif nStates==6:
        pass
        
    return aggregates
    

# def aggregate4sFits(phis, phiFits):
    # nTrials = len(phis)
    # Ga1s = [None for trial in range(nTrials)]
    # Ga2s = [None for trial in range(nTrials)]
    # e12s = [None for trial in range(nTrials)]
    # e21s = [None for trial in range(nTrials)]
    # Gd1s = [None for trial in range(nTrials)]
    # Gd2s = [None for trial in range(nTrials)]
    
    # if verbose > 1: # Report optimised parameters from each run
        # print('phi     |Ga1    |Ga2    |e12    |e21    |Gd1    |Gd2     ')
        # print('=========================================================')
    # for trial in range(nTrials):
        # #Models[trial].setLight(phis[trial])
        # v = phiFits[trial].valuesdict()
        # if verbose > 1: # Report optimised parameters from each run
            # print("{:.2g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}\t|{:.3g}".format(phis[trial],v['Ga1'],v['Ga2'],v['e12'],v['e21'],v['Gd1'],v['Gd2']))
        # Ga1s[trial] = v['Ga1']
        # Ga2s[trial] = v['Ga2']
        # e12s[trial] = v['e12']
        # e21s[trial] = v['e21']
        # Gd1s[trial] = v['Gd1']
        # Gd2s[trial] = v['Gd2']
    
    # print('Means (Medians): \tGd1 = {:.3g} ({:.3g}); Gd2 = {:.3g} ({:.3g})'.format(np.mean(Gd1s),np.median(Gd1s), np.mean(Gd2s),np.median(Gd2s)))
    # #print('Medians: \tGd1 = {:.3g}; Gd2 = {:.3g}'.format(np.median(Gd1s),np.median(Gd2s)))
    # aggregates = {'Ga1':Ga1s, 'Ga2':Ga2s, 'e12':e12s, 'e21':e21s, 'Gd1':Gd1s, 'Gd2':Gd2s}
    # return aggregates

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
##### New routines #####
    
def plotFit(I,t,onInd,offInd,phi,V,nStates,params,fitRates,index):
    
    RhO = models[str(nStates)]()
    #RhO.initStates(0.0)
    RhO.updateParams(params)
    RhO.phiFit = phi                    # Flux intensity at which the parameters were fit. 
    #RhO.setLight(0)                     # Re-initialise model to dark state
    
    ### Plot experimental curve
    totT = max(t)
    Ifig = plt.figure()
    gsPL = plt.GridSpec(4,1)
    axFit = Ifig.add_subplot(gsPL[:-1,:])
    plotLight(np.asarray([[t[onInd],t[offInd]]]), axFit) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    plt.xlim((0,totT))
    #plt.xlabel('$\mathrm{Time\ [ms]}$')
    plt.setp(axFit.get_xticklabels(), visible=False)
    plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    plt.plot(t,I,color='g',label='$\mathrm{Experimental\ Data}$')
    
    
    # IfitOn = calc4on(pOn, ton) # onPmin
    # IfitOff = calc4off(pOff, toff) # offPmin
#     tfit = np.append(ton,toff)
#     plot(tfit,Ifit,color='r')
    
    
    Idel, Ion, Ioff = I[:onInd+1], I[onInd:offInd+1], I[offInd:]
    tdel, ton, toff = t[:onInd+1], t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    
    # Delay phase
    RhO.setLight(RhO.phi_0)
    if RhO.useAnalyticSoln:
        soln = RhO.calcSoln(tdel, RhO.s_0)
    else:
        soln = odeint(RhO.solveStates, RhO.s_0, tdel, Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:], tdel[1:])
    
    
    
    # On phase
    RhO.setLight(phi) # Calculate transition rates for phi then override within resid4onPP
    if fitRates: # Override light-sensitive transition rates
        RhO.updateParams(params)
    
    RhO.s_on = soln[-1,:]
    if RhO.useAnalyticSoln:
        soln = RhO.calcSoln(ton, RhO.s_on)
    else:
        soln = odeint(RhO.solveStates, RhO.s_on, ton, Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:], ton[1:])
    #IfitOn = RhO.calcI(V, soln)
    
    if addTitles:
        if nStates == 3:
            plt.title('Three-state model fit to data (phi={:.3g}): [P={:.3g}; Gd={:.3g}; Gr={:.3g}] \n[k={:.3g}; p={:.3g}; phim={:.3g}; Gd={:.3g}; Gr0={:.3g}; Gr1={:.3g}]'.format(phi,RhO.P,RhO.Gd,RhO.Gr,RhO.k,RhO.p,RhO.phim,RhO.Gd,RhO.Gr0,RhO.Gr1))
        elif nStates == 4:
            plt.title('Four-state model fit to data (phi={:.3g}): \n[Ga1={:.3g}; Ga2={:.3g}; e12={:.3g}; e21={:.3g}; Gd1={:.3g}; Gd2={:.3g}]'.format(phi,RhO.Ga1,RhO.Ga2,RhO.e12,RhO.e21,RhO.Gd1,RhO.Gd2))
        elif nStates == 6:
            pass

    #IfitOn = calc4onPP(p4s,ton,RhO,V,phi) ##### p4s
    #plt.plot(ton+t[onInd],IfitOn,color='b')
    #IfitOff = calc4offPP(p4s,toff,RhO,V)
    #plt.plot(toff+t[offInd],IfitOff,color='g',linestyle=':')
    #IfitOff = calc4off(pOff,toff)
    #plt.plot(toff+t[offInd],IfitOff,color='b',linestyle=':',label='$\mathrm{Four-state\ model\ fit}$')#,linewidth=3)
    
    # Off phase
    RhO.setLight(0)
    #if fitRates: # Override light-sensitive transition rates
    #    RhO.updateParams(params)
        
    RhO.s_off = soln[-1,:]
    if RhO.useAnalyticSoln:
        soln = RhO.calcSoln(toff, RhO.s_off)
    else:
        soln = odeint(RhO.solveStates, RhO.s_off, toff, Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:],toff[1:])
    
    Ifit = RhO.calcI(V, RhO.states)
    
    # Plot model fit curve
    plt.plot(t,Ifit,color='b',label='$\mathrm{{Model\ fit\ ({}-states)}}$'.format(nStates)) #t[onInd:]


    # else: # fitParams
        # nPulses = 1
        # dt=t[1]-t[0]
        # delD = t[onInd]
        # onD = t[offInd]-t[onInd] # ton
        # offD = t[-1]-t[offInd]
        # padD = 0
        # Ifit, tfit, _ = runTrial(RhO,nPulses,V,phi,delD,onD,offD,padD,dt)
        # plt.plot(tfit,Ifit,color='b',label='$\mathrm{Four-state\ model\ fit}$')
        # # plt.plot(ton+t[onInd],IfitOn,color='g')
        # # plt.plot(toff+t[offInd],IfitOff,color='g')
    
    plt.legend(loc='best')
    
    axRes = Ifig.add_subplot(gsPL[-1,:],sharex=axFit)
    
    # axLag.set_aspect('auto')
    #plt.plot(t[onInd:],I[onInd:]-Ifit)
    plt.plot(t,I-Ifit)
    #plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff))
    #plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    #plt.plot(toff+t[offInd],Ioff-IfitOff)
    plt.ylabel('$\mathrm{Residuals}$')# % relative error')
    plt.xlabel('$\mathrm{Time\ [ms]}$')
    plt.axhline(y=0,linestyle=':',color='k')
    #plt.setp(axRes.get_xticklabels(), visible=False)
    plt.xlim((0,totT))
    plt.tight_layout()
    
    Ifig.savefig(fDir+'fit'+str(nStates)+'states'+str(index)+"."+saveFigFormat, format=saveFigFormat)

    
    print("Fit has been plotted for the {}-state model".format(nStates))# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
        

#fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=None,method=methods[3])
def fit3states(fluxSet,run,vInd,params,method=methods[3]): #,Ipmax #Iss=None, ### Modify to pass in p3s
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
    
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    Ions = [None for phiInd in range(nPhis)]
    Ioffs = [None for phiInd in range(nPhis)]
    tons = [None for phiInd in range(nPhis)]
    toffs = [None for phiInd in range(nPhis)]
    phis = []
    Is = []
    ts = []
    Vs = []
    
    # Trim off phase data
    #frac = 1
    #chop = int(round(len(Ioffs[0])*frac))
    
    for phiInd in range(nPhis):
        targetPC = fluxSet.trials[run][phiInd][vInd]
        I = targetPC.I
        t = targetPC.t
        onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
        offInd = targetPC.pulseInds[0,1]
        Ions[phiInd] = I[onInd:offInd+1]
        Ioffs[phiInd] = I[offInd:] #[I[offInd:] for I in Is]
        tons[phiInd] = t[onInd:offInd+1]-t[onInd]
        toffs[phiInd] = t[offInd:]-t[offInd] #[t[offInd:]-t[offInd] for t in ts]
        #args=(Ioffs[phiInd][:chop+1],toffs[phiInd][:chop+1])
        phi = targetPC.phi
        phis.append(phi)
                
        Is.append(I)
        ts.append(t)
        V = targetPC.V
        Vs.append(V)
    
    nTrials = nPhis
    
    
    #p3s = Parameters()
    
    #p3s.add('g', value=gmax, vary=False) # derived['gmax'] = True
    #p3s.add('Gr0', value=Gr0, vary=False)
    
    
        
    ### 3a. Fit exponential to off curve to find Gd
    ### Fit off curve
    pOffs = Parameters() # Create parameter dictionary
    #print(-Ioff[0])
    
    ### Original single exponential fit
    #pOffs.add('A', value=-Ioff[0])#vary=False) #-1 # Constrain fitting to start at the beginning of the experimental off curve
    #pOffs.add('Gd', value=0.1, min=0)
    #offPmin = minimize(resid3off,pOffs,args=(Ioff,toff),method=meth)
    
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        pOffs.add('Islow_'+str(phiInd), value=0.1, vary=True, min=0)
        pOffs.add('Ifast_'+str(phiInd), value=0.5, vary=True, min=0)
    
    ### Initialise to 0.5*Gd and 2*Gd
    #pOffs.add('A', value=-1) # Islow
    pOffs.add('Gd1', value=params['Gd'].value/5, min=params['Gd'].min, max=params['Gd'].max) #0)#, min=0))
    #pOffs.add('B', value=-1) # Ifast
    pOffs.add('Gd2', value=params['Gd'].value*5, min=params['Gd'].min, max=params['Gd'].max) #100)#, min=0))
    #pOffs.add('C', value=0)
    #offPmin = minimize(resid3off,pOffs,args=(Ioff[0:100],toff[0:100]),method=meth)
    
    def fit3off(p,t,trial):
        #v = p.valuesdict()
        #A = p['A'].value
        #B = p['B'].value
        #C = p['C'].value
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        Gd1 = p['Gd1'].value
        Gd2 = p['Gd2'].value
        return -(Islow * np.exp(-Gd1*t) + Ifast * np.exp(-Gd2*t))#+C)
    
    #err3off = lambda p,I,t: I - fit3off(p,t)
    def err3off(p,Ioffs,toffs):
        return np.r_[ [Ioffs[i]/Ioffs[i][0] - fit3off(p,toffs[i],i)/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    offPmin = minimize(err3off, pOffs, args=(Ioffs,toffs), method=method)
    #pOffs.add('Gd', value=max(pOffs['Gd1'].value,pOffs['Gd2'].value), min=0)
    #nf = pOffs['A'].value + pOffs['B'].value
    #pOffs.add('Gd', value=pOffs['A'].value*pOffs['Gd1'].value/nf+pOffs['B'].value*pOffs['Gd2'].value/nf, min=0)
    
    Gds = [None for phiInd in range(nPhis)]
    for phiInd in range(nPhis):
        v = pOffs.valuesdict()
        Islow = v['Islow_'+str(phiInd)]
        Ifast = v['Ifast_'+str(phiInd)]
        Gds[phiInd] = (Islow*v['Gd1'] + Ifast*v['Gd2'])/(Islow + Ifast)
    
    print(Gds)
    Gd = np.mean(Gds)
    
    fig = plt.figure()
    gs = plt.GridSpec(nTrials,1)
    for trial in range(nTrials):
        Islow = pOffs['Islow_'+str(trial)].value
        Ifast = pOffs['Ifast_'+str(trial)].value
        eq = 'I(t)=-({Islow:.3g}*exp(-{Gd1:.3g}*t) + {Ifast:.3g}*exp(-{Gd2:.3g}*t))'.format(Islow=Islow, Ifast=Ifast, Gd1=v['Gd1'], Gd2=v['Gd2'])
        ax = fig.add_subplot(gs[trial,:])
        ax.plot(toffs[trial], Ioffs[trial], 'g', linewidth=mp.rcParams['lines.linewidth']*3, label='Data: phi={phi:.3g}'.format(phi=phis[trial])) # Experimental data
        ax.plot(toffs[trial], fit3off(pOffs,toffs[trial],trial), 'b', label=eq) # Fits
        ax.plot(toffs[trial], -(Islow+Ifast)*np.exp(-Gd*toffs[trial]), 'r', label='-{I0:.3g}*exp(-{Gd:.3g}*t)'.format(I0=Islow+Ifast, Gd=Gd))
        plt.legend(loc=4) # Lower right
        if trial < nTrials-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.xlabel('')
        ax.set_ylim(-1,0.1)
    
    print(">>> Off-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOffs)) # offPmin?
    print(offPmin.message)
    print("Chi^2 (Off-phase): ",offPmin.chisqr)
    #p3s.add('Gd', value=Gd, vary=False, min=params['Gd'].min, max=params['Gd'].max) #pOffs['Gd'].value
    
    #print(pOffs['Gd'].value)
    #print(p3s['Gd'].value)
    
    ### 3b. Optimise to find P and Gr at a particular phi from the on curve
    ### Fit on curve
    #pOns = Parameters()
    #pOns.add('a0', value=-1, expr='-a2')
    #pOns.add('a1', value=1)
    #pOns.add('a2', value=1)
    
    pOns = Parameters() #params #deepcopy(params)
    #pOns['Gd'].set(value=Gd, vary=False)
    
    # Set parameters from Off-curve optimisation
    copyParam('Gd',params,pOns)
    #pOns.add('Gd', value=Gd, vary=False)
    pOns['Gd'].set(value=Gd, vary=False)
    
    copyParam('k',params,pOns)
    copyParam('p',params,pOns) #pOns.add('p',value=0.7,min=0.1,max=5)
    copyParam('phim',params,pOns) #pOns.add('phim',value=1e17,min=1e15,max=1e19) #1e19
    
    # Set parameters from general rhodopsin analysis routines
    copyParam('g',params,pOns) # pOns.add('g',value=gmax,vary=False)
    copyParam('Gr0',params,pOns)
    #pOns['Gr0'].set(vary=False) # Dark recovery rate
    copyParam('Gr1',params,pOns) # pOns.add('Gr',value=Gr0,vary=False)
    copyParam('E',params,pOns)
    if params['useIR'].value==True: # Redundant?
        copyParam('v0',params,pOns)
        copyParam('v1',params,pOns)
    
    
    
    RhO = models['3']()
    
    def fit3on(p,t,RhO,phi,V):
        RhO.updateParams(p)
        RhO.setLight(phi)
        #states = RhO.calcSoln(t, s0=[1,0,0]) # t starts at 0
        states = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
        return RhO.calcI(V, states)
        
    def err3on(p,Ions,tons,RhO,phis,Vs):
        return np.r_[ [Ions[i]/Ions[i][-1] - fit3on(p,tons[i],RhO,phis[i],Vs[i])/Ions[i][-1] for i in range(len(Ions))] ]
    
    onPmin = minimize(err3on, pOns, args=(Ions,tons,RhO,phis,Vs), method=method)
    
    
    
    
    # # Create dummy parameters for each phi
    # for phiInd in range(nPhis):
        # pOffs.add('a0_'+str(phiInd), value=-1, expr='-a2_{}'.format(phiInd))
        # pOffs.add('a1_'+str(phiInd), value=1)
        # pOffs.add('a2_'+str(phiInd), value=1)
    
    # pOns.add('tau_act', value=2, min=0)
    # pOns.add('tau_deact', value=10, min=0)
    
    # def fit3on(p,t,trial):
        # """Fit a biexponential curve to the on-phase to find lambdas"""
        # a0 = p['a0_'+str(trial)].value
        # a1 = p['a1_'+str(trial)].value
        # a2 = p['a2_'+str(trial)].value
        # #a0 = p['a0'].value
        # #a1 = p['a1'].value
        # #a2 = p['a2'].value
        # tau_act = p['tau_act'].value
        # tau_deact = p['tau_deact'].value
        # return -(a0 + a1*(1-np.exp(-t/tau_act)) + a2*np.exp(-t/tau_deact))
    
    # def err3on(p,Ions,tons):
        # return np.r_[ [Ions[i]/Ions[i][-1] - fit3on(p,tons[i],i)/Ions[i][-1] for i in range(len(Ions))] ]
    
    # onPmin = minimize(err3on, pOns, args=(Ions,tons), method=method)
    
    # # Calculate model parameters from fitting parameters
    # #tau_act = pOns['tau_act'].value
    # #tau_deact = pOns['tau_deact'].value
    # lam1 = 1/pOns['tau_deact'].value
    # lam2 = 1/pOns['tau_act'].value
    
    # #Gd = pOffs['Gd'].value
    
    # # findPeaks in experimental photocurrent
    # # if (V < E): # Find Minima
        # # Ip = min(I)
    # # else:       # Find Maxima
        # # Ip = max(I)
    # # Oss = Iplat/(gmax*(Vsaturate-E)) # f(V) here?
    # Iplats = [findPlateauCurrent(Ions[phiInd], p=0.1) for phiInd in range(nPhis)]
    # ###Ipmax = gmax*(V-E) # This is Ipmax (from saturate)
    # Oss = Iplat/Ipmax
    # alpha = (lam1+lam2)/Gd - 1
    # beta = alpha/((1/Oss)-1)
    
    # if (alpha**2 - 4*beta) < 0:
        # print('\n\nEarly Warning! No real solution exists!\n\n')
        # ### Rerun off-curve fitting with single exponential
    
    # x1 = (alpha + np.sqrt(alpha**2 - 4*beta))/2
    # Gr1 = x1 * Gd
    # P1 = lam1+lam2 - (Gr1+Gd)
    # print("Solution 1: P = {}; Gd = {}; Gr = {}".format(P1,Gd,Gr1))
    
    # x2 = (alpha - np.sqrt(alpha**2 - 4*beta))/2
    # Gr2 = x2 * Gd
    # P2 = lam1+lam2 - (Gr2+Gd)
    # print("Solution 2: P = {}; Gd = {}; Gr = {}".format(P2,Gd,Gr2))
    
    # if P1 > Gr1:
        # print("==> Solution 1 (assuming that P > Gr): {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(P1,Gd,Gr1))
        # P = P1
        # Gr = Gr1
    # else:
        # print("==> Solution 2 (assuming that P > Gr): {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(P2,Gd,Gr2))
        # P = P2
        # Gr = Gr2
    # #print("==> Original values: {{P={:.3}, Gd={:.3}, Gr={:.3}}}".format(  ,Gd,ChR.Gr))
    
    # if (Gr + Gd - 2*np.sqrt(Gr*Gd)) < P < (Gr + Gd + 2*np.sqrt(Gr*Gd)):
        # print('\n\nWarning! No real solution exists!\n\n')
    
    print("\n\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOns)) # onPmin?
    print(onPmin.message)
    print("Chi^2 (On-phase): ",onPmin.chisqr)
    
    #k = P/phi # To find the coefficient which scales phi to calculate the activation rate
    #p3s.add('k', value=k, vary=False)
    # N.B. The three-state model does not generalise well to other values of phi
    
    #GrPhi = Gr - Gr0 # To find Grlight for a particular phi
    #p3s.add('GrPhi', value=GrPhi, vary=False)

#     ci = lmfit.conf_interval(pOffs)
#     lmfit.printfuncs.report_ci(ci)
    

    ### 4. Optionally fit f(V) parameters with rectifier data
    # if 'rectifier' in dataset:
        # if hasattr(dataset['rectifier'], 'Iss'): # Use extracted values
            # print("Finish me!")
        # else: # Extract steady state values
            # print("Finish me!")
        # # Fit Curve of V vs Iss
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?
    
    # if verbose > 1:
        # for key, value in p3s.items() :
            # print (key, value)
    
    # RhO = selectModel(3)
    # #RhO.setParams(d3sp)
    # RhO.g = p3s['g'].value              # 16700     [pS]
    # RhO.k = p3s['k'].value              # 0.545e-14 [ms^-1 * photons^-1 * s * mm^2]
    # RhO.Gr0 = p3s['Gr0'].value          # 1/5000    [ms^-1] #Gr_dark
    # RhO.Gr1 = p3s['GrPhi'].value        # 1/165     [ms^-1] #Gr_light
    # RhO.Gd = p3s['Gd'].value            # 1/11      [ms^-1]
    # RhO.phiFit = phi                    # Flux intensity at which the parameters were fit. 
    # RhO.setLight(0)                     # Re-initialise model to dark state
    # RhO.reportParams()
    
    # ### Plot curves
    # totT = max(t)
    # Ifig = plt.figure()
    # gsPL = plt.GridSpec(4,1)
    # axFit = Ifig.add_subplot(gsPL[:-1,:])
    
    # #ax = Ifig.add_subplot(111)
    # # for p in range(0, nPulses):
        # # plt.axvspan(delD+(p*(onD+offD)),delD+(p*(onD+offD))+onD,facecolor='y',alpha=0.2)
    # plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    # plt.xlim((0,totT))
    # #plt.xlabel('$\mathrm{Time\ [ms]}$')
    # plt.setp(axFit.get_xticklabels(), visible=False)
    # plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    # plt.plot(t,I,color='b',label='$\mathrm{Experimental\ Data}$')
    # RhO.setLight(phi)
    # if addTitles:
        # plt.title('Three-state model fit to data (phi={:.3g}): [P={:.3g}; Gd={:.3g}; Gr={:.3g}]'.format(phi,RhO.P,RhO.Gd,RhO.Gr))
    # RhO.setLight(0)
    
    
    # if 1: # Use pseudo parameters
        # IfitOn = calc3on(pOns, ton) # onPmin
        # pOffs['A'].value = -IfitOn[-1] # Constrain off curve to start at the end of the fitted on curve
        # IfitOff = calc3off(pOffs, toff) # offPmin
        # ###pOffs.add('A',value=np.exp(pOffs['B'].value),vary=False)
    # elif 0: # Use parameters
        # I_RhO, t, states = runTrial(RhO, 1, V,phi,t[onInd],t[offInd]-t[onInd],totT-t[offInd],0,0.01) # Problem matching time arrays
        # #I_RhO, t, states = runTrial(RhO, 1, V,phi,round(t[onInd]),round(t[offInd]-t[onInd]),round(totT-t[offInd]),0,0.01)
        # IfitOn = I_RhO[onInd:offInd+1]
        # IfitOff = I_RhO[offInd:]
        # #Prot = protocols['custom']([phi], [-70], [[t[onInd],t[offInd]]], totT, 1, 0.1)
        # #Prot.runProtocol(RhO)
        # #Prot.plotProtocol()
        # #IfitOn = Prot.Is[0][0][0][onInd:offInd]
        # #IfitOff = Prot.Is[0][0][0][offInd:]
    # else: # Use analytic solution
        # #print(ton)
        # RhO.setLight(phi)
        # RhO.reportState()
        # onStates = RhO.calcSoln(ton)
        # IfitOn = RhO.calcI(V, onStates)
        # print(onStates)
        # RhO.setLight(0)
        # #print(toff)
        # IfitOff = RhO.calcI(V, RhO.calcSoln(toff, onStates[-1,:]))
    
# #     tfit = np.append(ton,toff)
# #     plot(tfit,Ifit,color='r')
    # plt.plot(ton+t[onInd],IfitOn,color='g',label='$\mathrm{Three-state\ model\ fit}$')
    # plt.plot(toff+t[offInd],IfitOff,color='g')
    # plt.legend(loc='best')
    

    # axRes = Ifig.add_subplot(gsPL[-1,:],sharex=axFit)
    
    # # axLag.set_aspect('auto')
                    
    # # Rfig = plt.figure() # Replace this with a subplot...
    # #print(IfitOn)
    # #print(Ion)
    # #print((IfitOn-Ion)/Ion)
    # plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    # plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    # #plt.plot(toff+t[offInd],Ioff-IfitOff)
    # plt.ylabel('$\mathrm{Residuals}$') # % relative error')
    # plt.xlabel('$\mathrm{Time\ [ms]}$')
    # #plt.setp(axRes.get_xticklabels(), visible=False)
    # plt.xlim((0,totT))
    # plt.tight_layout()
    
    #print(Is)
    #print(ts)
    #print(phis)
    #print(Vs)
    #print(pOns)
    RhO.setParams(pOns)
    
    nStates = 3
    for trial in range(nTrials):
        plotFit(Is[trial],ts[trial],onInd,offInd,phis[trial],Vs[trial],nStates,pOns,fitRates=False,index=trial)
        RhO.setLight(phis[trial])
        print('phi = {:.3g}: P = {:.3g}, Gd = {:.3g}, Gr = {:.3g}'.format(phis[trial], RhO.P, RhO.Gd, RhO.Gr))
    
    if verbose > 0:
        print("Parameters have been fit for the three-state model over a flux range of [{:.3g}, {:.3g}] [photons * s^-1 * mm^-2]".format(min(phis), max(phis)))
    
    return pOns #p3s # RhO


def fit4states(fluxSet,run,vInd,params,method=methods[3]): #fit4states(Is,ts,onInds,offInds,phis,V,Gr,gmax,params=None,method=methods[3]):#,Iss): # ,Ipeak
    
    
    # Specify initial values (and bounds) for dummy parameters for each phi
    #for phiInd in range(nPhis):
    #    pOffs.add('Islow_'+str(phiInd), value=0.1, vary=True, min=0)
    #    pOffs.add('Ifast_'+str(phiInd), value=0.5, vary=True, min=0)


    #Ion, Ioff = I[onInd:offInd+1], I[offInd:]
    #ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    Ions = [None for phiInd in range(nPhis)]
    Ioffs = [None for phiInd in range(nPhis)]
    tons = [None for phiInd in range(nPhis)]
    toffs = [None for phiInd in range(nPhis)]
    phis = []
    Is = []
    ts = []
    Vs = []
    
    # Trim off phase data
    #frac = 1
    #chop = int(round(len(Ioffs[0])*frac))
    
    for phiInd in range(nPhis):
        targetPC = fluxSet.trials[run][phiInd][vInd]
        I = targetPC.I
        t = targetPC.t
        onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
        offInd = targetPC.pulseInds[0,1]
        Ions[phiInd] = I[onInd:offInd+1]
        Ioffs[phiInd] = I[offInd:] #[I[offInd:] for I in Is]
        tons[phiInd] = t[onInd:offInd+1]-t[onInd]
        toffs[phiInd] = t[offInd:]-t[offInd] #[t[offInd:]-t[offInd] for t in ts]
        #args=(Ioffs[phiInd][:chop+1],toffs[phiInd][:chop+1])
        phi = targetPC.phi
        phis.append(phi)
                
        Is.append(I)
        ts.append(t)
        V = targetPC.V
        Vs.append(V)
    
    
    
    ### OFF PHASE
    ### 3a. OFF CURVE: Fit biexponential to off curve to find lambdas
    ### Fit off curve - if initial conditions can be calculated, it might be better to use an analytic solution relating directly to model parameters c.f. analysis_4state_off_new.m

    pOffs = Parameters() # Create parameter dictionary
    # pOffs.add_many(('Gd1',0.15,True,0.01,1,None), #('Gd1',0.11,False,0.01,1,None)
                    # ('Gd2',0.025,True,0.01,1,None), #('Gd2',0.023,False,0.01,1,None)
                    # ('e12d',0.01,True,0,1,None),
                    # ('e21d',0.01,True,0,1,None))
    copyParam('Gd1',params,pOffs)
    copyParam('Gd2',params,pOffs)
    copyParam('e12d',params,pOffs)
    copyParam('e21d',params,pOffs)
    
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        #pOffs.add('a0_'+str(phiInd), value=0, vary=False)
        pOffs.add('Islow_'+str(phiInd), value=0.1, vary=True, min=0)
        pOffs.add('Ifast_'+str(phiInd), value=0.5, vary=True, min=0)
    
    
    # lam1 + lam2 == Gd1 + Gd2 + e12d + e21d
    # lam1 * lam2 == Gd1*Gd2 + Gd1*e21d + Gd2*e12d
        
    calcB = lambda Gd1, Gd2, e12d, e21d: (Gd1 + Gd2 + e12d + e21d)/2
    calcC = lambda b, Gd1, Gd2, e12d, e21d: np.sqrt(b**2 - (Gd1*Gd2 + Gd1*e21d + Gd2*e12d))
    
    def lams(p): #(Gd1, Gd2, e12d, e21d):
        Gd1 = p['Gd1'].value
        Gd2 = p['Gd2'].value
        e12d = p['e12d'].value
        e21d = p['e21d'].value
        #v = p.valuesdict()
        b = calcB(Gd1, Gd2, e12d, e21d)
        c = calcC(b, Gd1, Gd2, e12d, e21d)
        return (b-c, b+c)
    
    def fit4off(p,t,trial):
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        lam1, lam2 = lams(p)
        return -(Islow*np.exp(-lam1*t) + Ifast*np.exp(-lam2*t))
        
    def err4off(p,Ioffs,toffs):
        """Normalise by the first element of the off-curve""" # [-1]
        return np.r_[ [Ioffs[i]/Ioffs[i][0] - fit4off(p,toffs[i],i)/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    #fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-lams(p)[0]*t) + p['a2'].value*np.exp(-lams(p)[1]*t))
    ##fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-p['lam1'].value*t) + p['a2'].value*np.exp(-p['lam2'].value*t))
    #errfunc = lambda p, Ioff, toff: Ioff - fitfunc(p,toff)
    
    offPmin = minimize(err4off, pOffs, args=(Ioffs,toffs), method=method)#, fit_kws={'maxfun':100000})
    
    print(">>> Off-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOffs))
        print("Error bars: ",offPmin.errorbars)
    #print('lambda1 = {:.5g}; lambda2 = {:.5g}'.format(pOff['lam1'].value, pOff['lam2'].value))
    print(offPmin.message)
    print("Chi^2: ",offPmin.chisqr)
    
    if verbose > 1:
        for phiInd in range(nPhis):
            print('phi_',str(phiInd),': Islow = ',pOffs['Islow_'+str(phiInd)].value,'; Ifast = ',pOffs['Ifast_'+str(phiInd)].value)
            #print('a2_'+str(phiInd),' = ',pOffs['a2_'+str(phiInd)].value)
    
    nTrials = nPhis
    fig = plt.figure()
    gs = plt.GridSpec(nTrials,1)
    #axf = axS.twinx()
    #axf.set_yscale('log')
    #axf.set_ylabel('$f\ \mathrm{[Hz]}$')
    for trial in range(nTrials):
        #a0 = pOffs['a0_'+str(trial)].value
        Islow = pOffs['Islow_'+str(trial)].value
        Ifast = pOffs['Ifast_'+str(trial)].value
        lam1,lam2 = lams(pOffs)
        #eq = '-({a0:.3g} + {a1:.3g}*exp(-{lam1:.3g}*t) + {a2:.3g}*exp(-{lam2:.3g}*t))'.format(a0=a0, a1=a1, a2=a2, lam1=lam1, lam2=lam2)
        eq = 'I(t)=-({Islow:.3g}*exp(-{lam1:.3g}*t) + {Ifast:.3g}*exp(-{lam2:.3g}*t))'.format(Islow=Islow, Ifast=Ifast, lam1=lam1, lam2=lam2)
        ax = fig.add_subplot(gs[trial,:])
        ax.plot(toffs[trial], Ioffs[trial], 'g', linewidth=mp.rcParams['lines.linewidth']*3, label='Data: phi={phi:.3g}'.format(phi=phis[trial])) # Experimental data
        ax.plot(toffs[trial], fit4off(pOffs,toffs[trial],trial), 'b', label=eq) # Fits
        plt.legend(loc=4) # Lower right
        if trial < nTrials-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.xlabel('')
        ax.set_ylim(-1,0.1)
        
    print('Gd1 = {}; Gd2 = {}; e12d = {}; e21d = {}'.format(pOffs['Gd1'].value, pOffs['Gd2'].value, pOffs['e12d'].value, pOffs['e21d'].value))
    
    pOffs['Gd1'].vary=False
    pOffs['Gd2'].vary=False
    pOffs['e12d'].vary=False
    pOffs['e21d'].vary=False
    
    
    
    
    
    
    ### ON PHASE
    
    pOns = Parameters() # deepcopy(params)
    
    # Set parameters from Off-curve optimisation
    copyParam('Gd1',pOffs,pOns)
    copyParam('Gd2',pOffs,pOns)
    copyParam('e12d',pOffs,pOns)
    copyParam('e21d',pOffs,pOns)
    
    # phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=pOns,method=method)
    copyParam('k1',params,pOns) #pOns.add('k1',value=3, min=0.01) #0.5 Ga1 = k1 * phi
    copyParam('k2',params,pOns) #pOns.add('k2',value=1.5, min=0.01) #0.2 Ga2 = k2 * phi
    copyParam('c1',params,pOns) #pOns.add('c1',value=0.05, min=0.01)
    copyParam('c2',params,pOns) #pOns.add('c2',value=0.01, min=0.01)
    copyParam('gam',params,pOns) #pOns.add('gam',value=0.05, min=0, max=1)
    # Place RhO,V,phi in onPmin?
    copyParam('p',params,pOns) #pOns.add('p',value=0.7,min=0.1,max=5)
    copyParam('q',params,pOns)
    copyParam('phim',params,pOns) #pOns.add('phim',value=1e17,min=1e15,max=1e19) #1e19
    
    # Set parameters from general rhodopsin analysis routines
    #Gr0,gmax
    copyParam('g',params,pOns) # pOns.add('g',value=gmax,vary=False)
    copyParam('Gr',params,pOns) # pOns.add('Gr',value=Gr0,vary=False)
    #copyParam('phi0',params,pOns) # pOns.add('phi0',value=5e18,min=1e14,max=1e20,vary=False) ################# Set this to be above the max flux??? 10**ceil(log10(max(phis)))
    copyParam('E',params,pOns)
    if params['useIR'].value==True:
        copyParam('v0',params,pOns)
        copyParam('v1',params,pOns)
    
    RhO = models['4']()
    
    # Normalise? e.g. /Ions[trial][-1] or /min(Ions[trial])
    def err4On(p,Ions,tons,RhO,Vs,phis):
        return np.r_[ [Ions[i]/Ions[i][-1] - calc4on(p,tons[i],RhO,Vs[i],phis[i])/Ions[i][-1] for i in range(len(Ions))]]
    
    ### Trim down ton? Take 10% of data or one point every ms? ==> [0::5]
    ### Instead try looping over a coarse grid of parameter values and saving RMSE for each combination c.f. analysis_4state_on_new.m
    
    if verbose > 1:
        print('Optimising',end='')
    onPmin = minimize(err4On, pOns, args=(Ions,tons,RhO,Vs,phis), method=method)
    
    print("\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOns))
        print("Error bars: ",onPmin.errorbars)
    if not onPmin.success:
        print(onPmin.success)
        print(onPmin.message)
        print(onPmin.lmdif_message)
    print("Chi^2: ",onPmin.chisqr)
    
    print('k1 = {}; k2 = {}; c1 = {}; c2 = {}'.format(pOns['k1'].value, pOns['k2'].value, pOns['c1'].value, pOns['c2'].value))
    print('gam = ',pOns['gam'].value)
    #print('phi0 = ',pOns['phi0'].value)
    print('phim = ',pOns['phim'].value)
    print('p = ',pOns['p'].value)
    print('q = ',pOns['q'].value)

    nStates = 4
    for trial in range(nTrials):
        plotFit(Is[trial],ts[trial],onInd,offInd,phis[trial],Vs[trial],nStates,pOns,fitRates=False,index=trial)


    if verbose > 0:
        print("Parameters have been fit for the four-state model")# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
    
    return pOns #p4s #RhO





def fit6states(fluxSet,dataSet,run,vInd,params,method=methods[3]): #shortPulseSet #fit4states(Is,ts,onInds,offInds,phis,V,Gr,gmax,params=None,method=methods[3]):#,Iss): # ,Ipeak
    
    
    # Specify initial values (and bounds) for dummy parameters for each phi
    #for phiInd in range(nPhis):
    #    pOffs.add('Islow_'+str(phiInd), value=0.1, vary=True, min=0)
    #    pOffs.add('Ifast_'+str(phiInd), value=0.5, vary=True, min=0)


    #Ion, Ioff = I[onInd:offInd+1], I[offInd:]
    #ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    Ions = [None for phiInd in range(nPhis)]
    Ioffs = [None for phiInd in range(nPhis)]
    tons = [None for phiInd in range(nPhis)]
    toffs = [None for phiInd in range(nPhis)]
    phis = []
    Is = []
    ts = []
    Vs = []
    
    Icycles = []
    nfs = [] # Normalisation factors: e.g. /Ions[trial][-1] or /min(Ions[trial])
    
    # Trim off phase data
    #frac = 1
    #chop = int(round(len(Ioffs[0])*frac))
    
    for phiInd in range(nPhis):
        targetPC = fluxSet.trials[run][phiInd][vInd]
        I = targetPC.I
        t = targetPC.t
        onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
        offInd = targetPC.pulseInds[0,1]
        Ions[phiInd] = I[onInd:offInd+1]
        Ioffs[phiInd] = I[offInd:] #[I[offInd:] for I in Is]
        tons[phiInd] = t[onInd:offInd+1]-t[onInd]
        toffs[phiInd] = t[offInd:]-t[offInd] #[t[offInd:]-t[offInd] for t in ts]
        #args=(Ioffs[phiInd][:chop+1],toffs[phiInd][:chop+1])
        phi = targetPC.phi
        phis.append(phi)
                
        Is.append(I)
        ts.append(t)
        V = targetPC.V
        Vs.append(V)
        
        Icycles.append(I[onInd:])
        nfs.append(I[offInd])
    
    
    ### OFF PHASE
    ### 3a. OFF CURVE: Fit biexponential to off curve to find lambdas
    ### Fit off curve - if initial conditions can be calculated, it might be better to use an analytic solution relating directly to model parameters c.f. analysis_4state_off_new.m
    
    
    pOffs = Parameters() # Create parameter dictionary
    copyParam('b1',params,pOffs) #Gd1
    copyParam('a4',params,pOffs) #Gd2
    copyParam('a30',params,pOffs) #Gf0
    copyParam('b20',params,pOffs) #Gb0
    
    ### Trim the first 10% of the off curve to allow I1 and I2 to empty?
    
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        #pOffs.add('a0_'+str(phiInd), value=0, vary=False)
        pOffs.add('Islow_'+str(phiInd), value=0.1, vary=True, min=0)
        pOffs.add('Ifast_'+str(phiInd), value=0.5, vary=True, min=0)
    
    
    # lam1 + lam2 == Gd1 + Gd2 + e12d + e21d
    # lam1 * lam2 == Gd1*Gd2 + Gd1*e21d + Gd2*e12d
        
    calcB = lambda b1, a4, a30, b20: (b1 + a4 + a30 + b20)/2 #Gd1, Gd2, e12d, e21d: (Gd1 + Gd2 + e12d + e21d)/2
    calcC = lambda b, b1, a4, a30, b20: np.sqrt(b**2 - (b1*a4 + b1*b20 + a4*a30)) #b, Gd1, Gd2, e12d, e21d: np.sqrt(b**2 - (Gd1*Gd2 + Gd1*e21d + Gd2*e12d))
    
    def lams(p): #(Gd1, Gd2, e12d, e21d):
        Gd1 = p['b1'].value
        Gd2 = p['a4'].value
        Gf0 = p['a30'].value
        Gb0 = p['b20'].value
        #v = p.valuesdict()
        b = calcB(Gd1, Gd2, Gf0, Gb0)
        c = calcC(b, Gd1, Gd2, Gf0, Gb0)
        return (b-c, b+c)
    
    def fit6off(p,t,trial):
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        lam1, lam2 = lams(p)
        return -(Islow*np.exp(-lam1*t) + Ifast*np.exp(-lam2*t))
        
    def err6off(p,Ioffs,toffs):
        """Normalise by the first element of the off-curve""" # [-1]
        return np.r_[ [Ioffs[i]/Ioffs[i][0] - fit6off(p,toffs[i],i)/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    #fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-lams(p)[0]*t) + p['a2'].value*np.exp(-lams(p)[1]*t))
    ##fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-p['lam1'].value*t) + p['a2'].value*np.exp(-p['lam2'].value*t))
    #errfunc = lambda p, Ioff, toff: Ioff - fitfunc(p,toff)
    
    offPmin = minimize(err6off, pOffs, args=(Ioffs,toffs), method=method)#, fit_kws={'maxfun':100000})
    
    print(">>> Off-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOffs))
        print("Error bars: ",offPmin.errorbars)
    #print('lambda1 = {:.5g}; lambda2 = {:.5g}'.format(pOff['lam1'].value, pOff['lam2'].value))
    print(offPmin.message)
    print("Chi^2: ",offPmin.chisqr)
    
    if verbose > 1:
        for phiInd in range(nPhis):
            print('phi_',str(phiInd),': Islow = ',pOffs['Islow_'+str(phiInd)].value,'; Ifast = ',pOffs['Ifast_'+str(phiInd)].value)
            #print('a2_'+str(phiInd),' = ',pOffs['a2_'+str(phiInd)].value)
    
    nTrials = nPhis
    fig = plt.figure()
    gs = plt.GridSpec(nTrials,1)
    #axf = axS.twinx()
    #axf.set_yscale('log')
    #axf.set_ylabel('$f\ \mathrm{[Hz]}$')
    for trial in range(nTrials):
        #a0 = pOffs['a0_'+str(trial)].value
        Islow = pOffs['Islow_'+str(trial)].value
        Ifast = pOffs['Ifast_'+str(trial)].value
        lam1,lam2 = lams(pOffs)
        #eq = '-({a0:.3g} + {a1:.3g}*exp(-{lam1:.3g}*t) + {a2:.3g}*exp(-{lam2:.3g}*t))'.format(a0=a0, a1=a1, a2=a2, lam1=lam1, lam2=lam2)
        eq = 'I(t)=-({Islow:.3g}*exp(-{lam1:.3g}*t) + {Ifast:.3g}*exp(-{lam2:.3g}*t))'.format(Islow=Islow, Ifast=Ifast, lam1=lam1, lam2=lam2)
        ax = fig.add_subplot(gs[trial,:])
        ax.plot(toffs[trial], Ioffs[trial], 'g', linewidth=mp.rcParams['lines.linewidth']*3, label='Data: phi={phi:.3g}'.format(phi=phis[trial])) # Experimental data
        ax.plot(toffs[trial], fit6off(pOffs,toffs[trial],trial), 'b', label=eq) # Fits
        plt.legend(loc=4) # Lower right
        if trial < nTrials-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.xlabel('')
        ax.set_ylim(-1,0.1)
        
    print('Gd1 = {}; Gd2 = {}; Gf0 = {}; Gb0 = {}'.format(pOffs['b1'].value, pOffs['a4'].value, pOffs['a30'].value, pOffs['b20'].value))
    
    pOffs['b1'].vary=False #Gd1
    pOffs['a4'].vary=False #Gd2
    pOffs['a30'].vary=False #Gf0
    pOffs['b20'].vary=False #Gb0
    
    
    ### Calculate Go (1/tau_opsin)
    print('Calculating opsin activation rate')
    # Assume that Gd1 > Gd2
    # Assume that Gd = Gd1 for short pulses
    
    def solveGo(tlag, Gd, Go0=1000, tol=1e-9):
        Go, Go_m1 = Go0, 0
        print(tlag, Gd, Go, Go_m1)
        while abs(Go_m1 - Go) > tol:
            Go_m1 = Go
            Go = ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
            #Go_m1, Go = Go, ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
            #print(Go, Go_m1)
        return Go
            
    if 'shortPulse' in dataSet: # Fit Gret
        from scipy.optimize import curve_fit
        # Fit tpeak = tpulse + tmaxatp0 * np.exp(-k*tpulse)
        #dataSet['shortPulse'].getProtPeaks()
        #tpeaks = dataSet['shortPulse'].IrunPeaks
        PD = dataSet['shortPulse']
        PCs = [PD.trials[p][0][0] for p in range(PD.nRuns)]
        #[pc.alignToTime() for pc in PCs]
        
        #tpeaks = np.asarray([PD.trials[p][0][0].tpeak for p in range(PD.nRuns)]) # - PD.trials[p][0][0].t[0]
        #tpulses = np.asarray([PD.trials[p][0][0].onDs[0] for p in range(PD.nRuns)])
        tpeaks = np.asarray([pc.tpeak for pc in PCs])
        tpulses = np.asarray([pc.onDs[0] for pc in PCs])
        
        devFunc = lambda tpulses, t0, k: tpulses + t0 * np.exp(-k*tpulses)
        p0 = (0,1)
        popt, pcov = curve_fit(devFunc, tpulses, tpeaks, p0=p0)
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        tsmooth = np.linspace(0,max(tpulses),101)
        ax.plot(tpulses,tpeaks,'x')
        ax.plot(tsmooth,devFunc(tsmooth,*popt))
        ax.plot(tsmooth,tsmooth,'--')
        ax.set_ylim([0,max(tpulses)+5])
        ax.set_xlim([0,max(tpulses)+5])
        
        #plt.tight_layout()
        #plt.axis('equal')
        
        #print(tpulses)
        #print(tpeaks)
        #print(*popt)
        #print(devFunc(tpeaks,*popt))
        
        # Solve iteratively Go = ((tlag*Gd) - np.log(Gd/Go))/tlag
        
        Gd1 = pOffs['b1'].value
        Go = solveGo(tlag=popt[0], Gd=Gd1, Go0=1000, tol=1e-9)
        print('t_lag = {}; Gd = {} --> Go = {}'.format(popt[0], Gd1, Go))
        # Gd2 = pOffs['a4'].value
        # Go2 = solveGo(tlag=popt[0], Gd=Gd2, Go0=1000, tol=1e-9)
        # print('Go2 = ', Go2)
        
    elif 'saturate' in dataSet:
        PD = dataSet['saturate']
        #PCs = [PD.trials[p][0][0] for p in range(PD.nRuns)]
        PC = PD.trials[0][0][0]
        tlag = PC.tlag ############################### Add to Photocurrent...
        
        Go = solveGo(tlag=tlag, Gd=Gd1, Go0=1000, tol=1e-9)
        #print('Go = ', Go)
        print('t_lag = {}; Gd = {} --> Go = {}'.format(tlag, Gd1, Go))
    
    else:
        Go = 1 # Default
        print('No data found to estimate Go: defaulting to Go = {}'.format(Go))
    
    
    ### ON PHASE
    
    pOns = Parameters() # deepcopy(params)
    
    # Set parameters from Off-curve optimisation
    copyParam('b1',pOffs,pOns) #Gd1
    copyParam('a4',pOffs,pOns) #Gd2
    copyParam('a30',pOffs,pOns) #Gf0
    copyParam('b20',pOffs,pOns) #Gb0
    
    # Set parameters from short pulse calculations
    pOns.add('a2', value=Go, vary=False, min=1e-9, max=1e9) #Go1
    pOns.add('b3', value=Go, vary=False, min=1e-9, max=1e9) #Go2
    #pOns.add('b3', value=Go2, vary=False, min=1e-9, max=1e9) #Go2
    #pOns.add('a2', value=Go, vary=True, min=Go/2, max=Go*2) #Go1
    #pOns.add('b3', value=Go, vary=True, min=Go/2, max=Go*2) #Go2
    
    # phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=pOns,method=method)
    copyParam('a10',params,pOns) #k1 #pOns.add('k1',value=3, min=0.01) #0.5 Ga1 = k1 * phi
    copyParam('b40',params,pOns) #k2 #pOns.add('k2',value=1.5, min=0.01) #0.2 Ga2 = k2 * phi
    copyParam('a31',params,pOns) #kf #pOns.add('c1',value=0.05, min=0.01)
    copyParam('b21',params,pOns) #kb #pOns.add('c2',value=0.01, min=0.01)
    copyParam('gam',params,pOns) #pOns.add('gam',value=0.05, min=0, max=1)
    # Place RhO,V,phi in onPmin?
    copyParam('p',params,pOns) #pOns.add('p',value=0.7,min=0.1,max=5)
    copyParam('q',params,pOns)
    copyParam('phim',params,pOns) #pOns.add('phim',value=1e17,min=1e15,max=1e19) #1e19
    
    # Set parameters from general rhodopsin analysis routines
    #Gr0,gmax
    copyParam('g',params,pOns) # pOns.add('g',value=gmax,vary=False)
    copyParam('a6',params,pOns) #Gr # pOns.add('Gr',value=Gr0,vary=False)
    #copyParam('phi0',params,pOns) # pOns.add('phi0',value=5e18,min=1e14,max=1e20,vary=False) ################# Set this to be above the max flux??? 10**ceil(log10(max(phis)))
    copyParam('E',params,pOns)
    if params['useIR'].value==True:
        copyParam('v0',params,pOns)
        copyParam('v1',params,pOns)
    
    RhO = models['6']()
    
    def calc6on(p,t,RhO,V,phi):
        """Simulate the on-phase from base parameters for the 6-state model"""
        if verbose > 1:
            print('.', end="") # sys.stdout.write('.')
        RhO.initStates(0)
        #RhO.setLight(phi) # Calculate transition rates for phi
        RhO.updateParams(p)
        RhO.setLight(phi) # Calculate transition rates for phi
        
        soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
        # soln,out = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian, full_output=True)
        # if out['message'] != 'Integration successful.':
            # #print(out)
            # print(RhO.reportParams())
        I_RhO = RhO.calcI(V, soln)
        return I_RhO
    
    # Normalise? e.g. /Ions[trial][-1] or /min(Ions[trial])
    def err6on(p,Ions,tons,RhO,Vs,phis):
        return np.r_[ [Ions[i]/Ions[i][-1] - calc6on(p,tons[i],RhO,Vs[i],phis[i])/Ions[i][-1] for i in range(len(Ions))]]
    
    ### Trim down ton? Take 10% of data or one point every ms? ==> [0::5]
    ### Instead try looping over a coarse grid of parameter values and saving RMSE for each combination c.f. analysis_4state_on_new.m
    
    if verbose > 1:
        print('Optimising',end='')
    onPmin = minimize(err6on, pOns, args=(Ions,tons,RhO,Vs,phis), method=method)
    
    print("\n>>> On-phase fitting summary: <<<")
    if verbose > 1:
        print(fit_report(pOns))
        print("Error bars: ",onPmin.errorbars)
    if not onPmin.success:
        print(onPmin.success)
        print(onPmin.message)
        print(onPmin.lmdif_message)
    print("Chi^2: ",onPmin.chisqr)
    
    fitParams = pOns
    
    #print('k1 = {}; k2 = {}; c1 = {}; c2 = {}'.format(pOns['k1'].value, pOns['k2'].value, pOns['c1'].value, pOns['c2'].value))
    print('k1 = {}; k2 = {}; kf = {}; kb = {}'.format(pOns['a10'].value, pOns['b40'].value, pOns['a31'].value, pOns['b21'].value))
    print('gam = ',pOns['gam'].value)
    #print('phi0 = ',pOns['phi0'].value)
    print('phim = ',pOns['phim'].value)
    print('p = ',pOns['p'].value)
    print('q = ',pOns['q'].value)


    def calcCycle(p,ton,toff,RhO,V,phi,fitRates=False):
        """Simulate the on-phase from base parameters for the 6-state model"""
        if verbose > 1:
            print('.', end="") # sys.stdout.write('.')
        
        #Idel, Ion, Ioff = I[:onInd+1], I[onInd:offInd+1], I[offInd:]
        #tdel, ton, toff = t[:onInd+1], t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
        
        RhO.initStates(0)
        #RhO.setLight(phi) # Calculate transition rates for phi
        RhO.updateParams(p)
                
        # # Delay phase
        # RhO.setLight(RhO.phi_0)
        # if RhO.useAnalyticSoln:
            # soln = RhO.calcSoln(tdel, RhO.s_0)
        # else:
            # soln = odeint(RhO.solveStates, RhO.s_0, tdel, Dfun=RhO.jacobian)
        # RhO.storeStates(soln[1:], tdel[1:])
        
        # On phase
        RhO.setLight(phi) # Calculate transition rates for phi
        if fitRates: # Override light-sensitive transition rates
            RhO.updateParams(params)
        RhO.s_on = RhO.states[-1,:] #soln[-1,:]
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(ton, RhO.s_on)
        else:
            soln = odeint(RhO.solveStates, RhO.s_on, ton, Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], ton[1:])
        
        # Off phase
        RhO.setLight(0)
        RhO.s_off = soln[-1,:]
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(toff, RhO.s_off)
        else:
            soln = odeint(RhO.solveStates, RhO.s_off, toff, Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], toff[1:])
        
        return RhO.calcI(V, RhO.states)
    
    
    def errCycle(p,Is,tons,toffs,nfs,RhO,Vs,phis):
        return np.r_[ [Is[i]/nfs[i] - calcCycle(p,tons[i],toffs[i],RhO,Vs[i],phis[i])/nfs[i] for i in range(len(Is))]]
    
    nStates = 6
    
    postOpt = True
    if postOpt: # Relax all parameters and reoptimise
        if verbose > 0:
            print('Performing post-fit optimisation!')
        postParams = Parameters()
        # for p,v in pOffs.items():
            # if p in params:
                # copyParam(p, pOffs, postParams)
                # v.vary = True
        for p in pOns:
            if p in params:
                copyParam(p, pOns, postParams)
                postParams[p].vary = True
        
        RhO = models[str(nStates)]()
        postPmin = minimize(errCycle, postParams, args=(Icycles,tons,toffs,nfs,RhO,Vs,phis), method=method)
        fitParams = postParams
    
    
    for trial in range(nTrials):
        plotFit(Is[trial],ts[trial],onInd,offInd,phis[trial],Vs[trial],nStates,fitParams,fitRates=False,index=trial)
    
    if verbose > 0:
        print("Parameters have been fit for the six-state model")# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
    
    return fitParams #pOns #p4s #RhO







def plotData(Is,ts,t_on,t_off,phis): ##### Replace with onInds and offInds...
    # Plot the photocurrents
    plt.figure()
    for i, phi in enumerate(phis):
        plt.plot(t,Is[i],label='$\phi={:.3g}$'.format(phi))
    plt.legend(loc='best')
    plt.xlabel('Time [ms]')
    plt.ylabel('Photocurrent [nA]')
    plt.axvspan(t_on,t_off,facecolor='y',alpha=0.2)

    

    
    

def fitRecovery(t_peaks, I_peaks, curveFunc, p0, eqString, fig=None):

    
    def expDecay(t, r, Imax): # Restrict so that a = -c to ensure (0,0) is passed through
        return Imax * np.exp(-r*t) - Imax


    def biExpDecay(t, a1, tau1, a2, tau2, I_ss):
        return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2) + I_ss

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


def calcIssfromfV(V,v0,v1,E):#,G): # Added E as another parameter to fit # ==> Fv() := fv()*(V-E)
    ##[s1s, s2s, s3s, s4s, s5s, s6s] = RhO.calcSteadyState(RhO.phiOn)
    ##psi = s3s + (RhO.gam * s4s) # Dimensionless
    
    #E = RhO.E
    if type(V) != np.ndarray:
        V = np.array(V)
    fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
    fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
    ##psi = RhO.calcPsi(RhO.steadyStates) ### This is not necessary for fitting!!!
    ##g_RhO = RhO.gbar * psi * fV # Conductance (pS * mu m^-2)
    ##I_ss = RhO.A * g_RhO * (V - E) # Photocurrent: (pS * mV)
    #I_ss = G * fV * (V-E)
    ##return I_ss * (1e-6) # 10^-12 * 10^-3 * 10^-6 (nA)
    return fV * (V - E)


def fitFV(Vs, Iss, curveFunc, p0, RhO, fig=None):#, eqString): =plt.gcf()
    """F(V):= f(V)(V-E)"""
    # RhO.fV(V) * (V-RhO.E)
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

        
def fitModels(dataSet, nStates=3, params=modelParams['3'], method=methods[3]): #, params, #fit3s=True, fit4s=False, fit6s=False):
    """Routine to fit as many models as possible and select between them according to some parsimony criterion"""
    
    if verbose > 0:
        t0 = time.perf_counter()
        print('Fitting parameters for the {}-state model'.format(nStates))
    ### Check contents of dataSet and produce report on model features which may be fit. 
    # e.g. if not 'rectifier': f(V)=1
    
    ### Trim data slightly to remove artefacts from light on/off transition ramps?
    
    ### Could use precalculated lookup tables to find the values of steady state O1 & O2 occupancies?
    
    
    # Determine the data type passed
    if isinstance(dataSet, PhotoCurrent): # Single photocurrent
        nRuns = 1
        nPhis = 1
        nVs = 1
        setPC = ProtocolData('custom',nRuns,[dataSet.phi],[dataSet.V])
        setPC.trials[0][0][0] = dataSet
        dataSet = {'custom':setPC}
    elif isinstance(dataSet['custom'], ProtocolData): # Set of photocurrents
        setPC = dataSet['custom']
        nRuns = setPC.nRuns
        nPhis = setPC.nPhis
        nVs = setPC.nVs
        #if (setPC.nPhis * setPC.nVs) > 1:
    elif isinstance(dataSet['custom'], PhotoCurrent): # Single photocurrent within ProtocolData
        #targetPC = dataSet['custom']
        nRuns = 1
        nPhis = 1
        nVs = 1
        #setPC = [[[dataSet['custom']]]]
        setPC = ProtocolData('custom',nRuns,[dataSet['custom'].phi],[dataSet['custom'].V])
        setPC.trials[0][0][0] = dataSet['custom']
    else:
        print(type(dataSet['custom']))
        print(dataSet['custom'])
        raise TypeError("dataSet['custom']")
    
    if nRuns == 1:
        run = 0
        
    if nVs == 1:
        vInd = 0
    else:
        try: 
            vInd = setPC.Vs.index(-70)
        except:
            #vInd = 0
            vInd = np.searchsorted(Vs, -70)
    
    if nPhis == 1:
        params['phim'].vary = False
        params['p'].vary = False
        # Fix other model specific parameters?
        if nStates == 4 or nStates == 6:
            params['q'].vary = False
    
    ### Extract the parameters relevant to all models - move inside loop for recovery protocol?
    # if params is not None: # Allow you to pass parameters outside of dataSet
        # params = params
    # elif 'params' in dataSet:
        # params = dataSet['params']
    # else: 
        # params = None
        
    
    # if 'params' in dataSet:
        # if verbose:
            # print('Parameters found in data set - overriding values passed.')
        # for p in dataSet['params']:
            # copyParam(p,dataSet['params'],params)
        
    nonOptParams = ['E', 'Gr0', 'v0', 'v1']
        
    # Now check for the following: E,[phi0,gam,A]; Gr0,gmax,Ipeak,Iss
    
    ### Reversal potential: E
    if 'E' in params:
        E = params['E'].value
    else:
        #global E # Set from global parameters
        E = dataSet['E']
    print('E = {}'.format(E))
    
    ### Peak recovery: Gr, Gr0, Gr_dark, a6 ### This currently fits to the first flux
    if 'Gr' in params:
        Gr0 = params['Gr'].value
    elif 'Gr0' in params:
        Gr0 = params['Gr0'].value
    elif 'Gr_dark' in params:
        Gr0 = params['Gr_dark'].value
    elif 'a6' in params:
        Gr0 = params['a6'].value
    else: ### 2. Fit exponential to peak recovery plots
        if hasattr(dataSet['recovery'], 'tau_r'):
            Gr0 = 1/dataSet['recovery'].tau_r # Gr,dark
        else:
            print("Extract the peaks and fit an exponential...")
            if not (hasattr(dataSet['recovery'], 'tpIPI') and hasattr(dataSet['recovery'], 'IpIPI')):
                # Extract peaks
                dataSet['recovery'].IpIPI = np.zeros(dataSet['recovery'].nRuns)
                dataSet['recovery'].tpIPI = np.zeros(dataSet['recovery'].nRuns)
                for r in range(dataSet['recovery'].nRuns):
                    ### Search only within the on phase of the second pulse
                    I_RhO = dataSet['recovery'].Is[r][0][vInd] # phiInd=0 and vInd=0 Run for each phi and V?
                    startInd = dataSet['recovery'].PulseInds[r][0][vInd][1,0]
                    endInd = dataSet['recovery'].PulseInds[r][0][vInd][1,1]
                    extOrder = int(1+endInd-startInd) #100#int(round(len(I_RhO)/5))
                    #peakInds = findPeaks(I_RhO[:endInd+extOrder+1],minmax,startInd,extOrder)
                    peakInds = findPeaks(I_RhO[:endInd+extOrder+1],startInd,extOrder)
                    if len(peakInds) > 0: # Collect data at the (second) peak
                        dataSet['recovery'].IpIPI[r] = I_RhO[peakInds[0]] #-1 peaks
                        dataSet['recovery'].tpIPI[r] = t[peakInds[0]] # tPeaks
            # Fit exponential
            popt, _, _ = fitPeaks(dataSet['recovery'].tpIPI, dataSet['recovery'].IpIPI, expDecay, p0IPI, '$I_{{peaks}} = {:.3}e^{{-t/{:g}}} {:+.3}$')
            Gr0 = 1/popt[1]
            ### calcGr0()
        params.add('Gr0', value=Gr0, vary=False)
    print('Gr0 = {}'.format(Gr0))
    
    ### Most extreme peak current: Ipmax        
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
        Vsat = setPC.trials[inds[0]][inds[1]][inds[2]].V
        #Ipmax = dataSet['custom'].Ipeak
        #Vsat = dataSet['custom'].V
    print('Ipmax = {}'.format(Ipmax))
    
    ### Maximum conductance: g        
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

    ############################################################# FINISH ME!!! ##########################################################
    ### Optionally fit f(V) (inward rectification) parameters with rectifier data: v0, v1
    # MUST MEASURE E AND FIT AFTER OTHER PARAMETERS
    # if 'v0' in params and 'v1' in params:
        # v0 = params['v0'].value
        # v1 = params['v1'].value
    # else:
    if 'rectifier' in dataSet:
        if hasattr(dataSet['rectifier'], 'Iss'): # Use extracted values
            Iss = dataSet['rectifier'].Iss
            Vs = dataSet['rectifier'].Vs
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
            params['v0'].value = popt[0]
            params['v1'].value = popt[1]
            params['E'].value = popt[2]
            
        # Fit Curve of V vs Iss
    else: # Assume f(V) = (V-E)
        #pass
        #params['v0'].value = 1e12
        #params['v1'].value = 1e12
        # Fix parameters
        params['useIR'].vary = False
        params['v0'].vary = False
        params['v1'].vary = False
    
        
    ### Should f(V) be incorporated into gmax (1.) and Oss (3b.) calculations?
    
    
    #Models = {'3':[[None for v in len(Vs)] for p in len(phis)]}
    
    ### Loop over phi and extrapolate for parameters - skip nRuns
    # for phiInd in range(nPhis):
        # for vInd in range(nVs):
            # targetPC = setPC.trials[0][phiInd][vInd] # Take the first run only
            #< Start fitting here...
    # targetPC = setPC.trials[0][0][0]
    # I = targetPC.I
    # t = targetPC.t
    # onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
    # offInd = targetPC.pulseInds[0,1]
    # V = targetPC.V
    # phi = targetPC.phi
    # Iss = targetPC.Iss # Iplat ############################# Change name to avoid clash with rectifier!    
    ###Iplat is only required for the 3-state fitting procedure

    useNew = True
    
    if useNew:
        from copy import deepcopy
        
        #fitParams = deepcopy(modelParams[str(nStates)])
        fitParams = params ###################################### Revise
        for p in fitParams: # Clear any units
            if fitParams[p].expr is not None:
                fitParams[p].expr=''
        
        fitParams['E'].value = E; fitParams['E'].vary=False
        if nStates == 3:
            fitParams['Gr0'].value = Gr0; fitParams['Gr0'].vary = False # Dark recovery rate
        elif nStates == 4: ### HACK!!!
            fitParams['Gr'].value = Gr0; fitParams['Gr'].vary=False
        elif nStates == 6:
            fitParams['a6'].value = Gr0; fitParams['a6'].vary = False # Dark recovery rate
        if params['useIR'].value == False:
            fitParams['useIR'].vary = False
            fitParams['v0'].vary = False
            fitParams['v1'].vary = False
        #fitParams['g'].value = gmax; fitParams['g'].vary=False
        #fitParams['g'].value = gmax; fitParams['g'].min = gmax*1/3; fitParams['g'].max = gmax*5/3 
        #fitParams['g'].value = gmax*1.25; fitParams['g'].min = gmax*2/3; fitParams['g'].max = gmax*5/3 
        
        if verbose > 0:
            print('\n\n*** Fitting a {}-state model to {} photocurrents ***'.format(nStates,nPhis))
        #phiFits[phiInd] = fitCurve(dataSet,nStates=nStates) #...
        if nStates==3:
            #phiFits[phiInd] = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=pOns,method=method)#,Iss)
            fittedParams = fit3states(setPC,run,vInd,fitParams,method)
        elif nStates==4:
            #phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=pOns,method=method)
            fittedParams = fit4states(setPC,run,vInd,fitParams,method)
        elif nStates==6:
            fittedParams = fit6states(setPC,dataSet,run,vInd,fitParams,method)
            # if 'shortPulse' in dataSet: # Fit Gret
                # from scipy.optimize import curve_fit
                # # Fit tpeak = tpulse + tmaxatp0 * np.exp(-k*tpulse)
                # #dataSet['shortPulse'].getProtPeaks()
                # #tpeaks = dataSet['shortPulse'].IrunPeaks
                # PD = dataSet['shortPulse']
                # PCs = [PD.trials[p][0][0] for p in range(PD.nRuns)]
                # #[pc.alignToTime() for pc in PCs]
                
                # #tpeaks = np.asarray([PD.trials[p][0][0].tpeak for p in range(PD.nRuns)]) # - PD.trials[p][0][0].t[0]
                # #tpulses = np.asarray([PD.trials[p][0][0].onDs[0] for p in range(PD.nRuns)])
                # tpeaks = np.asarray([pc.tpeak for pc in PCs])
                # tpulses = np.asarray([pc.onDs[0] for pc in PCs])
                
                # devFunc = lambda tpulses, t0, k: tpulses + t0 * np.exp(-k*tpulses)
                # p0 = (0,1)
                # popt, pcov = curve_fit(devFunc, tpulses, tpeaks, p0=p0)
                # fig = plt.figure()
                # ax = fig.add_subplot(111, aspect='equal')
                # tsmooth = np.linspace(0,max(tpulses),101)
                # ax.plot(tpulses,tpeaks,'x')
                # ax.plot(tsmooth,devFunc(tsmooth,*popt))
                # ax.plot(tsmooth,tsmooth,'--')
                # ax.set_ylim([0,max(tpulses)+5])
                # ax.set_xlim([0,max(tpulses)+5])
                
                # #plt.tight_layout()
                # #plt.axis('equal')
                
                # print(tpulses)
                # print(tpeaks)
                # print(*popt)
                # print(devFunc(tpeaks,*popt))
                
                # # Solve iteratively Gret = ((tmax*Gd) - np.log(Gd/Gret))/tmax
                
                # def solveGret(tmax, Gd, Gret0=1000, tol=1e-9):
                    # Gret, Gret_m1 = Gret0, 0
                    # print(tmax, Gd, Gret, Gret_m1)
                    # while abs(Gret_m1 - Gret) > tol:
                        # Gret_m1 = Gret
                        # Gret = ((tmax*Gd) - np.log(Gd/Gret_m1))/tmax
                        # #Gret_m1, Gret = Gret, ((tmax*Gd) - np.log(Gd/Gret_m1))/tmax
                        # print(Gret, Gret_m1)
                    # return Gret    
                
                # Gret = solveGret(tmax=popt[0], Gd=0.11, Gret0=500, tol=1e-9)
                # print('Gret = ', Gret)
                
                # fittedParams = {'blah':0.34213}
            
            #Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret)#,Iss)#...
        else:
            raise Exception('Invalid choice for nStates: {}!'.format(nStates))
    
    
    
    
    else:
        simultaneousFit = True
        
        if simultaneousFit:
            ### 3a. Fit biexponential to off curve to find lambdas
            ### Fit off curve - if initial conditions can be calculated, it might be better to use an analytic solution relating directly to model parameters c.f. analysis_4state_off_new.m

            pOffs = Parameters() # Create parameter dictionary
            pOffs.add_many(('Gd1',0.15,True,0.01,1,None), #('Gd1',0.11,False,0.01,1,None)
                            ('Gd2',0.025,True,0.01,1,None), #('Gd2',0.023,False,0.01,1,None)
                            ('e12d',0.01,True,0,1,None),
                            ('e21d',0.01,True,0,1,None))
            
            nTrials = nPhis
            for trial in range(nTrials):
                pOffs.add('a0_'+str(trial), value=0, vary=False)
                pOffs.add('a1_'+str(trial), value=0.1, vary=True, min=0)
                pOffs.add('a2_'+str(trial), value=0.5, vary=True, min=0)
            
            #Ion, Ioff = I[onInd:offInd+1], I[offInd:]
            #ton, toff = t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
            
            Ions = [None for trial in range(nTrials)]
            Ioffs = [None for trial in range(nTrials)]
            tons = [None for trial in range(nTrials)]
            toffs = [None for trial in range(nTrials)]
            
            for phiInd in range(nTrials):
                targetPC = setPC.trials[0][phiInd][0]
                I = targetPC.I
                t = targetPC.t
                onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
                offInd = targetPC.pulseInds[0,1]
                Ions[phiInd] = I[onInd:offInd+1]
                Ioffs[phiInd] = I[offInd:] #[I[offInd:] for I in Is]
                tons[phiInd] = t[onInd:offInd+1]-t[onInd]
                toffs[phiInd] = t[offInd:]-t[offInd] #[t[offInd:]-t[offInd] for t in ts]
            
            # lam1 + lam2 == Gd1 + Gd2 + e12d + e21d
            # lam1 * lam2 == Gd1*Gd2 + Gd1*e21d + Gd2*e12d
                
            calcB = lambda Gd1, Gd2, e12d, e21d: (Gd1 + Gd2 + e12d + e21d)/2
            calcC = lambda b, Gd1, Gd2, e12d, e21d: np.sqrt(b**2 - (Gd1*Gd2 + Gd1*e21d + Gd2*e12d))
            #lam1 = lambda : b - c
            #lam2 = lambda : b + c
            
            #sumLam = pOff['lam1'].value + pOff['lam2'].value
            #prodLam = pOff['lam1'].value * pOff['lam2'].value
            
            def lams(p): #(Gd1, Gd2, e12d, e21d):
                Gd1 = p['Gd1'].value
                Gd2 = p['Gd2'].value
                e12d = p['e12d'].value
                e21d = p['e21d'].value
                b = calcB(Gd1, Gd2, e12d, e21d)
                c = calcC(b, Gd1, Gd2, e12d, e21d)
                return (b-c, b+c)
            
            def fitfunc(p,t,trial):
                
                a0 = p['a0_'+str(trial)].value
                a1 = p['a1_'+str(trial)].value
                a2 = p['a2_'+str(trial)].value
                
                #Gd1 = p['Gd1'].value
                #Gd2 = p['Gd2'].value
                #e12d = p['e12d'].value
                #e21d = p['e21d'].value
                
                return -(a0 + a1*np.exp(-lams(p)[0]*t) + a2*np.exp(-lams(p)[1]*t))

                
            def errfunc(p,Ioffs,toffs):
                return np.r_[ [Ioffs[trial]/Ioffs[trial][-1] - fitfunc(p,toffs[trial],trial)/Ioffs[trial][-1] for trial in range(len(Ioffs))]]
                    
            
            #fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-lams(p)[0]*t) + p['a2'].value*np.exp(-lams(p)[1]*t))
            ##fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-p['lam1'].value*t) + p['a2'].value*np.exp(-p['lam2'].value*t))
            #errfunc = lambda p, Ioff, toff: Ioff - fitfunc(p,toff)
            
            frac = 1
            chop = int(round(len(Ioffs)*frac))
            
            offPmin = minimize(errfunc,pOffs,args=(Ioffs[:chop+1],toffs[:chop+1]),method=method)#, fit_kws={'maxfun':100000})
            print(">>> Off-phase fitting summary: <<<")
            if verbose > 1:
                print(fit_report(pOffs))
                print("Error bars: ",offPmin.errorbars)
            #print('lambda1 = {:.5g}; lambda2 = {:.5g}'.format(pOff['lam1'].value, pOff['lam2'].value))
            print(offPmin.message)
            print("Chi^2: ",offPmin.chisqr)
            
            for trial in range(nTrials):
                #print('a0_'+str(trial),' = ',pOffs['a0_'+str(trial)].value)
                print('a1_'+str(trial),' = ',pOffs['a1_'+str(trial)].value)
                print('a2_'+str(trial),' = ',pOffs['a2_'+str(trial)].value)
            
            
            # import matplotlib as mp
            # figWidth, figHeight = mp.rcParams['figure.figsize']
            # fig = plt.figure(figsize=(figWidth, 2*figHeight))
            fig = plt.figure()
            gs = plt.GridSpec(nTrials,1)
            #axf = axS.twinx()
            #axf.set_yscale('log')
            #axf.set_ylabel('$f\ \mathrm{[Hz]}$')
            for trial in range(nTrials):
                a0 = pOffs['a0_'+str(trial)].value
                a1 = pOffs['a1_'+str(trial)].value
                a2 = pOffs['a2_'+str(trial)].value
                lam1,lam2 = lams(pOffs)
                eq = '-({a0:.3g} + {a1:.3g}*exp(-{lam1:.3g}*t) + {a2:.3g}*exp(-{lam2:.3g}*t))'.format(a0=a0, a1=a1, a2=a2, lam1=lam1, lam2=lam2)
                ax = fig.add_subplot(gs[trial,:])
                ax.plot(toffs[trial], Ioffs[trial], 'g', linewidth=mp.rcParams['lines.linewidth']*3, label='Data') # Experimental data
                ax.plot(toffs[trial], fitfunc(pOffs,toffs[trial],trial), 'b', label=eq) # Fits
                plt.legend(loc=4) # Lower right
                if trial < nTrials-1:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.xlabel('')
                ax.set_ylim(-1,0.1)
                
            print('Gd1 = {}; Gd2 = {}; e12d = {}; e21d = {}'.format(pOffs['Gd1'].value, pOffs['Gd2'].value, pOffs['e12d'].value, pOffs['e21d'].value))
            
            pOffs['Gd1'].vary=False
            pOffs['Gd2'].vary=False
            pOffs['e12d'].vary=False
            pOffs['e21d'].vary=False
            
            pOns = Parameters()
            
            copyParam('Gd1',pOffs,pOns)
            copyParam('Gd2',pOffs,pOns)
            copyParam('e12d',pOffs,pOns)
            copyParam('e21d',pOffs,pOns)
            
            #Gr0,gmax
            pOns.add('g',value=gmax,vary=False)
            pOns.add('Gr',value=Gr0,vary=False)
            pOns.add('phi0',value=5e18,min=1e14,max=1e20,vary=False) ################# Set this to be above the max flux??? 10**ceil(log10(max(phis)))
            
            phiFits = [None for phi in range(nTrials)]
            Is = []
            ts = []
            phis = []
            Vs = []
            
            for phiInd in range(nTrials):
                targetPC = setPC.trials[0][phiInd][0]
                I = targetPC.I
                Is.append(I)
                t = targetPC.t
                ts.append(t)
                onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
                offInd = targetPC.pulseInds[0,1]
                V = targetPC.V
                Vs.append(V)
                phi = targetPC.phi
                phis.append(phi)
            
            if 0:
                #for phiInd, phi in enumerate(phis):
                for phiInd in range(nTrials):
                    # targetPC = setPC.trials[0][phiInd][0]
                    # I = targetPC.I
                    # t = targetPC.t
                    # onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
                    # offInd = targetPC.pulseInds[0,1]
                    # V = targetPC.V
                    # Vs.append(V)
                    # phi = targetPC.phi
                    # phis.append(phi)
                    
                    if verbose > 0:
                        print('\n\n*** Fitting phi: {:.3g} photocurrent ***'.format(phi))
                    #phiFits[phiInd] = fitCurve(dataSet,nStates=nStates) #...
                    if nStates==3:
                        phiFits[phiInd] = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=pOns,method=method)#,Iss)
                    elif nStates==4:
                        phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=pOns,method=method)
                    elif nStates==6:
                        Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret)#,Iss)#...
                    else:
                        raise Exception('Invalid choice for nStates: {}!'.format(nStates))
                    
            else:
                
                if nStates==4:
                    pOns.add('k1',value=3, min=0.01) #0.5 Ga1 = k1 * phi
                    pOns.add('k2',value=1.5, min=0.01) #0.2 Ga2 = k2 * phi
                    pOns.add('c1',value=0.05, min=0.01)
                    pOns.add('c2',value=0.01, min=0.01)
                    pOns.add('gam',value=0.05, min=0, max=1)
                    # Place RhO,V,phi in onPmin?
                    if 1:
                        pOns.add('p',value=0.7,min=0.1,max=5)
                        pOns.add('phim',value=1e17,min=1e15,max=1e19) #1e19
                    
                    RhO = models['4']()
                    
                    # Normalise? e.g. /Ions[trial][-1] or /min(Ions[trial])
                    def errfunc4On(p,Ions,tons,RhO,Vs,phis):
                        return np.r_[ [Ions[trial]/Ions[trial][-1] - calc4on(p,tons[trial],RhO,Vs[trial],phis[trial])/Ions[trial][-1] for trial in range(len(Ions))]]
                    
                    ### Trim down ton? Take 10% of data or one point every ms? ==> [0::5]
                    ### Instead try looping over a coarse grid of parameter values and saving RMSE for each combination c.f. analysis_4state_on_new.m
                    
                    if verbose > 1:
                        print('Optimising',end='')
                    onPmin = minimize(errfunc4On,pOns,args=(Ions,tons,RhO,Vs,phis),method=method)
                
                print("\n>>> On-phase fitting summary: <<<")
                if verbose > 1:
                    print(fit_report(pOns))
                    print("Error bars: ",onPmin.errorbars)
                if not onPmin.success:
                    print(onPmin.success)
                    print(onPmin.message)
                    print(onPmin.lmdif_message)
                print("Chi^2: ",onPmin.chisqr)
                
                print('k1 = {}; k2 = {}; c1 = {}; c2 = {}'.format(pOns['k1'].value, pOns['k2'].value, pOns['c1'].value, pOns['c2'].value))
                print('gam = ',pOns['gam'].value)
                print('phi0 = ',pOns['phi0'].value)
                print('phim = ',pOns['phim'].value)
                print('p = ',pOns['p'].value)

                for trial in range(nTrials):
                    plotFit(Is[trial],ts[trial],onInd,offInd,phis[trial],Vs[trial],nStates,pOns,fitRates=False,index=trial)
            
        else:
        
            nTrials = nPhis #len(phis)
            phiFits = [None for phi in range(nTrials)]
            phis = []
            #for phiInd, phi in enumerate(phis):
            for phiInd in range(nTrials):
                targetPC = setPC.trials[0][phiInd][0]
                I = targetPC.I
                t = targetPC.t
                onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
                offInd = targetPC.pulseInds[0,1]
                V = targetPC.V
                phi = targetPC.phi
                phis.append(phi)
                
                if verbose > 0:
                    print('\n\n*** Fitting phi: {:.3g} photocurrent ***'.format(phi))
                #phiFits[phiInd] = fitCurve(dataSet,nStates=nStates) #...
                if nStates==3:
                    phiFits[phiInd] = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=None,method=method)#,Iss)
                elif nStates==4:
                    phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=None,method=method)
                elif nStates==6:
                    Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret)#,Iss)#...
                else:
                    raise Exception('Invalid choice for nStates: {}!'.format(nStates))
            
            
            
            ### Fix light insensitive parameters
            # 3: Gd
            # 4: Gd1, Gd2, [Gr]
            # 6: a2, a4, b1, b3 [a6]
            paramsSet = Parameters()
            
            aggregates = aggregateFits(phis,phiFits,nStates)
            
            if nStates==3: #fit3s:    
                paramsSet.add('Gd', value=np.mean(aggregates['Gd']),vary=False)
            elif nStates==4: #fit4s:
                #aggregates = aggregate4sFits(phis,phiFits)
                paramsSet.add('Gd1',value=np.mean(aggregates['Gd1']),vary=False)
                paramsSet.add('Gd2',value=np.mean(aggregates['Gd2']),vary=False)
            elif nStates==6: #fit6s:
                pass
            else:
                raise Exception('Invalid choice for nStates: {}!'.format(nStates))
            
            ### Reoptimise
                
            for phiInd, phi in enumerate(phis):
                targetPC = setPC.trials[0][phiInd][0]
                I = targetPC.I
                t = targetPC.t
                onInd = targetPC.pulseInds[0,0] ### Consider multiple pulse scenarios
                offInd = targetPC.pulseInds[0,1]
                V = targetPC.V
                phi = targetPC.phi
                
                if verbose > 0:
                    print('\n*** Refitting phi: {:.3g} photocurrent ***'.format(phi))
                #phiFits[phiInd] = fitCurve(dataSet,nStates=nStates,params=paramsSet) #...
                if nStates==3:
                    phiFits[phiInd] = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,Iss,paramsSet,method=method)
                elif nStates==4:
                    phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,paramsSet,method=method)
                elif nStates==6:
                    Models['6'] = fit6states(I,t,onInd,offInd,phi,V,Gr0,gbar,Gret,paramsSet)#,Iss)#...
                else:
                    raise Exception('Invalid choice for nStates: {}!'.format(nStates))

            aggregates = aggregateFits(phis,phiFits,nStates)
    
    ### Fit functions to flux-dependent rates
    #########################################
    
    
    
    #RhO = models[nStates]()
    #RhO.setParams(fitParams)
    
    # Compare fit vs computational complexity...
    # RhO = selectModel(nStates)
    # Calculate chisqr for each model and select between them. 
    ###RhO = RhO3
    
    # for p,v in fittedParams.items():
        # v.vary = False
    
    # Run small signal analysis
    #runSSA(RhO)
    #characterise(RhO)    
    #import os
    #import pickle
    exportName = 'fitted{}sParams.pkl'.format(nStates)
    fh = open(os.path.join(dDir, exportName), "wb")
    pickle.dump(fittedParams, fh)
    fh.close()
    
    if verbose > 0:
        print('\n*** Parameters fit for the {}-state model in {:.3g}s ***\n'.format(nStates, time.perf_counter() - t0))
    
    return fittedParams #aggregates #fitCurve(dataSet, highestState)
    
    

    
    