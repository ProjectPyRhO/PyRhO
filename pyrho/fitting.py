from __future__ import print_function

#__all__ = ['fitModels', 'copyParam', 'getRecoveryPeaks', 'fitRecovery', 'fitfV']

from lmfit import minimize, Parameters, fit_report
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyrho.parameters import *
from pyrho.loadData import *
from pyrho.utilities import * # plotLight, round_sig, findPeaks, findPlateauCurrent
from pyrho.models import * # for fitPeaks
from pyrho.config import * #verbose, saveFigFormat, addTitles, fDir, dDir, eqSize
from pyrho import config
import os
import pickle
from copy import deepcopy
import warnings
from scipy.optimize import curve_fit


methods = ('leastsq', 'nelder', 'lbfgsb', 'powell', 'cg', 'cobyla', 'tnc', 'slsqp', 'differential_evolution')
defMethod = methods[3]

#methods=('leastsq', 'nelder', 'lbfgsb', 'powell', 'cg', 'newton', 'cobyla', 'tnc', 'trust-ncg', 'dogleg', 'slsqp', 'differential_evolution')
#'newton', 'trust-ncg', 'dogleg' : Require Jacobian

# Notes
# =====
# Empirically, 'powell', 'lbfgsb' and 'nelder' typically provide the best results
# 'newton', 'trust-ncg' and 'dogleg' require a Jacobian 
# http://scipy-lectures.github.io/advanced/mathematical_optimization/#choosing-a-method


# Use Akaike and Bayesian Information Criteria to compare model fits

    
##### Development ideas #####

def run(RhO, t):
    # Ga = RhO.Ga; Gd = RhO.Gd; Gr = RhO.Gr
    # if not RhO.useAnalyticSoln or 2*(Ga*Gd + Ga*Gr + Gd*Gr) > (Ga**2 + Gd**2 + Gr**2):
        # soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
    # else:
        # soln = RhO.calcSoln(t, RhO.states[-1,:])
        
    try:
        soln = RhO.calcSoln(t, RhO.states[-1,:])
    except: # Any exception e.g. NotImplementedError or ValueError
        soln = odeint(RhO.solveStates, RhO.states[-1,:], t, Dfun=RhO.jacobian)
    return soln


def plotData(Is, ts, t_on, t_off, phis): 
    # Replace with onInds and offInds...
    # Plot the photocurrents
    plt.figure()
    for i, phi in enumerate(phis):
        plt.plot(ts[i], Is[i], label='$\phi={:.3g}$'.format(phi))
    plt.legend(loc='best')
    plt.xlabel('Time [ms]')
    plt.ylabel('Photocurrent [nA]')
    plt.axvspan(t_on, t_off, facecolor='y', alpha=0.2)

    
# func(params, *args, **kws)

def errPhase(p,residual,Is,ts,RhO,Vs,phis):
    """     
    IN:
        p           : Model parameters
        residual    : Function to calculate the residual
        Is          : Set (list) of currents
        ts          : Set (list) of corresponding time arrays
        RhO         : Model object
        Vs          : List of Voltage clamp potentials
        phis        : List of flux values
    
    OUT:
        Concatenated set of residual values
    """
    #data = zip(Is, ts, nfs, Vs, phis) # Normalise? e.g. /Ions[trial][-1] or /min(Ions[trial])
    
    return np.r_[ [(Is[i] - residual(p,ts[i],RhO,Vs[i],phis[i]))/Is[i][-1] for i in range(len(Is))]]

'''    
def errSetOnPhase(p,Ions,tons,RhO,Vs,phis):
    return np.r_[ [Ions[i]/Ions[i][-1] - calcOnPhase(p,tons[i],RhO,Vs[i],phis[i])/Ions[i][-1] for i in range(len(Ions))]]
'''

    
##### Main fitting routines #####
    
def calcOnPhase(p,t,RhO,V,phi):
    """Simulate the on-phase from base parameters"""

    RhO.initStates(0)
    RhO.updateParams(p)
    RhO.setLight(phi) # Calculate transition rates for phi

    if verbose > 2:
        print('.', end="") # sys.stdout.write('.')    
        soln, out = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian, full_output=True)
        if out['message'] != 'Integration successful.':
            #print(out)
            print(RhO.reportParams())
    else:
        soln = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
        
    I_RhO = RhO.calcI(V, soln)
    return I_RhO


# Normalise? e.g. /Ions[trial][-1] or /min(Ions[trial])
def errOnPhase(p,Ions,tons,RhO,Vs,phis):
    return np.r_[ [(Ions[i] - calcOnPhase(p,tons[i],RhO,Vs[i],phis[i]))/Ions[i][-1] for i in range(len(Ions))]]

    
def reportFit(minResult, description, method):
    
    '''
    N           := number of data points
    N_{vars}    := number of variables
    \chi^2      := \Sum_i^N [Resid_i]^2
    \chi_v^2    := \chi^2 / (N - N_{vars})  # Reduced \chi^2
    AIC         := N ln(\chi^2 / N) + 2N_{vars}
    BIC         := N ln(\chi^2 / N) + ln(N) \cdot N_{vars}
    '''
    
    #Fitting parameters for the {}-state model
    print("\n--------------------------------------------------------------------------------")
    print("{} with the '{}' algorithm... ".format(description, method))
    print("--------------------------------------------------------------------------------\n")
    if hasattr(minResult, 'message'):
        print(minResult.message)
    
    if verbose > 1:
        print(fit_report(minResult))
        if verbose > 2:
            print(minResult.covar)
            print("Error bars: ", minResult.errorbars)
        
        if not minResult.success:
            print("Success: ", minResult.success)
            if method == 'leastsq':
                print("Integer error: ", minResult.ier)
                print(minResult.lmdif_message)
    else:
        print("Fit for {} variables over {} points ({} d.f.) with {} function evaluations".format(minResult.nvarys, minResult.ndata, minResult.nfree, minResult.nfev))
        print("Chi^2 (reduced): {}, ({})".format(minResult.chisqr, minResult.redchi))
        print("Akaike Info.:   {} \nBayesian Info.: {}".format(minResult.aic, minResult.bic))
        #print("Chi^2 \t rChi^2 \t AIC \t BIC")
        #print("{} \t {} \t {} \t {}".format(minResult.chisqr, minResult.redchi, minResult.aic, minResult.bic))
    
    
def copyParam(name, source, target):
    if not isinstance(name, (list, tuple)):
        names = [name]
    else:
        names = name
    for name in names:
        if name not in target:
            target.add(name, value=source[name].value, vary=source[name].vary, min=source[name].min, max=source[name].max, expr=source[name].expr)
        else:
            target[name].set(value=source[name].value, vary=source[name].vary, min=source[name].min, max=source[name].max, expr=source[name].expr)
    return #target
    

def plotOffPhaseFits(toffs, Ioffs, pOffs, phis, nStates, fitFunc, Exp1, Exp2, Gd=None):
    fig = plt.figure()
    #gs = plt.GridSpec(nTrials,1)
    ax = fig.add_subplot(111)
    lw = 2.5 * mpl.rcParams['lines.linewidth']
    
    nTrials = len(Ioffs)
    assert(len(toffs) == nTrials)
    assert(len(phis) == nTrials)
    
    for trial in range(nTrials):
        Islow = pOffs['Islow_'+str(trial)].value
        Ifast = pOffs['Ifast_'+str(trial)].value
        #ax = fig.add_subplot(gs[trial,:])
        
        ax.plot(toffs[trial], Ioffs[trial], color=colours[trial%len(colours)], linewidth=lw, markeredgecolor='None', label='Data: phi={phi:.3g}'.format(phi=phis[trial]))
        #ax.plot(toffs[trial], Ioffs[trial], 'g', linewidth=mpl.rcParams['lines.linewidth']*3, label='Data: phi={phi:.3g}'.format(phi=phis[trial])) # Experimental data
        
        eq = 'I(t)={Islow:.3g}*exp(-{Exp1:.3g}*t) + {Ifast:.3g}*exp(-{Exp2:.3g}*t)'.format(Islow=Islow, Ifast=Ifast, Exp1=Exp1, Exp2=Exp2)
        ax.plot(toffs[trial], fitFunc(pOffs,toffs[trial],trial), color='k', linestyle='--', label=eq) # Fits
        #ax.plot(toffs[trial], fit3off(pOffs,toffs[trial],trial), 'b', label=eq) # Fits
        
        if Gd is not None: # Plot single exponential decay too
            ax.plot(toffs[trial], (Islow+Ifast)*np.exp(-Gd*toffs[trial]), color=colours[trial%len(colours)], linestyle=':', label='I(t)={I0:.3g}*exp(-{Gd:.3g}*t)'.format(I0=Islow+Ifast, Gd=Gd)) # Removed - coefficients
            #ax.plot(toffs[trial], (Islow+Ifast)*np.exp(-Gd*toffs[trial]), 'r', label='I(t)={I0:.3g}*exp(-{Gd:.3g}*t)'.format(I0=Islow+Ifast, Gd=Gd)) # Removed - coefficients
        
        # if trial < nTrials-1:
            # plt.setp(ax.get_xticklabels(), visible=False)
            # plt.xlabel('')
        #ax.set_ylim(-1,0.1) ### Reconsider!!!
    
    plt.legend(loc='best') #loc=4 Lower right
    plt.xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right')
    plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    
    #ax.spines['left'].set_position('zero') # y-axis
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero') # x-axis
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    plt.tight_layout()
    
    fig.savefig(fDir+'OffPhaseFits'+str(nStates)+'states'+'.'+config.saveFigFormat, format=config.saveFigFormat)
    
    
def plotFit(PC, nStates, params, fitRates=False, index=None):
    
    RhO = models[str(nStates)]()
    RhO.updateParams(params)
    RhO.phiFit = PC.phi             # Flux intensity at which the parameters were fit. 

    phi = PC.phi
    V = PC.V
    ### Plot experimental curve
    #totT = t[-1] - t[0] #max(t)
    begT, endT = PC.begT, PC.endT #PC.t[0], PC.t[-1]
    I = PC.I
    t = PC.t
    
    Ifig = plt.figure()
    gsPL = plt.GridSpec(4,1)

    axFit = Ifig.add_subplot(gsPL[:-1,:])
    plotLight(PC.pulses, axFit)
    axFit.set_xlim((begT, endT))
    plt.setp(axFit.get_xticklabels(), visible=False)
    axFit.set_ylabel('$\mathrm{Photocurrent\ [nA]}$')
    axFit.plot(t, I, color='g', label='$\mathrm{Experimental\ Data}$')
    
    #axFit.spines['left'].set_position('zero') # y-axis
    axFit.spines['right'].set_color('none')
    axFit.spines['bottom'].set_position('zero') # x-axis
    axFit.spines['top'].set_color('none')
    axFit.spines['left'].set_smart_bounds(True)
    axFit.spines['bottom'].set_smart_bounds(True)
    axFit.xaxis.set_ticks_position('bottom')
    axFit.yaxis.set_ticks_position('left')
    
    
    ### Plot model-generated curve
    #onInd, offInd = PC.pulseInds[0,0], PC.pulseInds[0,1]
    #Idel, Ion, Ioff = I[:onInd+1], I[onInd:offInd+1], I[offInd:]
    #tdel, ton, toff = t[:onInd+1], t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    
    #TODO: Refactor to use calcCycle, runTrial or similar
    
    Idel, tdel  = PC.getDelayPhase()#;   tdel -= tdel[0]
    
    ## Delay phase
    RhO.setLight(RhO.phi_0)            
    if RhO.useAnalyticSoln:
        soln = RhO.calcSoln(tdel, RhO.s_0)
    else:
        soln = odeint(RhO.solveStates, RhO.s_0, tdel, Dfun=RhO.jacobian)
    RhO.storeStates(soln[1:], tdel[1:])
    
    
    for p in range(PC.nPulses):
        
        Ion, ton    = PC.getOnPhase(p)#;     ton -= ton[0]
        Ioff, toff  = PC.getOffPhase(p)#;    toff -= toff[0]
        
        ## On phase
        RhO.setLight(phi) # Calculate transition rates for phi
        if fitRates: # Override light-sensitive transition rates
            RhO.updateParams(params)
        
        RhO.s_on = soln[-1,:]
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(ton, RhO.s_on)
        else:
            soln = odeint(RhO.solveStates, RhO.s_on, ton, Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], ton[1:])
    
        if config.addTitles and p == 0:
            #TODO: Generate title with RhO.paramsList, RhO.photoRates, RhO.constRates
            if nStates == 3:
                plt.title('Three-state model fit to data (phi={:.3g}) [Ga={:.3g}; Gd={:.3g}; Gr={:.3g}] \n[k_a={:.3g}; p={:.3g}; k_r={:.3g}; q={:.3g}; phi_m={:.3g}; Gd={:.3g}; Gr0={:.3g}]'.format(phi,RhO.Ga,RhO.Gd,RhO.Gr,RhO.k_a,RhO.p,RhO.k_r,RhO.q,RhO.phi_m,RhO.Gd,RhO.Gr0))
            elif nStates == 4:
                plt.title('Four-state model fit to data (phi={:.3g}) \n[Ga1={:.3g}; Ga2={:.3g}; Gf={:.3g}; Gb={:.3g}; Gd1={:.3g}; Gd2={:.3g}]'.format(phi,RhO.Ga1,RhO.Ga2,RhO.Gf,RhO.Gb,RhO.Gd1,RhO.Gd2))
            elif nStates == 6:
                plt.title('Six-state model fit to data (phi={:.3g}) \n[Ga1={:.3g}; Ga2={:.3g}; Gf={:.3g}; Gb={:.3g}; Go1={:.3g}; Go2={:.3g}; Gd1={:.3g}; Gd2={:.3g}]'.format(phi,RhO.Ga1,RhO.Ga2,RhO.Gf,RhO.Gb,RhO.Go1,RhO.Go2,RhO.Gd1,RhO.Gd2))
        
        ## Off phase
        RhO.setLight(0)
        #if fitRates: # Override light-sensitive transition rates
        #    RhO.updateParams(params)
        
        RhO.s_off = soln[-1,:]
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(toff, RhO.s_off)
        else:
            soln = odeint(RhO.solveStates, RhO.s_off, toff, Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], toff[1:])
    
    Ifit = RhO.calcI(V, RhO.states)
    
    # Plot model fit curve
    axFit.plot(t, Ifit, color='b', label='$\mathrm{{Model\ fit\ ({}-states)}}$'.format(nStates)) #t[onInd:]
    axFit.legend(loc='best')
    
    ### Plot Residuals
    axRes = Ifig.add_subplot(gsPL[-1,:], sharex=axFit)
    plotLight(PC.pulses, axRes) #plotLight(np.asarray([[t[onInd],t[offInd]]]), axRes) #plt.axvspan(t[onInd],t[offInd],facecolor='y',alpha=0.2)
    
    # axLag.set_aspect('auto')
    #if minObj is not None:
    #    tcycle = PC.getCycle()[1]
    #    axRes.plot(tcycle, minObj.residual)
    #else:
    #plt.plot(t[onInd:],I[onInd:]-Ifit)
    PCofSS = bool(PC.ss_ is not None)
    if PCofSS:
        axRes.plot(t, (I-Ifit)*100/abs(PC.ss_))
        axRes.set_ylabel('$\mathrm{Error\ (\%\ I_{ss})}$')# relative error')
    else:
        axRes.plot(t, I-Ifit)
        axRes.set_ylabel('$\mathrm{Residuals}$')
    #plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff))
    #plt.plot(t[onInd:],np.append(Ion[:-1]-IfitOn[:-1],Ioff-IfitOff)*100/I[onInd:]) # Error relative to experimental curve
    #plt.plot(toff+t[offInd],Ioff-IfitOff)
    
    axRes.set_xlabel('$\mathrm{Time\ [ms]}$') #, position=(1,0), ha='right')
    
    #plt.setp(axRes.get_xticklabels(), visible=False)
    #plt.xlim((0,totT))
    
    plt.axhline(y=0, linestyle=':', color='k')
    #axRes.spines['left'].set_position('zero') # y-axis
    # axRes.spines['right'].set_color('none')
    # axRes.spines['bottom'].set_position('zero') # x-axis
    # axRes.spines['top'].set_color('none')
    # axRes.spines['left'].set_smart_bounds(True)
    # axRes.spines['bottom'].set_smart_bounds(True)
    # axRes.xaxis.set_ticks_position('bottom')
    # axRes.yaxis.set_ticks_position('left')
    
    plt.tight_layout()
    
    if index is None:
        index = ''
    Ifig.savefig(fDir+'fit'+str(nStates)+'states'+str(index)+'.'+config.saveFigFormat, format=config.saveFigFormat)

    if verbose > 1:
        print("Fit has been plotted for the {}-state model".format(nStates))# at a flux of {} [photons * s^-1 * mm^-2]".format(phi))
    
    return



def fit3states(fluxSet, run, vInd, params, method=defMethod, verbose=verbose): 
    """
    fluxSet := ProtocolData set (of Photocurrent objects) to fit
    run     := Index for the run within the ProtocolData set
    vInd    := Index for Voltage clamp value within the ProtocolData set
    params  := Parameters object of model parameters with initial values [and bounds, expressions]
    method  := Fitting algorithm for the optimiser to use
    """
    
    plotResult = bool(verbose > 1)
    
    nStates = 3
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    assert(0 < nPhis)
    assert(0 <= run < nRuns)
    assert(0 <= vInd < nVs)
    
    Ions = [None for phiInd in range(nPhis)]
    Ioffs = [None for phiInd in range(nPhis)]
    tons = [None for phiInd in range(nPhis)]
    toffs = [None for phiInd in range(nPhis)]
    phis = []
    Is = []
    ts = []
    Vs = []
    
    Icycles = []
    nfs = []
    
    #offTally = 0
    #soffs = [None for phiInd in range(nPhis)]
    #onTally = 0
    #sons = [None for phiInd in range(nPhis)]
    
    # Trim off phase data
    #frac = 1
    #chop = int(round(len(Ioffs[0])*frac))
    
    for phiInd in range(nPhis):
        targetPC = fluxSet.trials[run][phiInd][vInd]
        #targetPC.alignToTime()
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
        
        # if phiInd < nPhis-1:
            # soffs[phiInd] = slice(offTally, len(Ioffs[phiInd])+offTally)
            # sons[phiInd] = slice(onTally, len(Ions[phiInd])+onTally)
        # else:
            # soffs[phiInd] = slice(offTally, None)
            # sons[phiInd] = slice(onTally, None)
        # offTally += len(Ioffs[phiInd])
        # onTally += len(Ions[phiInd])
        
    
    nTrials = nPhis
    
    ### 3a. Fit exponential to off curve to find Gd
    ### Fit off curve
    iOffPs = Parameters() # Create parameter dictionary
    
    ### Original single exponential fit
    #pOffs.add('A', value=-Ioff[0])#vary=False) #-1 # Constrain fitting to start at the beginning of the experimental off curve
    #pOffs.add('Gd', value=0.1, min=0)
    #offPmin = minimize(resid3off,pOffs,args=(Ioff,toff),method=meth)
    
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        Iss = Ioffs[phiInd][0]
        if Iss < 0: # Excitatory
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, max=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, max=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
        else:
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, min=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, min=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
    
    iOffPs.add('Gd1', value=params['Gd'].value/5, min=params['Gd'].min, max=params['Gd'].max)
    iOffPs.add('Gd2', value=params['Gd'].value*5, min=params['Gd'].min, max=params['Gd'].max)
    
    def fit3off(p,t,trial):
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        Gd1 = p['Gd1'].value
        Gd2 = p['Gd2'].value
        return Islow * np.exp(-Gd1*t) + Ifast * np.exp(-Gd2*t)

    #err3off = lambda p,I,t: I - fit3off(p,t)
    def err3off(p,Ioffs,toffs):
        return np.r_[ [(Ioffs[i] - fit3off(p,toffs[i],i))/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    offPmin = minimize(err3off, iOffPs, args=(Ioffs,toffs), method=method)
    pOffs = offPmin.params
    
    #def err3off(p,Ioffs,toffs,soffs):
    #    return np.concatenate( [(Ioffs[s] - fit3off(p,toffs[s],i))/Ioffs[s][0] for i,s in enumerate(soffs)] )
    
    #offPmin = minimize(err3off, pOffs, args=(np.concatenate(Ioffs),np.concatenate(toffs),soffs), method=method)
    
    #pOffs.add('Gd', value=max(pOffs['Gd1'].value,pOffs['Gd2'].value), min=0)
    #nf = pOffs['A'].value + pOffs['B'].value
    #pOffs.add('Gd', value=pOffs['A'].value*pOffs['Gd1'].value/nf+pOffs['B'].value*pOffs['Gd2'].value/nf, min=0)
    
    Gds = [None for phiInd in range(nPhis)]
    for phiInd in range(nPhis):
        v = pOffs.valuesdict()
        Islow = v['Islow_'+str(phiInd)]
        Ifast = v['Ifast_'+str(phiInd)]
        Gds[phiInd] = (Islow*v['Gd1'] + Ifast*v['Gd2'])/(Islow + Ifast)
    
    Gd = np.mean(Gds)
    
    if plotResult:
        plotOffPhaseFits(toffs, Ioffs, pOffs, phis, nStates, fit3off, v['Gd1'], v['Gd2'], Gd=Gd)
        
    reportFit(offPmin, "Off-phase fit report for the 3-state model", method)
    print('Gd1 = {}; Gd2 = {} ==> Gd = {}'.format(pOffs['Gd1'].value, pOffs['Gd2'].value, Gd))
    
    
    ### Fit on curve
    
    iOnPs = Parameters()
    for p in ['Gd', 'k_a', 'k_r', 'p', 'q', 'phi_m', 'g0', 'Gr0', 'E', 'v0', 'v1']:
        copyParam(p, params, iOnPs)
    iOnPs['Gd'].set(value=Gd, vary=False)
    
    RhO = models['3']()
    
    def fit3on(p,t,RhO,phi,V):
        RhO.updateParams(p)
        RhO.setLight(phi)
        states = RhO.calcSoln(t, s0=[1,0,0]) # t starts at 0
        #states = odeint(RhO.solveStates, RhO.s_0, t, Dfun=RhO.jacobian)
        return RhO.calcI(V, states)
        
    def err3on(p,Ions,tons,RhO,phis,Vs):
        return np.r_[ [(Ions[i] - fit3on(p,tons[i],RhO,phis[i],Vs[i]))/Ions[i][-1] for i in range(len(Ions))] ]
    
    onPmin = minimize(err3on, iOnPs, args=(Ions,tons,RhO,phis,Vs), method=method)
    pOns = onPmin.params
    
    # if (Gr + Gd - 2*np.sqrt(Gr*Gd)) < Ga < (Gr + Gd + 2*np.sqrt(Gr*Gd)):
        # print('\n\nWarning! No real solution exists!\n\n')
    
    reportFit(onPmin, "On-phase fit report for the 3-state model", method)
    if verbose > 0:
        print('k_a = {}; p = {}; k_r = {}; q = {}; phi_m = {}'.format(pOns['k_a'].value, pOns['p'].value, 
                                                pOns['k_r'].value, pOns['q'].value, pOns['phi_m'].value))
    
    fitParams = pOns
    
    return fitParams



def fit4states(fluxSet, run, vInd, params, method=defMethod, verbose=verbose): 
    """
    fluxSet := ProtocolData set (of Photocurrent objects) to fit
    run     := Index for the run within the ProtocolData set
    vInd    := Index for Voltage clamp value within the ProtocolData set
    params  := Parameters object of model parameters with initial values [and bounds, expressions]
    method  := Fitting algorithm for the optimiser to use
    verbose := Text output (verbosity) level
    """
    
    plotResult = bool(verbose > 1)
    
    nStates = 4
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    assert(0 < nPhis)
    assert(0 <= run < nRuns)
    assert(0 <= vInd < nVs)
    
    Ions = [None for phiInd in range(nPhis)]
    Ioffs = [None for phiInd in range(nPhis)]
    tons = [None for phiInd in range(nPhis)]
    toffs = [None for phiInd in range(nPhis)]
    phis = []
    Is = []
    ts = []
    Vs = []
    
    Icycles = []
    nfs = []
    
    
    # Trim off phase data
    #frac = 1
    #chop = int(round(len(Ioffs[0])*frac))
    
    for phiInd in range(nPhis):
        targetPC = fluxSet.trials[run][phiInd][vInd]
        #targetPC.alignToTime()
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
    
    OffKeys = ['Gd1', 'Gd2', 'Gf0', 'Gb0']
    
    iOffPs = Parameters() # Create parameter dictionary
    for k in OffKeys:
        copyParam(k, params, iOffPs)
    
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        Iss = Ioffs[phiInd][0]
        if Iss < 0:
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, max=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, max=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
        else:
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, min=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, min=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
        
    # lam1 + lam2 == Gd1 + Gd2 + Gf0 + Gb0
    # lam1 * lam2 == Gd1*Gd2 + Gd1*Gb0 + Gd2*Gf0
    
    calcB = lambda Gd1, Gd2, Gf0, Gb0: (Gd1 + Gd2 + Gf0 + Gb0)/2
    calcC = lambda b, Gd1, Gd2, Gf0, Gb0: np.sqrt(b**2 - (Gd1*Gd2 + Gd1*Gb0 + Gd2*Gf0))
    
    def lams(p):
        Gd1 = p['Gd1'].value
        Gd2 = p['Gd2'].value
        Gf0 = p['Gf0'].value
        Gb0 = p['Gb0'].value
        b = calcB(Gd1, Gd2, Gf0, Gb0)
        c = calcC(b, Gd1, Gd2, Gf0, Gb0)
        return (b-c, b+c)
        
    def fit4off(p,t,trial):
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        lam1, lam2 = lams(p)
        return Islow*np.exp(-lam1*t) + Ifast*np.exp(-lam2*t)

    def err4off(p,Ioffs,toffs):
        """Normalise by the first element of the off-curve""" # [-1]
        return np.r_[ [(Ioffs[i] - fit4off(p,toffs[i],i))/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    #fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-lams(p)[0]*t) + p['a2'].value*np.exp(-lams(p)[1]*t))
    ##fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-p['lam1'].value*t) + p['a2'].value*np.exp(-p['lam2'].value*t))
    #errfunc = lambda p, Ioff, toff: Ioff - fitfunc(p,toff)
    
    offPmin = minimize(err4off, iOffPs, args=(Ioffs,toffs), method=method)#, fit_kws={'maxfun':100000})
    pOffs = offPmin.params
    
    reportFit(offPmin, "Off-phase fit report for the 4-state model", method)
    if verbose > 0:
        vd = pOffs.valuesdict()
        print('Gd1 = {Gd1}; Gd2 = {Gd2}; Gf0 = {Gf0}; Gb0 = {Gb0}'.format(**vd))
    
    if plotResult:
        lam1, lam2 = lams(pOffs)
        plotOffPhaseFits(toffs, Ioffs, pOffs, phis, nStates, fit4off, lam1, lam2, Gd=None)
        
    for k in OffKeys:
        pOffs[k].vary = False
    
    
    ### ON PHASE
    
    iOnPs = Parameters() # deepcopy(params)
    # Set parameters from Off-curve optimisation
    for k in OffKeys:
        copyParam(k, pOffs, iOnPs)
    
    for k in ['k1', 'k2', 'k_f', 'k_b', 'gam', 'p', 'q', 'phi_m', 'g0', 'Gr0', 'E', 'v0', 'v1']: #.extend(OffKeys):
        copyParam(k, params, iOnPs)
    
    RhO = models['4']()
    
    
    ### Trim down ton? Take 10% of data or one point every ms? ==> [0::5]
        
    if verbose > 2:
        print('Optimising ',end='')
    
    onPmin = minimize(errOnPhase, iOnPs, args=(Ions,tons,RhO,Vs,phis), method=method)
    pOns = onPmin.params
    
    reportFit(onPmin, "On-phase fit report for the 4-state model", method)
    if verbose > 0:
        print('k1 = {}; k2 = {}; k_f = {}; k_b = {}'.format(pOns['k1'].value, pOns['k2'].value, pOns['k_f'].value, pOns['k_b'].value))
        print('gam = {}; phi_m = {}; p = {}; q = {}'.format(pOns['gam'].value, pOns['phi_m'].value, pOns['p'].value, pOns['q'].value))
    
    fitParams = pOns
    
    return fitParams



def fit6states(fluxSet, quickSet, run, vInd, params, method=defMethod, verbose=verbose): 
    """
    fluxSet := ProtocolData set (of Photocurrent objects) to fit
    quickSet:= ProtocolData set (of Photocurrent objects) with short pulses to fit opsin activation rates
    run     := Index for the run within the ProtocolData set
    vInd    := Index for Voltage clamp value within the ProtocolData set
    params  := Parameters object of model parameters with initial values [and bounds, expressions]
    method  := Fitting algorithm for the optimiser to use
    verbose := Text output (verbosity) level
    """
    
    plotResult = bool(verbose > 1)
    
    nStates = 6
    
    ### Prepare the data
    nRuns = fluxSet.nRuns
    nPhis = fluxSet.nPhis
    nVs = fluxSet.nVs
    
    assert(0 < nPhis)
    assert(0 <= run < nRuns)
    assert(0 <= vInd < nVs)
    
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
        #targetPC.alignToTime()
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
        #nfs.append(targetPC.peak_)
    
    
    ### OFF PHASE
    ### 3a. OFF CURVE: Fit biexponential to off curve to find lambdas
        
    OffKeys = ['Gd1', 'Gd2', 'Gf0', 'Gb0']
    
    iOffPs = Parameters() # Create parameter dictionary
    for k in OffKeys:
        copyParam(k, params, iOffPs)

    ### Trim the first 10% of the off curve to allow I1 and I2 to empty?
    
    
    ### This is an approximation based on the 4-state model which ignores the effects of Go1 and Go2 after light off. 
    
    # lam1 + lam2 == Gd1 + Gd2 + Gf0 + Gb0
    # lam1 * lam2 == Gd1*Gd2 + Gd1*Gb0 + Gd2*Gf0
    
    calcB = lambda Gd1, Gd2, Gf0, Gb0: (Gd1 + Gd2 + Gf0 + Gb0)/2
    calcC = lambda b, Gd1, Gd2, Gf0, Gb0: np.sqrt(b**2 - (Gd1*Gd2 + Gd1*Gb0 + Gd2*Gf0))
    
    def lams(p):
        Gd1 = p['Gd1'].value
        Gd2 = p['Gd2'].value
        Gf0 = p['Gf0'].value
        Gb0 = p['Gb0'].value
        b = calcB(Gd1, Gd2, Gf0, Gb0)
        c = calcC(b, Gd1, Gd2, Gf0, Gb0)
        return (b-c, b+c)
        
    # Create dummy parameters for each phi
    for phiInd in range(nPhis):
        Iss = Ioffs[phiInd][0]
        if Iss < 0:
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, max=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, max=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
        else:
            iOffPs.add('Islow_'+str(phiInd), value=0.2*Iss, vary=True, min=0)
            iOffPs.add('Ifast_'+str(phiInd), value=0.8*Iss, vary=True, min=0, expr='{} - {}'.format(Iss, 'Islow_'+str(phiInd)))
    
    def fit6off(p,t,trial):
        Islow = p['Islow_'+str(trial)].value
        Ifast = p['Ifast_'+str(trial)].value
        lam1, lam2 = lams(p)
        return Islow*np.exp(-lam1*t) + Ifast*np.exp(-lam2*t)
    
    def err6off(p,Ioffs,toffs):
        """Normalise by the first element of the off-curve""" # [-1]
        return np.r_[ [(Ioffs[i] - fit6off(p,toffs[i],i))/Ioffs[i][0] for i in range(len(Ioffs))] ]
    
    #fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-lams(p)[0]*t) + p['a2'].value*np.exp(-lams(p)[1]*t))
    ##fitfunc = lambda p, t: -(p['a0'].value + p['a1'].value*np.exp(-p['lam1'].value*t) + p['a2'].value*np.exp(-p['lam2'].value*t))
    #errfunc = lambda p, Ioff, toff: Ioff - fitfunc(p,toff)
    
    offPmin = minimize(err6off, iOffPs, args=(Ioffs,toffs), method=method)#, fit_kws={'maxfun':100000})
    pOffs = offPmin.params
    
    reportFit(offPmin, "Off-phase fit report for the 6-state model", method)
    if verbose > 0:
        print('Gd1 = {}; Gd2 = {}; Gf0 = {}; Gb0 = {}'.format(pOffs['Gd1'].value, pOffs['Gd2'].value, 
                                                            pOffs['Gf0'].value, pOffs['Gb0'].value))
    
    if plotResult:
        lam1, lam2 = lams(pOffs)
        plotOffPhaseFits(toffs, Ioffs, pOffs, phis, nStates, fit6off, lam1, lam2, Gd=None)
        
    
    # Fix off-curve parameters
    for k in OffKeys:
        pOffs[k].vary = False
    
    
    ### Calculate Go (1/tau_opsin)
    print('\nCalculating opsin activation rate')
    # Assume that Gd1 > Gd2
    # Assume that Gd = Gd1 for short pulses
    
    def solveGo(tlag, Gd, Go0=1000, tol=1e-9):
        Go, Go_m1 = Go0, 0
        while abs(Go_m1 - Go) > tol:
            Go_m1 = Go
            Go = ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
            #Go_m1, Go = Go, ((tlag*Gd) - np.log(Gd/Go_m1))/tlag
        return Go
    
    #if 'shortPulse' in dataSet: # Fit Go
    if quickSet.nRuns > 1:
        from scipy.optimize import curve_fit
        # Fit tpeak = tpulse + tmaxatp0 * np.exp(-k*tpulse)
        #dataSet['shortPulse'].getProtPeaks()
        #tpeaks = dataSet['shortPulse'].IrunPeaks
        
        #PD = dataSet['shortPulse']
        PCs = [quickSet.trials[p][0][0] for p in range(quickSet.nRuns)] # Aligned to the pulse i.e. t_on = 0
        #[pc.alignToTime() for pc in PCs]
        
        #tpeaks = np.asarray([PD.trials[p][0][0].tpeak for p in range(PD.nRuns)]) # - PD.trials[p][0][0].t[0]
        #tpulses = np.asarray([PD.trials[p][0][0].onDs[0] for p in range(PD.nRuns)])
        tpeaks = np.asarray([pc.tpeak_ for pc in PCs])
        tpulses = np.asarray([pc.onDs[0] for pc in PCs])
        
        devFunc = lambda tpulses, t0, k: tpulses + t0 * np.exp(-k*tpulses)
        p0 = (0,1)
        popt, pcov = curve_fit(devFunc, tpulses, tpeaks, p0=p0)
        if plotResult:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal')
            nPoints = 10*int(round(max(tpulses))+1) # 101
            tsmooth = np.linspace(0, max(tpulses), nPoints)
            ax.plot(tpulses, tpeaks, 'x')
            ax.plot(tsmooth, devFunc(tsmooth, *popt))
            ax.plot(tsmooth, tsmooth,'--')
            ax.set_ylim([0,max(tpulses)]) #+5
            ax.set_xlim([0,max(tpulses)]) #+5
        #plt.tight_layout()
        #plt.axis('equal')        
        
        # Solve iteratively Go = ((tlag*Gd) - np.log(Gd/Go))/tlag
        Gd1 = pOffs['Gd1'].value
        Go = solveGo(tlag=popt[0], Gd=Gd1, Go0=1000, tol=1e-9)
        print('t_lag = {:.3g}; Gd = {:.3g} --> Go = {:.3g}'.format(popt[0], Gd1, Go))
        
    elif quickSet.nRuns == 1: #'delta' in dataSet:
        #PD = dataSet['delta']
        #PCs = [PD.trials[p][0][0] for p in range(PD.nRuns)]
        PC = quickSet.trials[0][0][0]
        tlag = PC.lag_ # := lags_[0] ############################### Add to Photocurrent...
        Go = solveGo(tlag=tlag, Gd=Gd1, Go0=1000, tol=1e-9)
        print('t_lag = {:.3g}; Gd = {:.3g} --> Go = {:.3g}'.format(tlag, Gd1, Go))
    
    else:
        Go = 1 # Default
        print('No data found to estimate Go: defaulting to Go = {}'.format(Go))
    
    
    ### ON PHASE
    
    iOnPs = Parameters() # deepcopy(params)
    
    # Set parameters from Off-curve optimisation
    for k in OffKeys:
        copyParam(k, pOffs, iOnPs)
    
    # Set parameters from general rhodopsin analysis routines
    for k in ['Go1', 'Go2', 'k1', 'k2', 'k_f', 'k_b', 'gam', 'p', 'q', 'phi_m', 'g0', 'Gr0', 'E', 'v0', 'v1']: #.extend(OffKeys):
        copyParam(k, params, iOnPs)
    
    # Set parameters from short pulse calculations
    iOnPs['Go1'].value = Go; iOnPs['Go1'].vary = False
    iOnPs['Go2'].value = Go; iOnPs['Go2'].vary = False
    
    RhO = models['6']()
    
    ### Trim down ton? Take 10% of data or one point every ms? ==> [0::5]
    
    if verbose > 2:
        print('Optimising ',end='')
    
    onPmin = minimize(errOnPhase, iOnPs, args=(Ions,tons,RhO,Vs,phis), method=method)
    pOns = onPmin.params
    
    reportFit(onPmin, "On-phase fit report for the 6-state model", method)
    if verbose > 0:
        print('k1 = {}; k2 = {}; k_f = {}; k_b = {}'.format(pOns['k1'].value, pOns['k2'].value, 
                                                        pOns['k_f'].value, pOns['k_b'].value))
        print('gam = {}; phi_m = {}; p = {}; q = {}'.format(pOns['gam'].value, pOns['phi_m'].value, 
                                                            pOns['p'].value, pOns['q'].value))
    
    fitParams = pOns
    
    return fitParams




#TODO: Tidy up and refactor getRecoveryPeaks and fitRecovery
def getRecoveryPeaks(recData, phiInd=None, vInd=None, usePeakTime=False):
    """Aggregate the times and currents of the photocurrent peaks for fitting Gr0
        
        usePeakTime := False {True, False} 
                        Use time of second peak (t_peak1) 
                        otherwise time of second pulse (t_on1)
    """
    
    if phiInd is None:
        phiMax, phiInd = getExt(recData.phis, 'max')
    
    if vInd is None:
        if recData.nVs == 1:
            vIndm70 = 0
        else:
            try: 
                vIndm70 = setPC.Vs.index(-70)
            except:
                vIndm70 = np.searchsorted(setPC.Vs, -70)
            #vIndm70 = getIndex(setPC.Vs, -70)
        vInd = vIndm70
    
    tpeaks1 = []
    Ipeaks1 = []
    
    ### Build array of second peaks
    for run in range(recData.nRuns):
        PC = recData.trials[run][phiInd][vInd]
        PC.alignToPulse(pulse=0, alignPoint=2) # End of the first pulse
        if usePeakTime:
            tpeaks1.append(recData.trials[run][phiInd][vInd].tpeaks_[1]) # Time of second peak
        else:
            tpeaks1.append(recData.trials[run][phiInd][vInd].pulses[1,0]) # Time of second pulse
        Ipeaks1.append(recData.trials[run][phiInd][vInd].peaks_[1])
    
    # Check for sorting...
    
    # Prepend t_off0 and Iss0
    run = 0 # Take comparators from the first run's first pulse
    tss0 = recData.trials[run][phiInd][vInd].pulses[0,1]
    Iss0 = recData.trials[run][phiInd][vInd].sss_[0]
    Ipeak0 = recData.trials[run][phiInd][vInd].peaks_[0]
    ### This would be correct except that O->C transitions confound the first ~500ms
    #t_peaks = np.r_[tss0, tpeaks1]
    #I_peaks = np.r_[Iss0, Ipeaks1]
    ### Therefore use the peaks only
    t_peaks = np.asarray(tpeaks1)
    I_peaks = np.asarray(Ipeaks1)
    
    return t_peaks, I_peaks, Ipeak0, Iss0    
    
    
def fitRecovery(t_peaks, I_peaks, params, Ipeak0, Iss0, ax=None, method=defMethod, verbose=verbose):
    
    
    if not params['Gr0'].vary:
        print('Gr0 fixed at {}'.format(params['Gr0'].value))
        return params
    
    plotResult = bool(verbose > 1)
    
    def errExpRec(p, t, I=None):
        # Restrict so that a = -c to ensure (0,0) is passed through?
        #model = p['a'].value * np.exp(-p['Gr0'].value*t) - p['Ipeak0'].value
        model = p['Ipeak0'].value - p['a'].value * np.exp(-p['Gr0'].value*t)
        if I is None:
            return model
        return I - model
    
    shift = t_peaks[0]
    # if np.isclose(shift, 0):
        # Iss0 = I_peaks[0]
    # else:
        # Iss0 = 0.5 * Ipeak0 ### Reconsider
    
    # if np.isclose(I_peaks[0], Iss0):
        # if np.isclose(shift, 0):
            # # Data has been prepended with first steady-state values
        # else:
            # shift = t_peaks[0]
    # else:
        
    # if not np.isclose(shift, 0):
        # Iss0 = 0.5 * Ipeak0 ### Reconsider
        # warnings.warn("Realigning peak times!")
    
    iRecPs = Parameters() # Create parameter dictionary
    copyParam('Gr0', params, iRecPs)
    ### a is now a dummy parameter
    #pRec.add('a', value=Iss0+Ipeak0, expr='{Iss0} + Ipeak0'.format(Iss0=Iss0)) # Iss = a - c
    #pRec.add('Ipeak0', value=-Ipeak0, vary=True) # Ipeak orig
    iRecPs.add('a', value=Ipeak0-Iss0) #, expr='Ipeak0 - {Iss0}'.format(Iss0=Iss0)) # Iss = a - c
    iRecPs.add('Ipeak0', value=Ipeak0, vary=False) # Ipeak orig
    #pRec.add('Ipeak0', value=Ipeak0-Iss0, vary=False) # Ipeak orig
    
    ### Shift is now handled in getRecoveryPeaks()
    #recMin = minimize(errExpRec, pRec, args=(t_peaks-shift, I_peaks), method=method)
    recMin = minimize(errExpRec, iRecPs, args=(t_peaks, I_peaks), method=method)
    #recMin = minimize(errExpRec, pRec, args=(t_peaks-shift, I_peaks-Iss0), method=method)
    
    pRec = recMin.params
    
    if verbose > 1:
        chosenFit = recMin.chisqr
        fits = {}
        for k in iRecPs: # Check all parameters have finite bounds for np.isfinite bug in differential_evolution
            if iRecPs[k].min is None or not np.isfinite(iRecPs[k].min):
                iRecPs[k].min = -1e15
            if iRecPs[k].max is None or not np.isfinite(iRecPs[k].max):
                iRecPs[k].max = 1e15
        
        for meth in methods:
            recMinAlt = minimize(errExpRec, iRecPs, args=(t_peaks, I_peaks), method=meth)
            fits[meth] = recMinAlt.chisqr
            #print(fit_report(recMin))
            if fits[meth] < chosenFit:
                print("Consider using the '{}' algorithm for a better fit (chisqr = {:.3}) ==> Gr0 = {:.3}".format(meth, fits[meth], recMinAlt.params['Gr0'].value))
                if verbose > 2:
                    print(fit_report(recMinAlt))
    
    # popt, pcov = curve_fit(curveFunc, t_peaks-shift, I_peaks, p0=p0) #Needs ball-park guesses (0.3, 125, 0.5)
    # peakEq = eqString.format(*[round_sig(p,3) for p in popt]) # *popt rounded to 3s.f.
    
    copyParam('Gr0', pRec, params)
    
    eqString = '$I_{{peak}} = {Ipeak0:+.3} - {a:.3}e^{{-{Gr0:g}\cdot t}}$'
    v = pRec.valuesdict()
    peakEq = eqString.format(a=round_sig(v['a'],3), Gr0=round_sig(v['Gr0'],3), Ipeak0=round_sig(v['Ipeak0'],3))
    
    if verbose > 1:
        print(peakEq)
    
    if ax is None:
        if plotResult:
            fig = plt.figure()
            ax = plt.subplot(111)
        else:
            return params
        
    #else:
    ax.scatter(t_peaks, I_peaks, color='r', marker='*')
    # Freeze axes
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax)
    nPoints = 10*int(round(xmax-xmin)+1)
    tsmooth = np.linspace(xmin, xmax, nPoints)#1+(xmax-xmin)*10) #totT
    Ismooth = errExpRec(pRec, tsmooth) #curveFunc(tsmooth,*popt)
    
    #ax.plot(tsmooth+shift, Ismooth, linestyle=':', color='r')#, linewidth=1.5*mpl.rcParams['lines.linewidth'])
    ax.plot(tsmooth, Ismooth, linestyle=':', color='r')
    ax.axhline(y=Iss0, linestyle=':', color='#aaaaaa')
    ax.axhline(y=Ipeak0, linestyle=':', color='#aaaaaa')
    
    x = 0.8
    y = 1.02*Ismooth[-1] #popt[2]
    ax.text(x*xmax, y, peakEq, ha='center', va='top')#, fontsize=config.eqSize) #, transform=ax.transAxes)
    
    # else:
        # fig = plt.figure()
        # ax1 = plt.subplot(211)
        # ax1.scatter(t_peaks, I_peaks)
        # tsmooth = np.linspace(xmin,xmax,1+(xmax-xmin)*10)
        # ax1.plot(tsmooth, errExpRec(pRec, tsmooth))
        # ax1.axhline(y=Iss0, linestyle=':')
        # ax1.axhline(y=Ipeak0, linestyle=':')
        # ax1.set_ylabel(r'$I_{peak} \mathrm{[nA]}$')
        # plt.setp(ax1.get_xticklabels(), visible=False)
        # ax2 = plt.subplot(212, sharex=ax1)
        # ax2.scatter(t_peaks, abs(I_peaks)/max(abs(I_peaks)))
        # ax2.plot(tsmooth, 1-(Iss0/Ipeak0)*np.exp(-Gr0*tsmooth))
        # ax2.set_ylabel(r'$\mathrm{Proportion of}\ I_{peak0}$')
        # ax2.set_xlabel(r'$\mathrm{Time [ms]}$')
        # ax2.axhline(y=1-Iss0/Ipeak0, linestyle=':')
        # ax2.text(x*xmax, 1-Iss0/Ipeak0, '$1-I_{ss0}/I_{peak0}$', ha='center', va='baseline', fontsize=eqSize)
        

    # if verbose > 1:
        # print("Parameters: {}".format(popt))
        # if type(pcov) in (tuple, list):
            # print("$\sigma$: {}".format(np.sqrt(pcov.diagonal())))
        # else:
            # print("Covariance: {}".format(pcov))
    # return popt, pcov, peakEq
    
    return params
    



    
#TODO; Refactor all fitting functions to do with fV and FV
def errfV(pfV, V, fVs=None):
    v = pfV.valuesdict()
    v0 = v['v0']
    v1 = v['v1']
    E = v['E']
    #if type(V) != np.ndarray:
    #    V = np.array(V)
    V = np.asarray(V)
    fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) #*(v1/(V-E)) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
    #zeroErrs = np.isclose(V, np.ones_like(V)*E)
    #fV[zeroErrs] = v1/v0
    fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
    if fVs is None:
        return fV
    return fVs - fV #calcfV(pfV, V)

# def errFV(pfV, V, FVs=None): #, v0, v1, E
    # V = np.asarray(V)
    # FV = errfV(pfV, V) * (V - pfV['E'].value)# * 1e-6
    # if FVs is None:
        # return FV
    # return FVs - FV
    
def errFV(pfV, V, FVs=None):
    """F(v) := f(v)*(v-E)"""
    v = pfV.valuesdict()
    v0 = v['v0']
    v1 = v['v1']
    E = v['E']
    ###if type(V) != np.ndarray:
    ###    V = np.array(V)
    V = np.asarray(V)
    FV = v1*(1-np.exp(-(V-E)/v0)) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
    #FV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
    ###zeroErrs = np.isclose(V, np.ones_like(V)*E)
    ###fV[zeroErrs] = v1/v0
    ###fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
    FV[np.isnan(FV)] = v1/v0 # Fix the error when dividing by zero
    #FV *= (V-E)
    if FVs is None:
        return FV #* 1e-6
    return FVs - FV #* 1e-6

# def calcFV(pfV, V): #, v0, v1, E
    # return calcfV(pfV, V) * (V - pfV['E'].value)
    
# def errFV(pfV, V, FVs): #, v0, v1, E
    # return FVs - (calcfV(pfV, V) * (V - pfV['E'].value))
    
def _calcfVnew(V, v0, E):
    if type(V) != np.ndarray:
        V = np.array(V)
    fV = (E+70)*(np.exp((E-V)/v0)-1)/((E-V)*(np.exp((E+70)/v0)-1)) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
    fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
    return fV #* (V - E)
    
def fitfV(Vs, Iss, params, relaxFact=2, method=defMethod, verbose=verbose):
    """Fitting function to find the parameters of the voltage dependence function"""
    
    # Use @staticmethod or @classmethod on RhodopsinModel.calcfV() and pass in parameters?

    
    ### Skip V == E
    #Prot.Vs = list(range(-100,80,5))
    #try:
    #    del Prot.Vs[Prot.Vs.index(0)]
    #except ValueError:
    #    pass
    
    plotResult = bool(verbose > 1)
    
    
    ifVPs = Parameters() # Create parameter dictionary
    
    # for p in ['E', 'v0', 'v1']:
        # copyParam(p, params, ifVPs)
        
    copyParam(['E', 'v0', 'v1'], params, ifVPs)
    
    Iss = np.asarray(Iss)
    #Vs = np.asarray(Vs)
    
    if 'g0' in params:
        g0 = params['g0'].value
    else:
        g0 = 25000
    pseudoV1 = calcV1(ifVPs['E'].value, ifVPs['v0'].value) * (g0 * 1e-6 * 0.5) # g0*f(phi)*v1 (assuming I in nA and f(phi)=0.5)
    
    Emeth = 'cf' #''#'model'#
    if Emeth == 'cf': # No bounds or expressions applied. Uses the Levenberg-Marquardt algorithm
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
        
        def calcRect(V, v0, v1, E): #, gpsi):
            if type(V) != np.ndarray:
                V = np.array(V)
            fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
            fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
            return fV * (V - E) # * gpsi
        
        #p0 = (pfV['v0'].value, pfV['v1'].value, pfV['E'].value)
        #calcV1(0,50)
        
        p0 = (ifVPs['v0'].value, pseudoV1, ifVPs['E'].value) # This routine is not actually finding v1 so the initial value needs reconsidering - it ncludes other scaling factors e.g. f_phi, g0
        #p0 = (pfV['v0'].value, 1, pfV['E'].value)
        popt, pcov = curve_fit(calcRect, Vs, Iss, p0) # (curveFunc, Vs, Iss, p0=p0)
        
        #ifVPs['v0'].value = popt[0] #; pfV['v0'].vary = False # Causes df error
        #ifVPs['v1'].value = popt[1] # Set for plotting
        #ifVPs['E'].value = popt[2]
        #pfV['E'].vary = False
        
        pfV = Parameters()
        copyParam(['E', 'v0', 'v1'], ifVPs, pfV)
        pfV['v0'].value = popt[0] #; pfV['v0'].vary = False # Causes df error
        pfV['v1'].value = popt[1] # Set for plotting
        pfV['E'].value = popt[2]
        
    
    elif Emeth == 'model':   ### Model based fitting: http://lmfit.github.io/lmfit-py/model.html
    
        def calcFV(Vs, E, v0, v1):
            Vs = np.asarray(Vs)
            FV = v1*(1-np.exp(-(Vs-E)/v0))#/((V-E)/v1)
            FV[np.isnan(FV)] = v1/v0
            return FV
        FVmod = Model(calcFV)
        FVmod.set_param_hint('E', value=ifVPs['E'].value, min=ifVPs['E'].min, max=ifVPs['E'].max, vary=ifVPs['E'].vary, expr=ifVPs['E'].expr)
        FVmod.set_param_hint('v0', value=ifVPs['v0'].value, min=ifVPs['v0'].min, max=ifVPs['v0'].max, vary=ifVPs['v0'].vary, expr=ifVPs['v0'].expr)
        #FVmod.set_param_hint('v1', value=pfV['v1'].value, min=pfV['v1'].min, max=pfV['v1'].max, vary=pfV['v1'].vary, expr=pfV['v1'].expr)
        FVmod.set_param_hint('v1', value=pseudoV1)
        modParams = FVmod.make_params()
        result = FVmod.fit(Iss, Vs=np.asarray(Vs), method=method) # , E=pfV['E'].value, v0=pfV['v0'].value, v1=pfV['v1'].value  method=method
        
        pfV = Parameters()
        for p in ['E', 'v0', 'v1']:
            copyParam(p, params, pfV)
            pfV[p].value = result.best_values[p]
        
        # pfV['E'].value = result.best_values['E']
        # pfV['v0'].value = result.best_values['v0']
        # pfV['v1'].value = result.best_values['v1'] # Set for plotting
        
    else:   # lmfit
    
        ### Set bounds on E based on zero-crossing in Iss array
        inds = np.argsort(Vs)
        IssSorted = Iss[inds]
        VsSorted = Vs[inds]
        signs = np.sign(IssSorted)
        # if np.any(signs == 0): # Test for zero current
        #max(Iss<0)
        #min(Iss>0)
        if signs[0] <= 0 and signs[-1] >= 0: # There is a change of sign or a zero
            for i,s in enumerate(signs):
                if s == 0: #np.isclose(s, 0) # Test for zero current
                    ifVPs['E'].value = VsSorted[i]
                    break
                elif s > 0:
                    ifVPs['E'].min = VsSorted[i-1]
                    ifVPs['E'].max = VsSorted[i]
                    print('Limits on E set: [{}, {}]'.format(ifVPs['E'].min, ifVPs['E'].max))
                    break
        
        # for v,i in zip(VsSorted,IssSorted):
            # if signs[]
        # pfV['E'].min = 
        
        #pfVfresh = deepcopy(pfV)
        
        ifVPs['v1'].value = pseudoV1
        ifVPs['v1'].min = None
        ifVPs['v1'].max = None
        
        if params['E'].vary:
            FVmin = minimize(errFV, ifVPs, args=(Vs, Iss), method=method) # kws={'FVs':Iss},
            print('1st stage v0: ', ifVPs['v0'].value)
            pfV = FVmin.params
            
            if verbose > 1:
                chosenFit = FVmin.chisqr
                fits = {}
                for meth in methods:
                    FVminAlt = minimize(errFV, ifVPs, args=(Vs, Iss), method=meth)
                    fits[meth] = FVminAlt.chisqr
                    if fits[meth] < chosenFit:
                        print("Consider using the '{}' algorithm for a better fit (chisqr = {:.3}) ==> E = {:.3}".format(meth, fits[meth], FVminAlt.params['E'].value))
                        if verbose > 2:
                            print(fit_report(FVminAlt))
            

    
    ###print('Estimate of f_phi = {}'.format(Ipeak/(pfV['g0'].value*pfV['v1'].value)))
    
    E = pfV['E'].value
    v0 = pfV['v0'].value
    
    setBounds(pfV['v0'], relaxFact)
    
    if plotResult:
        #nPoints = 10*int(round(endT-begT/self.dt))+1
        Vsmooth = np.linspace(min(Vs), max(Vs), 10*round(max(Vs)-min(Vs))+1) #1+(max(Vs)-min(Vs))/.1
        fig, ax1 = plt.subplots()
        ax1.plot(Vsmooth, errFV(pfV, Vsmooth), 'b', label='$I_{ss}$')
        ax1.scatter(Vs, Iss, c='b', marker='x')
        ax1.set_ylabel('$I_{ss}$ $\mathrm{[nA]}$', color='b') #$f(V) \cdot (V-E)$
        ax1.set_xlabel(r'$V_{clamp}\ \mathrm{[mV]}$')
    
    
    pfV['v1'].value = calcV1(E, v0)
    
    #pfV['E'].vary = False
    #pfV['E'].min = pfV['E'].value - 5
    #pfV['E'].max = pfV['E'].value + 5
    
    # if method != 'powell':  # Powell algorithm errors with only 1 d.f.: TypeError: zip argument #2 must support iteration
        # pfV['v1'].expr = '(70+E)/(exp((70+E)/v0)-1)'
    # else:
        # pfV['v1'].expr = '(70+E)/(exp((70+E)/v0)-1)'
        # pfV['E'].vary = True
        # pfV['E'].min = pfV['E'].value + 1e-9#* 2
        # pfV['E'].max = pfV['E'].value - 1e-9#/ 2
    
    pfV['v1'].expr = '(70+E)/(exp((70+E)/v0)-1)'
    pfV['E'].vary = True # Required since Powell algorithm errors with only 1 d.f.
    #setBounds(pfV['E'], relaxFact)
    pfV['E'].min = pfV['E'].value + 1e-9
    pfV['E'].max = pfV['E'].value - 1e-9
    
    
    try:
        vIndm70 = Vs.index(-70)
    except:
        cl = np.isclose(Vs, np.ones_like(Vs)*-70)
        vIndm70 = np.searchsorted(cl, True)
        #vIndm70 = np.searchsorted(Vs, -70)
    #vIndm70 = getIndex(Vs, -70)
    if verbose > 1:
        print('V=-70 at element {} ({})'.format(vIndm70, Vs[vIndm70]))
    
    gs = Iss / (np.asarray(Vs) - E) # 1e6 * 
    gm70 = Iss[vIndm70] / (-70 - E)# * -70 # 1e6 * 
    if verbose > 1:
        print('g(v=-70) = ', gm70)
    #gs[(Vs - E)==0] = None #(v1/v0)
    gNorm = gs / gm70 # Normalised conductance relative to V=-70
    
    if verbose > 2:
        print(np.c_[Vs,Iss,gs,gNorm]) #np.asarray(Vs)-E
    
    if params['v0'].vary or params['v1'].vary:
        fVmin = minimize(errfV, pfV, args=(Vs, gNorm), method=method)#, tol=1e-12)
        pfVfinal = fVmin.params
        
        if verbose > 1:
            chosenFit = fVmin.chisqr
            fits = {}
            for meth in methods:
                if method != 'powell':
                    pfV['E'].vary = False
                else:
                    pfV['E'].vary = True
                pfV['v1'].expr = '(70+E)/(exp((70+E)/v0)-1)'
                fVminAlt = minimize(errfV, pfV, args=(Vs, gNorm), method=meth)
                fits[meth] = fVminAlt.chisqr
                if fits[meth] < chosenFit:
                    print("Consider using the '{}' algorithm for a better fit (chisqr = {:.3}) ==> v0 = {:.3}, v1 = {:.3}".format(meth, fits[meth], fVminAlt.params['v0'].value, fVminAlt.params['v1'].value))
                    if verbose > 2:
                        print(fit_report(fVminAlt))
            
    pfVfinal['v0'].vary = False
    pfVfinal['v1'].vary = False
    
    # Corrections
    v0 = pfVfinal['v0'].value
    v1 = pfVfinal['v1'].value
    zeroErrs = np.isclose(Vs, np.ones_like(Vs)*E, rtol=1e-3)
    gNorm[zeroErrs] = v1/v0
    
    if plotResult:
        ax2 = ax1.twinx()
        eqString = r'$f(v) = \frac{{{v1:.3}}}{{v-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{v-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'
        fVstring = eqString.format(E=pfVfinal['E'].value, v0=pfVfinal['v0'].value, v1=pfVfinal['v1'].value)
        #vd = pfV.valuesdict()
        #fVstring = eqString.format(**vd)
        ax2.plot(Vsmooth, errfV(pfVfinal, Vsmooth), 'g', label=fVstring)
        ax2.scatter(Vs, gNorm, c='g', marker='+')
        ax2.set_ylabel('$f(v)$ $\mathrm{[1]}$', color='g')
        ax2.axvline(x=E, linestyle=':', color='k')
        ymin, ymax = ax2.get_ylim()
        revString = '$E = {}\ \mathrm{{[mV]}}$'.format(round_sig(E,3))
        ax2.text(E, 0.05*(ymax-ymin), revString, ha='center', va='center', fontsize=config.eqSize)
        # ax2.axvline(x=-70, linestyle=':', color='k')
        # ax2.axhline(y=1, linestyle=':', color='k')
        ymin, ymax = ax2.get_ylim();    ax2.set_ylim(ymin, ymax)
        xmin, xmax = ax2.get_xlim();    ax2.set_xlim(xmin, xmax)
        fVsmooth = errfV(pfVfinal, Vsmooth)
        m70ind = np.searchsorted(Vsmooth, -70)
        #ax2.vlines(x=-70, ymin=ymin, ymax=fVsmooth[m70ind], linestyles=':', colors='g') # ymax=errfV(pfV, -70)
        #ax2.hlines(y=1, xmin=Vsmooth[m70ind], xmax=xmax, linestyles=':', colors='g')
        ax2.vlines(x=-70, ymin=ymin, ymax=1, linestyles=':', colors='g') # ymax=errfV(pfV, -70)
        ax2.hlines(y=1, xmin=-70, xmax=xmax, linestyles=':', colors='g')
        plt.legend()
    
    # for p in ['E', 'v0', 'v1']:
        # copyParam(p, pfV, params)
    
    copyParam(['E', 'v0', 'v1'], pfVfinal, params)
    
    return params


def fitFV(Vs, Iss, p0, ax=None):
    """F(V):= f(V)*(V-E)"""
    # RhO.fV(V) * (V-RhO.E)
    if ax is None:
        ax = plt.gcf()
    #markerSize=40
    #eqString = r'$f(V) = \frac{{{v1:.3}}}{{V-{E:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{V-{E:+.2f}}}{{{v0:.3}}}}}\right)\right]$'

    def calcRect(V, v0, v1, E): #, gpsi):
        if type(V) != np.ndarray:
            V = np.array(V)
        fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
        fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
        return fV * (V - E) # * gpsi
    
    #psi = RhO.calcPsi(RhO.steadyStates)
    #sf = RhO.A * RhO.gbar * psi * 1e-6 # Six-state only
    #sf = RhO.g * psi * 1e-6 
    
    #sf = Iss[Vs.index(-70)]# * -70
    Iss = np.asarray(Iss)#/sf # np.asarray is not needed for the six-state model!!!
    #print(np.c_[Vs,Iss])
    
    #xfit = np.linspace(min(Vs), max(Vs), 1+(max(Vs)-min(Vs))/.1) #Prot.dt
    #yfit = calcRect(xfit, *p0FVnew)#*sf
    #ax.plot(xfit, yfit)

    #plt.figure()
    #plt.plot(xfit, yfit / sf)#((xfit - pFit[2]) * popt[3]))
    #plt.scatter(Vs, Iss, marker='x', s=markerSize)
    #plt.plot(xfit, calcRect(xfit, *p0FVnew))
    
    popt, pcov = curve_fit(calcRect, Vs, Iss, p0) # (curveFunc, Vs, Iss, p0=p0)
    
    #pFit = [round_sig(p,3) for p in popt]
    #print(pFit)
    print('Phase I curve fit: ', popt)
    #peakEq = eqString.format(v0=pFit[0], E=pFit[2], v1=pFit[1])
    
    v0 = popt[0]
    v1 = popt[1]
    E = popt[2]
    
    #vInd = np.searchsorted(Vs, (-70 - E))
    #sf = Iss[vInd]
    Im70 = Iss[Vs.index(-70)]# * -70
    gs = Iss / (Vs - E)
    
    #g0[(Vs - E)==0] = None #(v1/v0)
    gNorm = gs / (Im70 / (-70 - E))
    zeroErrs = np.isclose(Vs, np.ones_like(Vs)*E)
    gNorm[zeroErrs] = v1/v0
    
    def calcScale(V, v0, v1):
        if type(V) != np.ndarray:
            V = np.array(V)
        fV = (1-np.exp(-(V-E)/v0))/((V-E)/v1) # Dimensionless #fV = abs((1 - exp(-v/v0))/v1) # Prevent signs cancelling
        fV[np.isnan(fV)] = v1/v0 # Fix the error when dividing by zero
        return fV
    
    poptrel, pcov = curve_fit(calcScale, Vs, gNorm, p0=(v0, v1))
    print('Phase II curve fit: ', poptrel)
    
    if verbose > 1:
        print(np.c_[Vs,Iss,gs,gNorm])
    
    #popt[0] = poptrel[0]
    #popt[1] = poptrel[1]
    
    # Vrange = max(Vs) - min(Vs)
    # xfit = np.linspace(min(Vs), max(Vs), 1+Vrange/.1) #Prot.dt
    # yfit = calcRect(xfit, *popt)*sf
    
    #peakEq = eqString.format(*[round_sig(p,3) for p in popt])
    
    # ax.plot(xfit, yfit)#,label=peakEq)#,linestyle=':', color='#aaaaaa')
    # #col, = getLineProps(Prot, 0, 0, 0) #Prot, run, vInd, phiInd
    # #plt.plot(Vs,Iss,linestyle='',marker='x',color=col)
    # ax.scatter(Vs, Iss, marker='x', color=colours, s=markerSize)#,linestyle=''
    
    # x = 1 #0.8*max(Vs)
    # y = 1.2*yfit[-1]#max(IssVals[run][phiInd][:])
    # plt.text(-0.8*min(Vs),y,peakEq,ha='right',va='bottom',fontsize=eqSize)#,transform=ax.transAxes)
    
    return popt, poptrel


def getNormGs(Vs, Iss, E, v=-70):
    
    try:
        vIndm70 = Vs.index(-70)
    except:
        cl = np.isclose(Vs, np.ones_like(Vs)*-70)
        vIndm70 = np.searchsorted(cl, True)
        #vIndm70 = np.searchsorted(Vs, -70)
    #vIndm70 = getIndex(Vs, -70)
    if verbose > 1:
        print('V=-70 at element {} ({})'.format(vIndm70, Vs[vIndm70]))
    
    gs = Iss / (np.asarray(Vs) - E) # 1e6 * 
    gm70 = Iss[vIndm70] / (-70 - E)# * -70 # 1e6 * 
    if verbose > 1:
        print('g(v=-70) = ', gm70)
    #gs[(Vs - E)==0] = None #(v1/v0)
    gNorm = gs / gm70 # Normalised conductance relative to V=-70
    '''
    #vInd = np.searchsorted(Vs, (-70 - E))
    #sf = Iss[vInd]
    Im70 = Iss[Vs.index(-70)]# * -70
    gs = Iss / (Vs - E)
    
    #g0[(Vs - E)==0] = None #(v1/v0)
    gNorm = gs / (Im70 / (-70 - E))
    '''
    #zeroErrs = np.isclose(Vs, np.ones_like(Vs)*E)
    #gNorm[zeroErrs] = v1/v0
    return gNorm

    
    

def calcCycle(p, ton, toff, RhO, V, phi, fitRates=False): #,fitDelay=False
    """Simulate the on and off-phase from base parameters"""
    if verbose > 1:
        print('.', end="") # sys.stdout.write('.')
    
    #Idel, Ion, Ioff = I[:onInd+1], I[onInd:offInd+1], I[offInd:]
    #tdel, ton, toff = t[:onInd+1], t[onInd:offInd+1]-t[onInd], t[offInd:]-t[offInd]
    
    RhO.initStates(0)
    #RhO.setLight(phi) # Calculate transition rates for phi
    RhO.updateParams(p)
            
    if 0: # fitDelay: # Change to pass t array and pulse indices
        # Delay phase
        RhO.setLight(RhO.phi_0)
        if RhO.useAnalyticSoln:
            soln = RhO.calcSoln(tdel, RhO.s_0)
        else:
            soln = odeint(RhO.solveStates, RhO.s_0, tdel, Dfun=RhO.jacobian)
        RhO.storeStates(soln[1:], tdel[1:])
    
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
    return np.r_[ [(Is[i] - calcCycle(p,tons[i],toffs[i],RhO,Vs[i],phis[i]))/nfs[i] for i in range(len(Is))]]



def fitModels(dataSet, nStates=3, params=None, postFitOpt=True, relaxFact=2, method=defMethod, postFitOptMethod=None, verbose=verbose):
    '''
    #TODO """Routine to fit as many models as possible and select between them according to some parsimony criterion"""
    
    #fit3s=True, fit4s=False, fit6s=False
    ### Pass 'auto' to fit the highest model possible
    
    ### Pass 'all' then compare fit vs computational complexity...
    # Plot all model fits on the same data plots
    # RhO = selectModel(nStates)
    # Calculate chisqr for each model and select between them. 
    
    
    # Run small signal analysis
    #runSSA(RhO)
    #characterise(RhO)    
    
    
    if not isinstance(nStates, (list, tuple)):
        nStates = [nStates]
    else:
        nStates = nStates
    if not isinstance(params, (list, tuple)):
        params = [params]
    else:
        params = params
        
    for nSt in nStates:
        
    '''
    
    if not isinstance(nStates, (list, tuple)):
        nStates = [nStates]
    else:
        nStates = nStates
    if not isinstance(params, (list, tuple)):
        params = [params]
    else:
        params = params
    
    fitParams = [None for nSt in nStates]
    for i, nSt in enumerate(nStates):
        fitParams[i] = fitModel(dataSet, nStates=nStates[i], params=params[i], 
        postFitOpt=postFitOpt, relaxFact=relaxFact, method=method, postFitOptMethod=postFitOptMethod, verbose=verbose)

        
def fitModel(dataSet, nStates=3, params=None, postFitOpt=True, relaxFact=2, method=defMethod, postFitOptMethod=None, verbose=verbose):
    """Fit a model (with initial parameters) to a dataset of optogenetic photocurrents"""
    
    
    ### Define non-optimised parameters to exclude in post-fit optimisation
    nonOptParams = ['Gr0', 'E', 'v0', 'v1']
    
    if isinstance(nStates, str):
        nStates = int(nStates) # .lower()
    if nStates == 3 or nStates == 4 or nStates == 6: 
        pass
    else:
        print("Error in selecting model - please choose from 3, 4 or 6 states")
        raise NotImplementedError(nStates)
    
    if verbose > 0:
        t0 = wallTime()
        print("\n================================================================================")
        print("Fitting parameters for the {}-state model with the '{}' algorithm... ".format(nStates, method))
        print("================================================================================\n")
        
    ### Check contents of dataSet and produce report on model features which may be fit. 
    # e.g. if not 'rectifier': f(V)=1
    
    ### Trim data slightly to remove artefacts from light on/off transition ramps?
    
    ### Could use precalculated lookup tables to find the values of steady state O1 & O2 occupancies?
    
    if params is None:
        params = modelParams[str(nStates)]
    
    if isinstance(dataSet, dict):
        if 'step' in dataSet:
            fluxKey = 'step'
        elif 'custom' in dataSet:
            fluxKey = 'custom'
        else:
            raise KeyError("Flux set not found: Expected 'step' or 'custom'. ")
    else:
        fluxKey = 'step'
    
    
    # Determine the data type passed
    if isinstance(dataSet, PhotoCurrent): # Single photocurrent
        nRuns = 1
        nPhis = 1
        nVs = 1
        pc = copy.deepcopy(dataSet)
        setPC = ProtocolData(fluxKey, nRuns, [pc.phi], [pc.V])
        setPC.trials[0][0][0] = pc
        dataSet = {fluxKey:setPC}
    elif isinstance(dataSet, ProtocolData): # Set of photocurrents
        setPC = copy.deepcopy(dataSet)
        nRuns = setPC.nRuns
        nPhis = setPC.nPhis
        nVs = setPC.nVs
        dataSet = {fluxKey:setPC}
    elif isinstance(dataSet[fluxKey], PhotoCurrent): # Single photocurrent in dictionary
        nRuns = 1
        nPhis = 1
        nVs = 1
        pc = copy.deepcopy(dataSet[fluxKey])
        setPC = ProtocolData(fluxKey, nRuns, [pc.phi], [pc.V])
        setPC.trials[0][0][0] = pc
        dataSet = {fluxKey:setPC}
    elif isinstance(dataSet[fluxKey], ProtocolData): # Set of photocurrents in dictionary
        setPC = dataSet[fluxKey]
        nRuns = setPC.nRuns
        nPhis = setPC.nPhis
        nVs = setPC.nVs
        #if (setPC.nPhis * setPC.nVs) > 1:
    # elif isinstance(dataSet[fluxKey], PhotoCurrent): # Single photocurrent within ProtocolData
        # #targetPC = dataSet['custom']
        # nRuns = 1
        # nPhis = 1
        # nVs = 1
        # #setPC = [[[dataSet['custom']]]]
        # setPC = ProtocolData(fluxKey, nRuns, [dataSet[fluxKey].phi], [dataSet[fluxKey].V])
        # setPC.trials[0][0][0] = dataSet[fluxKey]
    else:
        print(type(dataSet[fluxKey]))
        print(dataSet[fluxKey])
        raise TypeError("dataSet[fluxKey]")
    
        
    if nRuns == 1:
        runInd = 0
        
    if nVs == 1:
        vIndm70 = 0
    else:
        try: 
            vIndm70 = setPC.Vs.index(-70)
        except:
            vIndm70 = np.searchsorted(setPC.Vs, -70)
        
        #vIndm70 = getIndex(setPC.Vs, -70)
    
    if nPhis == 1:
        params['phi_m'].vary = False
        params['p'].vary = False
        # Fix other model specific parameters?
        if 'q' in params: #nStates == 4 or nStates == 6:
            params['q'].vary = False
        nonOptParams.extend(['phi_m', 'p', 'q'])
        # if nStates == 4 or nStates == 6:
            # nonOptParams.append('Gf0')
            # nonOptParams.append('Gb0')
        ### Allow these to vary to *effectively* fit Gf and Gb
    
    PCs = [setPC.trials[runInd][phiInd][vIndm70] for phiInd in range(nPhis)]
    Vs = [pc.V for pc in PCs]
    phis = [pc.phi for pc in PCs]
    
        
    ### Extract the parameters relevant to all models - move inside loop for recovery protocol?
    
    ### Optionally fit f(V) (inward rectification) parameters with rectifier data: v0, v1
    # MUST MEASURE E AND FIT AFTER OTHER PARAMETERS

    if 'rectifier' in dataSet:
        rectKey = 'rectifier'
    elif setPC.nVs > 1:
        rectKey = fluxKey
    else:
        rectKey = None
        print("Only one voltage clamp value found [{}] - fixing parameters of f(v): ".format(setPC.Vs[0]), end='')
        
    if rectKey is not None:
        if verbose > 0:
            print('Rectifier protocol found: fitting E, v0 and v1 for f(v): ', end='')
        phiMax, phiIndMax = getExt(dataSet[rectKey].phis, 'max')
        IssSet, VsSet = dataSet[rectKey].getSteadyStates(run=0, phiInd=phiIndMax)
        if params['E'].vary or params['v0'].vary or params['v1'].vary:
            params = fitfV(VsSet, IssSet, params, relaxFact=relaxFact, verbose=verbose)

    params['E'].vary = False
    params['v0'].vary = False
    params['v1'].vary = False
    
    print('E = {} mV; v0 = {} mV**-1; v1 = {} mV**-1'.format(params['E'].value, params['v0'].value, params['v1'].value))
    
    
    ### Find most extreme peak current in the fluxSet: Ipmax 
    Ipmax, (rmax, pmax, vmax) = setPC.getIpmax(vIndm70)
    Vpmax = setPC.trials[rmax][pmax][vmax].V
    peakKey = fluxKey
    
    if 'delta' in dataSet:
        if isinstance(dataSet['delta'], ProtocolData):
            Ipsat, (rsat, psat, vsat) = dataSet['delta'].getIpmax()
            # try: 
                # vIndSat = dataSet['delta'].Vs.index(-70)
            # except:
                # vIndSat = np.searchsorted(dataSet['delta'].Vs, -70)
            Vsat = dataSet['delta'].trials[rsat][psat][vsat].V
        elif isinstance(dataSet['delta'], PhotoCurrent):
            Ipsat = dataSet['delta'].peak_
            Vsat = dataSet['delta'].V
        else:
            warnings.warn("Unknown data type for 'delta'!")
        
        ### Use peak from 'delta' instead of fluxSet
        if abs(Ipsat) > abs(Ipmax) and np.isclose(Vsat, -70): ##### Reconsider safeguard for choosing V=-70
            Ipmax = Ipsat
            Vpmax = Vsat
            peakKey = 'delta'
    
    print("Estimating g0 from '{}'; Ipmax = {:.3} nA: ".format(peakKey, Ipmax), end='')
    
    ### Maximum conductance: g0
    assert(Vpmax != params['E'].value)
    g0 = 1e6 * Ipmax / (Vpmax - params['E'].value)
    params['g0'].value = g0
    print('g0 = {} pS'.format(round_sig(g0, sig=3)))
    
    
    
    ### Peak recovery: Gr, Gr0, Gr_dark, a6 ### This currently fits to the first flux
    
    ### 2. Fit exponential to peak recovery plots
    if 'recovery' in dataSet and params['Gr0'].vary:
        if verbose > 0:
            print('Recovery protocol found, fitting dark recovery rate: ', end='')
        t_peaks, I_peaks, Ipeak0, Iss0 = getRecoveryPeaks(dataSet['recovery'])
        params = fitRecovery(t_peaks, I_peaks, params, Ipeak0, Iss0, ax=None, method=method, verbose=verbose)
    else:
        if verbose > 0:
            print('Recovery protocol not found, fixing initial value: ', end='')
    params['Gr0'].vary = False
    print('Gr0 = {} ms**-1'.format(params['Gr0'].value))
    
    
    if 'shortPulse' in dataSet:
        quickSet = dataSet['shortPulse'] # Override delta
    elif 'delta' in dataSet:
        quickSet = dataSet['delta']
    else:
        q = 0
        onD = dataSet[fluxKey].trials[q][0][0].onDs[0]
        qI = dataSet[fluxKey].trials[q][0][0]
        for run in range(1, setPC.nRuns): # Skip the first run
            if setPC.trials[q][0][0].onDs[0] < onD:
                q = run
                qI = setPC.trials[q][0][0]
        quickSet = ProtocolData(setPC.trials[q][0][0], nRuns=1, phis=[qI.phi], Vs=[qI.V])
    

    fitParams = params ###################################### Revise  ### Change units handling!!!
    
    for p in nonOptParams:
        params[p].vary = False
    
    if verbose > 0:
        if nPhis > 1:
            print('\nFitting over {} flux values [{:.3g}, {:.3g}] at {} mV (run {}) '.format(nPhis, min(phis), max(phis), setPC.trials[runInd][0][vIndm70].V, runInd), end='')
            print("{{nRuns={}, nPhis={}, nVs={}}}".format(nRuns, nPhis, nVs))
        else:
            print("\nOnly one flux value found [{}] at {} mV - fixing parameters of light-sensitive transitions. ".format(setPC.phis[0], setPC.trials[runInd][0][vIndm70].V))
    
    
    if nStates == 3:
        #phiFits[phiInd] = fit3states(I,t,onInd,offInd,phi,V,Gr0,gmax,Ipmax,params=pOns,method=method)#,Iss)
        fittedParams = fit3states(setPC, runInd, vIndm70, fitParams, method, verbose)
        constrainedParams = ['Gd']
    elif nStates == 4:
        #phiFits[phiInd] = fit4states(I,t,onInd,offInd,phi,V,Gr0,gmax,params=pOns,method=method)
        fittedParams = fit4states(setPC, runInd, vIndm70, fitParams, method, verbose)
        constrainedParams = ['Gd1', 'Gd2', 'Gf0', 'Gb0']
    elif nStates == 6:
        fittedParams = fit6states(setPC, quickSet, runInd, vIndm70, fitParams, method, verbose)
        constrainedParams = ['Gd1', 'Gd2', 'Gf0', 'Gb0', 'Go1', 'Go2']
        #constrainedParams = ['Go1', 'Go2', 'Gf0', 'Gb0']
        #nonOptParams.append(['Gd1', 'Gd2'])
    else:
        raise Exception('Invalid choice for nStates: {}!'.format(nStates))
    
    
    if postFitOpt: # Relax all parameters (except nonOptParams) and reoptimise
        PCs = [setPC.trials[runInd][phiInd][vIndm70] for phiInd in range(setPC.nPhis)]
        Icycles = [pc.getCycle()[0] for pc in PCs]
        nfs = [pc.I[pc.pulseInds[0,1]] for pc in PCs]
        tons = [pc.getOnPhase()[1] for pc in PCs]
        toffs = [pc.getOffPhase()[1] for pc in PCs]
        Vs = [pc.V for pc in PCs]
        phis = [pc.phi for pc in PCs]
        
        if postFitOptMethod is None:
            postFitOptMethod = method
        
        if verbose > 1:
            print("\n\nPerforming post-fit optimisation with the '{}' algorithm [relaxFact={}]!".format(postFitOptMethod, relaxFact))
        
        assert(relaxFact >= 1)
        
        for p in constrainedParams:
            setBounds(fittedParams[p], relaxFact)

        for p in fittedParams:
            if p not in nonOptParams:
                fittedParams[p].vary = True
            
        RhO = models[str(nStates)]()
        postPmin = minimize(errCycle, fittedParams, args=(Icycles,tons,toffs,nfs,RhO,Vs,phis), method=postFitOptMethod)
        optParams = postPmin.params
        
        if verbose > 0:
            reportFit(postPmin, "Post-fit optimisation report for the {}-state model".format(nStates), postFitOptMethod)
        
        # Create new Parameters object to ensure the default ordering
        orderedParams = Parameters()
        for p in params:
            copyParam(p, postPmin.params, orderedParams)
        
    else:
        # Create new Parameters object to ensure the default ordering
        orderedParams = Parameters()
        for p in params:
            copyParam(p, fittedParams, orderedParams)
    
    
    for trial in range(len(PCs)):
        plotFit(PCs[trial], nStates, orderedParams, fitRates=False, index=trial)#, postPmin, fitRates=False, index=trial)
    
    exportName = 'fitted{}sParams.pkl'.format(nStates)
    with open(os.path.join(dDir, exportName), "wb") as fh:
        pickle.dump(orderedParams, fh)
    
    # Plot set of curves
    plotFluxSetFits(fluxSet=setPC, nStates=nStates, params=orderedParams)
        
    if verbose > 0:
        print('')
        printParams(orderedParams)
        if verbose > 1:
            compareParams(params, orderedParams)
        print("\nParameters fit for the {}-state model in {:.3g}s".format(nStates, wallTime() - t0))
        print("--------------------------------------------------------------------------------\n")
    
    return orderedParams
    
    
def plotFluxSetFits(fluxSet, nStates, params, runInd=0, vInd=0):
    """Plot a (list of) model(s) to experimental data"""
    # Currently assumes square light pulses, all of the same duration
    # TODO: Generalise to handle multiple pulses of different durations
    from pyrho.config import colours
    
    setFig = plt.figure()
    setAx = setFig.add_subplot(111)
    lineWidth = 2.5 * mpl.rcParams['lines.linewidth']
    
    # Plot experimental data
    for phiInd in range(fluxSet.nPhis):
        PC = fluxSet.trials[runInd][phiInd][vInd]
        PC.alignToTime()
        phi = PC.phi
        setAx.plot(PC.t, PC.I, color=colours[phiInd%len(colours)], lw=lineWidth, markeredgecolor='None', 
        label="$\phi = {:.3g}\ \mathrm{{[ph. \cdot mm^{{-2}} \cdot s^{{-1}}]}}$".format(phi))
    
    if not isinstance(nStates, (list, tuple)):
        nStates = [nStates]
    
    if not isinstance(params, (list, tuple)):
        params = [params]
    
    for mInd, (s, pSet) in enumerate(zip(nStates, params)):
        RhO = models[str(s)]()
        RhO.updateParams(pSet)
        # see ProtocolData.plot()
        
        for phiInd in range(fluxSet.nPhis):
            PC = fluxSet.trials[runInd][phiInd][vInd]
            phi = PC.phi
            V = PC.V
            
            #TODO: Replace with calcCycle, runTrial or similar
            
            #I_RhO, t, soln = Sim.runTrial(RhO, phi, V, delD, cycles, self.dt, verbose) #self.totT, 
            
            ## Delay phase
            Idel, tdel  = PC.getDelayPhase()#;   tdel -= tdel[0]
            RhO.setLight(RhO.phi_0)            
            if RhO.useAnalyticSoln:
                soln = RhO.calcSoln(tdel, RhO.s_0)
            else:
                soln = odeint(RhO.solveStates, RhO.s_0, tdel, Dfun=RhO.jacobian)
            RhO.storeStates(soln[1:], tdel[1:])

            for p in range(PC.nPulses):
                ## On phase
                Ion, ton    = PC.getOnPhase(p)#;     ton -= ton[0]
                RhO.setLight(phi) # Calculate transition rates for phi
                RhO.s_on = soln[-1,:]
                if RhO.useAnalyticSoln:
                    soln = RhO.calcSoln(ton, RhO.s_on)
                else:
                    soln = odeint(RhO.solveStates, RhO.s_on, ton, Dfun=RhO.jacobian)
                RhO.storeStates(soln[1:], ton[1:])
                
                ## Off phase
                Ioff, toff  = PC.getOffPhase(p)#;    toff -= toff[0]
                RhO.setLight(0)
                RhO.s_off = soln[-1,:]
                if RhO.useAnalyticSoln:
                    soln = RhO.calcSoln(toff, RhO.s_off)
                else:
                    soln = odeint(RhO.solveStates, RhO.s_off, toff, Dfun=RhO.jacobian)
                RhO.storeStates(soln[1:], toff[1:])
                
            states, tfit = RhO.getStates()
            Ifit = RhO.calcI(V, states)
            
            # Plot model fit curve
            if len(nStates) > 1 and phiInd == fluxSet.nPhis-1:
                plt.plot(tfit, Ifit, color='k', ls=styles[mInd % len(styles)], label='$\mathrm{{Model\ fit\ ({}-states)}}$'.format(stateLabs[s]))
            else:
                plt.plot(tfit, Ifit, color='k', ls=styles[mInd % len(styles)])#c=colours[phiInd % len(colours)]) #,  label='$\phi={:.3g}$'.format(phi), label='$\mathrm{{Model\ fit\ ({}-states)}}$'.format(nStates))
            
    if config.addTitles: # and p == 0:
        if len(nStates) == 1:
            mString = '{}-state model'.format(stateLabs[s])
        else:
            states = ', '.join(str(stateLabs[st]) for st in nStates)
            mString = '{{{}}}-state models'.format(states)
        plt.title('{} fit to {} light intensities [{:.3g}, {:.3g}]'.format(mString, 
                            fluxSet.nPhis, min(fluxSet.phis), max(fluxSet.phis)))
                
    plt.legend(loc='best')
    plt.xlabel('$\mathrm{Time\ [ms]}$', position=(xLabelPos,0), ha='right') # (xpos,ypos) ypos is ignored
    plt.ylabel('$\mathrm{Photocurrent\ [nA]}$')
    
    # Assume square light pulses, all of the same duration
    for p in range(PC.nPulses):
        t_on, t_off = PC.pulses[p]
        plt.axvspan(t_on, t_off, facecolor='y', alpha=0.2)
    
    #setAx.spines['left'].set_position('zero') # y-axis
    setAx.spines['right'].set_color('none')
    setAx.spines['bottom'].set_position('zero') # x-axis
    setAx.spines['top'].set_color('none')
    setAx.spines['left'].set_smart_bounds(True)
    setAx.spines['bottom'].set_smart_bounds(True)
    setAx.xaxis.set_ticks_position('bottom')
    setAx.yaxis.set_ticks_position('left')
    plt.tight_layout()
    
    setFig.savefig(fDir+'fluxSetFit'+'-'.join(str(s) for s in nStates)+'states'+'.'+config.saveFigFormat, format=config.saveFigFormat)
    
    return


def setBounds(param, relaxFact=2):
    #s = np.sign(param.value)
    #param.min = param.value * relaxFact**-s
    #param.max = param.value * relaxFact**-s
    
    if param.value > 0:
        param.min = param.value / relaxFact #round_sig(param.value / relaxFact, sig=3)
        param.max = param.value * relaxFact #round_sig(param.value * relaxFact, sig=3)
    elif param.value < 0:
        param.min = param.value * relaxFact #round_sig(param.value * relaxFact, sig=3)
        param.max = param.value / relaxFact #round_sig(param.value / relaxFact, sig=3)
    else: #param.value == 0:
        param.min = -1e-9
        param.max = 1e-9
    