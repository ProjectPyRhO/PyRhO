
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp # For plotStates
from scipy.signal import * #for argrelextrema
from scipy.integrate import odeint
import warnings

from .parameters import *
#print(vars()['verbose'])
from .config import verbose, addTitles, saveFigFormat, fDir # For plotStates()
#from __init__ import verbose

# verbose = 0
#global verbose # global statement is so that subfunctions can assign to the global var

#global phi0
#global E
#global A

#addTitles = True

### Model Parameters ###

### Physical constants
h = 6.6260695729e-34   # Planck's constant (Js)
c = 2.99792458e8     # Speed of light    (m*s^-1)
#NA = 6.0221413e23    # Avogadro's Number (mol^-1)


'''
    == A few notes about colour ==

    Color   Wavelength(nm) Frequency(THz)
    Red     620-750        484-400
    Orange  590-620        508-484
    Yellow  570-590        526-508
    Green   495-570        606-526
    Blue    450-495        668-606
    Violet  380-450        789-668

    f is frequency (cycles per second)
    l (lambda) is wavelength (meters per cycle)
    e is energy (Joules)
    h (Plank's constant) = 6.6260695729 x 10^-34 Joule*seconds
                         = 6.6260695729 x 10^-34 m^2*kg/seconds
    c = 299792458 meters per second
    f = c/l
    l = c/f
    e = h*f
    e = c*h/l

    List of peak frequency responses for each type of 
    photoreceptor cell in the human eye:
        S cone: 437 nm
        M cone: 533 nm
        L cone: 564 nm
        rod:    550 nm in bright daylight, 498 nm when dark adapted. 
                Rods adapt to low light conditions by becoming more sensitive.
                Peak frequency response shifts to 498 nm.

'''

#import sys
#import os
#import traceback
#import optparse
#import time
#import logging


def lam2rgb(wav, gamma=0.8, output='norm'):
    """This converts a given wavelength of light to an 
    approximate RGB colour value with edge attenuation. 
    The wavelength must be given in nanometres in the 
    range from 380 nm - 750 nm (789 THz - 400 THz).
    
    Adapted from: http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    """
    
    
    wav = float(wav)
    if wav >= 380 and wav < 440:
        attenuation = 0.3 + 0.7 * (wav - 380) / (440 - 380)
        R = ((-(wav - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wav >= 440 and wav < 490:
        R = 0.0
        G = ((wav - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wav >= 490 and wav < 510:
        R = 0.0
        G = 1.0
        B = (-(wav - 510) / (510 - 490)) ** gamma
    elif wav >= 510 and wav < 580:
        R = ((wav - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wav >= 580 and wav < 645:
        R = 1.0
        G = (-(wav - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wav >= 645 and wav <= 750:
        attenuation = 0.3 + 0.7 * (750 - wav) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else: # Outside the visible spectrum
        R = 0.0
        G = 0.0
        B = 0.0
    
    if output == 'norm':
        return (R,G,B)
    elif output == 'hex':
        R *= 255
        R = max(0, min(round(R), 255))
        G *= 255
        G = max(0, min(round(G), 255))
        B *= 255
        B = max(0, min(round(B), 255))
        #return (int(R), int(G), int(B)) # int() truncates towards 0
        return "#{0:02x}{1:02x}{2:02x}".format(R,G,B), (R,G,B)

### Default parameter dictionaries

#d3sp = {'E':E, 'k':8.2, 'Gd':0.1, 'Gr0':1/5000, 'Gr1':0.016, 'phi0':1e15, 'phiSat':1e20, 'g':1.67e4}
# E         [mV]    Channel reversal potential
# k         [ms^-1] Quantum efficiency * number of photons absorbed by a RhO molecule per unit time
# eF        [ms^-1] Quantum efficiency * number of photons absorbed by a ChR2 molecule per unit time
# Gd        [ms^-1] @ 1mW mm^-2
# Gr_dark   [ms^-1] tau_r,dark = 5-10s p405 Nikolic et al. 2009
# Gr_light  [ms^-1] (Gr,light @ 1mW mm^-2)
# phi0 
# phiSat
# g         [nS]


### Output information
#verbose = 0


### Model functions ###

def irrad2flux(E,lam=470):   # E2phi
    """Converts irradiance [mW * mm^-2] and wavelength (default: 470) [nm] to flux [photons * s^-1 * mm^-2]"""
    Ep = 1e12 * h * c/lam    # Energy per photon [mJ] (using lambda in [nm])
    return E/Ep              # Photon flux (phi) scaled to [photons * s^-1 * mm^-2]


def flux2irrad(phi,lam=470):
    """Converts flux [photons * s^-1 * mm^-2] and wavelength (default: 470) [nm] to irradiance [mW * mm^-2]"""
    Ep = 1e12 * h * c/lam    # Energy per photon [mJ] (using lambda in [nm])
    return phi * Ep          # Irradiance (E) scaled to [mW * mm^-2]


def calcgbar(Ip,Vclamp,A):
    """Function to calculate a lower bound on the cell's maximum conductance from its peak current
    Ip      :=  Peak current [nA]
    Vclamp  :=  Clamp Voltage [mV]
    A       :=  Cell surface area [um^2]
    return gbar [pS/um^2]"""
    Gmax = Ip/Vclamp  # Maximum conductance for the whole cell
    gbar = Gmax/A     # Maximum conductance pS / um^2
    return gbar * (1e6) # 1e-12 S / (1e-6 m)^2 = (1e-6)*(1e-9 A / 1e-3 V)/(1e-6 m)^2 # Check the conversion factor

    
    
def times2cycles(times,totT):
    """Input    times:= [t_on,t_off], t_tot
       Output   cycles:= [onD,offD], delD"""
    times = np.array(times)
    #print("Times: ",times)
    nPulses = times.shape[0]
    assert(times.shape[1] <= 2)
    #delDs = [row[0] for row in times] # pulses[:,0]    # Delay Durations
    delD = times[0,0]
    onDs = [row[1]-row[0] for row in times] # pulses[:,1] - pulses[:,0]   # Pulse Durations
    offDs = np.append(times[1:,0],totT) - times[:,1]
    cycles = np.vstack((onDs,offDs)).transpose()
    return (cycles, delD)

def cycles2times(cycles,delD):
    """Function to convert pulse cycles to times. 
        Input    cycles:= [onD,offD], delD
        Output   times:= [t_on,t_off], totT"""
    
    ### Generalise to delDs c.f. varyIPI
        
    cycles = np.array(cycles)
    #print("Cycles: ",cycles)
    nPulses = cycles.shape[0]
    #print(nPulses)
    assert(cycles.shape[1] <= 2)
    times = np.zeros((nPulses, 2)) #[delD,delD+cycles[row,0] for row in pulses]
    lapsed = delD
    for p in range(nPulses):
        times[p,0] = lapsed
        times[p,1] = lapsed+cycles[p,0]
        lapsed += sum(cycles[p,:])
    return (times, lapsed)

def plotLight(times, ax=None, light='shade', lam=470, alpha=0.2): #=plt.gcf()
    """Function to plot light pulse(s)
    times   = [[t_on, t_off]...]
    ax      = Axes to plot on (default: gca()
    light   = Representation type: {'shade', 'borders', 'greyscale', 'hatch', 'spectral'}. Default: 'shade'
    lam     = Wavelength [nm] (default: 470)
    alpha   = Transparency (default: 0.2)"""
    
    if ax==None:
        ax=plt.gca()
    else:
        plt.sca(ax)
    nPulses = times.shape[0]
    if light == 'shade':
        for p in range(nPulses):
            plt.axvspan(times[p][0],times[p][1],facecolor='y',alpha=alpha)    
    elif light == 'borders': 
        for p in range(0, nPulses): ############################ Move to plotLight
            plt.axvline(x=times[p][0],linestyle='--',color='k')
            plt.axvline(x=times[p][1],linestyle='--',color='k')
    elif light == 'greyscale':
        # Set background to grey and illumination to white
        ax.set_axis_bgcolor('0.3')
        plt.axvspan(times[p][0],times[p][1],facecolor='w')#,alpha=alpha)  
    elif light == 'hatch':
        plt.axvspan(times[p][0],times[p][1],hatch='*')
    elif light == 'spectral':
        # Plot the colour corresponding to the wavelength
        if 380 <= lam <= 750:
            rgb = lam2rgb(lam)
            for p in range(0, nPulses): 
                plt.axvspan(times[p][0],times[p][1],facecolor=rgb,alpha=alpha)
        else: # Light is not in the visible spectrum - plot borders instead
            for p in range(0, nPulses): 
                plt.axvline(x=times[p][0],linestyle='--',color='k')
                plt.axvline(x=times[p][1],linestyle='--',color='k')
                plt.axvspan(times[p][0],times[p][1],hatch='/') #'*'
    elif light == 'None' or light == None:
        pass
    else:
        print(light)
        warnings.warn('Warning: Unrecognised light representation!')
    return
    
def round_sig(x, sig=3):
    """Round a float to n significant digits (default is 3). """
    if abs(x) == 0 or np.isinf(x) or np.isnan(x):
        return x
    else:
        return round(x, sig-int(np.floor(np.log10(abs(x))))-1)


def calcIssfromfV(V,v0,v1,E):#,G): # Added E as another parameter to fit
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
    
# def fitfV(Vs, Iss, curveFunc, p0, RhO, fig=None):#, eqString): =plt.gcf()
    # if fig==None:
        # fig=plt.gcf()
    # markerSize=40
    # eqString = r'$f(V) = \frac{{{:.3}}}{{V-{:+.2f}}} \cdot \left[1-\exp\left({{-\frac{{V-{:+.2f}}}{{{:.3}}}}}\right)\right]$'
    # psi = RhO.calcPsi(RhO.steadyStates)
    # #sf = RhO.A * RhO.gbar * psi * 1e-6 # Six-state only
    # sf = RhO.g * psi * 1e-6 
    # fVs = np.asarray(Iss)/sf # np.asarray is not needed for the six-state model!!!
    # popt, pcov = curve_fit(curveFunc, Vs, fVs, p0=p0) # (curveFunc, Vs, Iss, p0=p0)
    # pFit = [round_sig(p,3) for p in popt]
    # peakEq = eqString.format(pFit[0],pFit[2],pFit[2],pFit[1])
    
    # Vrange = max(Vs)-min(Vs)
    # xfit=np.linspace(min(Vs),max(Vs),Vrange/.1) #Prot.dt
    # yfit=curveFunc(xfit,*popt)*sf
    
    # #peakEq = eqString.format(*[round_sig(p,3) for p in popt])
    
    # fig.plot(xfit,yfit)#,label=peakEq)#,linestyle=':', color='#aaaaaa')
    # #col, = getLineProps(Prot, 0, 0, 0) #Prot, run, vInd, phiInd
    # #plt.plot(Vs,Iss,linestyle='',marker='x',color=col)
    # fig.scatter(Vs,Iss,marker='x',color=colours,s=markerSize)#,linestyle=''
    
    # # x = 1 #0.8*max(Vs)
    # # y = 1.2*yfit[-1]#max(IssVals[run][phiInd][:])
    # # plt.text(-0.8*min(Vs),y,peakEq,ha='right',va='bottom',fontsize=eqSize)#,transform=ax.transAxes)
    
    # if verbose:
        # print(peakEq)
    # return popt, pcov, peakEq


# def runTrial(RhO, nPulses, V,phiOn,delD,onD,offD,padD,dt): #dt; a1,a3,b2,b4,I_RhO
    # """Main routine for simulating a pulse train"""
    
    # if verbose >= 0:
        # print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse: [delD={:.4g}ms; onD={:.4g}ms; offD={:.4g}ms]".format(V,phiOn,delD,onD,offD+padD))
    
    # ### Delay phase (to allow the system to settle)
    # phi = 0
    # RhO.initStates(phi) # Reset state and time arrays from previous runs
    # start, end = 0.00, delD
    # t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
    # soln = odeint(RhO.solveStates, RhO.s_0, t_del, args=(None,), Dfun=RhO.jacobian) #delay
    # t = t_del
    # RhO.storeStates(soln,t)
    
    # for p in range(0, nPulses):
        
        # ### Light on phase
        # RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
        # start = end
        # end = start + onD
        # t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
        # onInd = len(RhO.t) - 1  # Start of on-phase
        # offInd = onInd + len(t) - 1 # Start of off-phase
        # RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
        # # Turn on light and set transition rates
        # phi = phiOn  # Light flux
        # RhO.setLight(phi)
        # if verbose > 1:
            # print("On-phase initial conditions:{}".format(RhO.s_on))
        # soln = odeint(RhO.solveStates, RhO.s_on, t, args=(None,), Dfun=RhO.jacobian)
        # RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
        
        # ### Light off phase
        # RhO.s_off = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
        # start = end
        # end = start + offD
        # if (p+1) == nPulses: # Add (or subtract) extra time after (during) the off phase
            # end += padD
        # t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # endpoint=True
        # # Turn off light and set transition rates
        # phi = 0  # Light flux
        # RhO.setLight(phi)
        # if verbose > 1:
            # print("Off-phase initial conditions:{}".format(RhO.s_off))
        # soln = odeint(RhO.solveStates, RhO.s_off, t, args=(None,), Dfun=RhO.jacobian)
        # RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
        
    # ### Calculate photocurrent
    # I_RhO = RhO.calcPhotocurrent(V, RhO.states)
    # states,t = RhO.getStates()
    
    
    
    # return I_RhO, t, states
    

# def runTrialPhi_t(RhO,V,phi_t,delD,stimD,totT,dt): #dt; a1,a3,b2,b4,I_RhO
    # """Main routine for simulating a pulse train"""
    # # Add interpolation of values for phi(t) to initialisation phi_t = interp1d(t,sin(w*t),kind='cubic')
    # #print("Simulating experiment at V = {:+}mV, phi = {:.3g}photons/s/mm^2, pulse = {}ms".format(V,phiOn,onD))
    
    # ### delD and stimD are only used for finding pulse indexes - could be removed along with separate delay phase?!!!
    # ### Combine this with the original runTrial
    
    # ### Delay phase (to allow the system to settle)
    # phi = 0
    # RhO.initStates(phi) # Reset state and time arrays from previous runs
    # start, end = 0.00, delD
    # t_del = np.linspace(start,end,((end-start)/dt)+1, endpoint=True) # Time vector
    # soln = odeint(RhO.solveStates, RhO.s_0, t_del, args=(None,), Dfun=RhO.jacobian) #delay # odeint(RhO.solveStates, RhO.s_0, t_del, Dfun=RhO.jacobian)
    # t = t_del
    # RhO.storeStates(soln,t)

    # ### Stimulation phase
    # RhO.s_on = soln[-1,:] # [soln[-1,0], soln[-1,1], soln[-1,2], soln[-1,3], soln[-1,4], soln[-1,5]]
    # start = end
    # end = totT #start + stimD
    # t = np.linspace(start,end,((end-start)/dt)+1, endpoint=True)
    # onInd = len(RhO.t) - 1  # Start of on-phase
    # #offInd = onInd + len(t) - 1 # Start of off-phase
    # offInd = onInd + int(round(stimD/dt)) # Consider rounding issues...
    # RhO.pulseInd = np.vstack((RhO.pulseInd,[onInd,offInd]))
    # # Turn on light and set transition rates
    # #phi = phiOn  # Light flux
    # #RhO.setLight(phi)
    # if verbose > 1:
        # print("On-phase initial conditions:{}".format(RhO.s_on))
        # #print('Pulse = [{}, {}]'.format(onInd,offInd))
    # #print(RhO.pulseInd)
    
    # #soln, out = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian, full_output=True)
    # #print(out)
    # soln = odeint(RhO.solveStates, RhO.s_on, t, args=(phi_t,), Dfun=RhO.jacobian)
    # RhO.storeStates(soln[1:],t[1:]) # Skip first values to prevent duplicating initial conditions and times
    # if verbose > 1:
        # print('Pulse = [{}, {}]'.format(RhO.t[onInd],RhO.t[offInd]))
    
    # ### Calculate photocurrent
    # I_RhO = RhO.calcPhotocurrent(V, RhO.states)
    # states,t = RhO.getStates()
    
    # # if verbose:
        # # print("Run={}, I_min={}, I_max={}, label={}".format(run,I_RhO.min(),I_RhO.max(),label))
    
    # return I_RhO, t, states

def expDecay(t, a, b, c):
    return a * np.exp(-t/b) + c


def biExpDecay(t, a1, tau1, a2, tau2, I_ss):
    return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2) + I_ss


# def fitPeaks(t_peaks, I_peaks, curveFunc, p0, eqString, fig=None):
    
    # shift = t_peaks[0] # ~ delD
# #     if protocol == 'varyIPI':
# #         plt.ylim(ax.get_ylim()) # Prevent automatic rescaling of y-axis
    # popt, pcov = curve_fit(curveFunc, t_peaks-shift, I_peaks, p0=p0) #Needs ball-park guesses (0.3, 125, 0.5)
    # peakEq = eqString.format(*[round_sig(p,3) for p in popt]) # *popt rounded to 3s.f.
    
    # if fig:
        # plt.figure(fig.number) # Select figure
# #     ext = 10 # Extend for ext ms either side
# #     xspan = t_peaks[-1] - t_peaks[0] + 2*ext 
# #     xfit=np.linspace(t_peaks[0]-ext-shift,t_peaks[-1]+ext-shift,xspan/dt)
        # plt.plot(t_peaks, I_peaks, linestyle='', color='r', marker='*')
        # xfit=np.linspace(-shift,Prot.totT-shift,Prot.totT/Prot.dt) #totT
        # yfit=curveFunc(xfit,*popt)
        
        # plt.plot(xfit+shift,yfit,linestyle=':',color='#aaaaaa')#,label="$v={:+} \mathrm{{mV}}$, $\phi={:.3g}$".format(V,phiOn))
        # #ylower = copysign(1.0,I_peaks.min())*ceil(abs((I_peaks.min()*10**ceil(abs(log10(abs(I_peaks.min())))))))/10**ceil(abs(log10(abs(I_peaks.min()))))
        # #yupper = copysign(1.0,I_peaks.max())*ceil(abs((I_peaks.max()*10**ceil(abs(log10(abs(I_peaks.max())))))))/10**ceil(abs(log10(abs(I_peaks.max()))))
    # #     if (len(Vs) == 1) and (len(phis) == 1) and (nRuns == 1):
    # #         x, y = 0.8, 0.9
    # #     else:
        # x = 0.8
        # y = yfit[-1] #popt[2]
        
        # plt.text(x*Prot.totT,y,peakEq,ha='center',va='bottom',fontsize=eqSize) #, transform=ax.transAxes)
    
    # print(peakEq)
    # if verbose > 1:
        # print("Parameters: {}".format(popt))
        # if type(pcov) in (tuple, list):
            # print("$\sigma$: {}".format(np.sqrt(pcov.diagonal())))
        # else:
            # print("Covariance: {}".format(pcov))
    # return popt, pcov, peakEq


# def getLineProps(Prot, run, vInd, phiInd):
    # if verbose > 1 and (Prot.nRuns>len(colours) or len(Prot.phis)>len(colours) or len(Prot.Vs)>len(colours)):
        # print("Warning: only {} line colours are available!".format(len(colours)))
    # if Prot.nRuns>1 and len(Prot.phis)>1 and len(Prot.Vs)>1:
        # print("Warning: Too many changing variables for one plot!")
    # if verbose:
        # print("Run={}/{}; phiInd={}/{}; vInd={}/{}".format(run,nRuns,phiInd,len(phis),vInd,len(Vs)))
    # if Prot.nRuns>1:
        # col = colours[run%len(colours)]
        # if len(Prot.phis)>1:
            # style=styles[phiInd%len(styles)]
        # elif len(Prot.Vs)>1:
            # style=styles[vInd%len(styles)]
        # else:
            # style = '-'
    # else:
        # if len(Prot.Vs)>1:
            # col = colours[vInd%len(colours)]
            # if len(Prot.phis)>1:
                # style = styles[phiInd%len(styles)]
            # else:
                # style = '-'
        # else:
            # if len(Prot.phis)>1:
                # col = colours[phiInd%len(colours)]
                # style = '-'
            # else:
                # col = 'b'   ### colours[0]
                # style = '-' ### styles[0]
    # return col, style




def findPeaks(I_phi,startInd=0,extOrder=5): # Revise this to take care of V<>E 
# findPeaks(I_phi,minmax,startInd=0,extOrder=5)
    
    if verbose > 1:
        print("findPeaks(): extOrder = ",extOrder)
    if np.isnan(I_phi).all() or len(I_phi)==0: # Not necessary
        warning.warn("Warning: No photocurrent for run[{}]phi[{}]v[{}]!".format(run,phiInd,vInd))
        #print("Warning: No photocurrent for run[{}]phi[{}]v[{}]!".format(run,phiInd,vInd))
        return []
    else:
        if abs(min(I_phi)) > abs(max(I_phi)): # Find Minima
            minmax = np.less #Ipmax = min(photocurrent.I_phi)
        else:       # Find Maxima
            minmax = np.greater #Ipmax = max(photocurrent.I_phi)
#         if (extOrder < 1):
#             if (V < E): # Find Minima
#                 Ip = min(I_phi)
#                 peakInds = [np.argmin(I_phi) + startInd]
#             else:       # Find Maxima
#                 Ip = max(I_phi)
#                 peakInds = [np.argmax(I_phi) + startInd]
#         else:
        # if (V < E): # Find Minima
            # peakInds = argrelextrema(I_phi[startInd:], np.less, order=extOrder)
        # else:       # Find Maxima
            # peakInds = argrelextrema(I_phi[startInd:], np.greater, order=extOrder)
        peakInds = argrelextrema(I_phi[startInd:], minmax, order=extOrder)
        peakInds = peakInds[0] + startInd # Remove tuple and transform into original co-ordinates
        peakInds = peakInds.astype(int)
#         Ip = I_phi[peakInds]

        ### Alternatively loop over pulses taking the min() or max() within each region

        if peakInds.any() and verbose > 1:
            print(peakInds)
            print("I[p-2]={}; I[p-1]={}; I[p]={}; I[p+1]={}; I[p+2]={};".format(*I_phi[peakInds-2:peakInds+3]))
        return peakInds #, Ip        

        
def findPlateauCurrent(I_phi,p=0.1): ### c.f. findPlateaus() in loadData.py
    """Trim the off-phase of I_phi. Optionally pass the proportion of data for the averaging window, p. """
    # Move a sliding window over the data and look for abs(dI/dt) < eps
    # Take the mean of the last p*t ms before t_off
    # Take the mean of the last n ms before t_off
    #onEndInd = np.searchsorted(t,onD+delD,side="left")
    windowInd = int(round(p*len(I_phi))) #np.searchsorted(t,t[onEndInd]-100,side="left") # Generalise
    I_ss = np.mean(I_phi[-windowInd:]) #:onEndInd+1
    # Fit an exponential from the t_peak to t_off
    # searchSlice = slice(peakInds[0], offInd+1)
    # popt = fitPeaks(t[peakInds[0]:offInd+1], I_phi[peakInds[0]:offInd+1], expDecay, p0inact, '$I_{{inact}} = {:.3}e^{{-t/{:g}}} {:+.3}$','')
    # Iss = popt[2]
    
    # Alternatively use a gradient-based approach
    #if abs(np.gradient(I_phi)) < tol:
    
    return I_ss


###### Model class definitions ######

def selectModel(nStates):
    """Model selection function"""
    if nStates == 3:
        return RhO_3states() #(E)
    elif nStates == 4:
        return RhO_4states() #(E,gam,phi0)
    elif nStates == 6:
        return RhO_6states() #(E,gam,phi0,A)
    else:
        print("Error in selecting model - please choose from 3, 4 or 6 states")
        raise NotImplementedError(nStates)
        

class OGmodel(object):
    "Common base class for all models"
    
    phi = 0.0  # Instantaneous Light flux
#     phi0 = 1e14 # Normalising light intensity parameter [photons * s^-1 * mm^-2] ~ 1/10 of threshold
    
    # Solver parameters
#     dt = 0.1 #0.01 #0.1 #0.001
    
    def __init__(self, nStates): #, rhoType='Rhodopsin'
        self.nStates = nStates
        #self.rhoType = rhoType # E.g. 'ChR2' or 'ArchT'
        
    def __str__(self):
        return self.rhoType+" {} state model (phi={:.3g})".format(self.nStates,self.phi) # Display transition rates?
        #self.__name__+

    def setParams(self, params):
        #for param, value in params.items():
        for p in params.keys():
            self.__dict__[p] = params[p].value #vars(self())[p]

    def exportParams(self, params):
        """Export parameters to lmfit dictionary"""
        for p in self.__dict__.keys():
            params[p].value = self.__dict__[p]
            
    def printParams(self):
        for p in self.__dict__.keys():
            print(p,' = ',self.__dict__[p])
    
    def storeStates(self,soln,t):#,pulseInds):
        self.states = np.append(self.states, soln, axis=0)
        self.t = np.append(self.t, t, axis=1)
        #self.pulseInd = np.append(self.pulseInd, pulseInds, axis=0)
    
    def getStates(self):
        """Returns states, t"""
        return self.states, self.t
    
    def initStates(self, phi, s0=None):
        """Clear state arrays and set transition rates"""
        self.states = np.empty([0,self.nStates])
        self.t = []
        self.pulseInd = np.empty([0,2],dtype=int) # Light on and off indexes for each pulse
        self.setLight(phi)
        #if s0 is not None: # Override default initial conditions
        #    self.s_0 = s0
    
    def calcfV(self, V): ############################################################### Finish this!
        if self.useInwardRect: 
            try:
                fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0)) # Dimensionless
            except ZeroDivisionError:
                if np.isscalar(V):
                    if (V == self.E):
                        fV = self.v1/self.v0
                else: #type(fV) in (tuple, list, array)
                    fV[np.isnan(fV)] = self.v1/self.v0 # Fix the error when dividing by zero
        else: 
            fV=1 ### Extend to vector
        return fV
                
    def plotStates(self,t,states,pulses,labels,phiOn=0,peaks=[],name=None): ### Check how to check for optional arguments
        if len(peaks) > 0:
            plotPieCharts = True
        else:
            plotPieCharts = False
        
        figWidth, figHeight = mp.rcParams['figure.figsize']
        fig = plt.figure(figsize=(figWidth, 1.5*figHeight))
        #fig = plt.figure()
        gs = plt.GridSpec(3,3)
        
        totT = t[-1]
        
        # Plot line graph of states
        axLine = fig.add_subplot(gs[0,:])
        plt.plot(t,states)
        sig, = plt.plot(t, np.sum(states,axis=1), color='k', linestyle='--')
        plt.setp(axLine.get_xticklabels(), visible=False)
        plt.ylabel('$\mathrm{State\ occupancy}$')# proportion')
        labelsIncSum = np.append(labels,'$\Sigma s_i$')
        plt.legend(labelsIncSum,loc=6)
        plt.xlim((0,totT)) #plt.xlim((0,delD+(nPulses*(onD+offD)))) # t_run
        plt.ylim((-0.1,1.1))
        if addTitles:
            plt.title('$\mathrm{State\ variables\ through\ time}$') #plt.title('State variables through time: $v={} \mathrm{{mV}},\ \phi={:.3g} \mathrm{{photons}} \cdot \mathrm{{s}}^{{-1}} \cdot \mathrm{{cm}}^{{-2}}$'.format(V,phiOn))
        plotLight(pulses, axLine)

        
        ### Plot stack plot of state variables
        axStack = fig.add_subplot(gs[1,:], sharex=axLine)
        plt.stackplot(t,states.T)
        plt.ylim((0,1))
        plt.xlim((0,totT)) #plt.xlim((0,delD+(nPulses*(onD+offD))))
        plotLight(pulses, axStack, 'borders')
        if addTitles:
            axStack.title.set_visible(False)
        plt.xlabel('$\mathrm{Time\ [ms]}$')
        plt.ylabel('$\mathrm{State\ occupancy}$')# proportion')
        
        if plotPieCharts:
            axS0 = fig.add_subplot(gs[2,0])
            initialStates = self.s_0 * 100
            #print(initialStates,labels)
            if verbose > 1:
                pct = {l:s for l,s in zip(labels,sizes)}
                print('Initial state occupancies (%):',sorted(pct.items(),key=lambda x: labels.index(x[0])))
            patches, texts, autotexts = plt.pie(initialStates, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False) #, explode=explode
            for lab in range(len(labels)):
                texts[lab].set_fontsize(mp.rcParams['ytick.labelsize'])
                autotexts[lab].set_fontsize(mp.rcParams['axes.labelsize'])
            plt.axis('equal')
            if addTitles:
                plt.title('$\mathrm{Initial\ state\ occupancies}$')
            #else:
            #    plt.title('$t_{0}$')
            
            if peaks: ### Plot peak state proportions
                pInd = peaks[0] # Plot the first peak
                axLine.axvline(x=t[pInd],linestyle=':',color='k')
                axStack.axvline(x=t[pInd],linestyle=':',color='k')
                axPeak = fig.add_subplot(gs[2,1])
                sizes = states[pInd,:] * 100
                #sizes = [s*100 for s in sizes]
                #explode = (0,0,0.1,0.1,0,0)
                if verbose > 1:
                    pct = {l:s for l,s in zip(labels,sizes)}
                    print('Peak state occupancies (%):',sorted(pct.items(),key=lambda x: labels.index(x[0])))
                patches, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False)#, explode=explode)
                for lab in range(len(labels)):
                    texts[lab].set_fontsize(mp.rcParams['ytick.labelsize'])
                    autotexts[lab].set_fontsize(mp.rcParams['axes.labelsize'])
                plt.axis('equal')
                if addTitles:
                    plt.title('$\mathrm{Simulated\ peak\ state\ occupancies}$')
                #else:
                #    plt.title('$t_{peak}$')
            
            if phiOn > 0: ### Plot steady state proportions
                axSS = fig.add_subplot(gs[2,2])
                steadyStates = self.calcSteadyState(phiOn) * 100 # Convert array of proportions to %
                #steadyStates = [s*100 for s in steadyStates]
                #explode = (0,0,0.1,0.1,0,0)
                if verbose > 1:
                    pct = {l:s for l,s in zip(labels,sizes)}
                    print('Steady state occupancies (%):',sorted(pct.items(),key=lambda x: labels.index(x[0])))
                patches, texts, autotexts = plt.pie(steadyStates, labels=labels, autopct='%1.1f%%', startangle=90, shadow=False) #, explode=explode
                for lab in range(len(labels)):
                    texts[lab].set_fontsize(mp.rcParams['ytick.labelsize'])
                    autotexts[lab].set_fontsize(mp.rcParams['axes.labelsize'])
                plt.axis('equal')
                if addTitles:
                    plt.title('$\mathrm{Analytic\ steady\ state\ occupancies}$')
                #else:
                #    plt.title('$t_{\inf}$')

            

        plt.tight_layout()
        
        if name is not None:
            plt.savefig(fDir+name+'.'+saveFigFormat, format=saveFigFormat)


class RhO_3states(OGmodel):
    
    # Class attributes
    nStates = 3
    hasAnalyticSoln = True
    labels = ['$C$','$O$','$D$']
    stateVars = ['C','O','D']
    connect = [[0,1,0],
               [0,0,1],
               [1,0,0]]
    conDir  = [[ 0,-1, 1],
               [ 1, 0,-1],
               [-1, 1, 0]]
    equations = """
                \begin{align*}
                \frac{dC}{dt} &= G_rD - \epsilon F C
                \frac{dO}{dt} &= \epsilon FC -G_{d}O
                \frac{dD}{dt} &= G_{d}O-G_{r}D
                C+O+D &= 1
                \epsilon F &= \phi\frac{\epsilon \sigma_{ret}}{w_{loss}}
                G_r &= G_{r,d} + \mathcal{H}(\phi) \cdot G_{r,l}
                I_{\phi} &= \bar{g} O \cdot (v-E)
                \end{align}
                """
                
    eqIss = """$I_{SS} = \bar{g} \cdot \frac{\epsilon F \cdot G_r}{G_d \cdot (G_r + \epsilon F) + \epsilon F \cdot G_r} \cdot (v-E) 
    = \bar{g} \cdot \frac{\tau_d}{\tau_d + \tau_r + \tau_\phi} \cdot (v-E)$"""
    
#     phi0 = 1e14 # Normalising light intensity parameter [photons * s^-1 * mm^-2] ~ 1/10 of threshold
    
    # Solver parameters
#     dt = 0.1 #0.01 #0.1 #0.001    
    
    def __init__(self, params=modelParams['3'], rhoType='Rhodopsin'): # E=0.0, #phi=0.0 RhO.setParams(d3sp)  # (self, E=0.0, rhoType='Rhodopsin')
        #self.equations = RhO_3states.equations
        self.rhoType = rhoType # E.g. 'ChR2' or 'ArchT'
        
        self.setParams(params)
        
        # Three-state model parameters
        #self.E = E                  # [mV]    Channel reversal potential
        #self.k = 8.2#1.64e-16 *1e18#8.2*(2e-14*1e-3) #0.5*1.2e-14/1.1    # [ms^-1] Quantum efficiency * number of photons absorbed by a RhO molecule per unit time
        #self.eF = self.set_eF(phi) # [ms^-1] Quantum efficiency * number of photons absorbed by a ChR2 molecule per unit time
        #self.Gd = 0.1 #0.117  #1/11              # [ms^-1] @ 1mW mm^-2
        #self.Gr0 = 1/5000 #Gr_dark = 1/5000       # [ms^-1] tau_r,dark = 5-10s p405 Nikolic et al. 2009
        #self.Gr1 = 0.016 #Gr_light = 0.016#1/165       # [ms^-1] (Gr,light @ 1mW mm^-2)
        #self.phi0 = 0.7e15#0.02/(2e-14*1e-3)
        #self.phiSat = 4e17#6.7e17 #13.6/(2e-14*1e-3)
        
        #self.g = 100e-3 * 1.67e5 *2   # [pS] g1 (100fS): Conductance of an individual ion channel * N (100000)
        
        self.useInwardRect = False # Implement inwardRect!!!
        # Initial conditions: Instantaneous photon flux density = 0.0
        self.s_0 = np.array([1,0,0])          # Initialise in dark 
        self.initStates(0.0)        # Initialise in dark
        if verbose:
            print("Nikolic et al. {}-state {} model initialised!".format(self.nStates,self.rhoType))
        
        if verbose > 1: 
            self.printParams()
    
    def reportParams(self): # Replace with local def __str__(self):
        report =  'Three-state model parameters\n'
        report += '============================\n'
        report += 'E    = {:12.3g}\n'.format(self.E)
        report += 'k    = {:12.3g}\n'.format(self.k)
        report += 'Gd   = {:12.3g}\n'.format(self.Gd)
        report += 'Gr0  = {:12.3g}\n'.format(self.Gr0)
        report += 'Gr1  = {:12.3g}\n'.format(self.Gr1)
        report += 'g    = {:12.3g}\n'.format(self.g)
        report += '----------------------------'
        #print('Three-state model parameters')
        #print('============================')
        #print('|  E        = {:12.3g} |'.format(self.E))
        #print('|  k        = {:12.3g} |'.format(self.k))
        #print('|  Gd       = {:12.3g} |'.format(self.Gd))
        #print('|  Gr_dark  = {:12.3g} |'.format(self.Gr_dark))
        #print('|  Gr_light = {:12.3g} |'.format(self.Gr_light))
        #print('|  g        = {:12.3g} |'.format(self.g))
        #print('----------------------------')
        return report
        
    def reportState(self):
        self.dispRates()
        #print('phi = {:.3g}'.format(self.phi))
        #print('V = {:.3g}'.format(self.V))
    
    def set_eF(self, phi):
        #return self.k * phi
        #return 0.5 * phi * (1.2e-8 * 1e-6) / 1.1  # eps * phi * sigma_ret (mm^-2) / wloss
        return self.k * phi #(1 - np.exp(-phi/self.phiSat))
    
    def set_Gr(self, phi):
        return self.Gr0 + (phi>0)*self.Gr1
        #return self.Gr0 + self.Gr1 * np.log(1 + phi/self.phi0) # self.Gr0 + self.Gr1
        #return self.Gr_dark + self.Gr_light * np.log(1 + phi/self.phi0) # self.Gr0 + self.Gr1
        # if phi>0:
            # return self.Gr_dark + self.Gr_light #*log(1 + phi/phi0)
        # else:
            # return self.Gr_dark
        # return 1/(taur_dark*exp(-log(1+phi/phi0))+taur_min) # Fig 6 Nikolic et al. 2009
        ### return Gr_dark + kr*(1-exp(-phi/phi0)) # a = Gr_max - Gr_dark
        
    #def set_Gd(self, phi):
    #    return 1/(taud_dark - a*log(1+phi/phi0)) # Fig 6 Nikolic et al. 2009
    
    def setLight(self, phi):
        """Set transition rates according to the instantaneous photon flux density"""
        #assert(phi >= 0)
        if phi < 0:
            phi = 0
        self.phi = phi
        self.eF = self.set_eF(phi)
        self.Gr = self.set_Gr(phi)
        if verbose > 1:
            self.dispRates()
    
    def dispRates(self):
        #print("Transition rates (phi={:.3g}): O <--[eF={:.3g}]-- C <--[Gr={:.3g}]-- D".format(self.phi,self.eF,self.Gr))
        #print(self.phi); print(type(self.phi))
        #print(self.eF); print(type(self.eF))
        #print(self.Gd); print(type(self.Gd))
        #print(self.Gr); print(type(self.Gr))
        print("Transition rates (phi={:.3g}): C --[eF={:.3g}]--> O --[Gd={:.3g}]--> D --[Gr={:.3g}]--> C".format(self.phi,self.eF,self.Gd,self.Gr))
    
    # def solveStates(self, s_0, t): # , phi e.g. f_phi(t) # http://stackoverflow.com/questions/5649313/a-possible-bug-in-odeint-interp1d-interplay
        # """Function describing the differential equations of the 3-state model to be solved by odeint"""
        # # Add interpolation of values for phi(t) to initialisation f_phi = interp1d(t,sin(w*t),kind='cubic')
        # # Then pass as an argument to integrator: odeint(func, y0, t, args=())
        # # self.setLight(f_phi(t))
        # C,O,D = s_0 # Split state vector into individual variables s1=s[0], s2=s[1], etc
        # f0 = -self.eF*C +             self.Gr*D   # C'
        # f1 =  self.eF*C - self.Gd*O               # O'
        # f2 =              self.Gd*O - self.Gr*D   # D'
        # #f2 = -(f0+f1)
        # return np.array([ f0, f1, f2 ])
    
    def solveStates(self, s_0, t, phi_t=None): # , phi e.g. f_phi(t) # solveStatesPhi_t http://stackoverflow.com/questions/5649313/a-possible-bug-in-odeint-interp1d-interplay
        """Function describing the differential equations of the 3-state model to be solved by odeint"""
        # Add interpolation of values for phi(t) to initialisation f_phi = interp1d(t,sin(w*t),kind='cubic')
        # Then pass as an argument to integrator: odeint(func, y0, t, args=())
        
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
            #print("dydt: phi_t({})={}".format(t,phi_t(t)))
        C,O,D = s_0 # Split state vector into individual variables s1=s[0], s2=s[1], etc
        f0 = -self.eF*C +             self.Gr*D   # C'
        f1 =  self.eF*C - self.Gd*O               # O'
        f2 =              self.Gd*O - self.Gr*D   # D'
        #f2 = -(f0+f1)
        return np.array([ f0, f1, f2 ])
    
    # def jacobian(self, s_0, t):
        # # Jacobian matrix used to improve precision / speed up ODE solver
        # # jac[i,j] = df[i]/dy[j]; where y'(t) = f(t,y)
        # return np.array([[-self.eF, 0, self.Gr],
                         # [self.eF, -self.Gd, 0],
                         # [0, self.Gd, -self.Gr]])
    
    def jacobian(self, s_0, t, phi_t=None): # jacobianPhi_t
        #print("Jac: phi_t({})={}".format(t,phi_t(t)))
        #self.setLight(phi_t(t))
        # Jacobian matrix used to improve precision / speed up ODE solver
        # jac[i,j] = df[i]/dy[j]; where y'(t) = f(t,y)
        return np.array([[-self.eF, 0, self.Gr],
                         [self.eF, -self.Gd, 0],
                         [0, self.Gd, -self.Gr]])
    
    # def hessian(self, s_0, t):
        # Hessian matrix for scipy.optimize.minimize (Only for Newton-CG, dogleg, trust-ncg.)
        # H(f)_ij(X) = D_iD_jf(X)
        # return np.array([[0, 0, 0],
        #                 [0, 0, 0],
        #                 [0, 0, 0]])
    
    def calcPhotocurrent(self, V, states):
        """Calculate the photocurrent from the cell membrane voltage and state matrix"""
        ### if useInwardRect: ... else: fV=1
        O = states[:,1] # time x C,O,D
        I_RhO = self.g*O*(V-self.E)
        return I_RhO * 1e-6 # pS * mV * 1e-6 = nA
    
    def calcPsi(self, states):
        return 1
    
    def calcfV(self, V):
        return 1
    
    def calcOn(self,t):
        """Calculate the on phase current for square light pulses from the analytic solution"""
        r = np.array([lam1, lam2])
        k = np.array([a1, a2])
        I = k * np.exp(-r*t)
        -(a0 + a1*(1-np.exp(-t/tau_act)) + a2*np.exp(-t/tau_deact))
        pass
    
    def calOff():
        """Calculate the off phase current for square light pulses from the analytic solution"""
        -(A*np.exp(-Gd*t))
        pass
    
    def calcSteadyState(self, phi):
        self.setLight(phi)
        denom3 = self.Gd * (self.Gr + self.eF) + self.eF * self.Gr
        Cs = self.Gd*self.Gr/denom3
        Os = self.eF*self.Gr/denom3
        Ds = self.eF*self.Gd/denom3
        self.steadyStates = np.array([Cs, Os, Ds])
        return np.array([Cs, Os, Ds])
    
    def calcSoln(self, t, s0=[1,0,0]):
        [C_0,O_0,D_0] = s0
        P = self.eF
        Gd = self.Gd
        Gr = self.Gr
        
        SP = P*Gd + P*Gr + Gd*Gr
        SQ = P**2 + Gd**2 + Gr**2
        lambda_1 = (Gd + Gr + P + (SQ-2*SP)**(1/2))/2
        lambda_2 = (Gd + Gr + P - (SQ-2*SP)**(1/2))/2
        Den_1 = (2*SP*(SQ-2*SP)**(1/2))
        Fac_1 = (C_0*Gd*P**2 + C_0*Gd**2*P + Gd*Gr**2*D_0 - Gd**2*Gr*D_0 - 2*Gd**2*Gr*O_0 - Gr*D_0*P**2 + Gr**2*D_0*P + Gd*O_0*P**2 - Gd**2*O_0*P + C_0*Gd*Gr*P - Gd*Gr*O_0*P - C_0*Gd*P*(SQ-2*SP)**(1/2) + Gd*Gr*D_0*(SQ-2*SP)**(1/2) + Gr*D_0*P*(SQ-2*SP)**(1/2) - Gd*O_0*P*(SQ-2*SP)**(1/2))
        Fac_2 = (C_0*Gd*P**2 + C_0*Gd**2*P + Gd*Gr**2*D_0 - Gd**2*Gr*D_0 - 2*Gd**2*Gr*O_0 - Gr*D_0*P**2 + Gr**2*D_0*P + Gd*O_0*P**2 - Gd**2*O_0*P + C_0*Gd*Gr*P - Gd*Gr*O_0*P + C_0*Gd*P*(SQ-2*SP)**(1/2) - Gd*Gr*D_0*(SQ-2*SP)**(1/2) - Gr*D_0*P*(SQ-2*SP)**(1/2) + Gd*O_0*P*(SQ-2*SP)**(1/2))
    
        C = (np.exp(-t*lambda_1)*np.exp(-t*lambda_2)*((P*np.exp(t*lambda_2)*Fac_1)/(2*(SP)) + (P*np.exp(t*lambda_1)*Fac_2)/(2*(SP)) + (P**2*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) - (P**2*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) - (Gd*P*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) - (Gr*P*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) + (Gd*P*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) + (Gr*P*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) + (2*Gd**2*Gr*P*np.exp(t*lambda_1)*np.exp(t*lambda_2)*(C_0 + D_0 + O_0))/(SP)))/(2*Gd*P)
        
        O = -(np.exp(-t*lambda_1)*np.exp(-t*lambda_2)*((np.exp(t*lambda_1)*Fac_2)/(2*(SP)) + (np.exp(t*lambda_2)*Fac_1)/(2*(SP)) - (Gd*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) + (Gr*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) - (P*np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) + (Gd*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) - (Gr*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) + (P*np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) - (2*Gd*Gr*P*np.exp(t*lambda_1)*np.exp(t*lambda_2)*(C_0 + D_0 + O_0))/(SP)))/(2*Gd)
        
        D = np.exp(-t*lambda_1)*np.exp(-t*lambda_2)*((np.exp(t*lambda_2)*Fac_1)/(2*(SP)*(SQ-2*SP)**(1/2)) - (np.exp(t*lambda_1)*Fac_2)/(2*(SP)*(SQ-2*SP)**(1/2)) + (Gd*P*np.exp(t*lambda_1)*np.exp(t*lambda_2)*(C_0 + D_0 + O_0))/(SP))
        
        return np.column_stack((C.T,O.T,D.T)) # Finish me!!!


class RhO_4states(OGmodel):
    
    # Class attributes
    nStates = 4
    hasAnalyticSoln = False
    labels = ['$C_1$','$O_1$','$O_2$','$C_2$']
    stateVars = ['C1','O1','O2','C2']
    connect = [[0,1,0,0],
               [1,0,1,0],
               [0,1,0,1],
               [1,0,1,0]]
    
    phi = 0.0  # Instantaneous Light flux
#     phi0 = 1e14 # Normalising light intensity parameter [photons * s^-1 * mm^-2] ~ 1/10 of threshold
    
    ###### gam = 0.1
    
    # Solver parameters
#     dt = 0.1 #0.01 #0.1 #0.001
    
    def __init__(self, params=modelParams['4'], rhoType='Rhodopsin'): #E=0.0, gam=0.05, phi0=1e14, #, A=31192)
        
        self.rhoType = rhoType # E.g. 'ChR2' or 'ArchT'
    
        self.setParams(params)
        # Four-state model parameters
        #self.E = E                  # [mV]      Channel reversal potential
        #self.gam = gam              # []        Ratio of single channel conductances: O2/O1
        #self.phi0 = phi0            # [photons * s^-1 * mm^-2] Normalising photon flux density ~ 1/10 of threshold
        ##self.A = A                  # [um^2]  Effective area of the cell / compartment
        #self.sigma_ret = 1.2e-8 * 1e-6 # Convert from m^2 to mm^2
        #self.w_loss = 1.1
        #self.k1 = 0.5 * self.sigma_ret / self.w_loss # Quantum efficiency * sigma_ret / w_loss
        #self.k2 = 0.15 * self.sigma_ret / self.w_loss # Quantum efficiency * sigma_ret / w_loss
        #self.Gr = 0.000400 # [ms^-1] ==> tau_r = 2.5s
        #self.Gd1 = 0.11 # [ms^-1]
        #self.Gd2 = 0.025 # [ms^-1]
        #self.c1 = 0.03
        #self.c2 = 0.0115
        #self.e12d = 0.01 # [ms^-1]
        #self.e21d = 0.015 # [ms^-1]
        #self.g = 100e-3 * 1.67e5   # [pS] g1 (100fS): Conductance of an individual ion channel * N (~150000)
        
        self.useInwardRect = False # Implement inwardRect!!!
        # Initial conditions
        self.s_0 = np.array([1,0,0,0])
        self.initStates(0.0) # phi
        if verbose:
            print("Nikolic et al. {}-state {} model initialised!".format(self.nStates,self.rhoType))
        
        if verbose > 1: 
            self.printParams()
    
    def reportParams(self): # Replace with local def __str__(self):
        report =  'Four-state model parameters\n'
        report += '===========================\n'
        report += 'E    = {:12.3g}\n'.format(self.E)
        report += 'gam  = {:12.3g}\n'.format(self.gam)
        report += 'phi0 = {:12.3g}\n'.format(self.phi0)
        report += 'k1   = {:12.3g}\n'.format(self.k1)
        report += 'k2   = {:12.3g}\n'.format(self.k2)
        report += 'Gr   = {:12.3g}\n'.format(self.Gr)
        report += 'Gd1  = {:12.3g}\n'.format(self.Gd1)
        report += 'Gd2  = {:12.3g}\n'.format(self.Gd2)
        report += 'c1   = {:12.3g}\n'.format(self.c1)
        report += 'c2   = {:12.3g}\n'.format(self.c2)
        report += 'e12d = {:12.3g}\n'.format(self.e12d)
        report += 'e21d = {:12.3g}\n'.format(self.e21d)
        report += 'g    = {:12.3g}\n'.format(self.g)
        report += '---------------------------'
        return report
    
    def set_Ga1(self, phi):
        # N.B. making Ga a function of time (as in Appendix 1) results in the Six-state model
        # Gai = ei * F * f(t,tChR) See App 1
        # Ga = 1/tauChR
        #e = 0.5
        #sigma_ret = 1.2e-8 * 1e-6 # Convert from m^2 to mm^2
        #w_loss = 1.1
        return self.k1 * phi/self.phi0 #e*phi*sigma_ret / w_loss
    
    def set_Ga2(self, phi):
        #e = 0.15
        #sigma_ret = 1.2e-8  * 1e-6 # Convert from m^2 to mm^2
        #w_loss = 1.1
        return self.k2 * phi/self.phi0 #e*phi*sigma_ret / w_loss
    
    def set_e12(self, phi):
        return self.e12d + self.c1*np.log(1+(phi/self.phi0)) # e12(phi=0) = e12d

    def set_e21(self, phi):
        return self.e21d + self.c2*np.log(1+(phi/self.phi0)) # e21(phi=0) = e21d

    def setLight(self, phi):
        """Set transition rates according to the instantaneous photon flux density"""
        #assert(phi >= 0)
        if phi < 0:
            phi = 0
        self.phi = phi
        self.Ga1 = self.set_Ga1(phi)
        self.Ga2 = self.set_Ga2(phi)
        self.e12 = self.set_e12(phi)
        self.e21 = self.set_e21(phi)
        if verbose > 1:
            self.dispRates()
            
    def dispRates(self):
        print("Transition rates (phi={:.3g}): C1 --[Ga1={:.3g}]--> O1 --[e12={:.3g}]--> O2".format(self.phi,self.Ga1,self.e12))
        print("Transition rates (phi={:.3g}): O1 <--[e21={:.3g}]-- O2 <--[Ga2={:.3g}]-- C2".format(self.phi,self.e21,self.Ga2))
    
    def solveStates(self, s_0, t, phi_t=None):
        """Function describing the differential equations of the 4-state model to be solved by odeint"""
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        C1,O1,O2,C2 = s_0 # Split state vector into individual variables s1=s[0], s2=s[1], etc
        f0 = -self.Ga1*C1   +       self.Gd1*O1     +                       self.Gr *C2 # C1'
        f1 =  self.Ga1*C1 - (self.Gd1+self.e12)*O1  +   self.e21*O2                     # O1'
        f2 =                        self.e12*O1 - (self.Gd2+self.e21)*O2 +  self.Ga2*C2 # O2'
        f3 =                                        self.Gd2*O2 - (self.Ga2+self.Gr)*C2 # C2'
        #f3 = -(f0+f1+f2)
        return np.array([ f0, f1, f2, f3 ])
    
    def jacobian(self, s_0, t, phi_t=None):
        # Jacobian matrix used to improve precision / speed up ODE solver
        # jac[i,j] = df[i]/dy[j]; where y'(t) = f(t,y)
        return np.array([[-self.Ga1, self.Gd1, 0, self.Gr],
                         [self.Ga1, -(self.Gd1+self.e12), self.e21, 0],
                         [0, self.e12, -(self.Gd2+self.e21), self.Ga2],
                         [0, 0, self.Gd2, -(self.Ga2+self.Gr)]])
    
    def calcPhotocurrent(self, V, states):
        O1 = states[:,1]
        O2 = states[:,2]
        I_RhO = self.g*(O1 + self.gam*O2)*(V-self.E)
        return I_RhO * 1e-6 # pS * mV * 1e-6 = nA
    
    def calcPsi(self, states):
        C1, O1, O2, C2 = states
        return O1 + self.gam * O2
    
    def calcfV(self, V):
        return 1
    
    def calcSteadyState(self, phi):
        self.setLight(phi)
        Ga1 = self.Ga1
        Ga2 = self.Ga2
        Gr = self.Gr
        Gd1 = self.Gd1
        Gd2 = self.Gd2
        e12 = self.e12
        e21 = self.e21
        denom4 = Ga1 * (e12 * (Gr + Gd2 + Ga2) + e21 * (Gr + Ga2) + Gd2 * Gr) + Gd1 * (e21 * (Gr + Ga2) + Gd2 * Gr) + e12 * Gd2 * Gr
        C1s = (Gd1 * (e21 * (Gr + Ga2) + Gd2 * Gr) + e12 * Gd2 * Gr) / denom4
        O1s = (Ga1 * (e21 * (Gr + Ga2) + Gd2 * Gr)) / denom4
        O2s = (e12 * Ga1 * (Gr + Ga2)) / denom4
        C2s = (e12 * Ga1 * Gd2) / denom4
        self.steadyStates = np.array([C1s, O1s, O2s, C2s])
        return np.array([C1s, O1s, O2s, C2s])


    
class RhO_6states(OGmodel):
    
    # Class attributes
    nStates = 6
    hasAnalyticSoln = False
    #labels = ['$s_1$','$s_2$','$s_3$','$s_4$','$s_5$','$s_6$']
    labels = ['$C_1$','$I_1$','$O_1$','$O_2$','$I_2$','$C_2$']
    # stateVars = ['s1','s2','s3','s4','s5','s6'] # 
    stateVars = ['C1','I1','O1','O2','I2','C2'] #['C1','I1','O1','O2','I2','C2']
    connect = [[0,1,0,0,0,0], # s_1 --> s_i=1...6
               [0,0,1,0,0,0], # s_2 -->
               [1,0,0,1,0,0],
               [0,0,1,0,0,1],
               [0,0,0,1,0,0],
               [1,0,0,0,1,0]]
    
#     phi = 0.0  # Instantaneous Light flux
#     phi0 = 1e14 # Normalising light intensity parameter [photons * s^-1 * mm^-2] ~ 1/10 of threshold
    
    # Solver parameters
#     dt = 0.1 #0.01 #0.1 #0.001
    
    def __init__(self, params=modelParams['6'], rhoType='Rhodopsin'): # E=0.0, gam=0.05, phi0=1e14, A=3.12e4# phi=0.0
        
        self.rhoType = rhoType # E.g. 'ChR2' or 'ArchT'
        
        self.setParams(params)
        #self.g = self.gbar * self.A # [pS]      Total conductance
        
        ### Model constants
        #self.E = E                  # [mV]    Channel reversal potential
        #self.gam = gam              # []      Ratio of single channel conductances: O2/O1
        #self.phi0 = phi0            # [photons * s^-1 * mm^-2] Normalising photon flux density ~ 1/10 of threshold
        #self.A = A                  # [um^2]  Effective area of the cell / compartment
        #self.v0 = 43  # (mV)
        #self.v1 = 4.1 # (mV) #-4.1 Changed to correct f(v) but preserve the form of the numerator
        #self.gbar = 2.4 # (pS * um^-2)
        # gam = 0.05 # (Dimensionless)
        # A = 31192#e-6 # (um^2)
        # E = 0 # 0:8 (mV) 0e-3
        #self.g = self.gbar * self.A # [pS]      Total conductance
        
        #self.a10 = 5
        #self.a2 = 1
        #self.a30 = 0.022
        #self.a31 = 0.0135
        #self.a4 = 0.025
        #self.a6 = 0.00033   # = 1/tau6, tau6 = 3s = 3000ms
        #self.b1 = 0.13
        #self.b20 = 0.011
        #self.b21 = 0.0048
        #self.b3 = 1
        #self.b40 = 1.1
        
        self.useInwardRect = True
        
        # Initial conditions
        self.s_0 = np.array([1,0,0,0,0,0])  # [s1_0=1, s2_0=0, s3_0=0, s4_0=0, s5_0=0, s6_0=0] # array not necessary
        self.initStates(0.0) # phi
        if verbose:
            print("Grossman et al. {}-state {} model initialised!".format(self.nStates,self.rhoType))
        
        if verbose > 1: 
            self.printParams()
    
    def set_a1(self,phi):
        return self.a10*(phi/self.phi0)

    def set_a3(self,phi):
        return self.a30 + self.a31*np.log(1+(phi/self.phi0))

    def set_b2(self,phi):
        return self.b20 + self.b21*np.log(1+(phi/self.phi0))

    def set_b4(self,phi):
        return self.b40*(phi/self.phi0)

    def setLight(self,phi):
        #assert(phi >= 0)
        if phi < 0:
            phi = 0
        self.phi = phi
        self.a1 = self.set_a1(phi)
        self.a3 = self.set_a3(phi)
        self.b2 = self.set_b2(phi)
        self.b4 = self.set_b4(phi)
        if verbose > 1:
            self.dispRates()
        
    def dispRates(self):
        print("Transition rates (phi={:.3g}): s3 --[a3]--> s4 = {}; s3 <--[b2]-- s4 = {}".format(self.phi,self.a3,self.b2))
        print("                  ^s2      s5^")
        print("                   \        /")
        print("                   [a1]   [b4]")
        print("                     \    /")
        print("                     s1  s6")
        print("Transition rates [a1] = {}; [b4] = {}".format(self.a1,self.b4))
    
    def solveStates(self, s_0, t, phi_t=None):
        """Function describing the differential equations of the 6-state model to be solved by odeint"""
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        s1,s2,s3,s4,s5,s6 = s_0 # Split state vector into individual variables s1=s[0], s2=s[1], etc
        d0 = -self.a1*s1 + self.b1*s3 + self.a6*s6           # s1'
        d1 =  self.a1*s1 - self.a2*s2                        # s2'
        d2 =  self.a2*s2 - (self.b1+self.a3)*s3 + self.b2*s4 # s3'
        d3 =  self.a3*s3 - (self.b2+self.a4)*s4 + self.b3*s5 # s4'
        d4 = -self.b3*s5 + self.b4*s6                        # s5'
        d5 =  self.a4*s4 - (self.b4+self.a6)*s6              # s6'
        #d5 = - (f0+f1+f2+f3+f4)
        #dr0 = -2*pi*f*sin(2*pi*f*t)*C(1+cos(2*pi*f*t))      # d/dt (A(1+cos(2*pi*f*t)))
        return np.array([ d0, d1, d2, d3, d4, d5 ])

    def jacobian(self, s_0, t, phi_t=None):
        return np.array([[-self.a1, 0, self.b1, 0, 0, self.a6],
                         [self.a1, -self.a2, 0, 0, 0, 0],
                         [0, self.a2, -(self.b1+self.a3), self.b2, 0, 0],
                         [0, 0, self.a3, -(self.b2+self.a4), self.b3, 0],
                         [0, 0, 0, 0, -self.b3, self.b4],
                         [0, 0, 0, self.a4, 0, -(self.b4+self.a6)]])
        
    def calcPhotocurrent(self, V, states): # Change to (V, psi)?
        """Takes Voltage [mV] and state variables s_3 and s_4 to calculate current [nA]
        By convention, positive ions entering the cell --> negative current (e.g. Na^+ influx). 
        Conversely, Positive ions leaving (or negative ions entering) the cell --> positive current (e.g. Cl^- in or K^+ out). """
        s3 = states[:,2] # O1
        s4 = states[:,3] # O2
        psi = s3 + (self.gam * s4) # Dimensionless
        ### if useInwardRect: ... else: fV=1
        try:
            fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0)) # Dimensionless
        except ZeroDivisionError:
            if np.isscalar(V):
                if (V == self.E):
                    fV = self.v1/self.v0
            else: #type(fV) in (tuple, list, array)
                fV[np.isnan(fV)] = self.v1/self.v0 # Fix the error when dividing by zero
        g_RhO = self.g * psi * fV # self.gbar * self.A # Conductance (pS * um^-2)
        I_RhO =  g_RhO * (V - self.E) # Photocurrent: (pS * mV)
        return I_RhO * (1e-6) # 10^-12 * 10^-3 * 10^-6 (nA)
    
##############################################################    
    def inwardRectifier(V):
        if self.v0 == 0:
            #print("f(V) undefined for v0 = 0")
            warnings.warn("f(V) undefined for v0 = 0")
        try:
            fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0)) # Dimensionless
        except ZeroDivisionError:
            pass
        fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0)) # Dimensionless
        if np.isscalar(fV):
            if (V == self.E):
                fV = self.v1/self.v0
        else: #type(fV) in (tuple, list, array)
            fV[np.isnan(fV)] = self.v1/self.v0 # Fix the error when dividing by zero

        return fV
    
    def equalFloats(x,y,eps=1e-6):
        diff = np.abs(x-y)
        if (x == y): # Catch infinities
            return True
        elif (x == 0 or y == 0 or diff < floatmin): # Either one number is 0 or both are very close
            return diff < (eps * floatmin)
        else: # Use relative error
            return diff / (np.abs(x) + np.abs(y)) < eps
#############################################################
        
    def calcSteadyState(self,phi): # implicitly depends on phi0
        self.setLight(phi)
        a1 = self.a1
        a2 = self.a2
        a3 = self.a3
        a4 = self.a4
        a6 = self.a6
        b1 = self.b1
        b2 = self.b2
        b3 = self.b3
        b4 = self.b4
        denom6 = (a1*a2*(a3*(b3*(b4+a4)+a4*b4)+b2*b3*b4)+b1*(a2*b2*b3*b4+a1*b2*b3*b4+a6*(a2*(b2*b3+a4*b3)+a1*(b2*b3+a4*b3)))+a6*(a1*(a2*(b2*b3+a4*b3+a3*b3)+a3*a4*b3)+a2*a3*a4*b3))
        s1s = (b1*(a2*b2*b3*b4 + a2*a6*(b2*b3 + a4*b3)) + a2*a3*a4*a6*b3)/denom6
        s2s = (b1*(a1*b2*b3*b4 + a1*a6*(b2*b3 + a4*b3)) + a1*a3*a4*a6*b3)/denom6
        s3s = (a1*a2*b2*b3*b4 + a1*a2*a6*(b2*b3 + a4*b3))/denom6
        s4s = (a1*a2*a3*b3*b4 + a1*a2*a3*a6*b3)/denom6
        s5s = (a1*a2*a3*a4*b4)/denom6
        s6s = (a1*a2*a3*a4*b3)/denom6
        self.steadyStates = np.array([s1s, s2s, s3s, s4s, s5s, s6s])
        return np.array([s1s, s2s, s3s, s4s, s5s, s6s])
    
    def calcPsi(self, states):
        s1,s2,s3,s4,s5,s6 = states
        return s3 + self.gam * s4
        

models = {'3': RhO_3states, '4': RhO_4states, '6': RhO_6states}