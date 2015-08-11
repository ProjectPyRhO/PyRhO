
TITLE Nikolic 4-state rhodopsin model

NEURON {
    POINT_PROCESS RhO4
    NONSPECIFIC_CURRENT i
    RANGE E, gam, v0, v1, g, gph, gv, i    : , gbar
    RANGE k1, k2, kf, kb, Gf0, Gb0, Gd1, Gd2, Gr0, p, q
    RANGE phiOn, phi, phim, delD, onD, offD, nPulses :phi0,
}


UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    :(photons)  =(1)
}


PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files

: Illumination
    phiOn   = 1e12 :ph/s/mm^2 :(mW/mm2)                             _______
    delD    = 25    (ms) <0, 1e9>   : delay before ON phase        |  ON   |  OFF
    onD     = 100   (ms) <0, 1e9>   : duration of ON phase  <-del->|<-ton->|<-toff->
    offD    = 50    (ms) <0, 1e9>   : duration of OFF phase _______|       |________
    nPulses = 1                     : num pulses to deliver        <-- one pulse -->

: Illumination constants   
:    Er0     = 0.1   : (mW/mm2)
:    phi0 = 1e14
    phim    = 1e16
:   lambda  = 470   : (nm)
:   flux0   = 2e16  : (photons/(sec*cm2))    

: Conductance
    E       = 0         (mV) : Channel reversal potential
    gam     = 0.05      (1)  : Ratio of open-state conductances
    g       = 75000     (pS) : gbar : defined in the main prog as HR_expression*Area

: Inward rectifier conductance    
    v0      = 43        (mV)
    v1      = 4.1       (mV)    

: State transition rate parameters (/ms)
    k1      = 0.05      (/ms)
    k2      = 0.015     (/ms)
    kf      = 0.03      (/ms)
    kb      = 0.0115    (/ms)
    Gf0     = 0.01      (/ms)
    Gb0     = 0.015     (/ms)
    p       = 0.7       (1)
    q       = 0.47      (1)
    Gd1     = 0.11      (/ms)
    Gd2     = 0.025     (/ms)
    Gr0     = 0.0004    (/ms)

}


ASSIGNED {

    phi     :(photons/s mm2)        : Instantaneous flux

    Ga1     (/ms)  : C1 -> O1
    Ga2     (/ms)  : C2 -> O2
    Gf     (/ms)  : O1 -> O2
    Gb     (/ms)  : O2 -> O1
    h1      (1)
    h2      (1)

    gph     (1) : Fractional conductance
    gv      (1) : Fractional conductance
    :g       (1) : Fractional conductance
    v       (mV) 
    :on            : Flag for light
    tally         : Number of pulses left to deliver
    i       (nA)

}


STATE { C1 O1 O2 C2 }


BREAKPOINT {
    SOLVE kin METHOD sparse
    gph     = O1+gam*O2        : Light dependency op=[0:1] 
    gv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : Voltage dependency (equal 1 at Vm=-70mV)
    :g       = gph*gv            : total relative conductance
    :i   = gbar*g*(v-E)*(1e-6)       : +ve ELECTRODE_CURRENT depolarises the cell :IRhO
    i = g*gph*gv*(v-E)*(1e-6)   : Photocurrent
    :iChR   = -i                : for plotting (by convention depolar. curr. are -ve)
}


INITIAL {   : Initialise variables
: State occupancy proportions - starts in C1
    C1  = 1     : Closed state 1 (Ground state)
    O1  = 0     : Open state 1 (High conductivity)
    O2  = 0     : Open state 2 (Low conductivity)
    C2  = 0     : Closed state 2

    phi = 0
    rates(phi) : necessary?
    i   = 0
    tally   = nPulses : Set the tally to the number of pulses required 
    if (tally > 0) {
        net_send(delD, 1)
    :    on = 0
        tally = tally - 1
    }
}


KINETIC kin {
    rates(phi)
    ~ C1 <-> O1 (Ga1, Gd1)
    ~ O1 <-> O2 (Gf,  Gb)
    ~ O2 <-> C2 (Gd2, Ga2)
    ~ C2 <-> C1 (Gr0, 0)
    CONSERVE C1 + O1 + O2 + C2 = 1
}


PROCEDURE rates(phi) { : Define equations for calculating transition rates

    :p1 = pow(phi,p)
    :p2 = pow(phi,q)
    :d1 = pow(phim,p)
    :d2 = pow(phim,q)
    
    if (phi > 0) {
        h1 = 1/(1 + pow(phim,p)/pow(phi,p)) : p1/(p1 + d1)
        h2 = 1/(1 + pow(phim,q)/pow(phi,q)) : p2/(p2 + d2)
    } else {
        h1 = 0
        h2 = 0
    }
    
        :Ga1 = k1*(phi/phi0)
    Ga1 = k1 * h1             :Ga1 = k1 * pow(phi,p)/(pow(phi,p) + pow(phim,p))
        :Ga2 = k2*(phi/phi0)
    Ga2 = k2 * h1             :Ga2 = k2 * pow(phi,p)/(pow(phi,p) + pow(phim,p))
        :Gf = Gf0 + kf*log(1+(phi/phi0))
    Gf = Gf0 + kf * h2      :Gf = Gf0 + kf * pow(phi,q)/(pow(phi,q) + pow(phim,q))
        :Gb = Gb0 + kb*log(1+(phi/phi0))
    Gb = Gb0 + kb * h2      :Gb = Gb0 + kb * pow(phi,q)/(pow(phi,q) + pow(phim,q))
    
}


: Add functions for calculating transition rates outside of simulations
:FUNCTION A1(phi) {
:    A1 = a10*(phi/phi0)
:}
: ...


NET_RECEIVE (w) {
    if (flag == 1) { : ignore any but self-events with flag == 1. This may not be necessary... see e.g. nrn/src/nrnoc/intfire1.mod
    :    if (on == 0) { : Turn the light on 
            phi = phiOn
    :        on = 1
            net_send(onD, 0) :net_send(onD, 1) : Schedule the next off phase
        } else { : Turn the light off
            phi = 0
    :        on = 0
            if (tally > 0) { : Schedule the next on phase
                net_send(offD, 1)
                tally = tally - 1
            }
        }
    :}
}
