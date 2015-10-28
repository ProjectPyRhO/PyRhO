TITLE Nikolic 4-state rhodopsin model

NEURON {
    POINT_PROCESS RhO4c
    NONSPECIFIC_CURRENT i
    RANGE E, gam, v0, v1, g, gph, gv, i
    RANGE k1, k2, kf, kb, Gf0, Gb0, Gd1, Gd2, Gr0, p, q
    RANGE phi, phim :, lambda
}


UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    :(photons)  =(1)
}


PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files

: Illumination constants   
    phim    = 1e16
:   lambda  = 470   : (nm)

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
    v       (mV) 
    i       (nA)
}


STATE { C1 O1 O2 C2 }


BREAKPOINT {
    SOLVE kin METHOD sparse
    gph     = O1+gam*O2        : Light dependency op=[0:1] 
    gv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : Voltage dependency (equal 1 at Vm=-70mV)
    i = g*gph*gv*(v-E)*(1e-6)   : Photocurrent
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