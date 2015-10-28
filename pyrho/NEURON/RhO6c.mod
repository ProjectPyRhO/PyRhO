TITLE Nikolic 6-state rhodopsin model

NEURON {
    : SUFFIX RhO6c : For density mechanisms only
    POINT_PROCESS RhO6c
    NONSPECIFIC_CURRENT i :IRhO     : This could be changed to update a specific ion
    RANGE E, gam, v0, v1, g, gph, gv, i
    RANGE k1, Go1, Gf0, kf, Gd2, Gr0, Gd1, Gb0, kb, Go2, k2, p, q
    RANGE phi, phim :, lambda
}


UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    :(photons) = (1)
}


PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files

: Illumination constants   
:   lambda  = 470   : (nm)
    phim    = 1e16  :(photons/s mm2)

: Conductance
    E       = 0     (mV) : Channel reversal potential
    gam     = 0.05  (1)  : Ratio of open-state conductances
    g       = 75000 (pS) : gbar : defined in the main prog as HR_expression*Area

: Inward rectifier conductance    : reversal potential for NpHr is about -400mV        e_rev = -400        :(mV)  
    v0      = 43    (mV)
    v1      = 4.1   (mV)    

: State transition rate parameters (/ms)
    k1     = 5         (/ms)
    Go1      = 1         (/ms)
    Gf0     = 0.022     (/ms)
    kf     = 0.0135    (/ms)
    Gd2      = 0.025     (/ms)
    Gr0      = 0.00033   (/ms)
    Gd1      = 0.13      (/ms)
    Gb0     = 0.011     (/ms)
    kb     = 0.0048    (/ms)
    Go2      = 1         (/ms)
    k2     = 1.1       (/ms)
    p       = 0.7       (1)
    q       = 0.47      (1)
}


ASSIGNED {
    phi     :(photons/s mm2)

    Ga1      (/ms)
    Gf      (/ms)
    Gb      (/ms)
    Ga2      (/ms)
    h1      (1)
    h2      (1)
    
    gph     (1) : Fractional conductance
    gv      (1) : Fractional conductance
    v       (mV) 
    i      (nA)
}


STATE { C1 I1 O1 O2 I2 C2 }


BREAKPOINT {
    SOLVE kin METHOD sparse
    gph     = O1+gam*O2        : light dependency op=[0:1] 
    gv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : voltage dependency (equal 1 at Vm=-70mV)
    i = g*gph*gv*(v-E)*(1e-6)
}


INITIAL {   : Initialise variables
: State occupancy proportions - starts in C1
    C1  = 1     : Closed state 1 (Ground state)
    I1  = 0     : Intermediate state 1
    O1  = 0     : Open state 1 (High conductivity)
    O2  = 0     : Open state 2 (Low conductivity)
    I2  = 0     : Intermediate state 2
    C2  = 0     : Closed state 2

    phi = 0
    rates(phi) : necessary?
    i   = 0
}


KINETIC kin {
    rates(phi)
    ~ C1 <-> I1 (Ga1, 0) :(a11, 0)
    ~ I1 <-> O1 (Go1, 0) :(a12, 0)
    ~ O1 <-> C1 (Gd1, 0)      
    ~ O1 <-> O2 (Gf, Gb) :(a2, b2)    
    ~ O2 <-> I2 (0, Go2) :(0, b31)
    ~ I2 <-> C2 (0, Ga2) :(0, b32)
    ~ O2 <-> C2 (Gd2, 0) :(a3, 0)    
    ~ C2 <-> C1 (Gr0, 0) :(a4, 0)
    CONSERVE C1 + I1 + O1 + O2 + I2 + C2 = 1
}


PROCEDURE rates(phi) { : Define equations for calculating transition rates

    if (phi>0) {
        h1 = 1/(1+pow(phim,p)/pow(phi,p))  : pow(phi,p)/(pow(phi,p) + pow(phim,p))
        h2 = 1/(1+pow(phim,q)/pow(phi,q))  : pow(phi,q)/(pow(phi,q) + pow(phim,q))
    } else {
        h1 = 0
        h2 = 0
    }
    
    Ga1 = k1 * h1
    Gf = Gf0 + kf * h2
    
    Gb = Gb0 + kb * h2
    Ga2 = k2 * h1
}


: Add functions for calculating transition rates outside of simulations
:FUNCTION A1(phi) {
:    A1 = a10*(phi/phi0)
:}
: ...