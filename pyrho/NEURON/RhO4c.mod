TITLE Four-state rhodopsin kinetic model

NEURON {
    POINT_PROCESS RhO4c
    NONSPECIFIC_CURRENT i
    RANGE i, E, gam, v0, v1, g0 :, fphi, fv
    RANGE k1, k2, k_f, k_b, Gf0, Gb0, Gd1, Gd2, Gr0, p, q
    RANGE phi, phi_m :, lambda
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
:   lambda  = 470   : (nm)
    phi_m   = 1e16           : Hill Constant

: Conductance
    g0      = 75000     (pS)    < 0, 1e9 >  : Biological conductance scaling factor
    gam     = 0.05      (1)     < 0, 1 >    : Ratio of open-state conductances gO2/gO1

: Inward rectifier conductance  fv(-70mV) := 1
    E       = 0         (mV)                : Channel reversal potential
    v0      = 43        (mV)                : Rectifier exponent scaling factor
    v1      = 4.1       (mV)                := (70+E)/(exp((70+E)/v0)-1) Since fv(-70) := 1

: State transition rate parameters (/ms)
    k1      = 0.05      (/ms)   < 0, 1e9 >  : Activation rate scaling factor (Ga1)
    k2      = 0.015     (/ms)   < 0, 1e9 >  : Activation rate scaling factor (Ga2)
    k_f     = 0.03      (/ms)   < 0, 1e9 >  : Forwards (O1 -> O2) activation rate scaling (Gf)
    k_b     = 0.0115    (/ms)   < 0, 1e9 >  : Backwards (O1 <- O2) activation rate scaling (Gb)
    Gf0     = 0.01      (/ms)   < 0, 1e9 >  : Forwards (O1 -> O2) dark rate (Gf)
    Gb0     = 0.015     (/ms)   < 0, 1e9 >  : Backwards (O1 <- O2) dark rate (Gb)
    p       = 0.7       (1)     < 0, 1000 > : Hill Coefficient Ga{1,2}
    q       = 0.47      (1)     < 0, 1000 > : Hill Coefficient G{f,b}
    Gd1     = 0.11      (/ms)   < 0, 1e9 >  : Deactivation rate (O1 -> C1)
    Gd2     = 0.025     (/ms)   < 0, 1e9 >  : Deactivation rate (O2 -> C2)
    Gr0     = 0.0004    (/ms)   < 0, 1e9 >  : Dark recovery rate (C1 <- C2)
}


ASSIGNED {
    phi     :(photons/s mm2)    : Instantaneous flux

    Ga1     (/ms)               : C1 -> O1
    Ga2     (/ms)               : C2 -> O2
    Gf      (/ms)               : O1 -> O2
    Gb      (/ms)               : O2 -> O1
    h1      (1)                 : Hill Equation (Ga1, Ga2)
    h2      (1)                 : Hill Equation (Gf, Gb)

    fphi    (1)                 : Photocycle fractional conductance
    fv      (1)                 : Rectifier fractional conductance
    v       (mV)                : Membrane voltage
    i       (nA)                : Photocurrent
}


STATE { C1 O1 O2 C2 }


BREAKPOINT {
    SOLVE kin METHOD sparse
    fphi    = O1+gam*O2                     : Light dependency {O1,O2}=[0:1] 
    fv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : Voltage dependency (equal 1 at Vm=-70mV)
    i       = g0*fphi*fv*(v-E)*(1e-6)       : Photocurrent
}


INITIAL {   : Initialise variables
: State occupancy proportions - starts in C1
    C1  = 1     : Closed state 1 (Ground state)
    O1  = 0     : Open state 1 (High conductivity)
    O2  = 0     : Open state 2 (Low conductivity)
    C2  = 0     : Closed state 2

    phi = 0
    rates(phi)
    setV1(E, v0)
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
    :d1 = pow(phi_m,p)
    :d2 = pow(phi_m,q)
    
    if (phi > 0) {
        h1 = 1/(1 + pow(phi_m,p)/pow(phi,p)) : p1/(p1 + d1)
        h2 = 1/(1 + pow(phi_m,q)/pow(phi,q)) : p2/(p2 + d2)
    } else {
        h1 = 0
        h2 = 0
    }
    
    Ga1 = k1 * h1           :Ga1 = k1 * pow(phi,p)/(pow(phi,p) + pow(phi_m,p))
    Ga2 = k2 * h1           :Ga2 = k2 * pow(phi,p)/(pow(phi,p) + pow(phi_m,p))
    Gf = Gf0 + k_f * h2     :Gf = Gf0 + k_f * pow(phi,q)/(pow(phi,q) + pow(phi_m,q))
    Gb = Gb0 + k_b * h2     :Gb = Gb0 + k_b * pow(phi,q)/(pow(phi,q) + pow(phi_m,q))
    
}


: Add functions for calculating transition rates outside of simulations

PROCEDURE setV1(E (mV), v0 (mV)) { :
: Calculate v1 from v0 and E to satisfy the definition: f(V=-70) := 1
    v1 = (70+E)/(exp((70+E)/v0)-1)
    : v1 = (-70-E)/(1-exp(-(-70-E)/v0))
}
