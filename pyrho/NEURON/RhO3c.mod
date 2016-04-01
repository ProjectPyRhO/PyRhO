TITLE Three-state rhodopsin kinetic model

NEURON  {
    POINT_PROCESS RhO3c
    NONSPECIFIC_CURRENT i
    RANGE i, E, v0, v1, g0      :, fphi, fv
    RANGE k_a, k_r, Gd, Gr0, p, q
    RANGE phi, phi_m
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
:   lambda  = 470       : (nm)              : Light wavelength
    phi_m   = 1e16      : (photons/sec mm2) < 0, 1e30 > : Hill Constant
  
: Conductance    
    g0      = 1         (pS)   < 0, 1e9 >   : Biological conductance scaling factor

: Inward rectifier conductance  fv(-70mV) := 1
    E       = 0         (mV)                : Channel reversal potential
    v0      = 43        (mV)                : Rectifier exponent scaling factor
    v1      = 4.1       (mV)                := (70+E)/(exp((70+E)/v0)-1) Since fv(-70) := 1
    
: State transition rate parameters (/ms)
    k_a     = 0.28      (/ms)   < 0, 1e9 >  : Activation rate scaling (Ga)
    p       = 0.4       (1)     < 0, 1000 > : Hill Coefficient (Ga)
    k_r     = 0.28      (/ms)   < 0, 1e9 >  : Recovery rate scaling (Gr)
    q       = 0.4       (1)     < 0, 1000 > : Hill Coefficient (Gr)
    Gd      = 0.0909    (/ms)   < 0, 1e9 >  : Deactivatrion Rate @ 1mW mm^-2
    Gr0     = 0.0002    (/ms)   < 0, 1e9 >  : tau_r,dark = 5-10s p405 Nikolic et al. 2009
}


ASSIGNED {
    phi     : (1)   : Instantaneous flux
    Ga      (/ms)   : Activation rate (C -> O)
    Gr      (/ms)   : Recovery Rate (D -> C)
    fphi    (1)     : Photocycle fractional conductance
    fv      (1)     : Rectifier fractional conductance
    v       (mV)    : Membrane voltage
    i       (nA)    : Photocurrent
}


STATE { C O D }


BREAKPOINT {
    SOLVE kin METHOD sparse
    fphi    = O                             : Light dependency on open state O=[0:1] 
    fv      = (1-exp(-(v-E)/v0))*v1/(v-E)   : Voltage dependency (equal 1 at Vm=-70mV)
    i       = g0*fphi*fv*(v-E)*(1e-6)       : Photocurrent
}


INITIAL {   : Initialise variables to fully dark-adapted state
: State occupancy proportions - starts in C
    C       = 1
    O       = 0
    D       = 0
    
    phi     = 0
    rates(phi)
    setV1(E, v0)
    i       = 0
}

:DERIVATIVE deriv {
:   rates(flux)
:   C' = (Gr * D) - (Ga * C)
:   O' = (Ga * C) - (Gd * D)
:   D' = (Gd * O) - (Gr * D)
:   D = 1 - (C + O)
:}

KINETIC kin {
    rates(phi)
    ~ C <-> O (Ga, 0)
    ~ O <-> D (Gd, 0)
    ~ D <-> C (Gr, 0)
    CONSERVE C + O + D = 1  : N.B. "NEURON's CVode method ignores conservation"
}


PROCEDURE rates(phi) {  : Define equations for calculating transition rates
    if (phi>0) {        : Safeguard against negative phi values
        Ga = k_a * 1/(1+pow(phi_m,p)/pow(phi,p))        : pow(phi,p)/(pow(phi,p) + pow(phi_m,p))
        Gr = Gr0 + k_r * 1/(1+pow(phi_m,q)/pow(phi,q))  : pow(phi,q)/(pow(phi,q) + pow(phi_m,q))
    } else {
        Ga = 0
        Gr = Gr0
    }
}


PROCEDURE setV1(E (mV), v0 (mV)) { :
: Calculate v1 from v0 and E to satisfy the definition: f(V=-70) := 1
    v1 = (70+E)/(exp((70+E)/v0)-1)
    : v1 = (-70-E)/(1-exp(-(-70-E)/v0))
}

: FUNCTION f_name(arg1 (units1), arg2 (units2), . . . ) (returned_units)
FUNCTION g(g0 (pS), fphi (1), fv (1)) (pS) { : Make the effective channel conductance available to hoc
    g = g0 * fphi * fv
}
