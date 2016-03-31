TITLE Nikolic 3-state rhodopsin model

NEURON  {
    POINT_PROCESS RhO3c
    NONSPECIFIC_CURRENT i
    RANGE i, E, v0, v1, g0      :, fphi, fv
    RANGE k_a, k_r, Gd, Gr0, p, q
    RANGE phi, phi_m            : phiOn, 
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
:   lambda  = 470       : (nm)
    phi_m   = 1e16      : (photons/sec mm2) < 0, 1e30 > : Hill Constant
  
: Conductance    
    E       = 0         (mV)                : Channel reversal potential
    g0      = 1         (pS)   < 0, 1e9 >   : Biological conductance scaling factor

: Inward rectifier conductance              fv = (1-exp(-(v-E)/v0))*v1/(v-E) := 1 at Vm=-70mV
    v0      = 43        (mV)
    v1      = 4.1       (mV)
    
: State transition rate parameters (/ms)
    k_a     = 0.28      (/ms)   < 0, 1e9 >  : Quantum efficiency * number of photons absorbed by a RhO molecule per unit time 0.5*1.2e-14 /1.1
    p       = 0.4       (1)                 : Hill Coefficient (Ga)
    k_r     = 0.28      (/ms)   < 0, 1e9 >  : Recovery rate scaling
    q       = 0.4       (1)                 : Hill Coefficient (Gr)
    Gd      = 0.0909    (/ms)   < 0, 1e9 >  : @ 1mW mm^-2
    Gr0     = 0.0002    (/ms)   < 0, 1e9 >  : tau_r,dark = 5-10s p405 Nikolic et al. 2009
}


ASSIGNED {
    phi     : (1)        : "flux" should now be "intensity". Er is normalised to be dimensionless
    Ga      (/ms)
    Gr      (/ms)
    fphi    (1)
    fv      (1)
    v       (mV) 
    i       (nA)
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
    rates(phi) : necessary?
    i       = 0
    setV1(E, v0)
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
