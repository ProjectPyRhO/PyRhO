TITLE Nikolic 6-state rhodopsin model

NEURON {
    : SUFFIX RhO6 : For density mechanisms only
    POINT_PROCESS RhO6
    :ELECTRODE_CURRENT i :IRhO
    NONSPECIFIC_CURRENT i :IRhO     : This could be changed to update a specific ion
    RANGE i, E, gam, v0, v1, g0 :, fphi, fv, i :IRhO :iChR   : , gbar
    RANGE k1, Go1, Gf0, kf, Gd2, Gr0, Gd1, Gb0, kb, Go2, k2, p, q
    RANGE phiOn, phi, phi_m, delD, onD, offD, nPulses :phi0,
}


UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    :(photons) = (1)
}


PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files

: Illumination
    phiOn   = 1e18  :(photons/s mm2):(mW/mm2): irradiance                _______
    delD    = 25    (ms)    <0, 1e9>    : delay before ON phase         |  ON   |  OFF
    onD     = 100   (ms)    <0, 1e9>    : duration of ON phase  <-delD->|<-onD->|<-offD->
    offD    = 50    (ms)    <0, 1e9>    : duration of OFF phase ________|       |________
    nPulses = 1     (1)     <0, 1e3>    : num pulses to deliver         <-- one pulse -->

: Illumination constants   
:    Er0     = 0.1   : (mW/mm2)
    :phi0    = 1e14  :(photons/s mm2)
:   lambda  = 470   : (nm)
    phi_m    = 1e16  :(photons/s mm2)    : Hill Constant

: Conductance
    E       = 0     (mV) : Channel reversal potential
    gam     = 0.05  (1)  : Ratio of open-state conductances
    g0      = 75000 (pS) : gbar : defined in the main prog as HR_expression*Area

: Inward rectifier conductance    : reversal potential for NpHr is about -400mV        e_rev = -400        :(mV)  
    v0      = 43    (mV)
    v1      = 4.1   (mV)    

: State transition rate parameters (/ms)
    k1      = 5         (/ms)
    Go1     = 1         (/ms)
    Gf0     = 0.022     (/ms)
    kf      = 0.0135    (/ms)
    Gd2     = 0.025     (/ms)
    Gr0     = 0.00033   (/ms)
    Gd1     = 0.13      (/ms)
    Gb0     = 0.011     (/ms)
    kb      = 0.0048    (/ms)
    Go2     = 1         (/ms)
    k2      = 1.1       (/ms)
    p       = 0.7       (1)     : Hill Coefficient Ga{1,2}
    q       = 0.47      (1)     : Hill Coefficient G{f,b}
    
:    e1      = 0.004 (/ms)     
:    b10     = 0.08  (/ms)   : Gd1 O1->C1
:    a2dark  = 0.012 (/ms)   
:    a2light = 0.010 (/ms)
:    b2dark  = 0.006 (/ms)   
:    b2light = 0.003 (/ms)

:    a30     = 0.08  (/ms)   : Gd2 O2->C2
:    e3      = 0.003 (/ms)   : Correct for illumination units
:    a40     = 0.00033   (/ms)   : 1/taudark, taudark=3s

}


ASSIGNED {

    phi     :(photons/s mm2)
:    flux        : "flux" should now be "intensity". Er is normalised to be dimensionless

    Ga1     (/ms)
    :a2      (/ms)
    Gf      (/ms)
    :a4      (/ms)
    :a6      (/ms)
    :b1      (/ms)
    Gb      (/ms)
    :b3      (/ms)
    Ga2      (/ms)
    h1      (1)
    h2      (1)
    
:    a11     (/ms) : ==a1 = k1*(phi/phi0)
:    a12     (/ms) : ==a2 (const)
:    b1      (/ms) : ==b1 (const)
:    a2      (/ms)      : ==a3 = a30 + a31*ln(1+phi/phi0)
:    b2      (/ms)      : ==b2 = b20 + b21*ln(1+phi/phi0)
:    a3      (/ms)      : ==a4 (const)
:    b31     (/ms)     : ==b3 (const)
:    b32     (/ms)      : ==b4 = b40*(phi/phi0)
:    a4      (/ms)      : ==a6 (const)
    
    fphi    (1) : Fractional conductance
    fv      (1) : Fractional conductance
    :g       (1) : Fractional conductance
    v       (mV) 
    :on            : Flag for light
    tally         : Number of pulses left to deliver
    i       (nA)
    :IRhO    (nA)
    :iChR
}


STATE { C1 I1 O1 O2 I2 C2 }


BREAKPOINT {
    SOLVE kin METHOD sparse
    fphi    = O1+gam*O2                     : light dependency op=[0:1] 
    fv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : voltage dependency (equal 1 at Vm=-70mV)
    :g   = fphi*fv            : total relative conductance
    :i   = gbar*g*(v-E)*(1e-6)       : +ve ELECTRODE_CURRENT depolarises the cell :IRhO
    i       = g0*fphi*fv*(v-E)*(1e-6)       : Photocurrent
    :iChR   = -i                : for plotting (by convention depolar. curr. are -ve)
}


INITIAL {   : Initialise variables
: State occupancy proportions - starts in C1
    C1  = 1     : Closed state 1 (Ground state)
    I1  = 0     : Intermediate state 1
    O1  = 0     : Open state 1 (High conductivity)
    O2  = 0     : Open state 2 (Low conductivity)
    I2  = 0     : Intermediate state 2
    C2  = 0     : Closed state 2

    :flux    = 0
    phi = 0
    rates(phi) : necessary?
    i   = 0
    :IRhO    = 0
    tally   = nPulses : Set the tally to the number of pulses required 
    if (tally > 0) {
        net_send(delD, 1)
    :    on = 0
        tally = tally - 1
    }
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
    CONSERVE C1 + I1 + O1 + O2 + I2 + C2 = 1 : Ignored when CVode is active. Add a normalisation function instead?
}


PROCEDURE rates(phi) { : Define equations for calculating transition rates

    if (phi>0) {
        h1 = 1/(1+pow(phi_m,p)/pow(phi,p))  : pow(phi,p)/(pow(phi,p) + pow(phi_m,p))
        h2 = 1/(1+pow(phi_m,q)/pow(phi,q))  : pow(phi,q)/(pow(phi,q) + pow(phi_m,q))
    } else {
        h1 = 0
        h2 = 0
    }

    
    :a1 = a10*(phi/phi0)
    Ga1 = k1 * h1
    : a2 (const)
    :a3 = a30 + a31*log(1+phi/phi0)
    Gf = Gf0 + kf * h2
    : a4 (const)
    : a6 (const)
    
    : b1 (const)
    :b2 = b20 + b21*log(1+phi/phi0)
    Gb = Gb0 + kb * h2
    : b3 (const)
    :b4 = b40*(phi/phi0)
    Ga2 = k2 * h1
    
    
:    a11 = e1*flux/Er0
:    a12 = 1
:    b1  = b10

:    a2  = a2dark + a2light*log(1+flux/Er0)  
:    b2  = b2dark + b2light*log(1+flux/Er0)  

:    a3  = a30
:    b31 = 1
:    b32 = e3*flux/Er0
:    a4  = a40
}


: Add functions for calculating transition rates outside of simulations
:FUNCTION A1(phi) {
:    A1 = a10*(phi/phi0)
:}
: ...


NET_RECEIVE (w) {
    if (flag == 1) { : ignore any but self-events with flag == 1. This may not be necessary... see e.g. nrn/src/nrnoc/intfire1.mod
    :    if (on == 0) { : Turn the light on 
            phi = phiOn :flux = Er 
    :        on = 1
            net_send(onD, 0) :net_send(onD, 1) : Schedule the next off phase
        } else { : Turn the light off
            phi = 0 :flux = 0
    :        on = 0
            if (tally > 0) { : Schedule the next on phase
                net_send(offD, 1)
                tally = tally - 1
            }
        }
    :}
}
