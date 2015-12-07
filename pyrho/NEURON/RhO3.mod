TITLE Nikolic 3-state rhodopsin model

NEURON  {
	POINT_PROCESS RhO3
	NONSPECIFIC_CURRENT i
	RANGE i, E, v0, v1, g0 :, fphi, fv
	RANGE k_a, k_r, Gd, Gr0, p, q
	RANGE phiOn, phi, phi_m, delD, onD, offD, nPulses
}


UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
    (pS) = (picosiemens)
    :(photons) = (1)
}


PARAMETER {	: Initialise parameters to defaults. These may be changed through hoc files

: Illumination
    phiOn   = 1e18      :(photons/sec mm2)  : Flux [Irradiance :(mW/mm^2)] 	     _______
	delD 	= 25	    (ms)    <0, 1e9>	: delay before ON phase 	        |  ON   |  OFF
	onD 	= 100	    (ms)    <0, 1e9>  	: duration of each ON phase	<-delD->|<-onD->|<-offD->
	offD 	= 50	    (ms)    <0, 1e9>  	: duration of each OFF phase________|       |________
	nPulses = 1		    (1)	    <0, 1e3>    : num pulses to deliver	            <-- one pulse -->

: Illumination constants   
:	lambda 	= 470	    :(nm)
    phi_m   = 1e16      :(photons/sec mm2)  : Hill Constant
  
: Conductance    
    E       = 0         (mV)                : Channel reversal potential
	g0    	= 1         (pS)	            : defined in the main prog as HR_expression*Area

: Inward rectifier conductance    	        fv = (1-exp(-(v-E)/v0))*v1/(v-E) := 1 at Vm=-70mV
	v0 		= 43	    (mV)
	v1 		= 4.1       (mV)
    
: State transition rate parameters (/ms)
	k_a		= 0.28      (/ms)               : Quantum efficiency * number of photons absorbed by a RhO molecule per unit time 0.5*1.2e-14 /1.1
    p       = 0.4       (1)                 : Hill Coefficient
    k_r		= 0.28      (/ms)
    q       = 0.4       (1)                 : Hill Coefficient
	Gd 		= 0.0909    (/ms) <0, 1e9>      : @ 1mW mm^-2
	Gr0     = 0.0002	(/ms) <0, 1e9>      : tau_r,dark = 5-10s p405 Nikolic et al. 2009
}


ASSIGNED {
	phi             : "flux" should now be "intensity". Er is normalised to be dimensionless
    Ga      (/ms)
    Gr      (/ms)

    fphi    (1)
	fv      (1)
	:g      (pS)
    v       (mV) 
    i       (nA)
	tally           : how many more pulses to deliver
}


STATE { C O D }


BREAKPOINT {
	SOLVE kin METHOD sparse
	fphi    = O		                        : Light dependency on open state O=[0:1] 
	fv 	    = (1-exp(-(v-E)/v0))*v1/(v-E)	: Voltage dependency (equal 1 at Vm=-70mV)
	i       = g0*fphi*fv*(v-E)*(1e-6)       : Photocurrent
}


INITIAL { 	: Initialise variables to fully dark-adapted state
: State occupancy proportions - starts in C
	C 	    = 1
	O 	    = 0
	D 	    = 0
    
    phi     = 0
    rates(phi) : necessary?
    i       = 0
    
	tally 	= nPulses
	if (tally > 0) {
		net_send(delD, 1)   : Schedule an on event at t+delD
		tally = tally - 1
	}
}

:DERIVATIVE deriv {
:	rates(flux)
:	C' = (Gr * D) - (Ga * C)
:	O' = (Ga * C) - (Gd * D)
:	D' = (Gd * O) - (Gr * D)
:}

KINETIC kin {
	rates(phi)
    ~ C <-> O (Ga, 0)
    ~ O <-> D (Gd, 0)
    ~ D <-> C (Gr, 0)
    CONSERVE C + O + D = 1
}


PROCEDURE rates(phi) {	        : Define equations for calculating transition rates
    if (phi>0) {
        :Ga = phi * k 			: k = quantum efficiency
        Ga = k_a * 1/(1+pow(phi_m,p)/pow(phi,p))    : pow(phi,p)/(pow(phi,p) + pow(phi_m,p))
        Gr = Gr0 + k_r * 1/(1+pow(phi_m,q)/pow(phi,q))    : pow(phi,q)/(pow(phi,q) + pow(phi_m,q))
    } else {
        Ga = 0
        Gr = Gr0
    }
}  


NET_RECEIVE (w) {
	if (flag == 1) {            : Switch light on
        phi = phiOn 
        net_send(onD, 0)        : Schedule an off event at t+onD
    } else {                    : Switch light off
        phi = 0   
        if (tally > 0) {        : More pulses to deliver
            net_send(offD, 1)   : Schedule an on event at t+offD
            tally = tally - 1
        }
    }
}
