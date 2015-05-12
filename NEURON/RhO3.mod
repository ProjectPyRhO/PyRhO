TITLE Nikolic 3-state rhodopsin model

NEURON  {
	POINT_PROCESS RhO3
	NONSPECIFIC_CURRENT i
	RANGE E, v0, v1, g, gph, gv, i             :gbar, gam
	RANGE k, Gd, Gr0, Gr1, p 
	RANGE phiOn, phi, phim, delD, onD, offD, nPulses      :phi0
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
    phiOn   = 1e12      :(photons/sec mm2) : Flux : Irradiance :(mW/mm^2) 	 _______
	delD 	= 25	    (ms) <0, 1e9>	: delay before ON phase 	        |  ON   |  OFF
	onD 	= 100	    (ms) <0, 1e9>  	: duration of each ON phase	<-delD->|<-onD->|<-offD->
	offD 	= 50	    (ms) <0, 1e9>  	: duration of each OFF phase________|       |________
	nPulses = 1		    (1)	 <0, 1e3>	: num pulses to deliver	            <-- one pulse -->

: Illumination constants   
:   phi0    = 1e14      (photons/sec mm2)     : Normalising Flux Constant
:	lambda 	= 470	    : (nm)
    phim    = 1e16
  
: Conductance    
    E       = 0         (mV)    : Channel reversal potential
	g    	= 1         (pS)	: defined in the main prog as HR_expression*Area

: Inward rectifier conductance    	gv = (1-exp(-(v-E)/v0))*v1/(v-E) := 1 at Vm=-70mV
	v0 		= 43	    (mV)
	v1 		= 4.1       (mV)
    
: State transition rate parameters (/ms)
	k 		= 0.28      (/ms):5.45e-15	(/ms):(mm2/photons) : Quantum efficiency * number of photons absorbed by a RhO molecule per unit time 0.5*1.2e-14 /1.1
    p       = 0.4       (1)
	Gd 		= 0.0909    (/ms) <0, 1e9> : @ 1mW mm^-2
	Gr0     = 0.0002	(/ms) <0, 1e9> : tau_r,dark = 5-10s p405 Nikolic et al. 2009
	Gr1     = 6.06e-3	(/ms) <0, 1e9> : (Gr,light @ 1mW mm^-2)
}


ASSIGNED {
    P   (/ms)
    Gr  (/ms)

    gph    
	gv
	:g    (pS)
    i    (nA)
	v    (mV) 
	phi		: "flux" should now be "intensity". Er is normalised to be dimensionless
	tally         : how many more to deliver
}


STATE { C O D }


BREAKPOINT {
	SOLVE kin METHOD sparse
	gph = O		                        : Light dependency on open state O=[0:1] 
	gv 	= (1-exp(-(v-E)/v0))*v1/(v-E)	: Voltage dependency (equal 1 at Vm=-70mV)
	i   = g*gph*gv*(v-E)*(1e-6)         : Photocurrent
}


INITIAL { 	: Initialise variables to fully dark-adapted state
: State occupancy proportions - starts in C
	C 	= 1
	O 	= 0
	D 	= 0
    
	phi = 0
    rates(phi) : necessary?
    i   = 0
    
	tally 	= nPulses
	if (tally > 0) {
		net_send(delD, 1)   : Schedule an on event at t+delD
		tally = tally - 1
	}
}

:DERIVATIVE deriv {
:	rates(flux)
:	C' = (Gr * D) - (P * C)
:	O' = (P * C) - (Gd * D)
:	D' = (Gd * O) - (Gr * D)
:}

KINETIC kin {
	rates(phi)
    ~ C <-> O (P, 0)
    ~ O <-> D (Gd, 0)
    ~ D <-> C (Gr, 0)
    CONSERVE C + O + D = 1
}


PROCEDURE rates(phi) {	: Define equations for calculating transition rates
    if (phi>0) {
        :P = phi * k 			: k = quantum efficiency
        P = k * 1/(1+pow(phim,p)/pow(phi,p))    : pow(phi,p)/(pow(phi,p) + pow(phim,p))
        Gr = Gr0 + Gr1
    } else {
        P = 0
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
