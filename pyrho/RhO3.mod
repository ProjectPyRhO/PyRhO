TITLE Nikolic 3-state rhodopsin model

NEURON  {
	SUFFIX RhO3
	POINT_PROCESS RhO
	NONSPECIFIC_CURRENT iRhO
	RANGE phi, g, gbar, gv, iRhO    : Er, flux
	RANGE k Gd Gr_dark Gr_light 
	RANGE del, ton, toff, num
}


UNITS {
	(nA) = (nanoamp)
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {	: Initialise parameters to defaults. These may be changed through hoc files

: Illumination
    phi     = 1e15 (photons/(sec*mm2)) : Flux
::	Er 		= 1	:(mW/mm2)			: irradiance 					_______
	del 	= 30	(ms)			: delay before ON phase 	   |  ON   |  OFF
	ton 	= 100	(ms) <0, 1e9>  	: duration of ON phase	<-del->|<-ton->|<-toff->
	toff 	= 50	(ms) <0, 1e9>  	: duration of OFF phase	_______|       |________
	num 	= 1		(1)				: num pulses to deliver	       <-- one pulse -->
    
: State transition rate parameters (/ms)
	k 			= 0.5*1.2e-14/1.1 	(mm2/photons) : Quantum efficiency * number of photons absorbed by a RhO molecule per unit time
	Gd 			= 1/11	    (/ms) <0, 1e9> : @ 1mW mm^-2
	Gr_dark 	= 1/5000	(/ms) <0, 1e9> : tau_r,dark = 5-10s p405 Nikolic et al. 2009
	Gr_light 	= 1/165		(/ms) <0, 1e9> : (Gr,light @ 1mW mm^-2)
       
: conductance    
:	gama 	= 0.05
:self.g = 100e-3 * 1.67e5    # [pS] g1 (100fS): Conductance of an individual ion channel * N (100000)
	gbar 	= 1 	: defined in the main prog as HR_expression*Area

: light constants   
    phi0    = 1e14  (photons/(sec*mm2))     : Normalising Flux Constant
::	Er0 	= 0.1	: (mW/mm2)
:	lambda 	= 470	: (nm)
:	flux0 	= 2e16	: (photons/(sec*cm2))

: reversal potential for NpHr is about -400mV        e_rev = -400        :(mV)  
	E 		= 0 	(mV)	: Channel reversal potential
	
: inward rectifier ChR2 conductance    	gv = (1-exp(-(v-E)/v0))*v1/(v-E) := 1 at Vm=-70mV
	v0 		= 43	(mV)
	v1 		= -4.1  (mV)     
}


ASSIGNED {
	a11   (/ms)
	a12   (/ms)  
	b1   (/ms)
	a2
	b2
	a3
	b31
	b32
	a4
	
	gv
	g    
	iRhO    (nA)
	v    (mV) 
	flux		: "flux" should now be "intensity". Er is normalised to be dimensionless
    on
	tally         : how many more to deliver
}


STATE { C O D }


BREAKPOINT {
	SOLVE kin METHOD sparse
:	gph 	= o1+gama*o2		: light dependency op=[0:1] 
	gv 		= (1-exp(-(v-E)/v0))*v1/(v-E)	: voltage dependency (equal 1 at Vm=-70mV)
:	g 		= gph*gv 			: total relative conductance
:	i 		= gbar*g			: +ve ELECTRODE_CURRENT depolarises the cell 
:	iChR	= -i				: for plotting (by convention depolar. curr. are -ve)
    g = gbar*gv*O               : total relative conductance
	iRhO = g*(v-E)
}


INITIAL { 	: Initialise variables to fully dark-adapted state
: State occupancy proportions - starts in C
	C 	= 1
	O 	= 0
	D 	= 0
	flux 	= 0
    rates(flux) : necessary?
	iRhO 	= 0
    
	tally 	= num
	if (tally > 0) {
		net_send(del, 1)
		on = 0
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
	rates(flux)
    ~ C <-> O (P, 0)
    ~ O <-> D (Gd, 0)
    ~ D <-> C (Gr, 0)
    CONSERVE C + O + D = 1
}


PROCEDURE rates() {	: Define equations for calculating transition rates
	P = flux * k 			: k = quantum efficiency
	Gr = Gr_dark + Gr_light
}  


NET_RECEIVE (w) {
	if (flag == 1) {        : ignore any but self-events with flag == 1
		if (on == 0) {      : turn it on   
			flux = Er 
			on = 1
			net_send(ton, 1)    : prepare to turn it off
		} else {            : turn it off
			flux = 0   
			on = 0
			if (tally > 0) {    : prepare to turn it on again
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}