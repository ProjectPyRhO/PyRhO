
### Fitting Hyper parameters

# Time window for steady-state averaging
# IPI curve initial parameters
# fV initial parameters
# Kinetics initial parameters
# Optimisation routine initial parameters
# ...


### Move this to models.py and add .setParams() method
from lmfit import Parameters, Parameter
from collections import OrderedDict, defaultdict
#from pyrho.utilities import irrad2flux, flux2irrad
from copy import deepcopy


### Hyperparameters

tFromOff = 50  # Time [ms] to start the sample window before the end of the pulse for Iss

# Optimisation initialisation values
#p0fV = (40,4,1)#25000)#,1)      # v0,v1,E,G
#p0FV = (40, 4, 1, 0.025)
#p0IPI = (0.5,4000,-1) # a*exp(-t/b)+c #(-1e-8,400,-1e-7)

# Used if plotKinetics
p0on = (-0.1,2,-1) # a*exp(-t/b)+c
p0inact = (-0.5,25,-0.5) # a*exp(-t/b)+c
p0off = (-0.1,7.5,-0.1,35,-0.1) # a1*exp(-t/tau1)+a2*exp(-t/tau2)+I_ss

### Add default kinetics parameters
### On phase
pOn = Parameters()
#pOn.add('a0', value=0, expr='-a2')
pOn.add('a1', value=1, min=1e-9)
pOn.add('a2', value=0.1, min=1e-9)
pOn.add('a0', value=0, expr='-a2')
pOn.add('tau_act', value=5, min=1e-9)
pOn.add('tau_deact', value=50, min=1e-9)

### Off phase
#Iss = pOn['a0'].value + pOn['a1'].value
pOffSing = Parameters() # copy.deepcopy(pOn)

# Single exponential
pOffSing.add('a0', value=0)#, expr='{}-a1-a2'.format(Iss))
pOffSing.add('a1', value=0, vary=True)
pOffSing.add('a2', value=-0, vary=False)
pOffSing.add('Gd1', value=0.1)#, min=1e-9)
pOffSing.add('Gd2', value=0, vary=False) #, expr='Gd1')#, min=1e-9)

# Double exponential
pOffDoub = Parameters()
pOffDoub.add('a0', value=0, vary=False)
pOffDoub.add('a1', value=0.1)
pOffDoub.add('a2', value=-0.1)#, expr='{}-a0-a1'.format(Iss))
pOffDoub.add('Gd1', value=0.1)#, min=1e-9)
pOffDoub.add('Gd2', value=0.01)#, vary=True) #, expr='Gd1')#, min=1e-9)





### Default model parameters

modelParams = OrderedDict([('3',Parameters()),('4',Parameters()),('6',Parameters())])
modelList = list(modelParams) # List of keys: list(modelParams.keys()) #This could be removed
stateLabs = {3:'Three', '3':'Three', 4:'Four', '4':'Four', 6:'Six', '6':'Six'}

modelFits = OrderedDict([   ('3', OrderedDict([('ChR2',Parameters()), ('NpHR',Parameters()), ('ArchT',Parameters())])), 
                            ('4', OrderedDict([('ChR2',Parameters())])), 
                            ('6', OrderedDict([('ChR2',Parameters())]))])

### Replace with defaultdict with default=key
modelLabels = OrderedDict([('E','E'), ('g0','g_0'), ('p','p'), ('k_a','k_a'), ('k_r','k_r'), ('phi_m','\phi_m'), ('Gd','G_d'), ('Gr0','G_{r0}'), ('v0','v_0'), ('v1','v_1'), 
                            ('gam','\gamma'), ('k1','k_1'), ('k2','k_2'), ('Gf0','G_{f0}'), ('Gb0','G_{b0}'), ('k_f','k_f'), ('k_b','k_b'), ('q','q'), 
                            ('Gd1','G_{d1}'), ('Gd2','G_{d2}'), ('Go1','G_{o1}'), ('Go2','G_{o2}'),
                            ('phi','\phi'), ('v','v')])

# Units for model parameters (for Brian)
# Replace with http://pythonhosted.org/NeuroTools/parameters.html
from brian2.units.allunits import psiemens, second, mole #*
from brian2.units.stdunits import *
pS = psiemens
sec = second
# Units used: ms, mm, mV, Hz,       # Nonstd: psiemens, second, mole
                                    #           nS, uS

modelUnits = OrderedDict([('g0',pS), ('gam',1), ('k_a',ms**-1), ('k_r',ms**-1), ('phi_m',mm**-2*second**-1), ('p',1), ('Gd',ms**-1), ('Gr0',ms**-1), 
                            ('k1',ms**-1), ('k2',ms**-1), ('Gf0',ms**-1), ('Gb0',ms**-1), ('k_f',ms**-1), ('k_b',ms**-1), ('q',1), 
                            ('Gd1',ms**-1), ('Gd2',ms**-1), ('Go1',ms**-1), ('Go2',ms**-1), ('E',mV), ('v0',mV), ('v1',mV), 
                            ('phi',mm**-2*second**-1), ('v',mV)])

#paramUnits
unitLabels = OrderedDict([('g0','pS'), ('gam',''), ('k_a','ms^-1'), ('k_r','ms^-1'), ('phi_m','ph./mm^2/s'), ('p',''), ('Gd','ms^-1'), ('Gr0','ms^-1'), #, ('k','ms^-1'), ('Gr1','ms^-1')
                            ('k1','ms^-1'), ('k2','ms^-1'), ('Gf0','ms^-1'), ('Gb0','ms^-1'), ('k_f','ms^-1'), ('k_b','ms^-1'), ('q',''), 
                            ('Gd1','ms^-1'), ('Gd2','ms^-1'), ('Go1','ms^-1'), ('Go2','ms^-1'), ('E','mV'), ('v0','mV'), ('v1','mV'),
                            ('phi','ph./mm^2/s'), ('v','mV')])

                            
####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80

#               (Name,    Value,  Vary, Min,  Max,  Expr=Units)                
modelFits['3']['ChR2'].add_many( # Depolarising: passively transports Na+, H+, K+ and Ca2+ down their electrochemical gradients
                ('g0',    1.57e5, True, 0.001,  1e6,  None),
                ('phi_m', 5e17,   True, 1e15,   1e19, None),
                ('k_a',   5,      True, 0.001,  1000, None),
                ('k_r',   0.1,    True, 0.001,  1000, None),
                ('p',     0.8,    True, 0.1,    5,    None),
                ('q',     0.25,   True, 0.1,    5,    None),
                ('Gd',    0.104,  True, 0.0001, 1,    None),
                ('Gr0',   0.0002, True, 0.0001, 0.1,  None),
                ('E',     0,      True, -1000,  1000, None),
                ('v0',    43,     True, -1e15,  1e15, None),
                ('v1',    17.1,   True, -1e15,  1e15, None))
                
modelFits['3']['NpHR'].add_many( # Hyperpolarising: pumps chloride ions into the cell
                ('g0',    1.57e5, True, 0.001,  1e6,  None),
                ('phi_m',  1.32e18,True, 1e15,   1e19, None),
                ('k_a',     0.01,   True, 0.001,  1000, None),
                ('k_r',     0.01,   True, 0.001,  1000, None),
                ('p',     0.793,  True, 0.1,    5,    None),
                ('q',     0.793,  True, 0.1,    5,    None),
                ('Gd',    0.1,  True, 0.0001, 1,    None),
                ('Gr0',   0.0002, True, 0.0001, 0.1,  None),
                ('E',     -400,   True, -1000,  1000, None),
                ('v0',    43,     True, -1e15,  1e15, None),
                ('v1',    17.1,   True, -1e15,  1e15, None))
                
modelFits['3']['ArchT'].add_many( # Hyperpolarising: actively extrudes Hydrogen ions
                ('g0',    1.57e5, True, 0.001,  1e6,  None),
                ('phi_m',  1.32e18,True, 1e15,   1e19, None),
                ('k_a',     0.01,   True, 0.001,  1000, None),
                ('k_r',     0.01,   True, 0.001,  1000, None),
                ('p',     0.793,  True, 0.1,    5,    None),
                ('q',     0.793,  True, 0.1,    5,    None),
                ('Gd',    0.1,  True, 0.0001, 1,    None),
                ('Gr0',   0.001, True, 0.0001, 0.1,  None),
                ('E',     0,   True, -1000,  1000, None),
                ('v0',    43,     True, -1e15,  1e15, None),
                ('v1',    17.1,   True, -1e15,  1e15, None))
                
    
modelFits['4']['ChR2'].add_many(
                ('g0',      1.14e5, True, 0.001,1e15,   None),
                ('gam',     0.00742,True, 0.0,  1,      None),
                ('phi_m',   2.33e17,True, 1e15, 1e19,   None),
                ('k1',      4.15,   True, 0.001,1e5,    None), #3
                ('k2',      0.868,  True, 0.001,1e5,    None), #1.5
                ('p',       0.833,  True, 0.1,  5,      None),
                ('Gf0',     0.0373, True, 0,    1e3,    None), #e12d
                ('k_f',     0.0581, True, 0.001,1e3,    None), #c1
                ('Gb0',     0.0161, True, 0,    1e3,    None), #e21d
                ('k_b',     0.063,  True, 0.001,1e3,    None), #c2
                ('q',       1.94,   True, 0.1,  5,      None),
                ('Gd1',     0.105,  True, 0.01, 1,      None),
                ('Gd2',     0.0138, True, 0.01, 1,      None),
                ('Gr0',     0.00033,True, 1e-6, 1,      None), #Gr #0.0004
                ('E',       0,      True, -1000,1000,   None),
                ('v0',      43,     True, -1e15,1e15,   None),
                ('v1',      17.1,   True, -1e15,1e15,   None))    
              

modelFits['6']['ChR2'].add_many(                
                ('g0',      2.52e4, True, 0.0,  1e15, None),
                ('gam',     0.0161, True, 0.0,  1,    None), # Max=1 if gO1 >= gO2
                ('phi_m',   3.54e17,True, 1e15, 1e19, None),
                ('k1',      13.4,   True, 0.0,  1000, None),
                ('k2',      2.71,   True, 0.0,  1000, None),
                ('p',       0.985,  True, 0.1,  5,    None),
                ('Gf0',     0.0389, True, 0.0,  1000, None),
                ('k_f',     0.103,  True, 0.0,  1000, None),
                ('Gb0',     0.0198, True, 0.0,  1000, None),
                ('k_b',     0.139,  True, 0.0,  1000, None),
                ('q',       1.58,   True, 0.1,  5,    None),
                ('Go1',     2,      True, 0.0,  1000, None),
                ('Go2',     0.0567, True, 0.0,  1000, None),
                ('Gd1',     0.112,  True, 0.0,  1000, None),
                ('Gd2',     0.0185, True, 0.0,  1000, None),
                ('Gr0',     0.00033,True, 0.0,  1000, None), #0.00163
                ('E',       0,      True, -1000,1000, None),
                ('v0',      43,     True, -1e15, 1e15,None),
                ('v1',      17.1,   True, -1e15, 1e15,None))


defaultOpsinType = 'ChR2'
rhoType = defaultOpsinType # Set this when selecting 
modelParams['3'] = modelFits['3'][defaultOpsinType]
modelParams['4'] = modelFits['4'][defaultOpsinType]
modelParams['6'] = modelFits['6'][defaultOpsinType]

unitPrefixes = {} ### Use a units library to convert between different prefixes


#Params = OrderedDict([('model', OrderedDict()), ('protocol', OrderedDict()), ('simulator', OrderedDict())])
# p, q: Hill coefficients
# phi_m: Hill constant

#from numpy import inf, nan, isfinite 
import numpy as np

### This needs serious refactoring! Create a superclass which hands off attribute/method calls to a Parameter(s)() attribute by default or looks in self for other properties
#class Parameters(OrderedDict):
class PyRhOparameters(Parameters):
    '''These classes are adapted from LMFIT since changes between 
    0.8.0 and 0.9.2 stopped to ability to set lists as values'''
    def __deepcopy__(self, memo):
        _pars = PyRhOparameters()

        # we're just about to add a lot of Parameter objects to the newly
        parameter_list = []
        for key, par in self.items():
            if isinstance(par, PyRhOparameter):
                param = PyRhOparameter(name=par.name,
                                        value=par.value,
                                        min=par.min,
                                        max=par.max)
                #param.vary = par.vary
                #param.stderr = par.stderr
                #param.correl = par.correl
                #param.init_value = par.init_value
                #param.expr = par.expr
                parameter_list.append(param)

        _pars.add_many(*parameter_list)

        return _pars
    
    def __setitem__(self, key, par):
        #if key not in self:
        #    if not valid_symbol_name(key):
        #        raise KeyError("'%s' is not a valid Parameters name" % key)
        if par is not None and not isinstance(par, (Parameter, PyRhOparameter)):
            raise ValueError("'%s' is not a Parameter" % par)
        OrderedDict.__setitem__(self, key, par)
        par.name = key
        #par._expr_eval = self._asteval
        #self._asteval.symtable[key] = par.value
    
    def add_many(self, *parlist):
        """
        Convenience function for adding a list of Parameters.
        Parameters
        ----------
        parlist : sequence
            A sequence of tuples, or a sequence of `Parameter` instances. If it
            is a sequence of tuples, then each tuple must contain at least the
            name. The order in each tuple is the following:
                name, value, vary, min, max, expr
        Example
        -------
        p = Parameters()
        # add a sequence of tuples
        p.add_many( (name1, val1, True, None, None, None),
                    (name2, val2, True,  0.0, None, None),
                    (name3, val3, False, None, None, None),
                    (name4, val4))
        # add a sequence of Parameter
        f = Parameter('name5', val5)
        g = Parameter('name6', val6)
        p.add_many(f, g)
        """
        for para in parlist:
            if isinstance(para, PyRhOparameter):
                self.__setitem__(para.name, para)
            else:
                param = PyRhOparameter(*para)
                self.__setitem__(param.name, param)
    
    def valuesdict(self):
        """
        Returns
        -------
        An ordered dictionary of name:value pairs for each Parameter.
        This is distinct from the Parameters itself, as it has values of
        the Parameter values, not the full Parameter object.
        """

        return OrderedDict(((p.name, p.value) for p in self.values()))
                
class PyRhOparameter(object):
    
    def __init__(self, name=None, value=None, min=-np.inf, max=np.inf, units=None, latex=None, descr=None): #, unitsLabel=None
        self.name = name
        self._val = value
        self._min = -np.inf
        self._max = np.inf
        self.min = min
        self.max = max

        self.units = units
        #self.label = label
        self.unitsLabel = str(self.units) #unitsLabel
        self.latex = latex
        self.descr = descr
        self.constant = True
        self._init_bounds()
    
    
    def set(self, value=None, vary=None, min=None, max=None, expr=None):
        """
        Set or update Parameter attributes (adapted from LMFIT).
        Parameters
        ----------
        value : float, optional
            Numerical Parameter value.
        vary : bool, optional
            Whether the Parameter is fixed during a fit.
        min : float, optional
            Lower bound for value(s). To remove a lower bound you must use -np.inf
        max : float, optional
            Upper bound for value(s). To remove an upper bound you must use np.inf
        """

        #self.__set_expression(expr)
        if value is not None:
            self._val = value
        if vary is not None:
            self.vary = vary
        if min is not None:
            self.min = min
        if max is not None:
            self.max = max
        
    '''    
    def set(self, value=None):
        if value is not None:
            self._val = value
    '''
    
    def _init_bounds(self):
        """make sure initial bounds are self-consistent"""
        # _val is None means - infinity.
        # _val is None means - infinity.
        if self._val is not None:
            if isinstance(self._val, str):
                return
            elif isinstance(self._val, (list, tuple)):
                self._clipList(self._val)
            else:
                if self.max is not None and self._val > self.max:
                    self._val = self.max
                if self.min is not None and self._val < self.min:
                    self._val = self.min
        elif self.min is not None: #and self._expr is None:
            self._val = self.min
        elif self.max is not None: #and self._expr is None:
            self._val = self.max
        #self.setup_bounds()
    
    def _clipList(self, values):
        for ind, val in enumerate(values):
            if isinstance(val, str):
                return
            elif isinstance(val, (list, tuple)): # Nested list e.g. cycles
                self._clipList(val)
            else:                    
                if self.max is not None and val > self.max:
                    values[ind] = self.max
                if self.min is not None and val < self.min:
                    values[ind] = self.min
        
    def get_max(self):
        return self._max

    def set_max(self, val):
        if val is None:
            val = np.inf
        self._max = val
        if self.min > self.max:
            self._min, self._max = self.max, self.min
        if np.isclose(self.min, self.max, atol=1e-13, rtol=1e-13):
            raise ValueError("Parameter '%s' has min == max" % self.name)

    def get_min(self):
        return self._min

    def set_min(self, val):
        if val is None:
            val = -np.inf
        self._min = val
        if self.min > self.max:
            self._min, self._max = self.max, self.min
        if np.isclose(self.min, self.max, atol=1e-13, rtol=1e-13):
            raise ValueError("Parameter '%s' has min == max" % self.name)

    min = property(get_min, set_min)
    max = property(get_max, set_max)
    

    
    def _getval(self):
        return self._val
    
    @property
    def value(self):
        return self._getval() # self._val # 
        
    @value.setter
    def value(self, val):
        self._val = val
        self._init_bounds()
        
    def __repr__(self):
        s = []
        if self.name is not None:
            s.append("'%s'" % self.name)
        sval = repr(self._getval())
        #if not self.vary and self._expr is None:
        #    sval = "value=%s (fixed)" % sval
        #elif self.stderr is not None:
        #    sval = "value=%s +/- %.3g" % (sval, self.stderr)
        s.append(sval)
        s.append("bounds=[%s:%s]" % (repr(self.min), repr(self.max)))
        #if self._expr is not None:
        #    s.append("expr='%s'" % self.expr)
        return "<Parameter %s>" % ', '.join(s)
    
    def __str__(self):
        return self.__repr__()
    
    def _repr_latex_(self):
        if self.latex is not None:
            s = self.latex
        else:
            s = self.name
        if self._val is not None:
            '''
            if isinstance(self._val, (list, tuple, np.ndarray)):
                #v = ",\,".join(str(self._val))
                #v = "\[" + v + "\]"
                v = str(self._val)
            else:
                v = self._val
            '''
            v = str(self._val)
        if self.units is not None:
            #u = "[{u}]".format(u=self.units._latex())
            u = self.units._latex()
            s = "\,".join([s, '=', v, u])
        return "$" + s + "$" # texIt(s)
    
    def latex(self):
        from IPython.display import Math
        return Math(self._repr_latex_())
#Params['g0'] = PyRhOparameter('g0', 2.5e4, psiemens, 'pS', 'g_0', 'Biological scaling factor for rhodopsin conductance')


####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80


### Protocols to be included in the next version:
### - Temperature (Q10)
### - pH (intracellular and extracellular)
### - Wavelength (lambda)

#protParams = OrderedDict([('step',Parameters()), ('delta',Parameters()), ('sinusoid',Parameters()), ('chirp',Parameters()), ('ramp',Parameters()), ('rectifier',Parameters()), ('shortPulse',Parameters()), ('recovery',Parameters()), ('custom',Parameters())])

protParams = OrderedDict([('step',PyRhOparameters()), ('delta',PyRhOparameters()), ('sinusoid',PyRhOparameters()), ('chirp',PyRhOparameters()), ('ramp',PyRhOparameters()), ('rectifier',PyRhOparameters()), ('shortPulse',PyRhOparameters()), ('recovery',PyRhOparameters()), ('custom',PyRhOparameters())])

protList = list(protParams) # List of keys #This could be removed

protParamLabels = OrderedDict([ ('phis', '\mathbf{\phi}'), 
                                ('Vs', '\mathbf{\mathrm{V}}'), 
                                ('delD', '\Delta t_{delay}'), 
                                ('onD', '\Delta t_{on}'), 
                                ('totT', 'T_{total}'), 
                                ('cycles', 'cycles'),
                                ('phi0', '\phi_0'), 
                                ('fs', '\mathbf{f}'), 
                                ('f0', 'f_0'), 
                                ('fT', 'f_T'), 
                                ('linear', 'linear'), 
                                ('startOn', '\phi_{t=0}>0'), 
                                #('phi_ton', '\phi_{t=0}'), 
                                ('pDs', '\mathbf{\Delta t_{on}}'), 
                                ('IPIs', '\mathbf{\Delta t_{off}}'), 
                                ('phi_ft', '\phi(t)') ])

protUnitLabels = defaultdict(lambda: '')
protUnitLabels['phis'] = 'ph./mm^2/s'
protUnitLabels['phi0'] = 'ph./mm^2/s'
#protUnitLabels['phi_ton'] = 'ph./mm^2/s' ### Revise!!!
protUnitLabels['Vs'] = 'mV'
protUnitLabels['delD'] = 'ms'
protUnitLabels['onD'] = 'ms'
protUnitLabels['cycles'] = 'ms'
protUnitLabels['pDs'] = 'ms'
protUnitLabels['IPIs'] = 'ms'
protUnitLabels['totT'] = 'ms'
protUnitLabels['fs'] = 'Hz'
protUnitLabels['f0'] = 'Hz'
protUnitLabels['fT'] = 'Hz'

protParamNotes = OrderedDict([ (prot, defaultdict(lambda: '')) for prot in protList]) 
for prot in protList:
    protParamNotes[prot]['phis'] = 'List of flux values'
    protParamNotes[prot]['Vs'] = 'List of voltage clamp values (if applied)'
    protParamNotes[prot]['delD'] = 'Delay duration before the first pulse' # cycle'
    protParamNotes[prot]['cycles'] = 'List of [on, off] durations for each pulse' # cycle'
    

#Exceptions
protParamNotes['custom']['phi_ft'] = 'Pulse generation function'
protParamNotes['sinusoid']['startOn'] = 'Start at maximum flux (else minimum)' # maximum of flux modulation 
protParamNotes['sinusoid']['phi0'] = 'Constant offset for modulation'
protParamNotes['sinusoid']['fs'] = 'List of modulation frequencies'
protParamNotes['chirp']['linear'] = 'Linear frequency sweep (else exponential)'
protParamNotes['chirp']['startOn'] = 'Start at maximum flux (else minimum)'
protParamNotes['chirp']['phi0'] = 'Constant offset for modulation'
protParamNotes['chirp']['f0'] = 'Starting frequency'
protParamNotes['chirp']['fT'] = 'Ending frequency'
protParamNotes['ramp']['phis'] = 'List of ending flux values'
#protParamNotes['ramp']['phi_ton'] = 'Starting flux value'
protParamNotes['ramp']['phi0'] = 'Constant offset for flux values'
protParamNotes['delta']['cycles'] = ''
protParamNotes['delta']['onD'] = 'On-phase duration'
protParamNotes['delta']['totT'] = 'Total simulation duration'
protParamNotes['shortPulse']['cycles'] = ''
protParamNotes['shortPulse']['pDs'] = 'List of pulse on-phase durations' #'List of cycle on-phase durations'
protParamNotes['shortPulse']['totT'] = 'Total simulation duration'
protParamNotes['recovery']['cycles'] = ''
protParamNotes['recovery']['onD'] = 'Pulse on-phase duration' #'Cycle on-phase duration'
protParamNotes['recovery']['IPIs'] = 'List of pulse off-phase durations' #'List of cycle off-phase durations'
protParamNotes['recovery']['totT'] = 'Total simulation duration'
    
#squarePulses = ['custom', 'delta', 'step', 'rectifier', 'shortPulse', 'recovery'] #{'custom': True, 'delta': True, 'step': True, 'rectifier': True, 'shortPulse': True, 'recovery': True}
#arbitraryPulses = ['custom', 'sinusoid', 'chirp', 'ramp'] #{'custom': True, 'sinusoid': True, 'chirp': True, 'ramp':True} # Move custom here

smallSignalAnalysis = ['delta', 'step', 'sinusoid'] 

protParams['custom'].add_many(('phis',  [1e16,1e17],        0,      None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), #'photons/s/mm^2'
                            ('Vs',      [-70,-20,10],       None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), #'mV'
                            ('delD',    25,                 0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), #'ms'
                            ('cycles',  [[150.,50.]],       0,      None,   ms, 'cycles',           'List of [on, off] durations for each pulse'))#, #'ms'#,

protParams['step'].add_many(('phis',    [1e16,1e17],        0,      None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), #'photons/s/mm^2'
                            ('Vs',      [-70,-40,-10,10,40,70], None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), #'mV'
                            ('delD',    25,                 0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), #'ms'
                            ('cycles',  [[150.,100.]],      0,      None,   ms, 'cycles',           'List of [on, off] durations for each pulse')) #'ms'
                            
protParams['sinusoid'].add_many(('phis',[1e12],             0,      None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), #'photons/s/mm^2'
                            ('phi0',    [0],                None,   None,   mole*mm**-2*second**-1, '\phi_0',        'Constant offset for flux'), #'photons/s/mm^2'
                            ('startOn', True,               False,  True,   1,  '\phi_{t=0}>0',         'Start at maximum flux (else minimum)'),
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}',  'List of voltage clamp values (if applied)'), #'mV'
                            ('fs',      [0.1,0.5,1,5,10],   0,      None,   Hz, '\mathbf{f}',           'List of modulation frequencies'), #'Hz' #50, 100, 500, 1000
                            ('delD',    25,                 0,      1e9,    ms, '\Delta t_{delay}',     'Delay duration before the first pulse'), #'ms'
                            ('cycles',  [[10000.,50.]],     0,      None,   ms, 'cycles',               'List of [on, off] durations for each pulse')) #'ms'

protParams['chirp'].add_many(('phis',   [1e12],             None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2'
                            ('phi0',    [0],                None,   None,   mole*mm**-2*second**-1, '\phi_0',        'Constant offset for flux'), # 'photons/s/mm^2'
                            ('linear',  True,               False,  True,   1,  'linear',           'Linear frequency sweep (else exponential)'), # False := exponential
                            ('startOn', False,              False,  True,   1,  '\phi_{t=0}>0',     'Start at maximum flux (else minimum)'),
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV'
                            ('delD',    100,                0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('cycles',  [[10000.,100.]],    0,      None,   ms, 'cycles',           'List of [on, off] durations for each pulse'), # 'ms'
                            ('f0',      0.1,                0,      None,   Hz, 'f_0',              'Starting frequency'), # 'Hz'
                            ('fT',      1000,               0,      None,   Hz, 'f_T',              'Ending frequency')) # 'Hz'

protParams['ramp'].add_many(('phis',    [1e16,1e17,1e18],   None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2' #1e12,1e13,1e14,1e15,
                            ('phi0',    0,                  None,   None,   mole*mm**-2*second**-1, '\phi_0',        'Constant offset for flux'), # 'photons/s/mm^2'
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV'
                            ('delD',    25,                 0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('cycles',  [[250.,25.]],       0,      None,   ms, 'cycles',           'List of [on, off] durations for each pulse')) # 'ms'#,

protParams['delta'].add_many(('phis',   [1e20],             None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2'
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV'
                            ('delD',    5,                  0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('onD',     1e-3,               0,      1e9,    ms, '\Delta t_{on}',    'On-phase duration'), # 'ms'
                            ('totT',    25.,                0,      None,   ms, 'T_{total}',        'Total simulation duration')) # 'ms'

protParams['rectifier'].add_many(('phis',[1e16],            None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2' # Change to 1e17?
                            ('Vs',      [-100,-70,-40,-10,20,50,80],None,None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV' #[-100,-80,-60,-40,-20,0,20,40,60,80]
                            ('delD',    50,                 0,      1e9,    ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('cycles',  [[250.,100.]],      None,   None,   ms, 'cycles',           'List of [on, off] durations for each pulse')) # 'ms' #,

protParams['shortPulse'].add_many(('phis',[1.5e15],         None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2' #1e12 
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV'
                            ('delD',    25,                 0,      None,   ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('pDs',     [1,2,3,5,8,10,20],  0,      None,   ms, '\mathbf{\Delta t_{on}}',   'List of pulse on-phase durations'), # 'ms' # [0.1, 0.2, 0.5, 1, 2, 5, 10]
                            ('totT',    100.,               0,      None,   ms, 'T_{total}',        'Total simulation duration')) # 'ms'

protParams['recovery'].add_many(('phis',[1e17],             None,   None,   mole*mm**-2*second**-1, '\mathbf{\phi}', 'List of flux values'), # 'photons/s/mm^2'
                            ('Vs',      [-70],              None,   None,   mV, '\mathbf{\mathrm{V}}', 'List of voltage clamp values (if applied)'), # 'mV'
                            ('delD',    100,                0,      None,   ms, '\Delta t_{delay}', 'Delay duration before the first pulse'), # 'ms'
                            ('onD',     100,                0,      None,   ms, '\Delta t_{on}',    'On-phase duration'), # 'ms'
                            ('IPIs',[500,1000,1500,2500,5000,7500,10000],None,None,ms,  '\mathbf{\Delta t_{off}}', 'List of pulse off-phase durations'), # 'ms' 
                            #('IPIs',[0.5,1,1.5,2.5,5,7.5,10],None,None,seconds), # 'ms' 
                            ('totT',    12000,              0,      None,   ms, 'T_{total}',        'Total simulation duration')) # 'ms'


simUnitLabels = defaultdict(lambda: '')
simUnitLabels['dt'] = 'ms'
simUnitLabels['v_init'] = 'mV'

simParamNotes = defaultdict(lambda: '')
simParamNotes['cell'] = 'List of hoc files'
simParamNotes['Vclamp'] = 'Use voltage clamp'
simParamNotes['Vcomp'] = 'Compartment to record from'
simParamNotes['expProb'] = 'Expresssion probability'
simParamNotes['v_init'] = 'Initialisation voltage'
simParamNotes['CVode'] = 'Use variable timestep integrator'
simParamNotes['dt'] = 'Numerical integration timestep'

#simParams = OrderedDict([('Python',Parameters()), ('NEURON',Parameters()), ('Brian',Parameters())])
simParams = OrderedDict([('Python',PyRhOparameters()), ('NEURON',PyRhOparameters()), ('Brian',PyRhOparameters())])
simList = list(simParams)


simParams['Python'].add_many(('dt', 0.1, 0, None)) #'ms'

# atol
simParams['NEURON'].add_many(('cell',   ['minimal.hoc'], None, None), #'morphology'
                             ('Vclamp', False,     False,  True), # Changed to False by default
                             ('Vcomp',  'soma',    None,   None),
                             ('expProb',1.0,       0.,     1.),
                             ('v_init', -65,       None,   None), # 'mV'
                             ('CVode',  False,     False,  True),
                             ('dt',     0.1,       0,   None)) # 'ms' #, 0.025

simParams['Brian'].add_many(('dt', 0.1, 0, None)) #, # 'ms'
                            

### Move somewhere else e.g. base.py
class PyRhOobject(object):
    """Common base class for all PyRhO objects"""
    
    # https://docs.python.org/3/reference/datamodel.html#special-method-names
    
    #def __new__(self):
    #    pass
        
    def __init__(self):
        pass
        
    def __del__(self):
        pass
    
    def __repr__(self):
        return str(self.__class__)
        
    def __str__(self):
        print("PyRhO object: ", self.__class__.__name__)
    
    def __call__(self):
        return
        
    def setParams(self, params):
        """Set all model parameters from a Parameters() object"""
        #for param, value in params.items():
        for p in params.keys():
            self.__dict__[p] = params[p].value #vars(self())[p]
    
    def updateParams(self, params):
        """Update model parameters which already exist"""
        pDict = params.valuesdict()
        count = 0
        for p, v in pDict.items():
            if p in self.__dict__: # Added to allow dummy variables in fitting parameters
                self.__dict__[p] = v #vars(self())[p]
                count += 1
        return count
    
    def getParams(self, params):
        """Export parameters to lmfit dictionary"""
        for p in self.__dict__.keys():
            params[p].value = self.__dict__[p]

    def exportParams(self, params):
        """Export parameters which are already in lmfit dictionary"""
        count = 0
        for p, v in self.__dict__.items():
            if p in params:
                params[p].value = v #self.__dict__[p]
                count += 1
        return count
            
    def printParams(self):
        for p in self.__dict__.keys():
            print(p, ' = ', self.__dict__[p])
    
    def printParamsWithLabels(self):
        for p in self.__dict__.keys():
            if p in unitLabels:
                print(p, ' = ', self.__dict__[p], ' [', unitLabels[p], ']')
            else:
                print(p, ' = ', self.__dict__[p])
    
    def printParamsWithUnits(self):
        for p in self.__dict__.keys():
            if p in modelUnits:
                print(p, ' = ', self.__dict__[p], ' * ', modelUnits[p])
            else:
                print(p, ' = ', self.__dict__[p])
    
    def getExt(self, var, ext='max'):
        if ext == 'max':
            mVal = max(self.__dict__[var])
        elif ext == 'min':
            mVal = min(self.__dict__[var])
        mInd = np.searchsorted(self.__dict__[var], mVal)
        return mVal, mInd
    
