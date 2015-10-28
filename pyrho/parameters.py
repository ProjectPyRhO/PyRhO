
#verbose = 1


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
#from .utilities import irrad2flux, flux2irrad
from copy import deepcopy


### Hyperparameters

#fittingParams = Parameters()
#methods=('leastsq','nelder','lbfgsb','powell','cg','newton','cobyla','tnc','trust-ncg','dogleg','slsqp')

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

#modelFits = OrderedDict([('ChR2', OrderedDict([('3',Parameters()),('4',Parameters()),('6',Parameters())]))])
modelFits = OrderedDict([('3', OrderedDict([('ChR2',Parameters()), ('ArchT',Parameters())])), ('4', OrderedDict([('ChR2',Parameters())])), ('6', OrderedDict([('ChR2',Parameters())]))])
#ChR2 = OrderedDict([('3',Parameters()),('4',Parameters()),('6',Parameters())])

modelLabels = OrderedDict([('E','E'), ('g','g'), ('k','k'), ('p','p'), ('phim','\phi_m'), ('Gd','G_d'), ('Gr0','G_{r0}'), ('Gr1','G_{r1}'), ('v0','v_0'), ('v1','v_1'),
                            ('gam','\gamma'), ('k1','k_1'), ('k2','k_2'), ('Gf0','G_{f0}'), ('Gb0','G_{b0}'), ('kf','k_f'), ('kb','k_b'), ('q','q'), 
                            ('Gd1','G_{d1}'), ('Gd2','G_{d2}'), ('Go1','G_{o1}'), ('Go2','G_{o2}')])

# Units for model parameters (for Brian)
# Replace with http://pythonhosted.org/NeuroTools/parameters.html
from brian2.units.allunits import *
from brian2.units.stdunits import *
modelUnits = OrderedDict([('g',psiemens), ('gam',1), ('phim',mm**-2*second**-1), ('k',ms**-1), ('p',1), ('Gd',ms**-1), ('Gr0',ms**-1), ('Gr1',ms**-1), 
                            ('k1',ms**-1), ('k2',ms**-1), ('Gf0',ms**-1), ('Gb0',ms**-1), ('kf',ms**-1), ('kb',ms**-1), ('q',1), 
                            ('Gd1',ms**-1), ('Gd2',ms**-1), ('Go1',ms**-1), ('Go2',ms**-1), ('E',mV), ('v0',mV), ('v1',mV)])

#paramUnits
unitLabels = OrderedDict([('g','pS'), ('gam',''), ('phim','ph./mm^2/s'), ('k','ms^-1'), ('p',''), ('Gd','ms^-1'), ('Gr0','ms^-1'), ('Gr1','ms^-1'), 
                            ('k1','ms^-1'), ('k2','ms^-1'), ('Gf0','ms^-1'), ('Gb0','ms^-1'), ('kf','ms^-1'), ('kb','ms^-1'), ('q',''), 
                            ('Gd1','ms^-1'), ('Gd2','ms^-1'), ('Go1','ms^-1'), ('Go2','ms^-1'), ('E','mV'), ('v0','mV'), ('v1','mV')])

#'$\mathrm{[ph. / mm^{2} / s]}$'
#'ph./mm^2/s'
                            
####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80

#d3sp = Parameters()
#                       (Name,    Value,  Vary, Min,  Max,  Expr=Units)
# modelParams['3'].add_many(
                # #('phi0',  1e15,   True, None, None, None),
                # #('phiSat',1e20,   True, None, None, None),
                # ('E',     0,      True, -1000,1000, 'mV'),
                # ('g',     1.67e4, True, 0.0,  1e15, 'pS'),
                # ('k',     5.5e-15,True, 0.0,  1e15, None),
                # ('p',     0.7,    True, 0.1,  5,    None),
                # ('phim',  1e17,   True, 1e15, 1e19, 'photons/s/mm^2'),
                # ('Gd',    0.0909, True, 0.0,  None, '/ms'),
                # ('Gr0',   1/5000, True, 0.0,  None, '/ms'),
                # ('Gr1',   1/165,  True, 0.0,  None, '/ms'),
                # ('useIR', False,  True, False,True, None),
                # ('v0',    43,     True, None, None, 'mV'),
                # ('v1',    4.1,    True, None, None, 'mV'))
                
# modelParams['3'].add_many(
                # ('g',     1.67e4, True, 0.001,  1e15, 'pS'),
                # ('phim',  1e17,   True, 1e15, 1e19, 'photons/s/mm^2'),
                # ('k',     25,     True, 0.001,  1000, '/ms'), #('k',     5.5e-15,True, 0.0,  1e15, None),
                # ('p',     0.7,    True, 0.1,  5,    None),
                # ('Gd',    0.0909, True, 0.0001,  1, '/ms'),
                # ('Gr0',   1/5000, True, 0.0001,  0.1, '/ms'),
                # ('Gr1',   1/165,  True, 0.0001,  1, '/ms'),
                # ('E',     0,      True, -1000,1000, 'mV'),
                # ('v0',    43,     True, -1e15, 1e15, 'mV'),
                # ('v1',    4.1,    True, -1e15, 1e15, 'mV'))
                
#modelFits['ChR2']['3'].add_many(                
modelFits['3']['ChR2'].add_many(
                ('g',     1.57e5, True, 0.001,  1e6,  None),
                ('phim',  1.32e18,True, 1e15,   1e19, None),
                ('k',     4.51,   True, 0.001,  1000, None),
                ('p',     0.793,  True, 0.1,    5,    None),
                ('Gd',    0.104,  True, 0.0001, 1,    None),
                ('Gr0',   0.0002, True, 0.0001, 0.1,  None),
                ('Gr1',   0.0386, True, 0.0001, 1,    None),
                ('E',     0,      True, -1000,  1000, None),
                ('v0',    43,     True, -1e15,  1e15, None),
                ('v1',    17.1,    True, -1e15, 1e15, None))
#modelLabels = OrderedDict([('3',OrderedDict()),('4',OrderedDict()),('6',OrderedDict())])
#modelLabels['3'] = {'E':'E', 'g':'g', 'k':'k', 'p':'p', 'phim':'\phi_m', 'Gd':'G_d', 'Gr0':'G_{r0}', 'Gr1':'G_{r1}', 'v0':'v_0', 'v1':'v_1'}
#modelLabels['3'] = [('E','E'), ('g','g'), ('k','k'), ('p','p'), ('phim','\phi_m'), ('Gd','G_d'), ('Gr0','G_{r0}'), ('Gr1','G_{r1}'), ('v0','v_0'), ('v1','v_1')]
#modelParams['3'] = deepcopy(modelFits['3']['ChR2'])
                

### Alternatively add another field
#d3sp['g'].units = 'pS'
#print(d3sp['g'].units)
#d4sp = Parameters()
# modelParams['4'].add_many(
                # ('phi0',  1e14,   True, None, None, 'photons/s/mm^2'),
                # ('E',     0,      True, -1000,1000, 'mV'),
                # ('gam',   0.05,   True, 0.0,  1e9,  None),
                # ('g',     1.67e4, True, 0.0,  1e15, 'pS'),
                # ('k1',    0.05,   True, 0.0,  1e3,  None), ### Add bounds checking?
                # ('k2',    0.015,  True, 0.0,  1e3,  None),
                # ('c1',    0.03,   True, 0.0,  None, None),
                # ('c2',    0.0115, True, 0.0,  None, None),
                # ('p',     0.7,    True, None, None, None),
                # ('phim',  5e17,   True, None, None, 'photons/s/mm^2'),
                # ('e12d',  0.01,   True, 0.0,  None, '/ms'),
                # ('e21d',  0.015,  True, 0.0,  None, '/ms'),
                # ('Gd1',   0.11,   True, 0.0,  None, '/ms'),
                # ('Gd2',   0.025,  True, 0.0,  None, '/ms'),
                # ('Gr',    0.0004, True, 0.0,  None, '/ms'),
                # ('useIR', False,  True, False,True, None),
                # ('v0',    43,     True, None, None, 'mV'),
                # ('v1',    4.1,    True, None, None, 'mV'))
              
# modelParams['4'].add_many(  #('phi0',    1e14,   True, 1e12, 1e21,   'photons/s/mm^2'),### Set this to be above the max flux??? 10**ceil(log10(max(phis)))
                # ('E',       0,      True, -1000,1000,   'mV'),
                # ('gam',     0.05,   True, 0.0,  1,      None),
                # ('g',       3.5e4,  True, 0.001,1e15,   'pS'),
                # ('k1',      1000,   True, 0.001,1e5,    '/ms'), #3
                # ('k2',      500,    True, 0.001,1e5,    '/ms'), #1.5
                # ('p',       0.7,    True, 0.1,  5,      None),
                # ('e12d',    0.01,   True, 0,    1e3,    '/ms'),
                # ('e21d',    0.01,   True, 0,    1e3,    '/ms'),
                # ('c1',      0.4,   True, 0.001, 1e3,    '/ms'),
                # ('c2',      0.2,   True, 0.001, 1e3,    '/ms'),
                # ('q',       0.47,    True, 0.1,  5,      None),
                # ('phim',    1e16,   True, 1e15, 1e19,   'photons/s/mm^2'),
                # ('Gd1',     0.15,   True, 0.01, 1,      '/ms'),
                # ('Gd2',     0.025,  True, 0.01, 1,      '/ms'),
                # ('Gr',      0.0004, True, 1e-6, 1,      '/ms'),
                # #('useIR',   False,  True, False,True,   None),
                # ('v0',      43,     True, -1e15, 1e15,   'mV'),
                # ('v1',      4.1,    True, -1e15, 1e15,   'mV'))


# modelParams['4'].add_many(  #('phi0',    1e14,   True, 1e12, 1e21,   'photons/s/mm^2'),### Set this to be above the max flux??? 10**ceil(log10(max(phis)))
                # ('g',       3.5e4,  True, 0.001,1e15,   'pS'),
                # ('gam',     0.05,   True, 0.0,  1,      None),
                # ('phim',    1e16,   True, 1e15, 1e19,   'photons/s/mm^2'),
                # ('k1',      1000,   True, 0.001,1e5,    '/ms'), #3
                # ('k2',      500,    True, 0.001,1e5,    '/ms'), #1.5
                # ('p',       0.7,    True, 0.1,  5,      None),
                # ('Gf0',    0.01,   True, 0,    1e3,    '/ms'), #e12d
                # ('kf',      0.4,   True, 0.001, 1e3,    '/ms'), #c1
                # ('Gb0',    0.01,   True, 0,    1e3,    '/ms'), #e21d
                # ('kb',      0.2,   True, 0.001, 1e3,    '/ms'), #c2
                # ('q',       0.47,    True, 0.1,  5,      None),
                # ('Gd1',     0.15,   True, 0.01, 1,      '/ms'),
                # ('Gd2',     0.025,  True, 0.01, 1,      '/ms'),
                # ('Gr0',      0.00033, True, 1e-6, 1,      '/ms'), #Gr #0.0004
                # ('E',       0,      True, -1000,1000,   'mV'),
                # ('v0',      43,     True, -1e15, 1e15,   'mV'),
                # ('v1',      4.1,    True, -1e15, 1e15,   'mV'))    
    
#modelFits['ChR2']['4'].add_many(
modelFits['4']['ChR2'].add_many(
                ('g',       1.14e5, True, 0.001,1e15,   None),
                ('gam',     0.00742,True, 0.0,  1,      None),
                ('phim',    2.33e17,True, 1e15, 1e19,   None),
                ('k1',      4.15,   True, 0.001,1e5,    None), #3
                ('k2',      0.868,  True, 0.001,1e5,    None), #1.5
                ('p',       0.833,  True, 0.1,  5,      None),
                ('Gf0',     0.0373, True, 0,    1e3,    None), #e12d
                ('kf',      0.0581, True, 0.001, 1e3,   None), #c1
                ('Gb0',     0.0161, True, 0,    1e3,    None), #e21d
                ('kb',      0.063,  True, 0.001, 1e3,   None), #c2
                ('q',       1.94,   True, 0.1,  5,      None),
                ('Gd1',     0.105,  True, 0.01, 1,      None),
                ('Gd2',     0.0138, True, 0.01, 1,      None),
                ('Gr0',     0.00033, True, 1e-6, 1,      None), #Gr #0.0004
                ('E',       0,      True, -1000,1000,   None),
                ('v0',      43,     True, -1e15, 1e15,  None),
                ('v1',      17.1,    True, -1e15, 1e15, None))    
              
#modelParams['4'] = deepcopy(modelFits['4']['ChR2'])
#d6sp = Parameters()
# modelParams['6'].add_many(
                # #('phi0',  1e14,   True, None, None, 'photons/s/mm^2'),
                # ('phim',    1e16,   True, 1e15, 1e19,   'photons/s/mm^2'),
                # ('E',     0,      True, -1000,1000, 'mV'),
                # ('gam',   0.05,   True, 0.0,  1,  None), #1e9
                # ('g',     75000,  True, 0.0,  1e15, 'pS'),
                # #('A',     31192,  True, 0.0,  1e15, 'um^2'),
                # #('gbar',  2.4,    True, 0.0,  1e15, 'pS/um^2'),
                # ('a10',   5e3,    True, 0.0,  None, '/ms'), #5
                # ('a2',    1,      True, 0.0,  None, '/ms'),
                # ('a30',   0.022,  True, 0.0,  None, '/ms'),
                # ('a31',   0.0135, True, 0.0,  None, '/ms'),
                # ('a4',    0.025,  True, 0.0,  None, '/ms'),
                # ('a6',    0.00033,True, 0.0,  None, '/ms'),
                # ('b1',    0.13,   True, 0.0,  None, '/ms'),
                # ('b20',   0.011,  True, 0.0,  None, '/ms'),
                # ('b21',   0.0048, True, 0.0,  None, '/ms'),
                # ('b3',    1,      True, 0.0,  None, '/ms'),
                # ('b40',   1.1e3,  True, 0.0,  None, '/ms'), #1.1
                # ('p',     0.7,    True, 0.1,  5,    None),
                # ('q',     0.47,   True, 0.1,  5,    None),
                # #('useIR', True,   True, False,True, None),
                # ('v0',    43,     True, -1e15, 1e15, 'mV'),
                # ('v1',    4.1,    True, -1e15, 1e15, 'mV'))

# modelParams['6'].add_many(                
                # ('g',     75000,  True, 0.0,  1e15, 'pS'),
                # ('gam',   0.05,   True, 0.0,  1,  None), # Max=1 if gO1 >= gO2
                # ('phim',    1e16,   True, 1e15, 1e19,   'photons/s/mm^2'),
                # ('k1',   5e3,    True, 0.0,  None, '/ms'), #a10 #5
                # ('k2',   1.1e3,  True, 0.0,  None, '/ms'), #b40 #1.1
                # ('p',     0.7,    True, 0.1,  5,    None),
                # ('Gf0',   0.022,  True, 0.0,  None, '/ms'), #a30
                # ('kf',   0.0135, True, 0.0,  None, '/ms'), #a31
                # ('Gb0',   0.011,  True, 0.0,  None, '/ms'), #b20
                # ('kb',   0.0048, True, 0.0,  None, '/ms'), #b21
                # ('q',     0.47,   True, 0.1,  5,    None),
                # ('Go1',    1,      True, 0.0,  None, '/ms'), #a2
                # ('Go2',    1,      True, 0.0,  None, '/ms'), #b3
                # ('Gd1',    0.13,   True, 0.0,  None, '/ms'), #b1
                # ('Gd2',    0.025,  True, 0.0,  None, '/ms'), #a4
                # ('Gr0',    0.00033,True, 0.0,  None, '/ms'), #Gr
                # ('E',     0,      True, -1000,1000, 'mV'),
                # ('v0',    43,     True, -1e15, 1e15, 'mV'),
                # ('v1',    17.1,    True, -1e15, 1e15, 'mV'))


#modelFits['ChR2']['6'].add_many(
modelFits['6']['ChR2'].add_many(                
                ('g',       2.52e4,  True, 0.0,  1e15, None),
                ('gam',     0.0161, True, 0.0,  1,    None), # Max=1 if gO1 >= gO2
                ('phim',    3.54e17,True, 1e15, 1e19, None),
                ('k1',      13.4,   True, 0.0,  1000, None),
                ('k2',      2.71,   True, 0.0,  1000, None),
                ('p',       0.985,  True, 0.1,  5,    None),
                ('Gf0',     0.0389, True, 0.0,  1000, None),
                ('kf',      0.103,  True, 0.0,  1000, None),
                ('Gb0',     0.0198, True, 0.0,  1000, None),
                ('kb',      0.139,  True, 0.0,  1000, None),
                ('q',       1.58,   True, 0.1,  5,    None),
                ('Go1',     2,      True, 0.0,  1000, None),
                ('Go2',     0.0567, True, 0.0,  1000, None),
                ('Gd1',     0.112,  True, 0.0,  1000, None),
                ('Gd2',     0.0185, True, 0.0,  1000, None),
                ('Gr0',     0.00033,True, 0.0,  1000, None), #0.00163
                ('E',       0,      True, -1000,1000, None),
                ('v0',      43,     True, -1e15, 1e15,None),
                ('v1',      17.1,    True, -1e15, 1e15,None))

#modelParams['6'] = deepcopy(modelFits['6']['ChR2'])
#modelParams = [d3sp,d4sp,d6sp]
#modelParamsDict = {'3':d3sp, '4':d4sp, '6':d6sp} # Rename these to be consistent with protParams

#defModelParams = modelFits['ChR2']
modelParams['3'] = modelFits['3']['ChR2']
modelParams['4'] = modelFits['4']['ChR2']
modelParams['6'] = modelFits['6']['ChR2']

unitPrefixes = {} ### Use a units library to convert between different prefixes

####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80

#simUnitLabels = OrderedDict([('dt','ms'), ('v_init','mV')])
#simUnitLabels = defaultdict(str)
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

simParams = OrderedDict([('Python',Parameters()), ('NEURON',Parameters()), ('Brian',Parameters())])
simList = list(simParams)

simParams['Python'].add_many(('dt',0.1,True,None,None,None)) #'ms'

simParams['NEURON'].add_many(('cell',   ['minimal.hoc'],True,None,None, None), #'morphology'
                             ('Vclamp', False,   True,   False,  True,   None), # Changed to False by default
                             ('Vcomp',  'soma', True,   None,   None,   None),
                             ('expProb',1.0,    True,   0.,     1.,     None),
                             ('v_init', -65,    True,   None,   None,   None), # 'mV'
                             ('CVode',  False,  True,   False,  True,   None),
                             ('dt',     0.1,    True,   None,   None,   None)) # 'ms' #, 0.025
                             #('nseg',3,True,1,1e9,None),
                             #('Vhold',-70,True,-200,200,'mV')) # Set by runTrial
# atol

simParams['Brian'].add_many(('dt', 0.1, True, None, None, None), # 'ms'
                            )#('method_choice', )
                             
                             
### Select simulation protocol
#protocols = ['custom', 'saturate', 'rectifier', 'shortPulse', 'recovery']
#protocol = protocols[2] #'recovery'#'shortPulse' # Set this interactively with radio buttons?
#Prot = selectProtocol(protocol)
#Prot = protocols['custom']([1e13,1e14,1e15], [-70,-40,-10,10,40], [[10.,160.]], 200., 1, 0.1)
#Prot = protocols['step']([1e14], [-70], [[25.,275.]], 300., 1, 0.1) # Unfinished
##Prot = protocols['sinusoid']([1e9,1e10,1e11,1e12,1e13,1e14,1e15], [1e9], [-70], np.logspace(0,2,num=15), [[25.,5000.]], 5050., 0.1)
#Prot = protocols['sinusoid']([1e11,1e12,1e13,1e14,1e15], [1e13], [-70], np.logspace(0,2,num=50), [[25.,5000.]], 5050., 0.1)
#Prot = protocols['ramp']() # Unfinished
#Prot = protocols['saturate']([irrad2flux(1000,470)], [-70], [[5.,5+1e-3]], 20., 1, 1e-3)
#Prot = protocols['rectifier']([irrad2flux(1,470),irrad2flux(10,470)], [-100,-80,-60,-40,-20,0,20,40,60,80], [[50.,300.]], 400., 1, 0.1)
#Prot = protocols['shortPulse']([1e12], [-70], 25, [1,2,3,5,8,10,20], 100, 0.1)
#Prot = protocols['recovery']([1e14], [-70], 100, 200, [500,1000,1500,2500,5000,7500,10000], 0.1)


### Protocols to be included in the next version:
### - Temperature (Q10)
### - pH (intracellular and extracellular)
### - Wavelength (lambda)

protParams = OrderedDict([('custom',Parameters()), ('step',Parameters()), ('sinusoid',Parameters()), ('chirp',Parameters()), ('ramp',Parameters()), ('saturate',Parameters()), ('rectifier',Parameters()), ('shortPulse',Parameters()), ('recovery',Parameters())])

protList = list(protParams) # List of keys #This could be removed

# Change these to include $ $
protParamLabels = OrderedDict([ ('phis', '\mathbf{\phi}'), 
                                ('Vs', '\mathbf{\mathrm{V}}'), 
                                ('delD', '\Delta t_{delay}'), 
                                ('onD', '\Delta t_{on}'), 
                                ('totT', 'T_{total}'), 
                                ('cycles', 'cycles'),
                                ('A0', 'A_0'), 
                                ('fs', '\mathbf{f}'), 
                                ('f0', 'f_0'), 
                                ('fT', 'f_T'), 
                                ('linear', 'linear'), 
                                ('startOn', '\phi_{t=0}>0'), 
                                ('phi_ton', '\phi_{t=0}'), 
                                ('pDs', '\mathbf{\Delta t_{on}}'), 
                                ('IPIs', '\mathbf{\Delta t_{off}}') ])

protUnitLabels = defaultdict(lambda: '')
protUnitLabels['phis'] = 'ph./mm^2/s'
protUnitLabels['A0'] = 'ph./mm^2/s'
protUnitLabels['phi_ton'] = 'ph./mm^2/s' ### Revise!!!
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
    protParamNotes[prot]['delD'] = 'Delay duration before the first pulse cycle'
    protParamNotes[prot]['cycles'] = 'List of [on, off] durations for each pulse cycle'
    

#Exceptions
protParamNotes['sinusoid']['startOn'] = 'Start at maximum of flux modulation (else minimum)'
protParamNotes['sinusoid']['A0'] = 'Constant offset for modulation'
protParamNotes['sinusoid']['fs'] = 'List of modulation frequencies'
protParamNotes['chirp']['linear'] = 'Use a linear frequency sweep (else exponential)'
protParamNotes['chirp']['startOn'] = 'Start at maximum of flux modulation (else minimum)'
protParamNotes['chirp']['A0'] = 'Constant offset for modulation'
protParamNotes['chirp']['f0'] = 'Starting frequency'
protParamNotes['chirp']['fT'] = 'Ending frequency'
protParamNotes['ramp']['phis'] = 'List of ending flux values'
protParamNotes['ramp']['phi_ton'] = 'Starting flux value'
protParamNotes['saturate']['cycles'] = ''
protParamNotes['saturate']['onD'] = 'On-phase duration'
protParamNotes['saturate']['totT'] = 'Total simulation duration'
protParamNotes['shortPulse']['cycles'] = ''
protParamNotes['shortPulse']['pDs'] = 'List of cycle on-phase durations'
protParamNotes['shortPulse']['totT'] = 'Total simulation duration'
protParamNotes['recovery']['cycles'] = ''
protParamNotes['recovery']['onD'] = 'Cycle on-phase duration'
protParamNotes['recovery']['IPIs'] = 'List of cycle off-phase durations'
protParamNotes['recovery']['totT'] = 'Total simulation duration'
    
#squarePulses = ['custom', 'saturate', 'step', 'rectifier', 'shortPulse', 'recovery'] #{'custom': True, 'saturate': True, 'step': True, 'rectifier': True, 'shortPulse': True, 'recovery': True}
#arbitraryPulses = ['custom', 'sinusoid', 'chirp', 'ramp'] #{'custom': True, 'sinusoid': True, 'chirp': True, 'ramp':True} # Move custom here
smallSignalAnalysis = ['saturate', 'step', 'sinusoid'] #{'sinusoid': True, 'step': True, 'saturate': True} 

#ProtParamsCustom = Parameters()
#phis=[1e14,1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[10.,160.]], totT=200., nRuns=1, dt=0.1
# protParams['custom'].add_many(('phis',[1e14,1e15,1e16,1e17],True,None,None,'photons/s/mm^2'),
                          # ('Vs',[-70,-40,-10,10,40],True,None,None,'mV'),
                          # ('pulses',[[10.,160.]],True,None,None,'ms'),
                          # ('totT', 200.,True,0,None,'ms'),#)#,
                          # #('nRuns', 1,True,1,None,None),
                          # ('dt',0.1,True,1e-9,10,'ms'))
                          
protParams['custom'].add_many(('phis',[1e15,1e16],True,None,None,None), #'photons/s/mm^2'
                            ('Vs',[-70,-20,10],True,None,None,None), #'mV'
                            ('delD', 25, True, 0, 1e9, None), #'ms'
                            ('cycles',[[150.,50.]],True,None,None,None)) #'ms'#,
                            #('totT', 200.,True,0,None,'ms'),#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))


#ProtParamsStep = Parameters()
#phis=[1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[50.,200.]], totT=300., nRuns=1, dt=0.1
# protParams['step'].add_many(('phis',[1e15,1e16,1e17],True,None,None,'photons/s/mm^2'),
                        # ('Vs',[-70,-40,-10,10,40],True,None,None,'mV'),
                        # ('pulses',[[50.,200.]],True,None,None,'ms'),
                        # ('totT', 300.,True,0,None,'ms'),#)#,
                        # #('nRuns', 1,True,1,None,None),
                        # ('dt',0.1,True,1e-9,10,'ms'))

protParams['step'].add_many(('phis',[1e15,1e16,1e17],True,None,None,None), #'photons/s/mm^2'
                            ('Vs',[-70,-40,-10,10,40],True,None,None,None), #'mV'
                            ('delD', 25, True, 0, 1e9, None), #'ms'
                            ('cycles',[[150.,100.]],True,None,None,None)) #'ms'#,
                            #('totT', 300.,True,0,None,'ms'),#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsSinusoid = Parameters()
#phis=[1e14], A0=[1e12], Vs=[-70], fs=np.logspace(-1,3,num=9), pulses=[[25.,275.]], totT=300., dt=0.1
# protParams['sinusoid'].add_many(('phis',[1e9],True,None,None,'photons/s/mm^2'),
                            # ('startOn',True,False,False,True,''),
                            # ('A0',[0],True,None,None,'photons/s/mm^2'),
                            # ('Vs',[-70],True,None,None,'mV'),
                            # ('fs',[0.1,0.5,1,5,10,50,100,500,1000],True,None,None,'Hz'),
                            # ('pulses',[[50.,1050.]],True,None,None,'ms'),
                            # ('totT', 1100.,True,0,None,'ms'),#)#,
                            # #('nRuns', 1,True,1,None,None),
                            # ('dt',0.1,True,1e-9,10,'ms'))
                            
protParams['sinusoid'].add_many(('phis',[1e9],True,None,None,None), #'photons/s/mm^2'
                            ('startOn',True,False,False,True,None),
                            ('A0',[0],True,None,None,None), #'photons/s/mm^2'
                            ('Vs',[-70],True,None,None,None), #'mV'
                            ('fs',[0.1,0.5,1,5,10,50,100,500,1000],True,None,None,None), #'Hz'
                            ('delD', 25, True, 0, 1e9, None), #'ms'
                            ('cycles',[[1000.,50.]],True,None,None,None)) #'ms' #,
                            #('totT', 1100.,True,0,None,'ms'),#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsChirp = Parameters()
# protParams['chirp'].add_many(('phis',[1e14],True,None,None,'photons/s/mm^2'),
                            # ('linear',True,True,False,True,''), # False := exponential
                            # ('startOn',False,True,False,True,''),
                            # ('A0',[0],True,None,None,'photons/s/mm^2'),
                            # ('Vs',[-70],True,None,None,'mV'),
                            # ('pulses',[[100.,1100.]],True,None,None,'ms'),
                            # ('totT', 1200.,True,0,None,'ms'),
                            # ('f0',0.1,True,None,None,'Hz'),
                            # ('fT',1000,True,None,None,'Hz'),#)#,
                            # ('dt',0.1,True,1e-9,10,'ms'))

protParams['chirp'].add_many(('phis',[1e14],True,None,None,None), # 'photons/s/mm^2'
                            ('linear',True,True,False,True,None), # False := exponential
                            ('startOn',False,True,False,True,None),
                            ('A0',[0],True,None,None,None), # 'photons/s/mm^2'
                            ('Vs',[-70],True,None,None,None), # 'mV'
                            ('delD', 100, True, 0, 1e9, None), # 'ms'
                            ('cycles',[[1000.,100.]],True,None,None,None), # 'ms'
                            #('totT', 1200.,True,0,None,'ms'),
                            ('f0',0.1,True,None,None,None), # 'Hz'
                            ('fT',1000,True,None,None,None)) # 'Hz'#,
                            #('dt',0.1,True,1e-9,10,'ms'))
                            
#ProtParamsRamp = Parameters()
# protParams['ramp'].add_many(('phis',[1e12,1e13,1e14,1e15,1e16,1e17,1e18],True,None,None,'photons/s/mm^2'),
                            # ('phi_ton',0,True,None,None,'photons/s/mm^2'),
                            # ('Vs',[-70],True,None,None,'mV'),
                            # ('pulses',[[25.,275.]],True,None,None,'ms'),
                            # ('totT', 300.,True,0,None,'ms'),#)#,
                            # #('nRuns', 1,True,1,None,None),
                            # ('dt',0.1,True,1e-9,10,'ms'))

protParams['ramp'].add_many(('phis',[1e12,1e13,1e14,1e15,1e16,1e17,1e18],True,None,None,None), # 'photons/s/mm^2'
                            ('phi_ton',0,True,None,None,None), # 'photons/s/mm^2'
                            ('Vs',[-70],True,None,None,None), # 'mV'
                            ('delD', 25, True, 0, 1e9, None), # 'ms'
                            ('cycles',[[250.,25.]],True,None,None,None)) # 'ms'#,
                            #('totT', 300.,True,0,None,'ms'),#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))
                            
#ProtParamsSaturate = Parameters()
#phis=[irrad2flux(1000,470)], Vs=[-70], pulses=[[5.,5+1e-3]], totT=20., nRuns=1, dt=1e-3
# protParams['saturate'].add_many(('phis',[1e20],True,None,None,'photons/s/mm^2'),
                            # ('Vs',[-70],True,None,None,'mV'),
                            # ('pulses',[[5.,5.+1e-3]],True,None,None,'ms'),
                            # ('totT', 20.,True,0,None,'ms'),#)#,
                            # #('nRuns', 1,True,1,None,None),
                            # ('dt',1e-3,True,1e-9,10,'ms'))
                            
protParams['saturate'].add_many(('phis',[1e20],True,None,None,None), # 'photons/s/mm^2'
                            ('Vs',[-70],True,None,None,None), # 'mV'
                            ('delD', 5, True, 0, 1e9, None), # 'ms'
                            ('onD', 1e-3, True, 0, 1e9, None), # 'ms'
                            #('pulses',[[5.,5.+1e-3]],True,None,None,'ms'),
                            ('totT', 25.,True,0,None,None)) # 'ms' #,#)#,
                            #('nRuns', 1,True,1,None,None),
                            #('dt',1e-3,True,1e-9,10,'ms'))

#ProtParamsInwardRect = Parameters()
#phis=[irrad2flux(1,470),irrad2flux(10,470)], Vs=[-100,-80,-60,-40,-20,0,20,40,60,80], pulses=[[50.,300.]], totT=400., nRuns=1, dt=0.1
# protParams['rectifier'].add_many(('phis',[1e16],True,None,None,'photons/s/mm^2'),
                            # ('Vs',[-100,-80,-60,-40,-20,0,20,40,60,80],True,None,None,'mV'),
                            # ('pulses',[[50.,300.]],True,None,None,'ms'),
                            # ('totT', 400.,True,0,None,'ms'),#)#,
                            # #('nRuns', 1,True,1,None,None),
                            # ('dt',0.1,True,1e-9,10,'ms'))

protParams['rectifier'].add_many(('phis',[1e16],True,None,None,None), # 'photons/s/mm^2'
                            ('Vs',[-100,-70,-40,-10,20,50,80],True,None,None,None), # 'mV' #[-100,-80,-60,-40,-20,0,20,40,60,80]
                            ('delD', 50, True, 0, 1e9, None), # 'ms'
                            ('cycles',[[250.,100.]],True,None,None,None)) # 'ms' #,
                            #('totT', 400.,True,0,None,'ms'),#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))
                            
#ProtParamsVaryPL = Parameters()
#phis=[1e12], Vs=[-70], delD=25, pDs=[1,2,3,5,8,10,20], totT=100, dt=0.1
protParams['shortPulse'].add_many(('phis',[1.5e15],True,None,None,None), # 'photons/s/mm^2' #1e12 # original
                            ('Vs',[-70],True,None,None,None), # 'mV'
                            ('delD',25,True,0,None,None), # 'ms'
                            ('pDs',[1,2,3,5,8,10,20],True,None,None,None), # 'ms' # [0.1, 0.2, 0.5, 1, 2, 5, 10]
                            ('totT', 100.,True,0,None,None)) # 'ms' #,#)#,
                            #('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsVaryIPI = Parameters()
#phis=[1e14], Vs=[-70], delD=100, onD=200, IPIs=[500,1000,1500,2500,5000,7500,10000], dt=0.1
protParams['recovery'].add_many(('phis',[1e14],True,None,None,None), # 'photons/s/mm^2'
                            ('Vs',[-70],True,None,None,None), # 'mV'
                            ('delD',100,True,0,None,None), # 'ms'
                            ('onD',100,True,0,None,None), # 'ms'
                            ('IPIs',[500,1000,1500,2500,5000,7500,10000],True,None,None,None), # 'ms' #)#,
                            ('totT', 12000, True, 0, None, None)) # 'ms' #,
                            #('dt',0.1,True,1e-9,10,'ms'))

# Automatically link fields?
#protocolParams = OrderedDict([('custom',ProtParamsCustom), ('step',ProtParamsStep), ('sinusoid',ProtParamsSinusoid), ('chirp',ProtParamsChirp), ('ramp',ProtParamsRamp), ('saturate',ProtParamsSaturate), ('rectifier',ProtParamsInwardRect), ('shortPulse',ProtParamsVaryPL), ('recovery',ProtParamsVaryIPI)])
#protocolParams = {'custom':ProtParamsCustom, 'step':ProtParamsStep, 'sinusoid':ProtParamsSinusoid, 'chirp':ProtParamsChirp, 'ramp':ProtParamsRamp, 'saturate':ProtParamsSaturate, 'rectifier':ProtParamsInwardRect, 'shortPulse':ProtParamsVaryPL, 'recovery':ProtParamsVaryIPI}

# Lists and Dicts
#protDict={'custom':'custom', 'step':'step', 'sinusoid':'sinusoid', 'chirp':'chirp', 'ramp':'ramp', 'saturate':'saturate', 'rectifier':'rectifier', 'shortPulse':'shortPulse', 'recovery':'recovery'}
#OrderedDict([('custom':'custom'), ('step':'step'), ('sinusoid':'sinusoid'), ('ramp':'ramp'), ('saturate':'saturate'), ('rectifier':'rectifier'), ('shortPulse':'shortPulse'), ('recovery':'recovery')])
#protIndDict={'custom':0, 'step':1, 'sinusoid':2, 'chirp':3, 'ramp':4, 'saturate':5, 'rectifier':6, 'shortPulse':7, 'recovery':8} # This must match tabs ordering
#{key:i for key in protInDict, i++}

#protList = ['custom', 'step', 'sinusoid', 'chirp', 'ramp', 'saturate', 'rectifier', 'shortPulse', 'recovery']



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
    
    def getExt(self, var, ext='max'):
        if ext == 'max':
            mVal = max(self.__dict__[var])
        elif ext == 'min':
            mVal = min(self.__dict__[var])
        mInd = np.searchsorted(self.__dict__[var], mVal)
        return mVal, mInd
    
