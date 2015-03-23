
#verbose = 1


### Fitting Hyper parameters

# Time window for steady-state averaging
# IPI curve initial parameters
# fV initial parameters
# Kinetics initial parameters
# Optimisation routine initial parameters
# ...




### Hyperparameters

# Optimisation initialisation values
p0fV = (40,4,1)#25000)#,1)      # v0,v1,E,G
p0IPI = (0.5,4000,-1) # a*exp(-t/b)+c #(-1e-8,400,-1e-7)

# Used if plotKinetics
p0on = (-0.1,2,-1) # a*exp(-t/b)+c
p0inact = (-0.5,25,-0.5) # a*exp(-t/b)+c
p0off = (-0.1,7.5,-0.1,35,-0.1) # a1*exp(-t/tau1)+a2*exp(-t/tau2)+I_ss

tFromOff = 50  # Time [ms] to start the sample window before the end of the pulse for Iss



### Default model parameters

### Move this to models.py and add .setParams() method
from lmfit import Parameters, Parameter
from collections import OrderedDict #, defaultdict

modelParams = OrderedDict([('3',Parameters()),('4',Parameters()),('6',Parameters())])
modelList = list(modelParams) # List of keys: list(modelParams.keys()) #This could be removed

####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80

#d3sp = Parameters()
#                       (Name,    Value,  Vary, Min,  Max,  Expr=Units)
modelParams['3'].add_many(
                #('phi0',  1e15,   True, None, None, None),
                #('phiSat',1e20,   True, None, None, None),
                ('E',     0,      True, -1000,1000, 'mV'),
                ('g',     1.67e4, True, 0.0,  1e15, 'pS'),
                ('k',     5.5e-15,True, 0.0,  1e15, None),
                ('Gd',    0.0909, True, 0.0,  None, '/ms'),
                ('Gr0',   1/5000, True, 0.0,  None, '/ms'),
                ('Gr1',   1/165,  True, 0.0,  None, '/ms'))
                

### Alternatively add another field
#d3sp['g'].units = 'pS'
#print(d3sp['g'].units)
#d4sp = Parameters()
modelParams['4'].add_many(
                ('phi0',  1e14,   True, None, None, 'photons/s/mm^2'),
                ('E',     0,      True, -1000,1000, 'mV'),
                ('gam',   0.05,   True, 0.0,  1e9,  None),
                ('g',     1.67e4, True, 0.0,  1e15, 'pS'),
                ('k1',    0.05,   True, 0.0,  1e3,  None), ### Add bounds checking?
                ('k2',    0.015,  True, 0.0,  1e3,  None),
                ('c1',    0.03,   True, 0.0,  None, None),
                ('c2',    0.0115, True, 0.0,  None, None),
                ('e12d',  0.01,   True, 0.0,  None, '/ms'),
                ('e21d',  0.015,  True, 0.0,  None, '/ms'),
                ('Gd1',   0.11,   True, 0.0,  None, '/ms'),
                ('Gd2',   0.025,  True, 0.0,  None, '/ms'),
                ('Gr',    0.0004, True, 0.0,  None, '/ms'))
              

#d6sp = Parameters()
modelParams['6'].add_many(
                ('phi0',  1e14,   True, None, None, 'photons/s/mm^2'),
                ('E',     0,      True, -1000,1000, 'mV'),
                ('gam',   0.05,   True, 0.0,  1e9,  None),
                ('g',     75000,  True, 0.0,  1e15, 'pS'),
                #('A',     31192,  True, 0.0,  1e15, 'um^2'),
                #('gbar',  2.4,    True, 0.0,  1e15, 'pS/um^2'),
                ('v0',    43,     True, None, None, 'mV'),
                ('v1',    4.1,    True, None, None, 'mV'),
                ('a10',   5,      True, 0.0,  None, '/ms'),
                ('a2',    1,      True, 0.0,  None, '/ms'),
                ('a30',   0.022,  True, 0.0,  None, '/ms'),
                ('a31',   0.0135, True, 0.0,  None, '/ms'),
                ('a4',    0.025,  True, 0.0,  None, '/ms'),
                ('a6',    0.00033,True, 0.0,  None, '/ms'),
                ('b1',    0.13,   True, 0.0,  None, '/ms'),
                ('b20',   0.011,  True, 0.0,  None, '/ms'),
                ('b21',   0.0048, True, 0.0,  None, '/ms'),
                ('b3',    1,      True, 0.0,  None, '/ms'),
                ('b40',   1.1,    True, 0.0,  None, '/ms'))



#modelParams = [d3sp,d4sp,d6sp]
#modelParamsDict = {'3':d3sp, '4':d4sp, '6':d6sp} # Rename these to be consistent with protParams


unitPrefixes = {} ### Use a units library to convert between different prefixes

####|###10####|###20####|###30####|###40####|###50####|###60####|###70####|###80



simParams = OrderedDict([('Python',Parameters()), ('NEURON',Parameters()), ('Brian',Parameters())])
simList = list(simParams)

simParams['Python'].add_many(('dt',0.1,True,None,None,'ms'))

simParams['NEURON'].add_many(('v_init',-65,True,None,None,'mV'),
                             ('CVODE',False,True,False,True,None),
                             ('dt',0.1,True,None,None,'ms'),
                             ('nseg',3,True,1,1e9,None),
                             ('expProb',1.0,True,0.,1.,None),
                             ('cell',['minimal.hoc'],True,None,None,None), #'morphology'
                             ('Vclamp',True,True,False,True,None))#,
                             #('Vhold',-70,True,-200,200,'mV')) # Set by runTrial
# atol

#simParams['Brian'].add_many()
                             
                             
### Select simulation protocol
#protocols = ['custom', 'saturate', 'inwardRect', 'varyPL', 'varyIPI']
#protocol = protocols[2] #'varyIPI'#'varyPL' # Set this interactively with radio buttons?
#Prot = selectProtocol(protocol)
#Prot = protocols['custom']([1e13,1e14,1e15], [-70,-40,-10,10,40], [[10.,160.]], 200., 1, 0.1)
#Prot = protocols['step']([1e14], [-70], [[25.,275.]], 300., 1, 0.1) # Unfinished
##Prot = protocols['sinusoid']([1e9,1e10,1e11,1e12,1e13,1e14,1e15], [1e9], [-70], np.logspace(0,2,num=15), [[25.,5000.]], 5050., 0.1)
#Prot = protocols['sinusoid']([1e11,1e12,1e13,1e14,1e15], [1e13], [-70], np.logspace(0,2,num=50), [[25.,5000.]], 5050., 0.1)
#Prot = protocols['ramp']() # Unfinished
#Prot = protocols['saturate']([irrad2flux(1000,470)], [-70], [[5.,5+1e-3]], 20., 1, 1e-3)
#Prot = protocols['inwardRect']([irrad2flux(1,470),irrad2flux(10,470)], [-100,-80,-60,-40,-20,0,20,40,60,80], [[50.,300.]], 400., 1, 0.1)
#Prot = protocols['varyPL']([1e12], [-70], 25, [1,2,3,5,8,10,20], 100, 0.1)
#Prot = protocols['varyIPI']([1e14], [-70], 100, 200, [500,1000,1500,2500,5000,7500,10000], 0.1)


### Protocols to be included in the next version:
### - Temperature (Q10)
### - pH (intracellular and extracellular)
### - Wavelength (lambda)

protParams = OrderedDict([('custom',Parameters()), ('step',Parameters()), ('sinusoid',Parameters()), ('chirp',Parameters()), ('ramp',Parameters()), ('saturate',Parameters()), ('inwardRect',Parameters()), ('varyPL',Parameters()), ('varyIPI',Parameters())])

protList = list(protParams) # List of keys #This could be removed

#squarePulses = ['custom', 'saturate', 'step', 'inwardRect', 'varyPL', 'varyIPI'] #{'custom': True, 'saturate': True, 'step': True, 'inwardRect': True, 'varyPL': True, 'varyIPI': True}
#arbitraryPulses = ['custom', 'sinusoid', 'chirp', 'ramp'] #{'custom': True, 'sinusoid': True, 'chirp': True, 'ramp':True} # Move custom here
smallSignalAnalysis = ['sinusoid', 'step', 'saturate'] #{'sinusoid': True, 'step': True, 'saturate': True} 

#ProtParamsCustom = Parameters()
#phis=[1e14,1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[10.,160.]], totT=200., nRuns=1, dt=0.1
protParams['custom'].add_many(('phis',[1e14,1e15,1e16,1e17],True,None,None,'photons/s/mm^2'),
                          ('Vs',[-70,-40,-10,10,40],True,None,None,'mV'),
                          ('pulses',[[10.,160.]],True,None,None,'ms'),
                          ('totT', 200.,True,0,None,'ms'),
                          #('nRuns', 1,True,1,None,None),
                          ('dt',0.1,True,1e-9,10,'ms'))


#ProtParamsStep = Parameters()
#phis=[1e15,1e16,1e17], Vs=[-70,-40,-10,10,40], pulses=[[50.,200.]], totT=300., nRuns=1, dt=0.1
protParams['step'].add_many(('phis',[1e15,1e16,1e17],True,None,None,'photons/s/mm^2'),
                        ('Vs',[-70,-40,-10,10,40],True,None,None,'mV'),
                        ('pulses',[[50.,200.]],True,None,None,'ms'),
                        ('totT', 300.,True,0,None,'ms'),
                        #('nRuns', 1,True,1,None,None),
                        ('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsSinusoid = Parameters()
#phis=[1e14], A0=[1e12], Vs=[-70], fs=np.logspace(-1,3,num=9), pulses=[[25.,275.]], totT=300., dt=0.1
protParams['sinusoid'].add_many(('phis',[1e14],True,None,None,'photons/s/mm^2'),
                            ('A0',[0],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('fs',[0.1,0.5,1,5,10,50,100,500,1000],True,None,None,'Hz'),
                            ('pulses',[[50.,1050.]],True,None,None,'ms'),
                            ('totT', 1100.,True,0,None,'ms'),
                            #('nRuns', 1,True,1,None,None),
                            ('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsChirp = Parameters()
protParams['chirp'].add_many(('phis',[1e14],True,None,None,'photons/s/mm^2'),
                            ('A0',[0],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('pulses',[[100.,1100.]],True,None,None,'ms'),
                            ('totT', 1200.,True,0,None,'ms'),
                            ('f0',0.1,True,None,None,'Hz'),
                            ('fT',1000,True,None,None,'Hz'),
                            ('dt',0.1,True,1e-9,10,'ms'))
                            
#ProtParamsRamp = Parameters()
protParams['ramp'].add_many(('phis',[1e12,1e13,1e14,1e15,1e16,1e17,1e18],True,None,None,'photons/s/mm^2'),
                            ('phi_ton',0,True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('pulses',[[25.,275.]],True,None,None,'ms'),
                            ('totT', 300.,True,0,None,'ms'),
                            #('nRuns', 1,True,1,None,None),
                            ('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsSaturate = Parameters()
#phis=[irrad2flux(1000,470)], Vs=[-70], pulses=[[5.,5+1e-3]], totT=20., nRuns=1, dt=1e-3
protParams['saturate'].add_many(('phis',[1e20],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('pulses',[[5.,5.+1e-3]],True,None,None,'ms'),
                            ('totT', 20.,True,0,None,'ms'),
                            #('nRuns', 1,True,1,None,None),
                            ('dt',1e-3,True,1e-9,10,'ms'))

#ProtParamsInwardRect = Parameters()
#phis=[irrad2flux(1,470),irrad2flux(10,470)], Vs=[-100,-80,-60,-40,-20,0,20,40,60,80], pulses=[[50.,300.]], totT=400., nRuns=1, dt=0.1
protParams['inwardRect'].add_many(('phis',[2.366e15,2.366e16],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-100,-80,-60,-40,-20,0,20,40,60,80],True,None,None,'mV'),
                            ('pulses',[[50.,300.]],True,None,None,'ms'),
                            ('totT', 400.,True,0,None,'ms'),
                            #('nRuns', 1,True,1,None,None),
                            ('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsVaryPL = Parameters()
#phis=[1e12], Vs=[-70], delD=25, pDs=[1,2,3,5,8,10,20], totT=100, dt=0.1
protParams['varyPL'].add_many(('phis',[1e12],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('delD',25,True,0,None,'ms'),
                            ('pDs',[1,2,3,5,8,10,20],True,None,None,'ms'),
                            ('totT', 100.,True,0,None,'ms'),
                            ('dt',0.1,True,1e-9,10,'ms'))

#ProtParamsVaryIPI = Parameters()
#phis=[1e14], Vs=[-70], delD=100, onD=200, IPIs=[500,1000,1500,2500,5000,7500,10000], dt=0.1
protParams['varyIPI'].add_many(('phis',[1e14],True,None,None,'photons/s/mm^2'),
                            ('Vs',[-70],True,None,None,'mV'),
                            ('delD',100,True,0,None,'ms'),
                            ('onD',100,True,0,None,'ms'),
                            ('IPIs',[500,1000,1500,2500,5000,7500,10000],True,None,None,'ms'),
                            ('dt',0.1,True,1e-9,10,'ms'))

# Automatically link fields?
#protocolParams = OrderedDict([('custom',ProtParamsCustom), ('step',ProtParamsStep), ('sinusoid',ProtParamsSinusoid), ('chirp',ProtParamsChirp), ('ramp',ProtParamsRamp), ('saturate',ProtParamsSaturate), ('inwardRect',ProtParamsInwardRect), ('varyPL',ProtParamsVaryPL), ('varyIPI',ProtParamsVaryIPI)])
#protocolParams = {'custom':ProtParamsCustom, 'step':ProtParamsStep, 'sinusoid':ProtParamsSinusoid, 'chirp':ProtParamsChirp, 'ramp':ProtParamsRamp, 'saturate':ProtParamsSaturate, 'inwardRect':ProtParamsInwardRect, 'varyPL':ProtParamsVaryPL, 'varyIPI':ProtParamsVaryIPI}

# Lists and Dicts
#protDict={'custom':'custom', 'step':'step', 'sinusoid':'sinusoid', 'chirp':'chirp', 'ramp':'ramp', 'saturate':'saturate', 'inwardRect':'inwardRect', 'varyPL':'varyPL', 'varyIPI':'varyIPI'}
#OrderedDict([('custom':'custom'), ('step':'step'), ('sinusoid':'sinusoid'), ('ramp':'ramp'), ('saturate':'saturate'), ('inwardRect':'inwardRect'), ('varyPL':'varyPL'), ('varyIPI':'varyIPI')])
#protIndDict={'custom':0, 'step':1, 'sinusoid':2, 'chirp':3, 'ramp':4, 'saturate':5, 'inwardRect':6, 'varyPL':7, 'varyIPI':8} # This must match tabs ordering
#{key:i for key in protInDict, i++}

#protList = ['custom', 'step', 'sinusoid', 'chirp', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI']





