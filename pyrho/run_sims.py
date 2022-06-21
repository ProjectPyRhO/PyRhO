# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 23:34:21 2022

@author: beers
"""

from pyrho import *
from pyrho.datasets import *

initParams = Parameters()
initParams.add_many(                
                # Name      Value   Vary  Min   Max   Expr
                ('g0',      2.52e4, True, 0.0,  1e15, None),
                ('gam',     0.0161, True, 0.0,  1,    None),  # Max=1 if gO1 >= gO2
                ('phi_m',   3.54e17,True, 1e15, 1e19, None),
                ('k1',      13.4,   True, 0.0,  1000, None),
                ('k2',      2.71,   True, 0.0,  1000, None),
                ('k3',      2.71,   True, 0.0,  1000, None), # New: find values
                ('p',       0.985,  True, 0.1,  5,    None),
                ('Gf0',     0.0389, True, 0.0,  1000, None),
                ('k_f',     0.103,  True, 0.0,  1000, None),
                #('Gb0',     0.0198, True, 0.0,  1000, None),
                ('k_b',     0.139,  True, 0.0,  1000, None),
                ('q',       1.58,   True, 0.1,  5,    None),
                ('Go1',     2,      True, 0.0,  1000, None),
                ('Go2',     0.0567, True, 0.0,  1000, None),
                ('Gd1',     0.112,  True, 0.0,  1000, None),
                ('Gd2',     0.0185, True, 0.0,  1000, None),
                ('Ga3',     250,    True, 0.0,  500, None),
                ('Gb',      40000,  True, 3400.0,  45000, None),
                ('E',       0,      True, -1000,1000, None),
                ('v0',      43,     True, -1e15, 1e15,None),
                ('v1',      17.1,   True, -1e15, 1e15,None))

saveData(initParams, 'initParams')

ChR2data = loadChR2()
fitParams = fitModels(ChR2data, nStates='6K', params=initParams, postFitOpt=True, relaxFact=2)