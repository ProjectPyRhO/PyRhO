from pyrho import *

RhO = models['6']()
Prot = protocols['step']()
Prot.phis = [1e16, 1e15, 1e14]
Sim = simulators['Python'](Prot, RhO)
Sim.run()
Sim.plot()
