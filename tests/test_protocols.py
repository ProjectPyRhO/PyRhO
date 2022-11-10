# import pytest
import numpy as np

# from pyrho import *
from pyrho import models, protocols, simulators


def test_recovery():
    RhO = models['6']()
    Prot = protocols['recovery']()
    #Prot.phis = [1e16, 1e15, 1e14]
    Sim = simulators['Python'](Prot, RhO)
    #Sim.Vclamp = True
    Sim.run()
    assert np.isclose(Prot.PD.params[0][0]['Gr0'].value, 0.00033)
    #Sim.plot()
