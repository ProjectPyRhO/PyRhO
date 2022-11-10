# import pytest
import numpy as np

# from pyrho import *
from pyrho import models, protocols, simulators


def test_recovery():
    rho = models['6']()
    prot = protocols['recovery']()
    #Prot.phis = [1e16, 1e15, 1e14]
    sim = simulators['Python'](prot, rho)
    #Sim.Vclamp = True
    sim.run()
    assert np.isclose(prot.PD.params[0][0]['Gr0'].value, 0.00033)
    #Sim.plot()
