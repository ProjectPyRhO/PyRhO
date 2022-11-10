import numpy as np

from pyrho import Parameters, fitModels
from pyrho.datasets import loadChR2

def test_fit_3_state_model():
    init_params = Parameters()
    init_params.add_many(
        # Name   Value   Vary    Min     Max     Expr
        ('g0',   1e5,    True,   0.001,  1e6,    None),
        ('phi_m',1e18,   True,   1e15,   1e19,   None),
        ('k_a',  5,      True,   0.001,  1000,   None),
        ('k_r',  0.1,    True,   0.001,  1000,   None),
        ('p',    0.8,    True,   0.1,    5,      None),
        ('q',    0.25,   True,   0.1,    5,      None),
        ('Gd',   0.1,    True,   0.0001, 1,      None),
        ('Gr0',  0.0002, True,   0.0001, 0.1,    None),
        ('E',    0,      True,   -1000,  1000,   None),
        ('v0',   43,     True,   -1e15,  1e15,   None),
        ('v1',   17.1,   True,   -1e15,  1e15,   None))
    data = loadChR2()
    fit_params, mini_objs = fitModels(data, nStates=3, params=init_params, postFitOpt=True, relaxFact=2)
    values = fit_params[0].valuesdict()
    # print(values, flush=True)
    # OrderedDict([('g0', 28551.460430437444), ('phi_m', 7.45862659417406e+17), ('k_a', 6.622531560387775), ('k_r', 0.08504416795822778),
    #             ('p', 0.8377064000308175), ('q', 0.27409731460831976), ('Gd', 0.05918146172064458), ('Gr0', 0.0002), ('E', 0), ('v0', 43), ('v1', 17.1)])
    assert np.isclose(values["g0"], 28551.460430437444)
    assert np.isclose(values["phi_m"], 7.45862659417406e+17)
    assert np.isclose(values["k_a"], 6.622531560387775)
    assert np.isclose(values["k_r"], 0.08504416795822778)
    assert np.isclose(values["p"], 0.8377064000308175)
    assert np.isclose(values["q"], 0.27409731460831976)
    assert np.isclose(values["Gd"], 0.05918146172064458)
    assert np.isclose(values["Gr0"], 0.0002)
    assert np.isclose(values["E"], 0)
    assert np.isclose(values["v0"], 43)
    assert np.isclose(values["v1"], 17.1)
