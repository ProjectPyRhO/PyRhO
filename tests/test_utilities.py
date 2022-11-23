# import pytest
import numpy as np
#from .. import pyrho
#from pyrho import utilities
from pyrho.utilities import (texIt, calcV1, getIndex, flux2irrad, irrad2flux,
                            times2cycles, cycles2times, round_sig)


def test_tex_it():
    assert texIt(r'\phi_\lambda') == r'$\phi_\lambda$'
    assert texIt(r'\phi_0') == r'$\phi_0$'


def test_calc_v1():
    assert np.isclose(calcV1(E=50, v0=-2.5), -120)
    # == (-70-50)/(1-np.exp(-(-70-50)/-2.5))


def test_get_index_list_None():
    assert getIndex([3.14, 7, None], None) == 2


def test_get_index_list_float():
    assert getIndex([5, 3.14, 7], 3.14) == 1


def test_get_index_array_float():
    assert getIndex(np.array([5, 3.14, 7]), 3.14) == 1


def test_flux2irrad():
    assert np.isclose(flux2irrad(1e16, 540), 3.678603)


def test_irrad2flux():
    assert np.isclose(irrad2flux(7.5, 650), 2.454132e16)


def test_times_to_cycles():
    times = [[25, 30], [50, 150]]
    t_end = 200
    cycles, Dt_delay = times2cycles(times, t_end)
    assert np.allclose(cycles, [[5,  20], [100,  50]])
    assert Dt_delay == 25


def test_cycles_to_times():
    cycles = [[50, 200], [150, 25]]
    Dt_delay = 50
    times, t_end = cycles2times(cycles, Dt_delay)
    assert np.allclose(times, [[50, 100], [300, 450]])
    assert t_end == 475


def test_round_sig_3():
    assert np.isclose(round_sig(8.348925e12, n=3), 8.35e12)
