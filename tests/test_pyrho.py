import numpy as np

import pyrho as pyr


def test_run():

    results = pyr.run(mods=6, prots='step', sims='Python')

    assert results['Python']['step']['6'].protocol == 'step'
    assert results['Python']['step']['6'].nRuns == 1
    n_runs = results['Python']['step']['6'].nRuns
    assert results['Python']['step']['6'].phis == [1e+17, 1e+16]
    assert results['Python']['step']['6'].nPhis == 2
    n_phis = results['Python']['step']['6'].nPhis
    assert results['Python']['step']['6'].Vs == [70, 40, 10, -10, -40, -70]
    assert results['Python']['step']['6'].nVs == 6
    n_Vs = results['Python']['step']['6'].nVs
    expected = [
        [
            [{'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': 0.28940366330873857, 'I_ss_': 0.11807373828840971},
            {'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': 0.21805856349284514, 'I_ss_': 0.08896566637421462},
            {'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': 0.07472090687056863, 'I_ss_': 0.030485366707663785},
            {'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': -0.09428461910084183, 'I_ss_': -0.03846716144867489},
            {'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': -0.5527999736167195, 'I_ss_': -0.2255367422250922},
            {'Dt_total': 275.0, 't_peak_': 2.300, 'I_peak_': -1.4739917108381357, 'I_ss_': -0.6013735607732101}],
            [{'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': 0.21553773644434285, 'I_ss_': 0.0830458094866926},
            {'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': 0.16240239895447683, 'I_ss_': 0.06257298098348528},
            {'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': 0.05564952063087659, 'I_ss_': 0.021441532998240993},
            {'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': -0.07021988992873568, 'I_ss_': -0.027055436775936154},
            {'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': -0.41170610509076594, 'I_ss_': -0.15862868067515787},
            {'Dt_total': 275.0, 't_peak_': 6.600, 'I_peak_': -1.0977775238209402, 'I_ss_': -0.42296919604863814}]
        ]]
    for run in range(n_runs):
        for phi in range(n_phis):
            for volt in range(n_Vs):
                pc = results['Python']['step']['6'].trials[run][phi][volt]
                features = expected[run][phi][volt]
                for k, v in features.items():
                    assert np.isclose(getattr(pc, k), v)
                
                

# for run in range(1):
#     for phi in range(2):
#         for volt in range(6):
#             pc = results['Python']['step']['6'].trials[run][phi][volt]
#             print(pc.Dt_total, pc.t_peak_, pc.I_peak_, pc.I_ss_)
