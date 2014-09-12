# -*- coding: utf-8 -*-
__author__ = 'mike'
from Structure.sample import Sample
from Functions import general
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import stats

matplotlib.rcParams['backend'] = 'TkAgg'


def hys_test():
    file = '/Users/mike/Google Drive/__code/RockPyV3/test_data/FeNi100-A-b-048-M001-140428-002.hys'
    files = ['test_data/FeNi100-A-b-048-M001-140428-002.hys',
             'test_data/FeNi100-A-c-105-M001-140428-001.hys',
             'test_data/FeNi100-A-d-140-M001-140428-001.hys',
             'test_data/FeNi100-A-d-200-M001-140428-001.hys',
             'test_data/FeNi100-A-d-322-M001-140428-001.hys']
    # file = ['/Users/mike/Google Drive/__PHD/__Projects/003 FeNi/04 data/FeNi10-B/VSM/FeNi10-B-a-054-GC-140623.hys']
    stats_all = []
    # for file in files:
    aux = []
    S = sample.Sample(name='FeNi100-A-d-322-M001')
    H1 = S.add_measurement(mtype='hys', mfile=file, machine='vsm')
    # test, stats = H1.approach_sat(branch='up_field', step=0.001, check=True)
    # aux.append(stats)
    # plt.plot(test[:, 0], test[:, 1])
    # test, stats = H1.approach_sat(branch='down_field', step=0.001, check=False)
    # plt.plot(H1.uf_interp[:,0],H1.uf_interp[:,1], label='interp1')
    # plt.plot(H1.up_field[:,0],H1.up_field[:,1], '+', label='data')
    # plt.plot(H1.df_interp[:,0],H1.df_interp[:,1])

    data = H1.up_field
    for i in np.arange(1.5 , 2, 0.1):
        plt.plot(data[:, 0] ** -i, data[:, 1], label = str(i))

    ''' approach to saturation test '''
    # test, stats = H1.approach_sat(branch='up_field', step=0.01, check=True, corrected=False)
    # test, stats = H1.approach_sat(branch='down_field', step=0.01, check=True, corrected=False)
    # plt.plot(test[:, 0], test[:, 1])
    plt.legend(loc='best')
    plt.xlim([0,6000])
    plt.show()

def initial_state_test():

    p11_nrm = '/Users/mike/Google Drive/__PHD/__Projects/001 Influence of Pressure on Paleointensities/04 data/LF6C/MUCSUSH-LF6C_P11_TRM(50µT)_xxx_NRM-2704.NRM'
    p11_af = '/Users/mike/Google Drive/__PHD/__Projects/001 Influence of Pressure on Paleointensities/04 data/LF6C/MUCSUSH-LF6C_P11_TRM(50µT)_xxx_AF-2706.AF'

    c14 = Sample(name='14c')

    sample = c14

    af = sample.add_measurement(mtype='af-demag', mfile=p11_af, machine='sushibar')
    af.add_treatment(ttype='pressure', options={'p_max': 0.3, 'p_seen': 0.3})
    af.add_initial_state(mtype='trm', machine='sushibar', mfile=p11_nrm)
    print af.initial_state
    # test = af.normalize(dtype='data', norm='initial_state')
    # print test.m#,af.data.m, af.initial_state.m, af.data.m / af.initial_state.m

if __name__ == '__main__':
    # hys_test()
    initial_state_test()