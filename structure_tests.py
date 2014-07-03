__author__ = 'mike'
from Structure import sample
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


if __name__ == '__main__':
    hys_test()