__author__ = 'mike'
from RockPyV3.Structure.sample import Sample

data = '../test_data/irm_test.irm'

S = Sample(name='irm_test')
M = S.add_measurement(mtype='irm', mfile=data, machine='vsm')

M.plt_irm_derivative(smoothing=1, norm=True)