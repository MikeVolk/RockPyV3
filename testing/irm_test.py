__author__ = 'mike'
from RockPyV3.Structure.sample import Sample
import RockPyV3.Plots.general as Rplt

data = '../test_data/irm_test.irm'

S = Sample(name='irm_test')
M = S.add_measurement(mtype='irm', mfile=data, machine='vsm')

# Rplt.IRM(S)
M.plt_irm_unmixing(norm=False)