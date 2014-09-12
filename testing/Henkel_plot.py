__author__ = 'mike'
from RockPyV3.Structure.sample import Sample
from RockPyV3.Plots.general import Henkel_Plot
data = '../test_data/FeNi06Eb006-GC-140624.irm'

S = Sample(name='irm_test')
M = S.add_measurement(mtype='irm', mfile=data, machine='vsm')
M = S.add_measurement(mtype='coe', mfile=data, machine='vsm')

Henkel_Plot(sample_list=S)