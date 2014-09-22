__author__ = 'mike'
from RockPyV3.Structure.sample import Sample

# test_data = '../test_data/MUCVSM-LF4C_IVe.forc'
test_data = '../test_data/FeNi06-F-d-24-M001-002.forc'

S = Sample(name='test')

F = S.add_measurement(mtype='forc', mfile=test_data, machine='VSM')
# F.calculate_forc(SF=4)
# F.plt_ha_hb_space()
# F.plt_residuals()
F.plt_drift_moment_sat()
F.plt_drift_field_sat()
# F.plt_forc()
# F.plt_forc_field_spacing()
# F.plt_ha_hb_space()
# F.plt_backfield()
# F.plt_hysteresis()
