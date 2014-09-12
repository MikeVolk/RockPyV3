__author__ = 'mike'
from RockPyV3.Structure.sample import Sample

# test_data = '/Users/mike/Google Drive/__code/RockPyV3/test_data/TEST.hys'
# test_data = '../test_data/TEST.forc'
test_data = '../test_data/FeNi06-F-d-24-M001-002.forc'

S = Sample(name='test')

F = S.add_measurement(mtype='forc', mfile=test_data, machine='VSM')
# print F.return_fitting_surface[4]
# F.calculate_forc(SF=2)
# print F.return_fitting_surface.shape
# F.plt_forc()
# F.plt_residuals()
# F.plt_drift_moment_sat()
# F.plt_drift_field_sat()
# F.plt_forcs()
# F.plt_forc_field_spacing()
# F.plt_ha_hb_space()
F.plt_backfield()
# F.plt_hysteresis()
