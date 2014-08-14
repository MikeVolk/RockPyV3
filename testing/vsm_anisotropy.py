__author__ = 'mike'
from Structure.sample import Sample
from Plots import general

x_up_files = ['/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.000',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.001',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.002',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.003',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.004',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.005',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.006',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.007',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.008',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.009',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.010',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.011',
              '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-up/LF4c-6c-IRM-aniso-X_up.hys.012',
]
x_left_files = ['/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.000',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.001',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.002',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.003',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.004',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.005',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.006',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.007',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.008',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.009',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.010',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.011',
                '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-left/LF4c-6c-IRM-aniso-X_left.hys.012',
]
x_front_files = ['/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.000',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.001',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.002',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.003',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.004',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.005',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.006',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.007',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.008',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.009',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.010',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.011',
                 '/Users/mike/Google Drive/__code/RockPyV3/test_data/X-front/LF4c-6c-IRM-aniso-X_front.hys.012',
]

sample = Sample(name='test_sample')

for i in x_up_files:
    hys = sample.add_measurement(mtype='hys', mfile=i, machine='vsm')
    hys.paramag_slope_correction()
for i in x_left_files:
    hys = sample.add_measurement(mtype='hys', mfile=i, machine='vsm')
    hys.paramag_slope_correction()
for i in x_front_files:
    hys = sample.add_measurement(mtype='hys', mfile=i, machine='vsm')
    hys.paramag_slope_correction()
general.Hysteresis(sample)