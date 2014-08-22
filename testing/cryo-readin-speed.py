import time


def speed_test_nl():
    from Structure.sample import Sample

    d_file = '/Users/mike/Google Drive/__code/RockPyV3/test_data/NLCRY-LF4C_P2-140801.TT'
    sample = Sample(name='1a')
    sample.add_measurement(mtype='palint', machine='cryo_nl', mfile=d_file)

def speed_test_sushi():
    from Structure.sample import Sample
    d_file = '/Users/mike/Google Drive/__code/RockPyV3/test_data/LF4c-6c_pdemag_P00.AF'
    sample = Sample(name='6c')
    af = sample.add_measurement(mtype='af-demag', machine='sushibar', mfile=d_file, mag_method='IRM')
    af.plot()

if __name__ == '__main__':
    import ReadIn
    start = time.time()
    # file = '/Users/mike/Google Drive/__code/RockPyV3/test_data/LF4c-6c_pdemag_P00.AF'
    # ReadIn.machines.sushibar_old(file, '6c')
    speed_test_sushi()

    elapsed = (time.time() - start)
    print elapsed