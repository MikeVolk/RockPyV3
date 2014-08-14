import time


def speed_test():
    from Structure.sample import Sample

    file = '/Users/mike/Google Drive/__code/RockPyV3/test_data/NLCRY-LF4C_P2-140801.TT'
    sample = Sample(name='1a')
    sample.add_measurement(mtype='palint', machine='cryo_nl', mfile=file)


if __name__ == '__main__':
    start = time.time()
    speed_test()
    elapsed = (time.time() - start)
    print elapsed