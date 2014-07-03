# coding=utf-8
__author__ = 'Mike'
import datetime
import time
import numpy as np
import csv
import helper
from pprint import pprint
import logging


def sushibar(file, sample=None):
    log = logging.getLogger('RockPy.READIN')
    log.info('IMPORTING\t automag file: << %s >>' % (file))
    header = {
        'sample': [0, str], 'site': [1, str], 'type': [2, str], 'run': [3, int], 'time': [4, str],
        'x': [5, float], 'y': [6, float], 'z': [7, float],
        'M': [8, float], 'Dc': [9, float], 'Ic': [10, float], 'Dg': [11, float], 'Ig': [12, float],
        'Ds': [13, float], 'Is': [14, float], 'a95': [15, float], 'sM': [16, float],
        'npos': [17, float], 'Dspin': [18, float], 'Ispin': [19, float],  # 'holder/sample': [20, float],
        #'cup/sample': [21, float],
        'bl diff/sample': [22, float],  #'steps/rev': [23, float],
        'par1': [24, float], 'par2': [25, float], 'par3': [26, float], 'par4': [27, float],
        'par5': [28, float], 'par6': [29, float], 'strat_level': [30, float], 'geoaz': [31, float],
        'hade': [32, float], 'dipdir': [33, float], 'dip': [34, float]
    }
    reader_object = csv.reader(open(file, 'rU'), delimiter='\t')
    aux = [i for i in reader_object][1:]
    for i in range(len(aux)):
        for j in range(len(aux[i])):
            if aux[i][j] == 'None':
                aux[i][j] = 0

    out = {column.lower(): np.array([header[column][1](i[header[column][0]]) for i in aux]) for column in header}
    out['time'] = [time.strptime(i, "%y-%m-%d %H:%M:%S") for i in out['time']]

    if sample:
        samples = list(set(out['sample']))
        if sample not in samples:
            log.error('UNKNOWN\t sample << %s >> can not be found in measurement file' % sample)
            return None
        else:
            if len(samples) > 1:
                sample_idx = np.where(out['sample'] == sample)[0]

                for key in out:
                    out[key] = np.array([out[key][i] for i in sample_idx])
                    # print out

    log.info(
        'RETURNING\t %i data types with %i data points for sample: << %s >>' % (len(out), len(out['sample']), sample))
    return out


def cryo_nl(file, sample=None):
    log = logging.getLogger('RockPy.READIN.CRYO_NL')
    log.info('IMPORTING\t cryomag file: << %s >>' % (file))
    header = {
        'sample': [0, str],
        # 'coreaz': [1, float], 'coredip': [2, float], 'bedaz': [3, float], 'beddip': [4, float],
        #'vol': [5, float], 'weight': [6, float],
        'step': [7, float],
        'type': [8, str], 'comment': [9, str], 'time': [10, str], 'mode': [11, str],
        'x': [12, float], 'y': [13, float], 'z': [14, float], 'M': [15, float],
        'sm': [16, float], 'a95': [17, float], 'dc': [18, float], 'ic': [19, float],
        # 'dg': [20, float], 'ig': [21, float],
        # 'ds': [22, float], 'is': [23, float],
    }
    out = {column.lower(): helper.get_data(file, header, column, header_skip=2) for column in header}

    if sample:
        samples = list(set(out['sample']))
        if sample not in samples:
            log.error('UNKNOWN\t sample << %s >> can not be found in measurement file' % sample)
            return None
        else:
            if len(samples) > 1:
                sample_idx = np.where((out['sample'] == sample) & (out['mode'] == 'results'))[0]
                acryl_idx = np.where((out['sample'] == 'acryl') & (out['mode'] == 'results'))[0]
                acryl = {}
                for key in out:
                    acryl[key] = np.array([out[key][i] for i in acryl_idx])
                    out[key] = np.array([out[key][i] for i in sample_idx])
                out['acryl'] = acryl
    out['time'] = [time.strptime(i, "%y-%m-%d %H:%M:%S") for i in out['time']]
    log.info(
        'RETURNING\t %i data types with %i data points for sample: << %s >>' % (len(out), len(out['sample']), sample))
    return out


def readMicroMagHeader(lines):
    sectionstart = False
    sectiontitle = None
    sectioncount = 0
    header = {}
    # section header: CAPITALS, no spaces
    lc = 0
    for l in lines:
        lc += 1
        sl = l.strip()  # take away any leading and trailing whitespaces
        if lc == 1 and not sl.startswith("MicroMag 2900/3900 Data File"):  # check first line
            print( "No valid MicroMag file. Header not found in first line.")
            return None

        if len(sl) == 0:  # empty line
            sectionstart = True
            continue  # go to next line
        if sectionstart:  # previous line was empty
            sectionstart = False
            if sl.isupper():  # we have a capitalized section header
                # log.INFO('reading header section %s' % sl)
                sectiontitle = sl
                header[sectiontitle] = {}  # make new dictionary to gather data fields
                sectioncount += 1
                continue  # go to next line
            else:  # no captitalized section header after empty line -> we are probably at the end of the header
                # verbous.INFO('reading header finished at line %d' % lc)
                break  # end of header
        if sectiontitle != None:  # we are within a section
            # split key value at fixed character position
            key = sl[:31].strip()
            value = sl[31:].strip(' "')
            if len(key) != 0:
                header[sectiontitle][key] = value
    header['meta'] = {}
    header['meta']['numberoflines'] = lc - 1  # store header length
    header['meta']['numberofsections'] = sectioncount
    return header


def vsm(file, sample=None):
    log = logging.getLogger('RockPy.READIN.vsm')
    log.info('IMPORTING\t VSM file: << %s >>' % (file))
    file = open(file, 'rU')
    reader_object = file.readlines()
    header = readMicroMagHeader(reader_object)  # get header
    out = [i for i in reader_object][header['meta']['numberoflines']:]  # without header
    aux = []
    out_data = []
    linecount = 0
    adj = False
    for i in out:
        if len(i) != 1:
            if i.strip() != '':
                d = i.strip('\n').split(',')
                try:
                    d = [float(i) for i in d]
                    aux.append(d)
                except:
                    log.debug('%s' % d)
                    if 'Adjusted' in d[0].split():
                        adj = True
                    pass
        else:
            out_data.append(aux)
            aux = []
    segments = np.array(out_data[0])
    if adj:
        header['adjusted'] = True
    out_data = np.array(out_data[1:])
    return segments, out_data, header


def mpms(files, sample=None):
    '''
    takes list of files, with files[0] being the cooling files[1] being the warming curve
    '''
    out = []
    for file in files:
        out_aux = {}
        log = logging.getLogger('RockPy.READIN.mpms')
        log.info('IMPORTING\t MPMS file: << %s >>' % (file))
        reader_object = csv.reader(open(file), delimiter=',')
        header = {'Time': [0, float], 'Comment': [1, str], 'Field (Oe)': [2, float], 'Temperature (K)': [3, float],
                  'Long Moment (emu)': [4, float], 'Long Scan Std Dev': [5, float], 'Long Offset (cm)': [6, float],
                  # 'Long Offset Std Dev': [7, float], 'Long Algorithm': [8, float], 'Long Reg Fit': [9, float],
                  #'Long Reg Factor': [10, float], 'Trans Moment (emu)': [11, float], 'Trans Scan Std Dev': [12, float],
                  #'Trans Offset (cm)': [13, float], 'Trans Offset Std Dev': [14, float],
                  # Trans Algorithm,Trans Reg Fit,Trans Reg Factor,Long Moment [w/o ABS] (emu),Long Scan Std Dev [w/o ABS],Long Offset [w/o ABS] (cm),
                  # Long Offset Std Dev [w/o ABS],Long Reg Fit [w/o ABS],Trans Moment [w/o ABS] (emu),Trans Scan Std Dev [w/o ABS],Trans Offset [w/o ABS] (cm),
                  # Trans Offset Std Dev [w/o ABS],Trans Reg Fit [w/o ABS],RSO Position (deg),Amplitude (cm),Frequency,Cycles to Average,Scans per Measurement,
                  'Delta Temp (K)':[34,float],#Error,EC Comp. Running,Using ABS,'}
        }

        aux = [i for i in reader_object][31:]
        out_aux = {column.lower(): np.array([header[column][1](i[header[column][0]]) for i in aux]) for column in header}
        out.append(out_aux)

    out = {key:[out[0][key], out[1][key]] for key in out[0]}
    out['time'] = [ map(time.localtime, i ) for i in out['time']]
    return out