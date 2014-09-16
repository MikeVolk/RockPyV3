# coding=utf-8
__author__ = 'Mike'
import datetime
import time
import numpy as np
import csv
import helper
from pprint import pprint
import logging


def sushibar_old(d_file, sample=None):
    log = logging.getLogger('RockPy.READIN')
    log.info('IMPORTING\t automag file: << %s >>' % d_file)
    header = {
        'sample': [0, str], 'site': [1, str], 'type': [2, str], 'run': [3, int], 'time': [4, str],
        'x': [5, float], 'y': [6, float], 'z': [7, float],
        'M': [8, float], 'Dc': [9, float], 'Ic': [10, float], 'Dg': [11, float], 'Ig': [12, float],
        'Ds': [13, float], 'Is': [14, float], 'a95': [15, float], 'sM': [16, float],
        'npos': [17, float], 'Dspin': [18, float], 'Ispin': [19, float],  # 'holder/sample': [20, float],
        # 'cup/sample': [21, float],
        'bl diff/sample': [22, float],  # 'steps/rev': [23, float],
        'par1': [24, float], 'par2': [25, float], 'par3': [26, float], 'par4': [27, float],
        'par5': [28, float], 'par6': [29, float], 'strat_level': [30, float], 'geoaz': [31, float],
        'hade': [32, float], 'dipdir': [33, float], 'dip': [34, float]
    }
    # print [i for i in header if header[i][1]==float]
    reader_object = csv.reader(open(d_file, 'rU'), delimiter='\t')
    aux = [i for i in reader_object][1:]
    for i in range(len(aux)):
        for j in range(len(aux[i])):
            if aux[i][j] == 'None':
                aux[i][j] = 0

    out = {column.lower(): np.array([header[column][1](i[header[column][0]]) for i in aux]) for column in header}
    out['time'] = [time.strptime(i[:19], "%Y-%m-%d %H:%M:%S") for i in out['time']]

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


def sushibar(file, sample, *args, **options):
    log = logging.getLogger('RockPy.READIN.sushibar')
    log.info('IMPORTING\t automag file: << %s >>' % file)

    floats = ['dspin', 'ispin', 'par1', 'dip', 'dipdir', 'geoaz', 'm', 'strat_level', 'a95', 'par5', 'par4', 'par3',
              'par2', 'sm', 'par6', 'dg', 'is', 'hade', 'dc', 'npos', 'bl diff/sample', 'y', 'x', 'ic', 'z', 'ds', 'ig']

    data_f = open(file)
    data = [i.strip('\n\r').split('\t') for i in data_f.readlines()]

    header = ['sample', 'site', 'type', 'run', 'time', 'x', 'y', 'z', 'M', 'Dc', 'Ic', 'Dg', 'Ig', 'Ds', 'Is',
              'a95', 'sM', 'npos', 'Dspin', 'Ispin', ' holder/sample', 'cup/sample', 'bl diff/sample', 'steps/rev',
              'par1', 'par2', 'par3', 'par4', 'par5', 'par6', 'strat_level', 'geoaz', 'hade', 'dipdir', 'dip']

    sample_data = np.array([i for i in data[2:] if i[0] == sample or sample in i[9]])

    if len(sample_data) == 0:
        log.error('UNKNOWN\t sample << %s >> can not be found in measurement file' % sample)
        return None

    out = {header[i].lower(): sample_data[:, i] for i in range(len(header))}

    for i in floats:
        try:
            out[i] = np.array(map(float, out[i]))
        except ValueError:
            for j in range(len(out[i])):
                if out[i][j] == 'None':
                    out[i][j] = '0'
                if out[i][j] == '-':
                    out[i][j] = np.nan
            try:
                out[i] = np.array(map(float, out[i]))
            except ValueError:
                continue

    out['run'] = map(int, out['run'])

    def time_conv(t):
        if t == '-':
            return None

        return time.strptime(t[:19], "%Y-%m-%d %H:%M:%S")

    out['time'] = map(time_conv, out['time'])
    log.info(
        'RETURNING\t %i data types with %i data points for sample: << %s >>' % (len(out), len(out['sample']), sample))
    return out


def jr6(file, sample=None):
    log = logging.getLogger('RockPy.READIN.jr6')
    log.info('IMPORTING\t Agico-JR6 file: << %s >>' % file)
    data_f = open(file)
    data = [i.strip('\n\r').split('\t') for i in data_f.readlines()]
    print data


def cryo_nl(file, sample=None):
    log = logging.getLogger('RockPy.READIN.CRYO_NL')
    log.info('IMPORTING\t cryomag file: << %s >>' % file)
    header = {
        'sample': [0, str],
        # 'coreaz': [1, float], 'coredip': [2, float], 'bedaz': [3, float], 'beddip': [4, float],
        # 'vol': [5, float], 'weight': [6, float],
        'step': [7, float],
        'type': [8, str], 'comment': [9, str], 'time': [10, str], 'mode': [11, str],
        'x': [12, float], 'y': [13, float], 'z': [14, float], 'M': [15, float],
        'sm': [16, float], 'a95': [17, float], 'dc': [18, float], 'ic': [19, float],
        # 'dg': [20, float], 'ig': [21, float],
        # 'ds': [22, float], 'is': [23, float],
    }
    header_no_mode = {
        'sample': [0, str],
        # 'coreaz': [1, float], 'coredip': [2, float], 'bedaz': [3, float], 'beddip': [4, float],
        # 'vol': [5, float], 'weight': [6, float],
        'step': [7, float],
        'type': [8, str], 'comment': [9, str], 'time': [10, str],
        'x': [11, float], 'y': [12, float], 'z': [13, float], 'M': [14, float],
        'sm': [15, float], 'a95': [16, float], 'dc': [17, float], 'ic': [18, float],
        # 'dg': [20, float], 'ig': [21, float],
        # 'ds': [22, float], 'is': [23, float],
    }

    data = helper.import_file(d_file=file, header_skip=1)
    file_header = data[0]

    if 'mode' not in file_header:
        header = header_no_mode

    data = np.array(data[1:], dtype=object)
    data = np.array([[value for value in line] for line in data]).T
    out = {column.lower(): map(header[column][1], data[header[column][0]]) for column in sorted(header)}

    if 'mode' not in file_header:
        out['mode'] = ['results' for i in out['sample']]

    if sample:
        samples = list(set(out['sample']))
        if sample not in samples:
            log.error('UNKNOWN\t sample << %s >> can not be found in measurement file' % sample)
            return None
        else:
            if len(samples) > 1:
                sample_idx = [i for i in range(len(out['sample']))
                              if out['sample'][i] == sample
                              if out['mode'][i] == 'results']
                holder_idx = [i for i in range(len(out['sample']))
                              if out['sample'][i] in ['holder', 'acryl']
                              if out['mode'][i] == 'results']
                holder = {}
                for key in out:
                    holder[key] = np.array([out[key][i] for i in holder_idx])
                    out[key] = np.array([out[key][i] for i in sample_idx])
                out['acryl'] = holder
    out['time'] = [time.strptime(i, "%y-%m-%d %H:%M:%S") for i in out['time']]
    log.info(
        'RETURNING\t %i data types with %i data points for sample: << %s >>' % (len(out), len(out['sample']), sample))

    print(out.keys())
    return out


def cryo_nl2(file, sample, *args, **options):
    log = logging.getLogger('RockPy.READIN.CRYO_NL')
    log.info('IMPORTING\t cryomag file: << %s >>' % file)

    floats = ['x', 'y', 'z', 'm', 'sm', 'a95', 'dc', 'ic', 'dg', 'ig', 'ds', 'is']
    data_f = open(file)
    data = [i.strip('\n\r').split('\t') for i in data_f.readlines()]

    header = ['name', 'coreaz', 'coredip', 'bedaz', 'beddip',
              'vol', 'weight', 'step', 'type', 'comment',
              'time', 'mode', 'x', 'y', 'z',
              'M', 'sM', 'a95', 'Dc', 'Ic', 'Dg', 'Ig', 'Ds', 'Is']
    sample_data = np.array([i for i in data[2:-1] if i[0] == sample or sample in i[9] if i[11] == 'results'])

    holder_data = np.array([i for i in data[2:-1] if i[0].lower() == 'acryl' if i[11] == 'results'])

    try:
        out = {header[i].lower(): sample_data[:, i] for i in range(len(header))}
        out['acryl'] = {header[i].lower(): holder_data[:, i] for i in range(len(header))}
    except IndexError:
        log.error('CANT find sample/holder')

    for i in floats:
        out[i] = map(float, out[i])
        out['acryl'][i] = map(float, out[i])
    out['step'] = map(int, out['step'])

    def time_conv(t):
        return time.strptime(t, "%y-%m-%d %H:%M:%S")

    out['time'] = map(time_conv, out['time'])
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


def vsm_forc(file, sample=None):
    log = logging.getLogger('RockPy.READIN.vsm_forc')
    log.info('IMPORTING\t VSM file: << %s >>' % file)

    file = open(file, 'rU')
    reader_object = file.readlines()
    header = readMicroMagHeader(reader_object)  # get header
    raw_out = [i for i in reader_object][header['meta']['numberoflines']:]  # without header

    # ### header part
    data_header = [i.split('\n')[0] for i in raw_out if
                   not i.startswith('+') and not i.startswith('-') and not i.split() == []][:-1]
    aux = [i for i in data_header[-3:]]
    h_len = len(aux[0]) / len(aux[-1].split()) + 1
    splits = np.array([[i[x:x + h_len] for x in range(0, len(i), h_len)] for i in aux]).T
    splits = ["".join(i) for i in splits]
    splits = [' '.join(j.split()) for j in splits]

    out = [i for i in raw_out if i.startswith('+') or i.startswith('-') or i.split() == []]
    out_data = []
    aux = []
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
            out_data.append(np.array(aux))
            aux = []
    out_data = np.array(out_data)
    out = {splits[i]: np.array([j[:, i] for j in out_data]) for i in range(len(splits))}
    log.info('RETURNING data << %s >> ' % (' - '.join(out.keys())))
    out.update(header)
    return out


def vsm(file, sample=None):
    log = logging.getLogger('RockPy.READIN.vsm')
    log.info('IMPORTING\t VSM file: << %s >>' % file)
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
        log.info('IMPORTING\t MPMS file: << %s >>' % file)
        reader_object = csv.reader(open(file), delimiter=',')
        header = {'Time': [0, float], 'Comment': [1, str], 'Field (Oe)': [2, float], 'Temperature (K)': [3, float],
                  'Long Moment (emu)': [4, float], 'Long Scan Std Dev': [5, float], 'Long Offset (cm)': [6, float],
                  # 'Long Offset Std Dev': [7, float], 'Long Algorithm': [8, float], 'Long Reg Fit': [9, float],
                  # 'Long Reg Factor': [10, float], 'Trans Moment (emu)': [11, float], 'Trans Scan Std Dev': [12, float],
                  # 'Trans Offset (cm)': [13, float], 'Trans Offset Std Dev': [14, float],
                  # Trans Algorithm,Trans Reg Fit,Trans Reg Factor,Long Moment [w/o ABS] (emu),Long Scan Std Dev [w/o ABS],Long Offset [w/o ABS] (cm),
                  # Long Offset Std Dev [w/o ABS],Long Reg Fit [w/o ABS],Trans Moment [w/o ABS] (emu),Trans Scan Std Dev [w/o ABS],Trans Offset [w/o ABS] (cm),
                  # Trans Offset Std Dev [w/o ABS],Trans Reg Fit [w/o ABS],RSO Position (deg),Amplitude (cm),Frequency,Cycles to Average,Scans per Measurement,
                  'Delta Temp (K)': [34, float],  # Error,EC Comp. Running,Using ABS,'}
        }

        aux = [i for i in reader_object][31:]
        out_aux = {column.lower(): np.array([header[column][1](i[header[column][0]]) for i in aux]) for column in
                   header}
        out.append(out_aux)

    out = {key: [out[0][key], out[1][key]] for key in out[0]}
    out['time'] = [map(time.localtime, i) for i in out['time']]
    return out


def simulation(files, sample=None):
    out = []
    return out


def vftb(file, *args, **options):
    '''
    '''
    log = logging.getLogger('RockPy.READIN.vftb')
    log.info('IMPORTING\t VFTB file: << %s >>' % file)
    reader_object = open(file)
    out = [i.strip('\r\n').split('\t') for i in reader_object.readlines()]

    mass = float(out[0][1].split()[1])
    out = np.array(out[4:])
    idx1 = [i for i in range(len(out)) if '' in out[i]]
    idx2 = [i for i in range(len(out)) if 'Set 2:' in out[i]]
    idx3 = [i for i in range(len(out)) if ' field / Oe' in out[i]]
    idx = idx1 + idx2 + idx3
    out = [['0.' if j == 'n/a' else j for j in i] for i in out]
    out = [out[i] for i in range(len(out)) if i not in idx]
    aux = []
    out_aux = []

    # out = [[np.nan if v is 'n/a' else v for v in i] for i in out]
    out = np.array([map(float, i) for i in out[4:]])
    header = {"field": [0, float], "moment": [1, float], "temp": [2, float], "time": [3, float], "std dev": [4, float],
              "suscep / emu / g / Oe": [5, float],
    }

    out = {column.lower(): np.array([header[column][1](i[header[column][0]]) for i in out]) for column in
           header}
    out['moment'] *= mass / 1E3
    out['field'] *= 0.0001
    return out
