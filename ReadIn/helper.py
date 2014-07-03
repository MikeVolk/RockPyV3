__author__ = 'Mike'
import logging
import csv
import numpy as np


def get_data(file, header, column_name, header_skip, delimiter='\t'):
    log = logging.getLogger('RockPy.READIN.get_data')
    '''
    help function
    '''
    reader_object = csv.reader(open(file), delimiter=delimiter)
    for i in range(header_skip):
        reader_object.next()
    out = [line[header[column_name][0]] for line in reader_object]
    for i in range(len(out)):
        if out[i] == 'None':
            out[i] = 0.
    try:
        out = map(header[column_name][1], out)
    except ValueError:
        log.error(column_name)
    return np.array(out[header_skip:])


def check_duplicates(T_xyz_list):
    log = logging.getLogger('RockPy.READIN.check_duplicates')

    if len(set(T_xyz_list[:, 0])) == len(T_xyz_list[:, 0]):
        log.debug('NO DUPLICATES FOUND in list')
        return T_xyz_list
    else:
        log.info('ATTENTION: DUPLICATES FOUND')
        log.info('ATTENTION: DELETING DUPLICATES - keeping last entry')
        T = T_xyz_list[:, 0]
        idx = [np.where(T == i)[0] for i in sorted(set(T))]
        out = np.array([T_xyz_list[i[-1]] for i in idx])
        return out


def get_type(dict, type):
    log = logging.getLogger('RockPy.READIN.get_type')
    if type in dict['type']:

        log.debug('FOUND:\t measurement steps << %s >>' % (type))
    else:
        log.error('CANT FIND\t << %s >> in dict' % (type))
        log.info('RETURNING:\t None')
        return None
    type_indx = np.where(dict['type'] == type)[0]
    out = np.array([[dict['step'][i],
                     dict['x'][i], dict['y'][i], dict['z'][i], dict['m'][i]]
                    for i in type_indx])
    out = check_duplicates(out)
    return out


def claculate_difference(T_xyzm_list1, T_xyzm_list2):
    # print T_xyzm_list1[0]
    # print T_xyzm_list2[0]
    out = np.array([[T_xyzm_list1[i][0],
                     T_xyzm_list1[i][1] - T_xyzm_list2[i][1],  #X
                     T_xyzm_list1[i][2] - T_xyzm_list2[i][2],  #Y
                     T_xyzm_list1[i][3] - T_xyzm_list2[i][3],  #Z
                     np.linalg.norm([T_xyzm_list1[i][1] - T_xyzm_list2[i][1],
                                     T_xyzm_list1[i][2] - T_xyzm_list2[i][2],
                                     T_xyzm_list1[i][3] - T_xyzm_list2[i][3]])
                    ]
                    for i in range(len(T_xyzm_list1))])
    return out