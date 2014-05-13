__author__ = 'Mike'
import logging


def get_data(file, header, column_name, header_skip):
    log = logging.getLogger('RockPy.READIN')
    '''
    help function
    '''
    reader_object = csv.reader(open(file, 'rU'), delimiter='\t')
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