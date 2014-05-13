__author__ = 'Mike'
import logging


def create_logger(name):
    log = logging.getLogger(name=name)
    log.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler('RPV3.log')
    fh.setFormatter(formatter)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(fh)
    log.addHandler(ch)

    return ch, fh