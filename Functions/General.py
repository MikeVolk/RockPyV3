__author__ = 'Mike'
import logging
import numpy as np
import matplotlib.pyplot as plt
from math import degrees, radians
from math import sin, cos, tan


def create_logger(name):
    log = logging.getLogger(name=name)
    log.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s \t %(levelname)s - %(message)s')
    # fh = logging.FileHandler('RPV3.log')
    # fh.setFormatter(formatter)
    # ch = logging.FileHandler('RPV3.log')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    # log.addHandler(fh)
    log.addHandler(ch)

    return ch#, fh


def differentiate(data_list, diff=1, smoothing=1, norm=False, check=False):
    """
    caclulates a smoothed (if smoothing > 1) derivative of the data

    :param data_list: ndarray - 2D array with xy[:,0] = x_data, xy[:,1] = y_data
    :param diff: int - the derivative to be calculated f' -> diff=1; f'' -> diff=2
          default = 1
    :param smoothing: int - smoothing factor
    :param norm:
    :param check:
    :return:
    """
    log = logging.getLogger('RockPy.FUNCTIONS.general.diff')
    log.info('DIFFERENTIATING\t data << %i derivative - smoothing: %i >>' % (diff, smoothing))
    data_list = np.array(data_list)
    # getting X, Y data
    X = data_list[:, 0]
    Y = data_list[:, 1]

    # ''' test with X^3'''
    # X = np.linspace(-10,10,1000)
    # Y = X**3
    # Y2 = 3*X**2
    # Y3 = 6*X

    # # derivative
    for i in range(diff):
        deri = [[X[i], (Y[i + smoothing] - Y[i - smoothing]) / (X[i + smoothing] - X[i - smoothing])] for i in
                range(smoothing, len(Y) - smoothing)]
        deri = np.array(deri)
        X = deri[:, 0]
        Y = deri[:, 1]
    MAX = max(abs(deri[:, 1]))

    if norm:
        deri[:, 1] /= MAX
    if check:
        if norm:
            plt.plot(data_list[:, 0], data_list[:, 1] / max(data_list[:, 1]))
        if not norm:
            plt.plot(data_list[:, 0], data_list[:, 1])
        plt.plot(deri[:, 0], deri[:, 1])
        plt.show()

    return deri


def rotate(xyz, axis='x', degree=0, *args):
    """


    :rtype : object
    :param x:
    :param y:
    :param z:
    :param axis:
    :param degree:
    :return:
    """
    a = radians(degree)

    RX = [[1, 0, 0],
          [0, cos(a), -sin(a)],
          [0, sin(a), cos(a)]]

    RY = [[cos(a), 0, sin(a)],
          [0, 1, 0],
          [-sin(a), 0, cos(a)]]

    RZ = [[cos(a), -sin(a), 0],
          [sin(a), cos(a), 0],
          [0, 0, 1]]

    if axis.lower() == 'x':
        out = np.dot(xyz, RX)
    if axis.lower() == 'y':
        out = np.dot(xyz, RY)
    if axis.lower() == 'z':
        out = np.dot(xyz, RZ)

    return out