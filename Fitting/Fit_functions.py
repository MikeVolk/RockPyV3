import numpy as np
from scipy import tanh, cosh, sinh
__author__ = 'Mike'


def langevin(X, p): # p[0] = y0, p[1]= C, p[2]=xc
    X = np.float128(X)
    p = np.float128(p)
    diff = X - p[2]
    diff = np.array(diff)
    COSH = cosh(diff) / sinh(diff)
    Back = 1 / diff
    Y = p[0] + p[1] * (COSH - Back)
    Y = np.float64(Y)
    return Y


def func_tanh(X, p): # p[0] = y0, p[1]= C, p[2]=xc
    '''

    :param X: list
    :param p: list
              p[0] = y0, p[1]= C, p[2]=xc
    Y = p[0] + p[1] * tanh( (X - p[2]) / p[3] )
    :return: list


    '''
    X = np.float128(X)
    X = np.array(X)
    p = np.float128(p)

    if p[3]==0:
        p[3]=1
    Y = p[0] + p[1] * tanh( (X - p[2]) / p[3] )

    if len(p) > 8:
        Y += p[4] + p[5] * tanh( (X - p[6]) / p[7])

    Y = np.float64(Y)

    return Y

def curvature(X, p):
    X = np.array(X)
    Y = p[0] + p[1] * X + p[2] * X**2
    return Y


def gauss_ramon(x,m,s,k,p):
    pass