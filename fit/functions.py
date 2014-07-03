__author__ = 'mike'
import numpy as np

def hyperbolic_basis(x, parameters):
    '''
    calculation for fit of hyperbolic basis function according to von Dobeneck 1995
    typically only works on saturated samples

    :param x: values for fit
    :param parameters: lmfit(parameters)
    :return:


    ms = saturation magnetizetion
    mrs = remanence saturation
    '''


    func1 = a1 * np.tanh()