__author__ = 'mike'
import numpy as np
from scipy import special


def hyperbolic_basis(parameters, x):
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


def func_tanh(parameters, x):
    '''
    '''
    amp = parameters['amp'].value
    shift = parameters['shift'].value
    out = -amp * np.tanh(x + shift)
    return out


def log_normal(parameters, x):
    '''

    :math:`E_i` represents the mean value of a correspond- ing log-Gaussian distribution, and therefore, the
    logarithmic field value of the maximum gradient.
    :math:`G_i` describes the standard deviation or half-width of the distribution

    math::

       M_{\text{fit}} = B * G * \sqrt{\frac{\pi}{2} \erf (\frac{(x-E)}{G}) - erf(\frac{-E}{G})
    '''
    B = parameters['B'].value
    E = parameters['E'].value
    G = parameters['G'].value

    pi = np.pi

    out = B * G * ( np.sqrt(pi) / 2 ) * ( special.erf((x - E) / G) - special.erf(-E / G) )
    return out