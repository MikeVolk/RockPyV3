__author__ = 'Mike'
import numpy as np
from scipy.special import erfc


def normal(x, parameters, dfunc='pdf', *args, **kwargs):
    '''
    parameters:
        amp: amplitude of function
        mean: mean of function - called 'center' in parameters
        standard deviation:
        variance:
        skewness:
        kurtosis:

    :option:
        output:
            pdf [standard] gives probability density function
            cdf

    :math:
       e^(-(x-mu)^2/(2 sigma^2))/(sqrt(2 pi) sigma)

    '''

    amp = parameters['amp'].value
    center = parameters['center'].value
    sig = parameters['sig'].value
    print amp, center, sig

    if dfunc == 'pdf':
        out_pdf = np.exp(-(x - center) ** 2 / 2 * sig ** 2) / np.sqrt(2 * np.pi) * sig
        out_pdf *= amp / max(out_pdf)
        return out_pdf
    if dfunc == 'cdf':
        out_cdf = amp * 0.5 * erfc((center - x) / np.sqrt(2) * sig)
        return out_cdf


def normal_skew(x, parameters, dfunc='pdf', *args, **kwargs):
    '''
    parameters:
        amp: amplitude of function
        mean: mean of function - called 'center' in parameters
        standard deviation:
        variance:
        skewness:
        kurtosis:

    :option:
        output:
            pdf [standard] gives probability density function
            cdf

    :math:
        (e^(-(x-mu)^2/(2 sigma^2)) erfc(-(alpha (x-mu))/(sqrt(2) sigma)))/(sqrt(2 pi) sigma)
    '''

    amp = parameters['amp'].value
    center = parameters['center'].value
    sig = parameters['sig'].value
    skew = parameters['skew'].value
    print amp, center, sig, skew

    if dfunc == 'pdf':
        out_pdf = np.exp(-(x - center) ** 2 / 2 * sig ** 2) * erfc(-(skew * (x - center))/ (np.sqrt(2)*sig)) / (np.sqrt(2 * np.pi) * center)
        out_pdf *= amp / max(out_pdf)
        return out_pdf
    if dfunc == 'cdf':
        out_cdf = amp * 0.5 * erfc((center - x) / np.sqrt(2) * sig)
        return out_cdf

