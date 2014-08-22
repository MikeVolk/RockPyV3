# -*- coding: utf-8 -*-
__author__ = 'Mike'
import numpy as np
import math
from scipy.special import erf, erfc
import matplotlib.pyplot as plt
from scipy import pi, sqrt, exp
import scipy.integrate as sp_int


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
    center = parameters['mu'].value
    sig = parameters['sig'].value

    if dfunc == 'pdf':
        out_pdf = np.exp(-(x - center) ** 2 / 2 * sig ** 2) / np.sqrt(2 * np.pi) * sig
        out_pdf *= amp / max(out_pdf)
        return out_pdf
    if dfunc == 'cdf':
        out_cdf = amp * 0.5 * erfc((center - x) / np.sqrt(2) * sig)
        return out_cdf


def normal_skew(x, parameters, dfunc='pdf', reverse = False, check=False, *args, **kwargs):
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
    # general.create_logger('RockPy.FIT.SKEW_NORMAL')
    # log = logging.getLogger('RockPy.MEASUREMENT.HYSTERESIS')

    amp = parameters['amp'].value
    center = parameters['mu'].value
    sig = parameters['sig'].value
    skew = parameters['skew'].value

    norm_pdf = (1 / (sig * np.sqrt(2 * math.pi))) * np.exp(-(np.power((x - center), 2) / (2 * np.power(sig, 2))))
    norm_cdf = (0.5 * (1 + erf((skew * ((x - center) / sig)) / (np.sqrt(2)))))

    out_pdf = 2 * amp * norm_pdf * norm_cdf

    if check:
        plt.plot(x, out_pdf)
        plt.show()

    if dfunc == 'pdf':
        return out_pdf

    if dfunc == 'cdf':
        # log.info('INTEGRATING < NOTE: numerical integration: check consistency >')
        y_int = sp_int.cumtrapz(out_pdf, x, initial = 0)
        if reverse:
            y_int = max(y_int) - y_int

        return y_int

def linear(x, parameters):
    print parameters.keys()
    slope = parameters['slope'].value
    intercept = parameters['intercept'].value

    return intercept + slope *x

def GGD(x, parameters):
    '''
    :parameter: center ist die Median der Verteilung
    :parameter: sig ist ahnlich wie die Standardabweichung und ist genau gleich der Standardabweichung im Speziallfall einer Gaußschen Verteilung (die man mit skew=1 und p=2 bekommt).
    :parameter: skew ist ein Parameter, der die Asymmetrie kontrolliert.
                Mit skew=1 ist die Funktion symmetrisch. Positive werte fuhren zu eine langsamere Abnahme nach links. Negative Werte bedeuten eine langsamere Abnahme nach rechts.
    :parameter: kurt ist ein Parameter, der die Rechteckigkeit der Funktion kontrolliert. kurt=2 ergibt die "normale" Form der Gaußschen Verteilung. 0<kurt<2 fuhrt zu spitzigere Verteilungen. p>2 ergibt "Kofferformige" Verteilungen die mit p->Unendlich zu einen Rechteck konvergieren.

    skew =1 und kurt = 2 ergibt die Gaußsche Verteilung wo center den Mittelwert ist und sig die Standardabweichung
    kurt = 2 ergibt asymmetrische "Gauß-Artige" Verteilungen

    In manchen Fallen kann p=2 gesetzt werden ohne dieses Parameter weiter zu optimieren.
    '''

    amp = parameters['amp'].value
    center = parameters['center'].value
    sig = parameters['sig'].value
    skew = parameters['skew'].value
    kurt = parameters['kurt'].value

    # print '-------------'
    # print 'amp, center, sig, skew, kurt'
    # print amp, center, sig, skew, kurt

    c1 = special.gamma(1 + 1 / kurt)

    out = amp * abs(skew * np.exp(skew * (x - center) / sig) + np.exp((x - center) / (skew * sig)) / skew) \
          * np.exp(
        -abs(np.log((np.exp(skew * (x - center) / sig) + np.exp((x - center) / (skew * sig))) / 2)) ** kurt / 2) / (
          2 ** (1 + 1 / kurt) * sig * c1 * ( np.exp(skew * (x - center) / sig) + np.exp((x - center) / (skew * sig))))

    return out