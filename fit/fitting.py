__author__ = 'Mike'
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math as math
import scipy.special as sp
import distributions
from lmfit import Parameters, minimize, printfuncs
import matplotlib.pyplot as plt
class Fit(object):
    def __init__(self, x_data, y_data, paramters = None):
        self.x_data = x_data
        self.y_data = y_data
        
        if not paramters:
            self.parameters = Parameters()
        else:
            self.parameters = paramters


    def check(self):
        pass


class HyperbolicBase(Fit):

    def __init__(self, x_data, y_data, paramters = None):
        super(HyperbolicBase, self).__init__(x_data, y_data, paramters = None)


def normal_skewed(x, y, parameters=None, dfunc='pdf', check=False, *args, **kwargs):
    # verbous.INFO('FITTING SKEWED-GAUSSIAN')

    if parameters == None:
        parameters = Parameters()
        parameters = Parameters()
        parameters.add('amp', value=max(y))
        max_idx = np.argmin(abs(y - max(y)))
        parameters.add('mu', value=x[max_idx])
        parameters.add('sig', value=15.)
        parameters.add('skew', value=-5.)

    if y[0] > y[-1]:
        reverse = True
    else:
        reverse = False
    errfunc = lambda p, x, y: distributions.normal_skew(x, parameters, reverse=reverse, dfunc=dfunc) - y
    fitout = minimize(errfunc, parameters, args=(x, y))

    # verbous.INFO('FITTED PARAMETERS:')

    if check:
        printfuncs.report_fit(fitout.params)
        x_fit = np.linspace(0, 600, 1000)
        y_fit = distributions.normal_skew(x_fit, parameters, reverse=reverse, dfunc=dfunc)
        plt.plot(x, y)

        plt.plot(x_fit, y_fit)
        plt.show()

    return parameters