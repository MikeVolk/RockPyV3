__author__ = 'Mike'
import logging
import numpy as np
import matplotlib.pyplot as plt
import math as math
import scipy.special as sp
import distributions, functions
from lmfit import Parameters, minimize, printfuncs
import matplotlib.pyplot as plt
from RockPyV3.Functions import general
from scipy.optimize import curve_fit


class Fit(object):
    def __init__(self, x_data, y_data, paramters=None, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.FITTING')

        self.x = x_data
        self.y = y_data
        
        if not paramters:
            self.parameters = Parameters()
        else:
            self.parameters = paramters


    def check(self):
        pass


class HyperbolicBase(Fit):

    def __init__(self, x_data, y_data, paramters = None):
        super(HyperbolicBase, self).__init__(x_data, y_data, paramters = None)


class Tanh(Fit):
    general.create_logger('RockPy.FITTING.tanh')

    def __init__(self, x_data, y_data, paramters=None):
        self.log = logging.getLogger('RockPy.FITTING.TANH')
        super(Tanh, self).__init__(x_data, y_data, paramters=None, log=self.log)

        self.fit_x = np.arange(min(self.x), max(self.x), 0.01)

        amp = np.max(abs(self.y))
        zero_idx = np.argmin(abs(self.y))

        self.parameters.add('amp', value=amp)
        self.parameters.add('shift', value=self.x[zero_idx])

        ''' Fitting '''
        errfunc = lambda p, x, y: functions.func_tanh(self.parameters, x) - y

        fitout = minimize(errfunc, self.parameters, args=(self.x, self.y))

        self.fit_y = functions.func_tanh(self.parameters, self.fit_x)

    def calc(self, x_value):
        out = functions.func_tanh(self.parameters, x_value)
        return out


class Log_Normal(Fit):
    general.create_logger('RockPy.FITTING.log_normal')

    def __init__(self, x_data, y_data, paramters=None):
        self.log = logging.getLogger('RockPy.FITTING.TANH')
        super(Log_Normal, self).__init__(x_data, y_data, paramters=None, log=self.log)
        self.log.info('FITTING\t log normal distribution to data')
        self.fit_x = np.linspace(min(self.x), max(self.x), 1000)
        self.parameters.add('B', value=0.01)
        self.parameters.add('E', value=1.0)
        self.parameters.add('G', value=0.5)

        ''' Fitting '''
        errfunc = lambda p, x, y: functions.log_normal(self.parameters, x) - y
        fitout = minimize(errfunc, self.parameters, args=(self.x, self.y))

        self.fit_y = functions.log_normal(self.parameters, self.fit_x)

    def calc(self, x_value):
        out = functions.log_normal(self.parameters, x_value)
        return out


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