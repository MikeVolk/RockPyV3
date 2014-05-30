__author__ = 'Mike'

import distributions
from lmfit import minimize, Parameters
from scipy.stats import norm
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import copy


def main():
    def normal():
        data = norm.rvs(10.0, 2.5, size=5000)
        aux = np.histogram(data, bins=100)
        x = aux[1][:-1]
        y = aux[0]

        ''' Guessing '''
        params = Parameters()
        params.add('amp', value=max(y), vary = True)
        params.add('center', value=x[np.where(y == max(y))], vary = True)
        params.add('sig', value=(2 * (x[np.argmin(abs(y) - max(y))] - x[np.argmin(abs(y) / max(y) - 0.5)])))
        p_init = copy.deepcopy(params)

        ''' Fitting '''
        errfunc = lambda p, x, y: distributions.normal(x, params) - y
        fitout = minimize(errfunc, params, args=(x, y))

        print '-----------------------------------'
        print ' Fit results'
        print '-----------------------------------'
        for i in params:
            # if params[i].value != p_init[i].value:
            print '%5s\t: %.3f \t %.3f' % (i, params[i].value, p_init[i].value)

        ''' Plotting '''
        x_fit = np.linspace(min(x), max(x), 100)
        y_fit_pdf = distributions.normal(x_fit, params, dfunc='pdf')
        y_fit_cdf = distributions.normal(x_fit, params, dfunc='cdf')
        y_fit = []
        diff = []
        int = []
        for i in range(len(y_fit_pdf)):
            y_fit.append(sum(y_fit_pdf[0:i]))
            if i > 1:
                int.append(integrate.simps(y_fit_pdf[0:i], x_fit[0:i]))
            else:
                int.append(0)
        # plt.plot(x, y)
        print abs(y_fit/sum(y_fit)-int/sum(int))
        d2 = y_fit - (y_fit/sum(y_fit) - int/sum(int))
        # plt.plot(x_fit, y_fit_cdf/sum(y_fit_cdf), '-')
        plt.plot(x_fit, y_fit/sum(y_fit), '--')
        plt.plot(x_fit, int/sum(int), ':')
        plt.plot(x_fit, d2/sum(d2), '-', color = 'grey')
        plt.xlim(15,17)
        # plt.plot(x_fit, y_fit/max(y_fit)-y_fit_cdf/max(y_fit_cdf))
        # plt.plot(x_fit, int/max(int)-y_fit_cdf/max(y_fit_cdf), '--')
        plt.show()

    normal()


if __name__ == '__main__':
    main()
