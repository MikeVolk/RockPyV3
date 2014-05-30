# -*- coding: utf-8 -*-
#  from RockPy.helper import skewed_gauss
from matplotlib import pyplot as plt, pyplot
from numpy import amax, amin, arange  #, range
from pylab import *
from scipy import optimize, stats, special
from scipy import cosh
import numpy as np
import Fit_functions
from lmfit import Parameters
from pprint import pprint
import warnings
# from RockPy_OLD import PalInt


def AraiLine(X, Y, *args):
    slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
    x = np.array(X)
    line = slope * x + intercept
    return line, x, slope, r_value, std_err



    #                 DataDict[Sample]['PalInt'][Pressure]['Stats'] = {
    #                     'B_anc': -LabField * specimen_b,
    #                     'Sigma B_anc' : specimen_b_sigma * LabField,
    #                     'Slope': specimen_b,
    #                     'SigmaSlope': specimen_b_sigma,
    #                     'beta': specimen_b_beta,
    #                     'f': specimen_f,
    #                     'g': specimen_g,
    #                     'MAD': specimen_int_MAD,
    #                     'fvds': specimen_fvds,
    #                     'q': specimen_q,
    #                     'NRM/TRM' : NRMTRMratio,
    #                 }


def AraiPlot(DataDict, Sample):
    xi = [v for v in DataDict['pTRM'][:, 1] / DataDict['TH'][0, 1]]
    #mpTRM = DataDict['pTRM'][-1,1]/DataDict['TH'][0,1]
    y = [v for v in DataDict['TH'][:, 1] / DataDict['TH'][0, 1]]
    slope, intercept, r_value, p_value, std_err = stats.linregress(xi, y)
    x = np.array(xi)
    line = slope * x + intercept
    return line, x, slope, r_value, std_err


def PalIntTHSteps(DataDict, Sample, Quantity='M', InitialGuess=[20., 25., -1e-6]):
    #Get Data
    DataX = np.array(PalInt.ExtractPalIntData(DataDict, Sample, 'TH', Quantity)[:, 0])
    DataY = np.array(PalInt.ExtractPalIntData(DataDict, Sample, 'TH', Quantity)[:, 1])
    DataY = DataY / DataY[0]
    # Fit the function
    fitfunc = lambda p, x: 0.5 * (tanh((x / -p[0]) + p[1]) + 1) + p[2] * x ** 2  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y  # Distance to the target function

    FinalFit, success = optimize.leastsq(errfunc, InitialGuess[:], args=(DataX, DataY))

    print '----------------------------------------------'
    print 'Fitting TH-Step Data'
    print 'initial Parameters: %.3f, %.3f, %.3e' % (InitialGuess[0], InitialGuess[1], InitialGuess[2])
    print 'fitted Parameters:  %.3f, %.3f, %.3e  ' % (FinalFit[0], FinalFit[1], FinalFit[2])
    print '----------------------------------------------'

    XValues = linspace(DataX.min(), DataX.max(), 1000)  #DataX#
    Data = fitfunc(FinalFit, XValues)
    FittedData = np.array([(XValues[i], Data[i]) for i in range(len(XValues))])
    FinalFit = list(FinalFit)
    FinalFit.append(0.5)
    return FittedData, FinalFit


def PalIntpTRMSteps(DataDict, Sample, Quantity='M', InitialGuess=[0.5, 20., 25., 0.]):
    if not 'pTRM' in DataDict[Sample]['PalInt'][0].keys():
        DataDict = PalInt.GetPalIntDataDict(DataDict)
    #Get Data
    DataX = np.array(PalInt.ExtractPalIntData(DataDict, Sample, 'TH', Quantity)[:, 0])
    THmax = max(np.array(PalInt.ExtractPalIntData(DataDict, Sample, 'TH', Quantity)[:, 1]))
    DataY = np.array(PalInt.ExtractPalIntData(DataDict, Sample, 'pTRM', Quantity)[:, 1])
    DataY = DataY / THmax
    #     for i in range(len(DataX)):
    #             print DataX[i] , ', ' ,DataY[i]
    # Fit the function
    fitfunc = lambda p, x: p[3] * ( tanh(x / p[0] - p[1]) + 1) - p[2] * x ** 2  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y  # Distance to the target function
    FinalFit, success = optimize.leastsq(errfunc, InitialGuess[:], args=(DataX, DataY))
    print '----------------------------------------------'
    print 'Fitting pTRM-gained Data'
    print 'initial Parameters: %.3f, %.3f, %.3f, %.3e' % (
        InitialGuess[3], InitialGuess[0], InitialGuess[1], InitialGuess[2])
    print 'fitted Parameters:  %.3f, %.3f, %.3f, %.3e  ' % (FinalFit[3], FinalFit[0], FinalFit[1], FinalFit[2])
    print '----------------------------------------------'

    XValues = linspace(DataX.min(), DataX.max(), 1000)  #DataX#linspace(DataX.min(), DataX.max(), 100)
    Data = fitfunc(FinalFit, XValues)
    #### replace 0.5 from TH fit with value from pTRM Fit
    InitialGuess.pop(3)
    print InitialGuess
    InitialGuess.append(FinalFit[3])
    print InitialGuess
    DataTHFit = fitfunc(InitialGuess, XValues)
    FittedData = np.array([(XValues[i], Data[i]) for i in range(len(XValues))])
    THFittedData = np.array([(XValues[i], DataTHFit[i]) for i in range(len(XValues))])
    return FittedData, THFittedData, FinalFit


def fit_langevin(XY, **kwargs):
    check = kwargs.get('check', False)
    XY = XY[:-5]

    X = XY[:, 0]
    Y = XY[:, 1]

    Y /= max(Y)

    ymax = max(Y)
    ymin = min(Y)
    C = (ymax - ymin) / 2

    max_idx = np.argmin(abs(Y - 0.5))

    y0 = Y[max_idx] + 0.1
    xc = X[max_idx] + 25

    #print Y

    p0 = [y0, C, xc]  # Inital guess is a normal distribution

    errfunc = lambda p, X, Y: Fit_functions.langevin(X, p) - Y  # Distance to the target function

    p1, success = optimize.leastsq(errfunc, p0[:], args=(X, Y))

    X_fit = arange(min(X), max(X), 0.01)
    Y_fit = Fit_functions.langevin(X_fit, p1)

    plt.plot(X, Y, 'b')
    #plt.plot(xc, y0, 'o')
    plt.plot(X_fit, Y_fit, 'r')
    plt.show()


def fit_tanh(XY, **kwargs):
    '''
    This function accepts a numpy.list with list[:,0] = X values and
    list[:,1] = Y-values. A tanh function is then fitted to the data.
    The functions initial guess is calculated as follows.

    .. math::

        f(x) = c_1 + c_2 * \\tanh( \\frac{X-c_3}{c_4} )

    :param XY:
    :param options:
    'check' : makes a plot with the fit and the fitted data
    :return:
    '''
    check = kwargs.get('check', False)
    # XY = XY[:-5]
    X = XY[:, 0]
    Y = XY[:, 1]
    Y /= max(Y)

    ymax = max(Y)
    ymin = min(Y)
    C = (ymax - ymin) / 2

    max_idx = np.argmin(abs(Y - 0.5 * max(Y)))
    idx75 = np.argmin(abs(Y - 0.75 * max(Y)))

    y0 = Y[max_idx] + 0.1
    xc = X[max_idx] + 10
    x75 = X[idx75] + 10
    x75 -= xc

    y0 *= 1
    p0 = [y0, C, xc, x75]  # Inital guess is a normal distribution

    errfunc = lambda p, X, Y: Fit_functions.func_tanh(X, p) - Y  # Distance to the target function
    p1, success = optimize.leastsq(errfunc, p0[:], args=(X, Y))

    X_fit = arange(min(X), max(X), 0.01)
    Y_fit = Fit_functions.func_tanh(X_fit, p1)

    if check:
        print 'initial guess'
        print p0[:4]
        #print p0[4:8]
        print 'fitted guess'
        print p1[:4]
        #print p1[4:8]


        plt.plot(X, Y, 'b')
        plt.plot(xc, y0 / 2, 'o')
        plt.plot(X_fit, Y_fit, 'r')
        plt.show()

    out = np.array([[X_fit[i], Y_fit[i]] for i in range(len(Y_fit))])

    return out


def linear(XY):
    '''
    output: out: a list np.array with x/y values of the fitted curve for plotting
            r_value: r-value of fit
            slope: slope of fit
            intercept: intercept of fit
    '''
    XY = np.array(XY)
    X = XY[:, 0]
    Y = XY[:, 1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
    X_fit = np.linspace(X[0], X[-1], 100)
    Y_fit = intercept + slope * X_fit

    out = [[X_fit[i], Y_fit[i]] for i in range(len(Y_fit))]
    out = np.array(out)
    return out, r_value, slope, intercept, std_err


def curvature(XY, **kwargs):
    '''
    This function accepts a numpy.list with list[:,0] = X values and
    list[:,1] = Y-values. A tanh function is then fitted to the data.
    The functions initial guess is calculated as follows.

    .. math::

        f(x) = c_1 + c_2 * \\tanh( \\frac{X-c_3}{c_4} )

    :param XY:
    :param options:
    'check' : makes a plot with the fit and the fitted data
    :return:
    '''

    check = kwargs.get('check', False)

    X = XY[:, 0]

    Y = XY[:, 1]

    X /= max(Y)
    Y /= max(Y)

    const = 0.

    slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)

    p0 = [intercept, slope, const]  # Inital guess is a normal distribution

    errfunc = lambda p, X, Y: Fit_functions.curvature(X, p) - Y  # Distance to the target function
    p1, success = optimize.leastsq(errfunc, p0[:], args=(X, Y))

    if check:
        #print 'initial guess'
        #print p0
        #print 'fitted guess'
        #print p1
        #
        #print 'intensity: %.2f, fitted %.2f' %(-p0[1]*35, -p1[1]*35)
        #print 'fitted %.2f' %(( -p0[1] + p1[2] * p1[1]) * 35)
        X_fit = arange(min(X), max(X), (max(X) - min(X)) / 1000)

        Y_fit = Fit_functions.curvature(X_fit, p1)

        plt.plot(X, Y, 'ob')

        plt.plot([1, 0], [0, 1], '-b')

        plt.plot(X_fit, Y_fit, 'g')
        plt.plot(X_fit, p0[0] + p0[1] * X_fit, '--r')

        #plt.plot(X_fit, p1[0] + p1[1]*X_fit, '--g')
        #plt.plot(X_fit, p1[0] + p1[2]*X_fit**2, '--g')

        #plt.plot([0,1], [0,0], '-k')
        plt.show()

    out = [[X_fit[i], Y_fit[i]] for i in range(len(Y_fit))]
    out = np.array(out)
    return out


def gauss(XY, **kwargs):
    check = kwargs.get('check')
    X = XY[:, 0]
    Y = XY[:, 1]

    YMax = max(np.fabs(Y))
    max_idx = np.argmin(Y - YMax)
    Y /= Y.sum()

    # Fit a guassian
    p0 = [X[max_idx], 18, 0]  # Inital guess is a normal distribution

    errfunc = lambda p, x, y: skewed_gauss(x, p) - y  # Distance to the target function
    p1, success = optimize.leastsq(errfunc, p0[:], args=(X, Y))

    X_fit = np.linspace(min(X), max(X), 1000)
    Y_fit = skewed_gauss(X_fit, p1)

    if check:
        print p0
        print p1
        plt.clf()
        plt.plot(X, Y / max(Y), 'ob')
        plt.plot(X_fit, Y_fit / max(Y_fit), 'g')
        plt.show()

    out = [[X_fit[i], Y_fit[i]] for i in range(len(X_fit))]
    out = np.array(out)

    return out

def func_tanh(X, parameters): # p[0] = y0, p[1]= C, p[2]=xc
    '''

    :param X: list
    :param p: list
              p[0] = y0, p[1]= C, p[2]=xc
    Y = p[0] + p[1] * tanh( (X - p[2]) / p[3] )
    :return: list


    '''
    amp = parameters['amp'].value
    center = parameters['center'].value
    sig = parameters['sig'].value
    p2 = parameters['p2'].value

    X = np.float128(X)
    X = np.array(X)

    Y = 0.5* (abs(amp) + amp * tanh( (X - center) / sig )+ p2*X**2)


    Y = np.float64(Y)

    return Y

def Gauss_Hermite(x, parameters):
    amp = parameters['amp'].value
    center = parameters['center'].value
    sig = parameters['sig'].value
    c1 = -np.sqrt(3)
    c2 = -np.sqrt(6)
    c3 = 2 / np.sqrt(3)
    c4 = np.sqrt(6) / 3
    c5 = np.sqrt(6) / 4
    skew = parameters['skew'].value
    kurt = parameters['kurt'].value

    gauss_hermite = amp * np.exp(-.5 * ((x - center) / sig) ** 2) * (1 + skew * (
        c1 * ((x - center) / sig + c3 * ((x - center) / sig) ** 3) + kurt * (
            c5 + c2 * ((x - center) / sig) ** 2 + c4 * ((x - center) / sig) ** 4)))

    return gauss_hermite


def Gauss_Ramon(x, parameters):
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
      * np.exp(-abs(np.log((np.exp(skew * (x - center) / sig) + np.exp((x - center) / (skew * sig))) / 2)) ** kurt / 2) / ( 2 ** (1 + 1 / kurt) * sig * c1 * ( np.exp(skew * (x - center) / sig) + np.exp((x - center) / (skew * sig))))

    return out

### DISTRIBUTIONS ###

def Normal(x, parameters):
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


    return out

