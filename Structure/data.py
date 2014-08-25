__author__ = 'mike'
import numpy as np
import copy
import logging


class tensor():
    '''

    '''

    def __init__(self, tensor):
        '''

        '''
        self.tensor = tensor
        self.get_eigensystem()


    def get_eigensystem(self, norm=False):
        '''
        :Parameters:
           norm : boolean
              if True :math:`\\sum_i \\tau_i =1`


        '''
        evals, evecs = np.linalg.eigh(self.tensor)
        idx = np.argsort(-evals)
        self.evals = evals[idx]
        evecs = evecs[:, idx]
        norm_evals = [i / sum(evals) for i in evals]
        self.evecs = np.array(evecs)
        self.evals_norm = norm_evals


class data(object):
    def __init__(self, variable, measurement, time=None):
        '''
        Generic 3D / 1D data containe with rudimentary functions concerning paleomagnetic data
        '''
        self.log = logging.getLogger('RockPy.DATA.data%s' % str(measurement.shape))
        self.log.debug('CREATING data structure: dimension << %s >>' % str(measurement.shape))
        self.variable = variable
        self.measurement = measurement

        if time is None:
            time = np.empty(len(variable))

        self.time = time

    ''' data properties '''

    @property
    def x(self):
        if self.measurement.shape[0] == self.measurement.size:
            return self.measurement
        else:
            return self.measurement[:, 0]

    @property
    def y(self):
        if self.measurement.shape[0] == self.measurement.size:
            return self.measurement
        else:
            return self.measurement[:, 1]

    @property
    def z(self):
        if self.measurement.shape[0] == self.measurement.size:
            return self.measurement
        else:
            return self.measurement[:, 2]

    @property
    def m(self):
        if self.measurement.shape[0] == self.measurement.size:
            return self.measurement
        else:
            out = np.array(map(np.linalg.norm, self.measurement))
            return out

    @property
    def d(self):
        if self.measurement.shape[0] == self.measurement.size:
            self.log.error('DATA is only 1d, << dec >> can not be calculated')
        else:
            out = np.arctan2(self.y, self.x)
            out = np.degrees(out)

            for i in range(len(out)):
                if out[i] < 0:
                    out[i] += 360
                if out[i] > 360:
                    out[i] -= 360
            return out

    @property
    def i(self):
        if self.measurement.shape[0] == self.measurement.size:
            self.log.error('DATA is only 1d, << inc >> can not be calculated')
        else:
            out = np.arcsin(self.z / self.m)
            out = np.degrees(out)
            return out

    @property
    def len(self):
        return len(self.measurement)

    # ## functions

    def max(self, component='m'):
        out = getattr(self, component, None)
        if out:
            out = np.fabs(out)
            out = np.max(out)
            return out
        else:
            self.log.error('COMPONENT << %s >> not found' % component)

    def min(self, component='m'):
        out = getattr(self, component, None)
        if out:
            out = np.fabs(out)
            out = np.min(out)
            return out
        else:
            self.log.error('COMPONENT << %s >> not found' % component)

    def diff(self, strength=1):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        aux = [[self_copy.variable[i],
             np.array((self_copy.measurement[i - strength] - self_copy.measurement[i + strength]) / (
             self_copy.variable[i - strength] - self_copy.variable[i + strength]))]
            for i in range(strength, len(self_copy.variable) - strength)]

        self_copy.variable = np.array([ i[0] for i in aux])
        self_copy.measurement = np.array([ i[1] for i in aux])

        return self_copy

    ### calculations
    def __sub__(self, other):  #
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))
        self_copy.measurement = np.array(
            [self_copy.measurement[i] - other.measurement[i] for i in range(len(self_copy.measurement))
             if self_copy.variable[i] == other.variable[i]])
        return self_copy

    def __add__(self, other):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        self_copy.measurement += other.measurement
        return self_copy

