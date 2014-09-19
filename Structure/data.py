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
    def __init__(self, variable, measurement, std=None, time=None, ):
        '''
        Generic 3D / 1D data containe with rudimentary functions concerning paleomagnetic data
        '''
        self.log = logging.getLogger('RockPy.DATA.data%s' % str(measurement.shape))
        self.variable = np.array(variable)
        # print np.c_[measurement, np.array([np.linalg.norm(i) for i in measurement])]

        if len(measurement.T) >= 4:
            # check for x,y,z,m data
            self.log.info('DIMENSION of << m >> is larger than 3(x,y,z) -> reducing to 3')
            measurement = measurement[:, :3]

        self.log.debug('CREATING data structure: dimension << %s >>' % str(measurement.shape))
        self.measurement = np.c_[measurement, np.array([np.linalg.norm(i) for i in measurement])]

        self.std_dev = std

        if time is None:
            try:
                time = np.empty(len(self.variable))
            except:
                time = None

        self.time = time

    def __repr__(self):
        return '<Structure.data.data object len,dim: (%i,%i)> ' % (self.len, self.dim)

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
            return self.measurement[:, 3]

    def recalc_m(self):
        self.log.info('RECALCULATING << m >> data')
        self.measurement[:, 3] = np.array([np.linalg.norm(i) for i in self.measurement[:, :3]])

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

    @property
    def dim(self):
        return len(self.measurement.T)

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

    def derivative(self, strength=1):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        aux = [[self_copy.variable[i],
                np.array((self_copy.measurement[i - strength] - self_copy.measurement[i + strength]) / (
                    self_copy.variable[i - strength] - self_copy.variable[i + strength]))]
               for i in range(strength, len(self_copy.variable) - strength)]

        self_copy.variable = np.array([i[0] for i in aux])
        self_copy.measurement = np.array([i[1] for i in aux])
        self_copy.recalc_m()
        return self_copy


    # ## calculations
    def __sub__(self, other):  #
        '''
        vector subtraction
        '''
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))
        # self_copy = copy.deepcopy(self)

        if isinstance(other, data):
            self_copy.measurement = np.array(
                [self_copy.measurement[i] - other.measurement[i] for i in range(len(self_copy.measurement))
                 if self_copy.variable[i] == other.variable[i]
                ])
            self_copy.recalc_m()

        if isinstance(other, np.ndarray):
            try:
                self_copy.measurement -= other
            except:
                print 'nope'
        return self_copy

    def __add__(self, other):

        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))
        if isinstance(other, data):
            self_copy.measurement += other.measurement
            self_copy.recalc_m()

        if isinstance(other, np.ndarray) or isinstance(other, float) or isinstance(other, int):
            self_copy.measurement += other

        return self_copy

    def __mul__(self, other):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        if isinstance(other, data):
            if other.len != 1:
                self_copy.measurement = np.array(
                    [self_copy.measurement[i] * other.measurement[i] for i in range(len(self_copy.measurement))
                     if self_copy.variable[i] == other.variable[i]])
                self.recalc_m()
            else:
                self_copy.measurement *= other.measurement
                print self_copy.m
        if isinstance(other, np.ndarray) or isinstance(other, float) or isinstance(other, int):
            self_copy.measurement *= other

        return self_copy

    def __div__(self, other):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        if isinstance(other, data):
            if other.len != 1:
                self_copy.measurement = np.array(
                    [self_copy.measurement[i] / other.measurement[i] for i in range(len(self_copy.measurement))
                     if self_copy.variable[i] == other.variable[i]])
            else:
                self_copy.measurement /= other.measurement
                print self_copy.m
        if isinstance(other, np.ndarray) or isinstance(other, float) or isinstance(other, int):
            self_copy.measurement /= other

        return self_copy

    def equal_var(self, other):
        '''
        returns data object that has only the same variables as the other data_obj
        '''

        #todo interpolation of data
        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))

        idx = [i for i in range(len(self_copy.variable)) if
               self_copy.variable[i] not in other.variable]  # indices of self, not in common with other

        if idx:
            self.log.info('FOUND different variables deleting << %i >> non-equal' % len(idx))
            for i in idx[::-1]:
                self_copy.variable = np.delete(self_copy.variable, i)
                self_copy.measurement = np.delete(self_copy.measurement, i, axis=0)
                self_copy.time = np.delete(self_copy.time, i)

        return self_copy


    def retrieve(self, idx):
        return self.variable[idx], self.measurement[idx], self.time[idx]

    def append(self, idx):
        pass

    def slice_range(self, low_lim=None, up_lim=None):
        '''
        Generates a copy of data with data within the specified range of variables (inclusive)

        :param low_lim:
        :param up_lim:
        :rtype : data_obj
        '''

        self_copy = self.__class__(copy.deepcopy(self.variable),
                                   copy.deepcopy(self.measurement),
                                   copy.deepcopy(self.time))
        if low_lim is None:
            low_lim = min(self_copy.variable)
        if up_lim is None:
            up_lim = max(self_copy.variable)

        idx = [i for i in range(len(self_copy.variable)) if self_copy.variable[i] < low_lim if
               self_copy.variable[i] > up_lim]
        if idx:
            self.log.debug('FOUND variables not within specified range deleting << %i >> non-equal' % len(idx))
            for i in idx:
                self_copy.variable = np.delete(self_copy.variable, i)
                self_copy.measurement = np.delete(self_copy.measurement, i)
                self_copy.time = np.delete(self_copy.time, i)

        return self_copy

    @property
    def self_copy(self):
        self_copy = self.__class__(copy.deepcopy(self.variable),
                           copy.deepcopy(self.measurement[:,:3]),
                           copy.deepcopy(self.time))
        return self_copy

    def normalize(self, value = 1):
        out = self / value
        out.recalc_m()
        return out

    def normalize_initial_m(self):
        norm_factor = abs(self.m[0])
        out = self / norm_factor
        out.recalc_m()
        return out
