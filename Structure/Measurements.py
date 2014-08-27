# coding=utf-8
__author__ = 'Mike'

from RockPyV3.Functions import general
from RockPyV3.ReadIn import machines, helper
from RockPyV3.fit import fitting, distributions, functions

import RockPyV3.Plots.general as RPplt
from RockPyV3.Plots import hysteresis, viscotity
from RockPyV3.Paleointensity import statistics
import treatments, data
import logging
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d, splrep
from scipy.interpolate import UnivariateSpline
from scipy import stats, interpolate
import matplotlib.pyplot as plt
import csv
import time

# data [variable, x,y,z,m, (time)] - no dataclass
class Measurement(object):
    def __init__(self, sample_obj, mtype, mfile, machine, log=None, **options):

        if not log:
            self.log = logging.getLogger('RockPy.MEASUREMENT')
        else:
            self.log = logging.getLogger(log)

        implemented = ['af-demag', 'af', 'parm-spectra',
                       'hys', 'irm', 'coe',
                       'palint', 'thellier', 'pseudo-thellier',
                       'zfc', 'fc',
                       'visc',
        ]
        self.normalization = {}

        self.raw_data = None
        ''' initial state '''
        self.is_raw_data = None
        self.initial_state = None
        self.initial_state_data = None

        if mtype.lower() in implemented:
            self.log.debug('FOUND\t measurement type: << %s >>' % mtype.lower())
            self.mtype = mtype.lower()
            self.sample_obj = sample_obj
            self.mfile = mfile

            if machine:
                self.machine = machine.lower()
            else:
                self.machine = None

            if self.machine and self.mfile:
                self.import_data()
            else:
                self.log.warning('NO machine or mfile passed -> no raw_data will be generated')
        else:
            self.log.error('UNKNOWN\t measurement type: << %s >>' % mtype)
            return None

        self.treatment = None

    def import_data(self, rtn_raw_data=None, **options):
        implemented = {'sushibar': {'af-demag': machines.sushibar,
                                    'af': machines.sushibar,
                                    'parm-spectra': machines.sushibar,
                                    'nrm': machines.sushibar,  # externally applied magnetization
                                    'trm': machines.sushibar,  # externally applied magnetization
                                    'irm': machines.sushibar,  # externally applied magnetization
                                    'arm': machines.sushibar,  #externally applied magnetization
        },
                       'vsm': {'hys': machines.vsm,
                               'irm': machines.vsm,
                               'coe': machines.vsm,
                               'visc': machines.vsm},
                       'cryo_nl': {'palint': machines.cryo_nl2},
                       'mpms': {'zfc': machines.mpms,
                                'fc': machines.mpms, },
                       'simulation': {'palint': machines.simulation,
                       },
                       'vftb': {'hys': machines.vftb,
                                'coe': machines.vftb}}

        self.log.info(' IMPORTING << %s , %s >> data' % (self.machine, self.mtype))

        machine = options.get('machine', self.machine)
        mtype = options.get('mtype', self.mtype)
        mfile = options.get('mfile', self.mfile)

        if machine in implemented:
            if mtype in implemented[machine]:
                raw_data = implemented[machine][mtype](mfile, self.sample_obj.name)
                if raw_data is None:
                    self.log.error('IMPORTING\t did not transfer data - CHECK sample name and data file')
                    return
                else:
                    if rtn_raw_data:
                        self.log.info(' RETURNING raw_data for << %s , %s >> data' % (machine, mtype))
                        return raw_data
                    else:
                        self.raw_data = raw_data
            else:
                self.log.error('IMPORTING UNKNOWN\t measurement type << %s >>' % self.mtype)
        else:
            self.log.error('UNKNOWN\t machine << %s >>' % self.machine)

    def add_initial_state(self,
                          mtype, mfile, machine,  # standard
                          **options):
        initial_state = options.get('initial_variable', 0.0)
        self.log.info(' ADDING  initial state to measurement << %s >> data' % self.mtype)
        self.is_raw_data = self.import_data(machine=machine, mfile=mfile, mtype=mtype, rtn_raw_data=True)
        components = ['x', 'y', 'z', 'm']
        self.initial_state = np.array([self.is_raw_data[i] for i in components]).T
        self.initial_state = np.c_[initial_state, self.initial_state]
        self.initial_state_data = data.data(0, self.initial_state[:,:3], self.is_raw_data['time'], self.is_raw_data['sm'])
        self.__dict__.update({mtype:self.initial_state})


    def add_treatment(self, ttype, options=None):
        self.ttype = ttype.lower()
        self.log.info('ADDING\t treatment to measurement << %s >>' % self.ttype)

        implemented = {
            'pressure': treatments.Pressure,
        }

        if ttype.lower() in implemented:
            self.treatment = implemented[ttype.lower()](self.sample_obj, self, options)
        else:
            self.log.error('UNKNOWN\t treatment type << %s >> is not know or not implemented' % ttype)

    def interpolate(self, what, x_new=None):
        self.log.info('INTERPOLATIN << %s >> using interp1d (Scipy)' % what)
        xy = np.sort(getattr(self, what), axis=0)
        mtype = getattr(self, 'mtype')

        if mtype == 'coe':
            xy = np.array([[-i[0], i[1]] for i in xy])
            xy = xy[::-1]
        x = xy[:, 0]
        y = xy[:, 1]

        if x_new is None:
            x_new = np.linspace(min(x), max(x), 5000)

        fc = splrep(x, y, s=0)
        y_new = interpolate.splev(x_new, fc, der=0)

        out = np.array([[x_new[i], y_new[i]] for i in range(len(x_new))])

        if mtype == 'coe':
            out = np.array([[-i[0], i[1]] for i in out])
        self.log.debug('INTERPOLATION READY')
        return out


    def normalize(self, dtype, norm, value=1, **options):
        implemented = {
            'mass': self.sample_obj.mass_kg,
            'value': value,
            'nrm': getattr(self, 'nrm', None),
            'trm': getattr(self, 'trm', None),
            'irm': getattr(self, 'irm', None),
            'arm': getattr(self, 'arm', None),
            'th': getattr(self, 'th', None),
            'pt': getattr(self, 'pt', None),
            'ptrm': getattr(self, 'ptrm', None),
            'initial_state': getattr(self, 'initial_state_data', None),
            'i_s': getattr(self, 'initial_state', None),
        }
        norm_factor = implemented[norm.lower()]
        print norm_factor
        data = getattr(self, dtype)
        # try:
        # self.log.info(' Normalizing datatype << %s >> by << %s [ %.3e ]>>' % (dtype, norm, norm_factor))
        out = data / norm_factor
        return out
        # except:
        #     self.log.error('CANT normalize by data-type << %s >>' % dtype)
        #     self.log.warning('RETURNING NON NORMALIZED data-type << %s >>' % dtype)
        #     return data


# data_object [af_field, x,y,z, sm, time]
class Af_Demag(Measurement):
    def __init__(self, sample_obj,
                 mtype, mfile, machine,  # standard
                 mag_method,  # measurement specific
                 **options):
        log = 'RockPy.MEASUREMENT.af-demag'
        super(Af_Demag, self).__init__(sample_obj, mtype, mfile, machine, log)

        self.mag_method = mag_method

        try:
            self.raw_data.pop('sample')
            self.__dict__.update(self.raw_data)
            self.fields = self.__dict__.pop('par1')
        except AttributeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

        except TypeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None


    @property
    def data(self):
        measurement = np.c_[self.x, self.y, self.z]
        out = data.data(variable=self.fields, measurement=measurement, std_dev=self.sm, time=self.time)
        return out

    # old data format
    # def data(self):
    # out = np.vstack((self.fields, self.x, self.y, self.z, self.m))
    #     return out.T

    def plot(self):
        RPplt.Af_Demag(self.sample_obj)


# data [af_field, x,y,z,m]
class pARM_spectra(Measurement):
    # general.create_logger('RockPy.MEASUREMENT.PARM-SPECTRA')

    def __init__(self, sample_obj,
                 mtype, mfile, machine,
                 **options):
        log = 'RockPy.MEASUREMENT.parm-spectra'
        Measurement.__init__(self, sample_obj, mtype, mfile, machine, log)
        try:
            self.raw_data.pop('sample')
            self.__dict__.update(self.raw_data)
            self.ac_field = self.__dict__.pop('par1')
            self.dc_field = self.__dict__.pop('par2')
            self.u_window_limit = self.__dict__.pop('par3')
            self.l_window_limit = self.__dict__.pop('par4')
            self.window_size = self.__dict__.pop('par5')
            self.windows = np.vstack((self.l_window_limit, self.u_window_limit  ))
            self.fields = np.mean(self.windows, axis=0)

        except AttributeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

        except TypeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

    @property
    def data(self):
        measurement = np.c_[self.x, self.y, self.z]
        out = data.data(self.fields, measurement)
        return out

    # def data(self):
    #     out = np.vstack((self.fields, self.x, self.y, self.z, self.m))
    #     return out.T

    def subtract_af3(self):
        self.log.info('SUBTRACTING\t AF3 data of pARM measurement')
        self.x -= self.x[0]
        self.y -= self.y[0]
        self.z -= self.z[0]
        self.m = np.array([np.sqrt(self.x[i] ** 2 + self.y[i] ** 2 + self.z[i] ** 2) for i in range(len(self.x))])

        keys = ['fields', 'is', 'site', 'strat_level', 'mtype', 'par6', 'ic', 'steps/rev', 'ig', 'window_size',
                'hade', 'ac_field', 'a95', 'type', 'bl diff/sample', 'run',
                'u_window_limit', 'dspin', 'dg', 'geoaz', 'dc', 'l_window_limit', 'ispin', 'cup/sample', 'ds',
                'windows', 'm', 'dipdir', 'npos', 'sm', 'time', 'y', 'x', 'z', 'dip', ' holder/sample',
                'dc_field']

        self.log.debug('DELETING\t AF3 data of pARM measurement')
        for k in keys:
            self.__dict__[k] = self.__dict__[k][1:]


#data #todo find Coe data structure
class Coe(Measurement):
    def __init__(self, sample_obj, mtype, mfile, machine, mag_method):
        log = 'RockPy.MEASUREMENT.COE'
        Measurement.__init__(self, sample_obj, mtype, mfile, machine, log)

        # todo homogenize input
        # todo write test() function

        ''' VSM '''

        if machine.lower() == 'vsm':
            self.info = self.raw_data[-1]
            self.measurement_settings = self.info['SCRIPT']

            if self.measurement_settings['Include IRM?'] == 'Yes':
                self.fields = self.raw_data[1][1][:, 0]
                self.rem = self.raw_data[1][1][:, 1]
                if self.measurement_settings['Include direct moment?'] == 'Yes':
                    self.dmom = self.raw_data[1][1][:, 2]
                    self.direct_moment = np.column_stack((self.fields, self.dmom))

            else:
                self.fields = self.raw_data[1][0][:, 0]
                self.rem = self.raw_data[1][0][:, 1]
                if self.measurement_settings['Include direct moment?'] == 'Yes':
                    self.dmom = self.raw_data[1][0][:, 2]
                    self.direct_moment = np.column_stack((self.fields, self.dmom))

        ''' VFTB '''

        if machine.lower() == 'vftb':
            self.info = {}
            self.measurement_settings = {'Include IRM?': 'No',
                                         'Include direct moment?': 'No'}
            self.fields = self.raw_data['field']
            self.rem = self.raw_data['moment']
            self.moments = np.zeros(len(self.fields))

        self.remanence = np.column_stack((self.fields, self.rem))
        self.remanence_interpolated = self.interpolate('remanence')

        if self.measurement_settings['Include direct moment?'] == 'Yes':
            self.direct_moment_interpolated = self.interpolate('remanence')
        else:
            self.direct_moment_interpolated = []

        ''' Bcr calculation '''
        # ### taking the 0 crossing in the interpolated data
        bcri_idx = np.argmin(abs(self.remanence_interpolated[:, 1]))
        self.bcri = abs(self.remanence_interpolated[bcri_idx][0])  # in T
        # ### calculation using log normal fitting
        self.bcr = self.calculate_bcr()  # in  T
        # print self.remanence_interpolated

    def calculate_bcr(self, check=False):
        '''
        Calculate using log normal gaussian distributions according to. Only one function will be fitted.

        Leonhardt, R. (2006). Analyzing rock magnetic measurements: The RockMagAnalyzer 1.0 software. Computers \& Geosciences, 32(9), 1420–1431. doi:10.1016/j.cageo.2006.01.006

        '''
        log_fields = np.log10(np.fabs(self.remanence[:, 0]))
        # tanh_fit=fitting.Tanh(log_fields, self.remanence[:,1])

        log_norm_fit = fitting.Log_Normal(log_fields, self.remanence[:, 1] - min(self.remanence[:, 1]))
        idx = np.argmin(abs(log_norm_fit.fit_y + min(self.remanence[:, 1])))
        out = 10 ** log_norm_fit.fit_x[idx]

        if check:
            plt.axhline(color='k')
            plt.axvline(color='k')
            plt.plot(10 ** log_fields * 1000, self.rem, '.-', label='data')
            # plt.plot(10**(tanh_fit.fit_x) *1000,tanh_fit.fit_y, label ='tanh fit')
            plt.plot(10 ** log_norm_fit.fit_x * 1000, log_norm_fit.fit_y + min(self.remanence[:, 1]),
                     label='log_norm fit')
            plt.xlabel('log(Field)')
            plt.ylabel('Moment')
            plt.title('calculate_bcr check')
            plt.legend(loc='best')
            plt.show()

        return out


#data #todo find Irm data structure
class Irm(Measurement):
    general.create_logger('RockPy.MEASUREMENT.IRM')

    def __init__(self, sample_obj, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.IRM')
        Measurement.__init__(self, sample_obj, mtype, mfile, machine, self.log)

        self.fields = self.raw_data[1][0][:, 0]
        self.rem = self.raw_data[1][0][:, 1]
        self.dmom = self.raw_data[1][0][:, 2]

        self.remanence = np.column_stack((self.fields, self.rem))
        self.direct_moment = np.column_stack((self.fields, self.dmom))


#data #todo find Hysteresis data structure
class Hysteresis(Measurement):
    '''
    A subclass of the measurement class. It devides the raw_data give into an **down_field** and **up_field** branch.
    It will also look or a **virgibn** branch or an :math:`M_{si}` branch.



    '''
    general.create_logger('RockPy.MEASUREMENT.HYSTERESIS')

    def __init__(self, sample_obj, mtype, mfile, machine, mag_method):
        """


        :param sample_obj:
        :param mtype:
        :param mfile:
        :param machine:
        :param mag_method:
        :rtype : hysteresis_object


        """

        self.log = logging.getLogger('RockPy.MEASUREMENT.HYSTERESIS')
        Measurement.__init__(self, sample_obj, mtype, mfile, machine, self.log)

        if machine.lower() == 'vsm':
            self.info = self.raw_data[2]
            self.calibration_factor = float(self.info['SETTINGS']['Calibration factor'])
            adj = self.info.get('adjusted', False)
            if adj:
                self.fields = np.array([[j[2] for j in i] for i in self.raw_data[1]])
                self.moments = np.array([[j[3] for j in i] for i in self.raw_data[1]])
            else:
                self.fields = np.array([[j[0] for j in i] for i in self.raw_data[1]])
                self.moments = np.array([[j[1] for j in i] for i in self.raw_data[1]])

        if machine.lower() == 'vftb':
            self.info = {}
            self.calibration_factor = 1
            self.fields = self.raw_data.get('field')
            self.moments = self.raw_data.get('moment')


        # ## initialize
        self.virgin = None
        self.msi = None
        self.up_field = []
        self.down_field = []

        if machine.lower() == 'vsm':
            for i in range(len(self.fields)):
                if self.fields[i][0] > self.fields[i][-1]:
                    self.down_field = np.column_stack((self.fields[i], self.moments[i]))
                    self.log.debug('FOUND\t << down_field branch >> stored as measurement.down_field [field,moment]')
                else:
                    if abs(self.fields[i][0]) < self.fields[i][len(self.fields[i]) / 2]:
                        self.log.debug('FOUND\t << virgin branch >> stored as measurement.virgin [field,moment]')
                        self.virgin = np.column_stack((self.fields[i], self.moments[i]))
                    else:
                        self.log.debug('FOUND\t << up_field branch >> stored as measurement.up_field [field,moment]')
                        self.up_field = np.column_stack((self.fields[i], self.moments[i]))[::-1]

        if machine.lower() == 'vftb':
            idx = np.where(np.diff(self.fields) < 0)[0]
            self.virgin = np.column_stack((self.fields[0:idx[0] + 1], self.moments[0:idx[0] + 1]))

            self.down_field = np.column_stack((self.fields[idx[0]:idx[-1] + 2], self.moments[idx[0]:idx[-1] + 2]))
            self.up_field = np.column_stack((self.fields[idx[-1] + 1:], self.moments[idx[-1] + 1:]))[::-1]

        self.log.debug('CALCULATING\t << irreversible >> stored in measurement.irrev [field, moment]')

        # indices = [np.argmin(abs(i - self.down_field[:, 0])) for i in self.up_field[:, 0]]

        self.uf_interp = self.interpolate('up_field')
        self.df_interp = self.interpolate('down_field')
        self.down_field_interpolated = self.df_interp
        self.up_field_interpolated = self.uf_interp

        self.irr = np.array(
            [[self.uf_interp[i, 0], (self.df_interp[i, 1] - self.uf_interp[i, 1]) / 2] for i in
             range(len(self.uf_interp))])

        self.rev = np.array(
            [[self.uf_interp[i, 0], (self.df_interp[i, 1] + self.uf_interp[i, 1]) / 2] for i in
             range(len(self.uf_interp))])

        self.mrs = self.calculate_mrs()
        self.ms = self.calculate_ms()


        # ## check for msi branch ###
        # ## msi is found, if the first virgin value is 0.9*mrs ###

        if self.virgin is not None:

            self.mrs_msi = self.virgin[0, 1] / self.mrs[0]

            if self.virgin[0, 1] >= 0.9 * self.mrs[0]:  # todo change back to 0.9
                self.log.debug('FOUND\t << Msi branch >> saved as measurement.msi [field,moment]')
                self.msi = self.virgin
            else:
                self.log.error('NO\t << Msi branch >> FOUND\t virgin(0)/MRS\t value: %.2f' % self.mrs_msi)

            self.log.debug('MSI/MRS\t value: %.2f' % self.mrs_msi)

        ''' CORRECTIONS '''
        self.paramag_corrected = False

    # # Energies ##

    def E_delta_t(self):
        '''
        Method calculates the :math:`E^{\Delta}_t` value for the hysteresis.
        It uses scipy.integrate.simps for calculation of the area under the down_field branch for positive fields and later subtracts the area under the Msi curve.
        '''
        if self.msi is not None:
            df_positive = np.array([i for i in self.down_field if i[0] > 0])[::-1]
            df_energy = scipy.integrate.simps(df_positive[:, 1], x=df_positive[:, 0])
            virgin_energy = scipy.integrate.simps(self.msi[:, 1], x=self.msi[:, 0])
            out = 2 * (df_energy - virgin_energy)
            return out

        else:
            self.log.error('UNABLE\t to find Msi branch')
            return 0.0


    def E_hys(self, check=False):
        '''
        '''
        # shifting hysteresis by Ms
        df = np.array([[i[0], i[1] - min(self.down_field[:, 1])] for i in self.down_field])[::-1]
        uf = np.array([[i[0], i[1] - min(self.up_field[:, 1])] for i in self.up_field])
        df_energy = scipy.integrate.simps(df[:, 1], x=df[:, 0])
        uf_energy = scipy.integrate.simps(uf[:, 1], x=uf[:, 0])

        out = df_energy - uf_energy

        if check:
            plt.plot(df[:, 0], df[:, 1])
            plt.plot(uf[:, 0], uf[:, 1])
            plt.plot([-1, 1], [0, 0])
            plt.show()

        return out


    def approach_sat(self, branch='up_field', step=0.01, check=False, corrected=False):
        '''
        This method is supposed to calculate the approach to saturation :math:`M_s` value. It does not work yet.

        :param branch: str
        :param step: float
        :param check: bool
        :param corrected: bool
        '''
        self.log.warning('THIS METHOD IS EXPERIMENTAL DO NOT USE FOR DATA ANALYSIS')

        data = getattr(self, branch)
        x_new = np.arange(min(data[:, 0]), max(data[:, 0]), step)
        MB = self.interpolate(branch, x_new)
        M = MB[:, 1]
        B = MB[:, 0]
        M_dashdash2 = [(2 * M[i] - M[i - 1] - M[i + 1]) / step ** 2 for i in range(1, len(B) - 1)]
        M_dashdash = general.differentiate(MB, diff=2, smoothing=2)
        out2 = np.vstack((np.log10(B[1:-1]), np.log10(np.fabs(M_dashdash2)))).T
        out = np.vstack((np.log10(M_dashdash[:, 0]), np.log10(np.fabs(M_dashdash[:, 1])))).T

        aux = np.array([i for i in out if i[0] > -0.7 if i[0] <= 0. if not np.isnan(i[1])])
        # aux2 = np.array([i for i in out2 if i[0] > -1 if i[0] <= 0. if not np.isnan(i[1])])
        self.log.info('CALCULATING SLOPE\t in range -1:0')
        self.log.warning('SLOPE\t could be wrong, check with check option!')

        slope, intercept, r_value, p_value, std_err = stats.linregress(aux[:, 0], aux[:, 1])
        beta = slope + 2
        alpha = np.e ** intercept / (beta ** 2 - beta)
        if alpha > 0 or beta > 0:
            self.log.error('!!!\t alpha and/or beta not below zero.')
        if check:
            # import matplotlib.pyplot as plt

            x = np.arange(min(aux[:, 0]), max(aux[:, 0]), step)
            y = intercept + x * slope
            plt.plot(out[:, 0], out[:, 1], label='diff')
            plt.plot(out2[:, 0], out2[:, 1], label='approx Fabian[2006]')
            plt.xlabel('$\log(B)$')
            plt.ylabel('$\log(|M^{\'\'}|)$')
            plt.plot(x, y)
            plt.xlim([plt.axis()[0], plt.axis()[1]])
            plt.fill_between(x, plt.axis()[2], y, color='#555555', alpha=0.1)
            plt.fill_between([plt.axis()[0], x[0], max(x)], [y[0], y[0], y[-1]], [y[-1], y[-1], y[-1]], color='#555555',
                             alpha=0.1)
            plt.legend(loc='best')
            plt.show()

        self.log.info('RETURNING:\t alpha: %.3f, beta: %.3f' % (alpha, beta))

        if corrected:
            correction = alpha * MB[:, 0] ** beta
            MB[:, 1] -= correction
            plt.plot(MB[:, 0], MB[:, 1])
            plt.show()
            out = MB

        return out, [alpha, beta, slope, intercept, std_err]


    def fit_irrev(self):
        """
        This method will fit a function to the irreversible part of the hysteresis. It does not work yet.


        """
        params = fitting.normal_skewed(self.irr[:, 0], self.irr[:, 1])
        self.field_fit = np.linspace(min(self.irr[:, 0]), max(self.irr[:, 0]), 500)
        self.irrev_fit = distributions.normal_skew(self.field_fit, parameters=params, check=False)
        self.irrev_fit /= max(self.irrev_fit)


    def correction_holder(self):
        """
        This method is supposed to calculate the approach to saturation :math:`M_s` value. It does not work yet.


        """
        pass


    # # calculations ##


    def calculate_bc(self, poly=5, check=False):
        '''
        Calulates Bc from hysteresis loop. Calculates seperately for up_field and down_field branch.
        The data gets fitted by a polynomial of degree **poly** (default = 5).

        '''
        # find where M=0 or closest to 0
        df_idex = np.argmin(abs(self.down_field[:, 1]))
        uf_idex = np.argmin(abs(self.up_field[:, 1]))

        # generate data with +- 10 points around Bc
        df_data = np.array([self.down_field[i] for i in range(df_idex - 6, df_idex + 6)])
        uf_data = np.array([self.up_field[i] for i in range(uf_idex - 6, uf_idex + 6)])[::-1]

        df_x_new = np.arange(min(df_data[:, 0]), max(df_data[:, 0]), 0.001)
        uf_x_new = np.arange(min(uf_data[:, 0]), max(uf_data[:, 0]), 0.001)

        df_poly = np.poly1d(np.polyfit(df_data[:, 0], df_data[:, 1], poly))
        uf_poly = np.poly1d(np.polyfit(uf_data[:, 0], uf_data[:, 1], poly))
        df_new = map(df_poly, df_x_new)
        uf_new = map(uf_poly, uf_x_new)

        df_bc = abs(df_poly(0))
        uf_bc = abs(uf_poly(0))

        out = [np.mean([df_bc, uf_bc]), np.std([df_bc, uf_bc]), df_bc, uf_bc]

        if check:
            plt.plot(df_data[:, 0], df_data[:, 1], '.')
            plt.plot(df_x_new, df_new)
            plt.plot(uf_data[:, 0], uf_data[:, 1], '.')
            plt.plot(uf_x_new, uf_new)
            plt.axhline(color='k')
            plt.axvline(color='k')
            plt.title('BC calculation check')
            plt.xlabel('Field')
            plt.ylabel('Moment')
            plt.show()

        self.log.info('CALCUATING\t bc for up/downfield using %i polinomial' % poly)
        self.log.info('RETURNING\t mean, std, bc_downfield, bc_upfield')

        self.bc_mean = np.mean([df_bc, uf_bc])
        self.bc_std = np.std([df_bc, uf_bc])
        self.bc_downfield = df_bc
        self.bc_upfield = uf_bc
        return out


    def calculate_brh(self):
        MRS = self.mrs[0]
        idx = np.argmin(abs(self.up_field_interpolated[:, 1] - MRS))
        out = self.up_field_interpolated[idx]
        return out

    def calculate_mrs(self, check=None):
        '''
        '''
        ''' interpoalted '''
        dfi_idx = np.argmin(self.df_interp[:, 0] ** 2)
        ufi_idx = np.argmin(self.uf_interp[:, 0] ** 2)
        mrs_dfi = self.df_interp[dfi_idx, 1]
        mrs_ufi = self.uf_interp[ufi_idx, 1]
        mrs_i = (abs(mrs_dfi) + abs(mrs_ufi) ) / 2
        self.log.info(
            'CALCULATING\t Mrs from interpolated data: df: %.1e, uf %.1e, mean:%.2e' % (mrs_dfi, mrs_ufi, mrs_i))

        ''' non interpolated '''
        df_idx = np.argmin(self.down_field[:, 0] ** 2)
        uf_idx = np.argmin(self.up_field[:, 0] ** 2)
        mrs_df = self.down_field[df_idx, 1]
        mrs_uf = self.up_field[uf_idx, 1]
        mrs = (abs(mrs_df) + abs(mrs_uf) ) / 2
        self.log.info(
            'CALCULATING\t Mrs from NON interpolated data: df: %.1e, uf %.1e, mean: %.2e' % (mrs_df, mrs_uf, mrs))

        out = [mrs_i, mrs_dfi, mrs_ufi, mrs, mrs_df, mrs_uf]
        self.log.info('RETURNING\t mrs_i, mrs_dfi, mrs_ufi, mrs, mrs_df, mrs_uf')

        if check:
            plt.axhline(0, color='black')
            plt.axvline(0, color='black')
            plt.plot(self.down_field[:, 0], self.down_field[:, 1], 'b.-')
            plt.plot(self.up_field[:, 0], self.up_field[:, 1], 'b.-')
            plt.plot(self.down_field[df_idx, 0], self.down_field[df_idx, 1], 'ro')
            plt.plot(self.up_field[uf_idx, 0], self.up_field[uf_idx, 1], 'ro')
            plt.plot(self.down_field_interpolated[:, 0], self.down_field_interpolated[:, 1], 'g--')
            plt.plot(self.up_field_interpolated[:, 0], self.up_field_interpolated[:, 1], 'g--')
            plt.plot(self.down_field_interpolated[dfi_idx, 0], self.down_field_interpolated[dfi_idx, 1], 'go')
            plt.plot(self.up_field_interpolated[ufi_idx, 0], self.up_field_interpolated[ufi_idx, 1], 'go')
            plt.show()

        return out


    def calculate_ms(self, check=None):
        df_max = max(self.down_field[:, 1])
        uf_max = max(self.up_field[:, 1])

        df_min = min(self.down_field[:, 1])
        uf_min = min(self.up_field[:, 1])

        ms = np.mean([abs(df_max), abs(uf_max), abs(df_min), abs(uf_min)])
        ms_sigm = np.std([abs(df_max), abs(uf_max), abs(df_min), abs(uf_min)])

        self.log.info(
            'CALCULATING\t Ms: df:(%.2e,%.2e), uf:(%.2e,%.2e), mean: %.2e' % (
                df_max, df_min, uf_max, uf_min, float(ms)))

        if df_max / df_min >= 0.1 or df_min / df_max >= 0.1:
            print 'test'
        out = [ms, ms_sigm]

        if check:
            plt.plot(self.down_field[:, 0], self.down_field[:, 1], label='df+')
            plt.plot(-self.down_field[:, 0], -self.down_field[:, 1], label='df-')
            plt.plot(self.up_field[:, 0], self.up_field[:, 1], label='uf+')
            plt.plot(-self.up_field[:, 0], -self.up_field[:, 1], label='uf-')
            plt.legend(loc='best')
            plt.ylim([0.95 * df_max, 1.005 * df_max])
            plt.xlim([0.7 * max(self.down_field[:, 0]), max(self.down_field[:, 0])])
            plt.ylabel('Field [T]')
            plt.xlabel('moment')
            plt.show()

        return out

    def calculate_sigma_hys(self):
        '''
        Calculation accoding to Fabian2003

        math::

           \sigma_{\text{hys}} = log ( \frac{E_{\text{hys}}}{4 M_s B_c}
        '''
        # todo

        pass


    def paramag_slope_correction(self, percentage_of_field=80.):
        p = percentage_of_field / 100
        df_start = np.array([i for i in self.down_field if i[0] >= max(self.down_field[:, 0]) * p])  # 1
        df_end = np.array([i for i in self.down_field if i[0] <= min(self.down_field[:, 0]) * p])  # 2
        uf_start = np.array([i for i in self.up_field if i[0] >= max(self.up_field[:, 0]) * p])  # 3
        uf_end = np.array([i for i in self.up_field if i[0] <= min(self.up_field[:, 0]) * p])  # 4

        slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(df_start[:, 0], df_start[:, 1])
        slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(df_end[:, 0], df_end[:, 1])
        slope_3, intercept_3, r_value_3, p_value_3, std_err_3 = stats.linregress(uf_start[:, 0], uf_start[:, 1])
        slope_4, intercept_4, r_value_4, p_value_4, std_err_4 = stats.linregress(uf_end[:, 0], uf_end[:, 1])
        self.log.info('PARAMAG-CORRECTION - CALCULATING: high field slope: << %.2e, %.2e, %.2e, %.2e >>' % (
            slope_1, slope_2, slope_3, slope_4))
        slope = np.mean([slope_1, slope_2, slope_3, slope_4])

        self.down_field[:, 1] -= self.down_field[:, 0] * slope
        self.up_field[:, 1] -= self.up_field[:, 0] * slope

        if self.virgin is not None:
            self.virgin[:, 1] -= self.virgin[:, 0] * slope
        if self.msi:
            self.msi[:, 1] -= self.msi[:, 0] * slope

        self.uf_interp = self.interpolate('up_field')
        self.df_interp = self.interpolate('down_field')
        self.down_field_interpolated = self.df_interp
        self.up_field_interpolated = self.uf_interp

        self.irr = np.array(
            [[self.uf_interp[i, 0], (self.df_interp[i, 1] - self.uf_interp[i, 1]) / 2] for i in
             range(len(self.uf_interp))])

        self.rev = np.array(
            [[self.uf_interp[i, 0], (self.df_interp[i, 1] + self.uf_interp[i, 1]) / 2] for i in
             range(len(self.uf_interp))])

        self.mrs = self.calculate_mrs()
        self.ms = self.calculate_ms()

        self.paramag_sus = slope
        self.paramag_corrected = True

    def plot_Fabian2003(self, norm='mass', out='show', folder=None, name=None):
        RPplt.Hys_Fabian2003(self.sample_obj, norm=norm, log=None, value=None, plot=out, measurement=self).show(out=out,
                                                                                                                folder=folder
        )

    def plot(self, norm='mass', out='show', virgin=False, folder=None, name='output.pdf', figure=None):
        factor = {'mass': self.sample_obj.mass(),
                  'max': max(self.down_field[:, 1])}
        norm_factor = factor[norm]

        if figure is None:
            fig = plt.figure()
        else:
            fig = figure

        ax = fig.add_subplot(111)

        hysteresis.plot_hys(hys_obj=self, ax=ax, norm_factor=norm_factor, out='rtn', folder=folder, name=name)

        if self.virgin is not None and virgin:
            hysteresis.plot_virgin(hys_obj=self, ax=ax, norm_factor=norm_factor, out='rtn', folder=folder, name=name)
        if out == 'rtn':
            return fig
        if out == 'show':
            plt.show()
        if out == 'save':
            if folder is not None:
                plt.savefig(folder + self.samples[0].name + '_' + name, dpi=300)


class Viscosity(Measurement):
    def __init__(self, sample_obj, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.visc')
        Measurement.__init__(self, sample_obj, mtype, mfile, machine, self.log)

        if machine.lower() == 'vsm':
            self.info = self.raw_data[-1]
            self.measurement_settings = self.info['SCRIPT']

            self.times = self.raw_data[1][0][1:, 0]  # Time in seconds after start
            self.moments = self.raw_data[1][0][1:, 1]  # Magnetic Moment
            self.fields = self.raw_data[1][0][1:, 2]  # Fields in Tesla

            self.log_times = np.log10(self.times)

        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = stats.linregress(self.log_times,
                                                                                                self.moments)

    def return_fit(self):
        # todo refactor
        print '%s \t %.3e \t %.3e \t %.3e' % (self.sample_obj.name, self.intercept, self.slope, self.std_err)

    def plot(self, norm='max', out='show', virgin=False, folder=None, name='output.pdf', figure=None, **kwargs):

        plt_field = kwargs.get('plt_field', False)

        factor = {'mass': self.sample_obj.mass(),
                  'max': max(self.moments)}

        self.norm = norm
        norm_factor = factor[norm]

        params = {'backend': 'ps',
                  'text.latex.preamble': [r"\usepackage{upgreek}",
                                          r"\usepackage[nice]{units}"],
                  'axes.labelsize': 12,
                  'text.fontsize': 12,
                  'legend.fontsize': 8,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10,
                  'text.usetex': False,
                  'axes.unicode_minus': True}
        plt.rcParams.update(params)

        if figure is None:
            fig = plt.figure(figsize=(11.69, 8.27))
            plt.subplots_adjust(left=0.08, bottom=None, right=0.94, top=0.87,
                                wspace=None, hspace=None)
        else:
            fig = figure

        fig.suptitle('%s IRM viscous decay' % self.sample_obj.name, fontsize=16)

        ax = fig.add_subplot(111)
        ax2 = ax.twiny()

        ax.set_ylim([min(self.moments) / norm_factor, max(self.moments) / norm_factor])
        ax.set_xlim([min(self.log_times), max(self.log_times)])

        viscotity.plot_log_visc(visc_obj=self, ax=ax, norm_factor=norm_factor, out='rtn', folder=folder, name=name)
        viscotity.plot_visc(visc_obj=self, ax=ax2, norm_factor=norm_factor, out='rtn', folder=folder, name=name,
                            plt_opt={'color': 'k', 'linestyle': '--', 'marker': ''},
        )
        viscotity.add_formula_text(self, ax=ax)
        viscotity.add_label(self, ax, log=True)
        viscotity.add_label(self, ax2, log=False)

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()

        # ## plot linear fit
        x_fit = np.linspace(min(self.log_times), max(self.log_times), 2)
        y_fit = (self.intercept + self.slope * x_fit) / norm_factor
        ax.plot(x_fit, y_fit, '--', color='#808080')

        if plt_field:
            ax3 = ax.twinx()
            ax3.plot(self.log_times, self.fields * 1000000, 'r', label='residual field - log(t)')
            ax3.set_ylabel('Field [$\\mu T $]')
            lines3, labels3 = ax3.get_legend_handles_labels()
            ax.legend(lines + lines2 + lines3, labels + labels2 + labels3, loc='best').draw_frame(False)
        else:
            ax.legend(lines + lines2, labels + labels2, loc='best').draw_frame(False)

        if out == 'rtn':
            return fig
        if out == 'show':
            plt.show()
        if out == 'save':
            if folder is not None:
                plt.savefig(folder + self.sample_obj.name + '_' + name, dpi=300)


#data #todo find Thellier data structure
class Thellier(Measurement):
    '''

    '''
    # general.create_logger('RockPy.MEASUREMENT.thellier-thellier')

    def __init__(self, sample_obj,
                 mtype, mfile, machine,
                 mag_method='',
                 lab_field=35.2,
                 **options):
        log = 'RockPy.MEASUREMENT.thellier-thellier'

        Measurement.__init__(self, sample_obj, mtype, mfile, machine, log)

        self.correction = ''
        self.components = {'temp': 0, 'x': 1, 'y': 2, 'z': 3, 'm': 4, 'time': 5}

        if mfile:
            self.holder = self.raw_data['acryl']
            self.trm = helper.get_type(self.raw_data, 'TRM')
            self.nrm = helper.get_type(self.raw_data, 'NRM')

            # get_type gives - T, x, y, z, m
            self.th = helper.get_type(self.raw_data, 'TH')
            test = data.data(self.th[:, 0], self.th[:, 1], self.th[:, 5])

            ''' STDEVS for TH, PTRM, SUM '''
            # initializing #todo reaad stdev from file
            self.th_stdev = None
            self.ptrm_stdev = None
            self.sum_stdev = None

            types = list(set(self.raw_data['type']))

            if self.nrm is not None:
                self.th = np.vstack((self.nrm, self.th))

            # check if trm if not replace with NRM for P0 steps
            if not 'NRM' in types and not 'TRM' in types:
                self.trm = np.array([self.th[0]])
                self.nrm = np.array([self.th[0]])

            if self.trm is None:
                self.trm = self.nrm

            if self.nrm is None:
                self.nrm = self.trm

            if self.trm[:, 0] == 0:
                self.trm[:, 0] = 20

            self.pt = np.vstack((self.th[0], helper.get_type(self.raw_data, 'PT')))

            self.temps = [i for i in self.th[:, 0] if i in self.pt[:, 0]]

            self.ck = helper.get_type(self.raw_data, 'CK')
            self.ac = helper.get_type(self.raw_data, 'AC')
            self.tr = helper.get_type(self.raw_data, 'TR')

            self.ptrm = self.calculate_ptrm()
            self.sum = self.calculate_sum()
            self.difference = self.calculate_difference()

        self.lab_field = lab_field

    @property
    def trm_data(self):
        out = data.data(self.trm[:,0], self.trm[:,1:4], self.trm[:,-1])
        return out
    @property
    def nrm_data(self):
        out = data.data(self.nrm[:,0], self.nrm[:,1:4], self.nrm[:,-1])
        return out
    @property
    def th_data(self):
        out = data.data(self.th[:,0], self.th[:,1:4], self.th[:,-1])
        return out
    @property
    def pt_data(self):
        out = data.data(self.pt[:,0], self.pt[:,1:4], self.pt[:,-1])
        return out
    @property
    def ptrm_data(self):
        # print self.th_data.retrieve(0)
        out = self.pt_data - self.th_data
        return out
    @property
    def sum_data(self):
        out = self.pt_data + self.th_data
        return out
    @property
    def difference_data(self):
        out = self.ptrm_data + self.th_data
        return out


    def get_data(self, step, dtype):
        if step not in self.__dict__.keys():
            self.log.error('CANT FIND\t << %s >> in data' % step)
        types = {'temps': 0,
                 'x': 1, 'y': 2, 'z': 3,
                 'm': 4}
        out = getattr(self, step)[:, types[dtype]]
        return out

    def _get_th_ptrm_data(self, t_min=20, t_max=700, **options):
        """
        Returns th, ptrm data with [x,y,z] if option 'm'=True is given, it outputs [x,y,z,m]

        :param t_min: float
        :param t_max: float
        :return: ndarray
        """
        m = options.get('m')

        if m:
            ran = 5
        else:
            ran = 4
        xy_data = np.array([[i[1:ran], j[1:ran]] for i in self.th for j in self.ptrm
                            if i[0] == j[0]
                            if t_min <= i[0] <= t_max
                            if t_min <= j[0] <= t_max])
        x = xy_data[:, 0]
        y = xy_data[:, 1]
        return x, y


    def _get_th_ptrm_stdev_data(self, t_min=20, t_max=700):
        """
        Returns th, ptrm data

        :param t_min: float
        :param t_max: float
        :return: ndarray
        """
        xy_data = np.array([[i[1:5], j[1:5]] for i in self.th_stdev for j in self.ptrm_stdev
                            if i[0] == j[0]
                            if t_min <= i[0] <= t_max
                            if t_min <= j[0] <= t_max])

        x = xy_data[:, 0]
        y = xy_data[:, 1]
        return x, y

    def _get_ck_data(self):
        '''
        Helper function, returns the preceding th steps to each ck step

        :returns: list [ck_ij, th_i, ptrm_j, th_j]
           where ck_ij = the ptrm check to the ith temperature after heating to the jth temperature
        '''
        out = []

        for ck in self.ck:
            th_j = [0, 0, 0, 0, 0]
            for th in self.th:
                if ck[-1] - th[-1] > 0:  # if the time diff >0 => ck past th step
                    if th_j[-1] < th[-1]:
                        th_j = th
                if ck[0] == th[0]:
                    th_i = th
            for ptrm in self.ptrm:
                if ptrm[0] == th_j[0]:
                    ptrm_j = ptrm
            for pt in self.pt:
                if pt[0] == th_i[0]:
                    pt_i = pt
                    # print ptrm
            d_ck = ck[1:4] - th_j[1:4]
            d_ck_m = np.linalg.norm(d_ck)
            d_ck = np.array([ck[0], d_ck[0], d_ck[1], d_ck[2], d_ck_m, ck[-1]])

            out.append([d_ck, th_i, ptrm_j, th_j])

        # for i in out:
        # print i[0][0], i[1][0], i[2][0], i[3][0]
        # print i[0][4], i[1][4], i[2][4], i[3][4]
        return out

    def _get_ac_data(self):
        '''
        Helper function, returns the preceding th steps to each ck step

        :returns: list [ck_ij, th_i, ptrm_j, th_j]
           where ck_ij = the ptrm check to the ith temperature after heating to the jth temperature
        '''
        out = []

        for ac in self.ac:
            th_j = [0, 0, 0, 0, 0]
            for th in self.th:
                if ac[-1] > th[-1]:
                    if th_j[-1] < th[-1]:
                        th_j = th
                if ac[0] == th[0]:
                    th_i = th
            for ptrm in self.ptrm:
                if ptrm[0] == th_j[0]:
                    ptrm_j = ptrm
            for pt in self.pt:
                if pt[0] == th_j[0]:
                    pt_i = pt

            d_ac = pt_i[1:4] - ac[1:4]
            d_ac_m = np.linalg.norm(d_ac)
            d_ac = np.array([ac[0], d_ac[0], d_ac[1], d_ac[2], d_ac_m, ac[-1]])

            out.append([d_ac, th_i, ptrm_j, th_j])
        # for i in out:
        # print i[0][0], i[1][0], i[2][0], i[3][0]
        # print i[0][4], i[1][4], i[2][4], i[3][4]
        return out

    def flip_ptrm(self):
        self.ptrm[:, 1:3] = - self.ptrm[:, 1:3]
        self.pt[:, 1:3] = self.th[:, 1:3] + self.ptrm[:, 1:3]

    def calculate_ptrm(self):
        self.log.info('CALCULATING\t PTRM')
        aux = np.array([[i, j]
                        for i in self.th for j in self.pt
                        if i[0] == j[0]])
        th = aux[:, 0][:, 1:4]
        pt = aux[:, 1][:, 1:4]
        ptrm = th - pt
        times = aux[:, 1][:, -1]
        temps = aux[:, 0][:, 0]
        ptrm_m = [np.linalg.norm(i) for i in ptrm]

        out = np.c_[temps, ptrm, ptrm_m, times]

        return out

    def calculate_sum(self):

        self.log.info('CALCULATING\t SUM')
        out = np.array([[i[0],
                         abs(j[1]) + abs(i[1]),  # X
                         abs(j[2]) + abs(i[2]),  # Y
                         abs(j[3]) + abs(i[3])  # Z
                        ]
                        for i in self.th for j in self.ptrm
                        if i[0] == j[0]])

        aux = out[:, 1:]
        aux = [np.linalg.norm(i) for i in aux]
        out = np.c_[out, aux]
        return out

    def calculate_difference(self):

        self.log.info('CALCULATING\t DIFFERENCE')
        out = np.array([[i[0],
                         j[1] - i[1],  # X
                         j[2] - i[2],  # Y
                         j[3] - i[3],  # Z
                         np.linalg.norm([j[1] + i[1], j[2] + i[2], j[3] + i[3]])
                        ]
                        for i in self.th for j in self.ptrm
                        if i[0] == j[0]])
        return out

    def correction_holder(self):
        holder_th = helper.get_type(self.holder, 'TH')
        self.log.info('CALCULATING\t DIFFERENCE')
        self.log.warning('DELETING HOLDER DATA not fully implemented. Some steps might be missing')
        self.th = np.array([[i[0],
                             i[1] - i[1],  # X
                             i[2] - j[2],  # Y
                             i[3] - j[3],  # Z
                             np.linalg.norm([j[1] - i[1], j[2] - i[2], j[3] - i[3]])
                            ]
                            for i in self.th for j in holder_th
                            if i[0] == j[0]])
        pass

    def correction_last_step(self):
        self.log.info('CALCULATING\t DIFFERENCE')
        d = {'th': [], 'pt': [], 'nrm': [], 'trm': [], 'ac': [], 'tr': [], 'ck': []}

        last_step = self.th[-1] * 0.95

        self.th[:, 1:4] -= last_step[1:4]
        th_m = [np.linalg.norm(i[1:4]) for i in self.th]
        self.th[:, 4] = th_m

        self.pt[:, 1:4] -= last_step[1:4]
        pt_m = [np.linalg.norm(i[1:4]) for i in self.pt]
        self.pt[:, 4] = pt_m

        self.ck[:, 1:4] -= last_step[1:4]
        ck_m = [np.linalg.norm(i[1:4]) for i in self.ck]
        self.ck[:, 4] = ck_m

        self.ac[:, 1:4] -= last_step[1:4]
        ac_m = [np.linalg.norm(i[1:4]) for i in self.ac]
        self.ac[:, 4] = ac_m

        self.tr[:, 1:4] -= last_step[1:4]
        tr_m = [np.linalg.norm(i[1:4]) for i in self.tr]
        self.tr[:, 4] = tr_m

        self.nrm[:, 1:4] -= last_step[1:4]
        nrm_m = [np.linalg.norm(i[1:4]) for i in self.nrm]
        self.nrm[:, 4] = nrm_m

        self.trm[:, 1:4] -= last_step[1:4]
        trm_m = [np.linalg.norm(i[1:4]) for i in self.trm]
        self.trm[:, 4] = trm_m

        self.sum = self.calculate_sum()
        self.ptrm = self.calculate_ptrm()
        self.differnence = self.calculate_difference()

    def optimize(self, component='m'):
        '''
        finds combination of t_min and t_max that maximizes q/gap_max factor
        '''
        idx = self.components[component] - 1
        opt = np.array(
            [[i, j, self.q(t_min=i, t_max=j)[idx] / self.gap_max(t_min=i, t_max=j)] for i in self.temps for j in
             self.temps if i < j])
        opt_max = np.nanmax(opt[:, 2])
        opt_max_idx = np.where(opt[:, 2] == opt_max)[0]
        return opt[opt_max_idx]

    def intensity(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        out = self.lab_field * abs(slopes)
        return out

    def std_error(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        out = self.lab_field * sigmas
        return out


    ''' Statistic Methods '''

    def calculate_arai_fit(self, component='m', t_min=0, t_max=700, **options):
        '''
        maybe not appropriate, because it only minimizes the y residual
        '''
        self.log.info('CALCULATING\t << %s >> arai line fit << t_min=%.1f , t_max=%.1f >>' % (component, t_min, t_max))
        idx = self.components[component]
        xy = np.array([[i[idx], j[idx]] for i in self.th for j in self.ptrm
                       if i[0] == j[0]
                       if i[0] >= t_min if i[0] <= t_max
                       if j[0] >= t_min if j[0] <= t_max])

        # self.slope, self.intercept, self.r_value, self.p_value, self.std_err = stats.linregress(xy[:, 0], xy[:, 1])
        # return self.slope, self.intercept, self.std_err, self.r_value, self.p_value

    def MAD(self, step='th'):
        step = step.lower()
        xyz = getattr(self, step)
        x = xyz[:, 1]
        y = xyz[:, 2]
        z = xyz[:, 3]
        MAD = statistics.MAD(x, y, z)
        return MAD


    def calculate_slope_data(self, t_min=20, t_max=700, comp = 'm', **options):
        self.log.info('NEW CALCULATING\t sloep fit << t_min=%.1f , t_max=%.1f >>' % (t_min, t_max))
        #getting only same th - ptrm data
        ptrm = self.ptrm_data.equal_var(self.th_data) # data object with equal var of th_data
        th = self.th_data.equal_var(self.ptrm_data) # data object with equal var of ptrm_data
        # for comp in ['x', 'y', 'z', 'm']:
        x = getattr(ptrm.slice_range(low_lim=t_min, up_lim=t_max), comp)
        y = getattr(th.slice_range(low_lim=t_min, up_lim=t_max), comp)

        x_mean = np.mean(x)
        y_mean = np.mean(y)

        x_diff = x - x_mean
        y_diff = y - y_mean

        ''' square differences '''
        x_diff_sq = x_diff ** 2
        y_diff_sq = y_diff ** 2

        ''' sum squared differences '''
        x_sum_diff_sq = np.sum(x_diff_sq)
        y_sum_diff_sq = np.sum(y_diff_sq)

        mixed = x_diff * y_diff
        mixed_sum = np.sum(x_diff * y_diff)
        ''' calculate slopes '''
        N = len(x)

        slopes = np.sqrt(y_sum_diff_sq / x_sum_diff_sq) * np.sign(mixed_sum)

        sigmas = np.sqrt((2 * y_sum_diff_sq - 2 * slopes * mixed_sum) / ( (N - 2) * x_sum_diff_sq))

        y_intercept = y_mean + abs(slopes * x_mean)
        x_intercept = - y_intercept / slopes
        return slopes, sigmas, y_intercept, x_intercept

    def calculate_slope(self, t_min=20, t_max=700, **options):
        self.log.info('CALCULATING\t arai line fit << t_min=%.1f , t_max=%.1f >>' % (t_min, t_max))

        y, x = self._get_th_ptrm_data(t_min=t_min, t_max=t_max, m=True)

        x = np.fabs(x)
        y = np.fabs(y)

        ''' calculate averages '''

        x_mean = np.mean(x, axis=0)
        y_mean = np.mean(y, axis=0)

        ''' calculate differences '''
        x_diff = x - x_mean
        y_diff = y - y_mean

        ''' square differences '''
        x_diff_sq = x_diff ** 2
        y_diff_sq = y_diff ** 2

        ''' sum squared differences '''
        x_sum_diff_sq = np.sum(x_diff_sq, axis=0)
        y_sum_diff_sq = np.sum(y_diff_sq, axis=0)

        mixed = x_diff * y_diff
        mixed_sum = np.sum(x_diff * y_diff, axis=0)
        ''' calculate slopes '''
        N = len(x)

        slopes = np.sqrt(y_sum_diff_sq / x_sum_diff_sq) * np.sign(mixed_sum)
        sigmas = np.sqrt((2 * y_sum_diff_sq - 2 * slopes * mixed_sum) / ( (N - 2) * x_sum_diff_sq))

        y_intercept = y_mean + abs(slopes * x_mean)
        x_intercept = - y_intercept / slopes

        return slopes, sigmas, y_intercept, x_intercept

    def slope(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        return -np.fabs(slopes)

    def sigma(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        return sigmas

    def y_intercept(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        return y_intercept

    def x_intercept(self, t_min=20, t_max=700, **options):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        return x_intercept

    def N(self, t_min=20, t_max=700, debug=False, **options):
        idx = 0
        xy = np.array([[i[idx], j[idx]] for i in self.th for j in self.ptrm
                       if i[0] == j[0]
                       if i[0] >= t_min if i[0] <= t_max
                       if j[0] >= t_min if j[0] <= t_max])
        out = len(xy)
        return out

    def x_dash(self, t_min=20, t_max=700, **options):
        # todo comments
        '''
        :math:`x_0 and :math:`y_0` the x and y points on the Arai plot projected on to the best-ﬁt line. These are used to
        calculate the NRM fraction and the length of the best-ﬁt line among other parameters. There are
        multiple ways of calculating :math:`x_0 and :math:`y_0`, below is one example.

        ..math:

          x_i' = \frac{1}{2} \left( x_i + \frac{y_i - Y_{int}}{b}

        '''

        y, x = self._get_th_ptrm_data(t_min=t_min, t_max=t_max)

        x_m = [np.linalg.norm(i) for i in x]
        x = np.c_[x, x_m]

        y_m = [np.linalg.norm(i) for i in y]
        y = np.c_[y, y_m]

        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)

        x_dash = 0.5 * (x + ((y - y_intercept) / slopes))
        return x_dash

    def y_dash(self, t_min=20, t_max=700, **options):
        # todo comments
        '''
        :math:`x_0 and :math:`y_0` the x and y points on the Arai plot projected on to the best-ﬁt line. These are used to
        calculate the NRM fraction and the length of the best-ﬁt line among other parameters. There are
        multiple ways of calculating :math:`x_0 and :math:`y_0`, below is one example.

        ..math:

          x_i' = \frac{1}{2} \left( x_i + \frac{y_i - Y_{int}}{b}

        '''
        y, x = self._get_th_ptrm_data(t_min=t_min, t_max=t_max)

        x_m = [np.linalg.norm(i) for i in x]
        x = np.c_[x, x_m]

        y_m = [np.linalg.norm(i) for i in y]
        y = np.c_[y, y_m]

        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)

        y_dash = 0.5 * ( y + slopes * x + y_intercept)

        return y_dash

    def delta_x_dash(self, t_min=20, t_max=700, **options):
        '''
        ∆x0 and ∆y0 are TRM and NRM lengths of the best-ﬁt line on the Arai plot, respectively (Figure 1).
        '''
        x_dash = self.x_dash(t_min=t_min, t_max=t_max)
        out = abs(np.max(x_dash, axis=0) - np.min(x_dash, axis=0))
        return out

    def delta_y_dash(self, t_min=20, t_max=700, **options):
        '''
        ∆x0 and ∆y0 are TRM and NRM lengths of the best-ﬁt line on the Arai plot, respectively (Figure 1).
        '''
        y_dash = self.y_dash(t_min=t_min, t_max=t_max)
        out = abs(np.max(y_dash, axis=0) - np.min(y_dash, axis=0))
        return out

    def scatter(self, t_min=20, t_max=700, **options):
        # todo change to new xyzm version
        '''
        The “scatter” parameter :math:´beta´: the standard error of the slope σ (assuming uncertainty in both the pTRM and NRM
        data) over the absolute value of the best-fit slope |b| (Coe et al. 1978).

        :param quantity:
        :param t_min:
        :param t_max:
        '''
        self.log.debug('CALCULATING\t scatter parameter')
        if t_min != 20 or t_max != 700:
            slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(component=component, t_min=t_min,
                                                                            t_max=t_max)
        else:
            slope, sigma, intercept = self.slope, self.sigma, self.intercept
        scatter = sigma / abs(slope)
        return scatter

    def f(self, component='m', t_min=20, t_max=700, **options):
        """
        :Parameters:

        :Returns:

        The remanence fraction, f, was defined by Coe et al. (1978) as:

        .. math::

           f =  \\frac{\\Delta y^T}{y_0}

        where :math:`\Delta y^T` is the length of the NRM/TRM segment used in the slope calculation.
        """
        self.log.debug('CALCULATING\t f parameter')
        delta_y_dash = self.delta_y_dash(t_min=t_min, t_max=t_max)
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        f = delta_y_dash / abs(y_intercept)
        return f

    def f_VDS(self, component='m', t_min=20, t_max=700, **options):
        delta_y = self.delta_y_dash(t_min=t_min, t_max=t_max)
        VDS = self.VDS(t_min=t_min, t_max=t_max)
        out = delta_y / self.VDS(t_min=t_min, t_max=t_max)
        return out

    def VDS(self, t_min=20, t_max=700, **options):
        """
        An equal area projection may be the most useful way to present demagnetization data from a specimen with several
        strongly overlapping remanence components. In order to represent the vector nature of paleomagnetic data, it is
        necessary to plot intensity information. Intensity can be plotted versus demagnetization step in an intensity
        decay curve. However, if there are several components with different directions, the intensity
        decay curve cannot be used to determine, say, the blocking temperature spectrum or mdf, because it is the vector
        sum of the two components. It is therefore advantageous to consider the decay curve of the vector difference
        sum (VDS) of Gee et al. (1993). The VDS “straightens out” the various components by summing up the vector
        differences at each demagnetization step, so the total magnetization is plotted, as opposed to the resultant.

        :return: float
        """
        NRM_max = np.linalg.norm(self.th[-1][1:4])
        NRM_sum = np.sum(self.VD(t_min=t_min, t_max=t_max))
        out = NRM_max + NRM_sum

        return out

    def VD(self, t_min=20, t_max=700, **options):
        NRM_diff = [np.linalg.norm(self.th[i + 1][1:4] - self.th[i][1:4])
                    for i in range(0, len(self.th) - 1)
                    if self.th[i][0] >= t_min
                    if self.th[i + 1][0] <= t_max]
        return NRM_diff

    def FRAC(self, t_min=20, t_max=700, **options):
        '''
        NRM fraction used for the best-fit on an Arai diagram determined entirely by vector difference sum calculation (Shaar and Tauxe, 2013).
        '''
        NRM_sum = np.sum(np.fabs(self.VD(t_min=t_min, t_max=t_max)))
        out = NRM_sum / self.VDS(t_min=t_min, t_max=t_max)
        return out

    def beta(self, t_min=20, t_max=700):
        '''
        :param component: str
           magnetic component (x,y,z,m)
        :param t_min: float
           min temperature for slope calculation
        :param t_max: float
           max temperature for slope calculation

        :math`\beta` is a measure of the relative data scatter around the best-fit line and is the ratio of the
        standard error of the slope to the absolute value of the slope (Coe et al., 1978)

        .. math:

           \beta = \frac{\sigma_b}{|b|}

        '''
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        out = sigmas / abs(slopes)

        return out

    def g(self, t_min=20, t_max=700):
        """
        Gap factor: A measure of the gap between the points in the chosen segment of the Arai plot and the least-squares
        line. ‘g’ approaches (n-2)/(n-1) (close to unity) as the points are evenly distributed.

        :Parameters:
           quantity : str [defaulf = 'M']
              useful for calculations in Z-direction
           t_min : int [default = 20]
           t_max : int [default = 700]
        :Return:
           gap : float
              returns gap-factor for chosen line segment
        """
        y_dash = self.y_dash(t_min=t_min, t_max=t_max)
        delta_y_dash = self.delta_y_dash(t_min=t_min, t_max=t_max)
        y_dash_diff = [(y_dash[i + 1] - y_dash[i]) ** 2 for i in range(len(y_dash) - 1)]
        y_sum_dash_diff_sq = np.sum(y_dash_diff, axis=0)

        out = 1 - y_sum_dash_diff_sq / delta_y_dash ** 2

        # todo check if calculation wrong do to rounding error g=0.883487222
        return out

    def q(self, t_min=20, t_max=700):
        self.log.debug('CALCULATING\t quality parameter')

        beta = self.beta(t_min=t_min, t_max=t_max)
        f = self.f(t_min=t_min, t_max=t_max)
        gap = self.g(t_min=t_min, t_max=t_max)

        quality = (f * gap) / beta
        return quality

    def gap_max(self, t_min=20, t_max=700):
        '''
        1D
        :return: float

        '''
        vd = self.VD(t_min=t_min, t_max=t_max)
        max_vd = np.max(vd)
        sum_vd = np.sum(vd)
        return max_vd / sum_vd

    def w(self, t_min=20, t_max=700):
        xy_data = np.array([[i, j] for i in self.th for j in self.ptrm
                            if i[0] == j[0]
                            if t_min <= i[0] <= t_max
                            if t_min <= j[0] <= t_max])

        y = xy_data[:, 0]
        y = y[:, 1:5]

        x = xy_data[:, 1]
        x = np.fabs(x[:, 1:5])

        ''' calculate averages '''

        x_mean = np.mean(x, axis=0)
        y_mean = np.mean(y, axis=0)

        ''' calculate differences '''
        x_diff = x - x_mean
        y_diff = y - y_mean

        ''' square differences '''
        x_diff_sq = x_diff ** 2
        y_diff_sq = y_diff ** 2

        ''' sum squared differences '''
        x_sum_diff_sq = np.sum(x_diff_sq, axis=0)
        y_sum_diff_sq = np.sum(y_diff_sq, axis=0)

        mixed_sum = np.sum(x_diff * y_diff, axis=0)

        s = 2 + 2 * mixed_sum / np.sqrt(x_sum_diff_sq * y_sum_diff_sq)
        s = np.sqrt(s)
        out = self.f(t_min=t_min, t_max=t_max) * self.g(t_min=t_min, t_max=t_max) / s  # # calculation with s
        out = self.q(t_min=t_min, t_max=t_max) / np.sqrt(self.N(t_min=t_min, t_max=t_max) - 2)
        return out

    def area(self, t_min=20, t_max=700):
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        # print slopes

        y, x = self._get_th_ptrm_data(t_min=t_min, t_max=t_max)

        x = np.c_[x, map(np.linalg.norm, x)]
        y = np.c_[y, map(np.linalg.norm, y)]

        y2 = np.fabs(y) - ( x * slopes + y_intercept)
        # y2 = y

        area = scipy.integrate.simps(y2, x=x, axis=0)
        sum_distance = np.sum(y2, axis=0)

        plt.plot(x[:, 0], y2[:, 0], label='x')
        plt.plot(x[:, 1], y2[:, 1], label='y')
        plt.plot(x[:, 2], y2[:, 2], label='z')
        plt.plot(x[:, 3], y2[:, 3], label='m')
        plt.axhline(color='k')
        plt.legend(loc='best')
        plt.show()

    def print_statistics_table(self, component='m', t_min=20, t_max=700, lab_field=35, header=True, **options):
        csv = options.get('csv')
        csv_hdr = options.get('csv_header')

        idx = self.components[component] - 1
        slopes, sigmas, y_intercept, x_intercept = self.calculate_slope(t_min=t_min, t_max=t_max)
        f = self.f(t_min=t_min, t_max=t_max)
        f_VDS = self.f_VDS(t_min=t_min, t_max=t_max)
        VDS = self.VDS(t_min=t_min, t_max=t_max)
        FRAC = self.FRAC(t_min=t_min, t_max=t_max)
        beta = self.beta(t_min=t_min, t_max=t_max)
        g = self.g(t_min=t_min, t_max=t_max)
        q = self.q(t_min=t_min, t_max=t_max)
        gap_max = self.gap_max(t_min=t_min, t_max=t_max)
        w = self.w(t_min=t_min, t_max=t_max)

        out_header = '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s' % (
            'intensity', 'stdev', 'slope', 'sigma', 'y_intercept', 'x_intercept', 'f', 'f_VDS', 'VDS', 'FRAC', 'beta',
            'g',
            'q', 'gap_max', 'w')
        csv_header = ['intensity', 'stdev', 'slope', 'sigma', 'y_intercept', 'x_intercept', 'f', 'f_VDS', 'VDS', 'FRAC',
                      'beta', 'g', 'q', 'gap_max', 'w']

        data = [-slopes[idx] * lab_field, sigmas[idx] * lab_field, slopes[idx], sigmas[idx], y_intercept[idx],
                x_intercept[idx], f[idx], f_VDS[idx], VDS, FRAC, beta[idx], g[idx], q[idx], gap_max, w[idx]]

        out = '%.2f\t %.2f\t %.2f\t %.2f\t %.3e\t %.3e\t %.2f\t %.2f\t %.3e\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f' % (
            -slopes[idx] * lab_field, sigmas[idx] * lab_field, slopes[idx], sigmas[idx], y_intercept[idx],
            x_intercept[idx], f[idx], f_VDS[idx], VDS, FRAC, beta[idx], g[idx], q[idx], gap_max, w[idx])

        if header: print out_header
        if csv_hdr:
            return csv_header
        if csv:
            return data
        else:
            return data

    ''' NON STATISTIC FUNCTIONS '''

    def __get_M(self, step='th'):
        implemented = {'th': self.th,
                       'pt': self.pt,
                       'ck': self.ck,
                       'ac': self.ac,
                       'tr': self.tr,
                       'ptrm': self.ptrm,
                       'sum': self.sum}
        OUT = [np.sqrt(i[1] ** 2 + i[2] ** 2 + i[3] ** 2) for i in implemented[step]]
        return np.array(OUT, dtype=np.ndarray)

    def __get_D(self, step='th'):
        """
        :Parameter:
           step : str [default = 'th']
                The paleomagnetic step
        :Return:

        """
        from math import degrees

        if step not in ['th', 'pt', 'ptrm', 'ac', 'ck', 'tr', 'sum']:
            print 'No such step: %s' % step
            return

        implemented = {'th': self.th,
                       'pt': self.pt,
                       'ck': self.ck,
                       'ac': self.ac,
                       'tr': self.tr,
                       'ptrm': self.ptrm,
                       'sum': self.sum}

        aux = [np.arctan2(i[2], i[1]) for i in implemented[step]]
        D = map(degrees, aux)
        D = np.array(D)

        for i in range(len(D)):
            if D[i] < 0:
                D[i] += 360
            if D[i] > 360:
                D[i] -= 360
        return D

    def __get_I(self, step='th'):
        """
        Calculates the Inclination from a given step.

        :Parameter:
           step : str [default = 'th']
                The paleomagnetic step
        :Return:
           I : inclination Data

        Inclination is calculated with,

        .. math::

           I = \\tan^{-1} \\left( \\sqrt{\\frac{z}{x^2 + y^2} } \\right)
        """
        from math import degrees

        implemented = {'th': self.th,
                       'pt': self.pt,
                       'ck': self.ck,
                       'ac': self.ac,
                       'tr': self.tr,
                       'ptrm': self.ptrm,
                       'sum': self.sum}

        aux = [np.arctan2(i[3], np.sqrt(i[2] ** 2 + i[1] ** 2)) for i in implemented[step]]
        I = map(degrees, aux)
        I = np.array(I)

        return I

    def fit_degree2(self, component='m'):
        from scipy.optimize import fmin

        idx = self.components[component] - 1
        slope = self.slope()[idx]
        # parametric function, x is the independent variable
        # and c are the parameters.
        # it's a polynomial of degree 2
        fp = lambda c, x: c[0] + slope * x + c[2] * x ** 2

        # error function to minimize
        e = lambda p, x, y: (abs((fp(p, x) - y))).sum()

        y, x = self._get_th_ptrm_data()
        y = y[:, idx]
        mx = max(y)
        x = x[:, idx]
        x /= mx
        y /= mx
        # fitting the data with fmin
        p0 = [y[0], self.slope()[idx], 0.3]  # initial parameter value
        p = fmin(e, p0, args=(x, y))

        print 'initial parameters: ', p0
        print 'estimater parameters: ', p

        xx = np.linspace(0, max(x), 100)
        plt.plot(x, y, 'bo',
                 xx, fp(p, xx), 'r',
        )
        plt.show()

    def export_tdt(self, folder=None, name='_output'):

        """
        NOT FINISHED
        :param folder:
        :param name:
        """
        if folder is None:
            from os.path import expanduser

            folder = expanduser("~") + '/Desktop/'

        ADD = {'TH': 0, 'PT': 0.1, 'CK': 0.2, 'TR': 0.3, 'AC': 0.4, 'NRM': 0}
        output = csv.writer(open(str(folder) + '%s' % self.sample_obj.name + name + '.tdt', 'wb'), delimiter='\t')

        for i in range(len(self.raw_data['time'])):
            if not self.raw_data['type'][i] in ['TRM']:
                sample_t = self.raw_data['step'][i] + ADD[self.raw_data['type'][i]]
                sample_m = self.raw_data['m'][i] / self.sample_obj.volume(length_unit='m')
                sample_d = self.raw_data['dc'][i]
                sample_i = self.raw_data['ic'][i]
                out = [self.sample_obj.name, sample_t, sample_m, sample_d, sample_i]
                output.writerow(out)

    def print_table(self, step='th', **options):
        if step in ['th', 'sum', 'ptrm', 'pt', 'ck', 'ac', 'tr']:
            for i in getattr(self, step):
                print i[0], '\t', i[1], '\t', i[2], '\t', i[3]

    def dunlop(self, component='m'):
        RPplt.Dunlop(self, component=component)


class Zfc_Fc(Measurement):
    def __init__(self, sample_obj, mtype, mfile, machine, mag_method='IRM'):
        log = 'RockPy.MEASUREMENT.zfc-fc'
        super(Zfc_Fc, self).__init__(sample_obj, mtype, mfile, machine, log)
        self.mag_method = mag_method
        self.cooling_data = {key: self.raw_data[key][0] for key in self.raw_data.keys()}
        self.warming_data = {key: self.raw_data[key][1] for key in self.raw_data.keys()}

        self.cooling = np.array([[self.cooling_data['temperature (k)'][i], self.cooling_data['long moment (emu)'][i],
                                  self.cooling_data['field (oe)'][i]] for i in
                                 range(len(self.cooling_data['temperature (k)']))])

        self.warming = np.array([[self.warming_data['temperature (k)'][i], self.warming_data['long moment (emu)'][i],
                                  self.warming_data['field (oe)'][i]] for i in
                                 range(len(self.warming_data['temperature (k)']))])[::-1]

        self.field = np.mean(self.raw_data['field (oe)']) / 10000  # conversion to T
        self.c_max = max(self.cooling[10:, 1])
        self.w_max = max(self.warming[10:, 1])
        self.c_min = min(self.cooling[10:, 1])
        self.w_min = min(self.warming[10:, 1])

    def plot(self, direction='cooling', option='show', norm=None, norm2=None):
        '''
        plots data
        '''
        implemented = ['cooling', 'warming']
        line = {'cooling': '-', 'warming': ':'}

        ''' Y-LABEL '''

        y_label = {'max': 'Magnetic Moment norm a.u.',
                   'warming': 'Magnetic Moment norm a.u.',
                   'cooling': 'Magnetic Moment norm a.u.',
        }

        if type(norm) == float or type(norm) == int:
            y_label = 'M(T)/M(T=%i)' % norm
            if type(norm2) == float or type(norm2) == int:
                y_label += '-M(T=%i)' % norm2
        else:
            if type(norm2) == float or type(norm2) == int:
                y_label = 'M(T)-M(T=%i)' % norm2
            else:
                if norm is not None:
                    y_label = y_label[norm]
                else:
                    y_label = 'Am^2'
        if direction in implemented:
            self.log.info('PLOTTING\t %s curve' % direction)
            xy = self.get_data(direction=direction, norm=norm, norm2=norm2)

            plot = plt.plot(xy[10:, 0], xy[10:, 1],
                            label=self.sample_obj.name + ' ' + direction + '\t %.1fT' % self.field,
                            linestyle=line[direction])
            if option == 'show':
                plt.ylabel(y_label)
                plt.xlabel('Temperature')
                plt.legend(loc='best')
                plt.ylim([plt.axis()[2], plt.axis()[3] * 1.1])
                plt.show()
            if option == 'return':
                return plot

    def get_derivative(self, direction, diff=1, smoothing=2, norm=False):
        """
        :param direction:
        :param diff:
        :param smoothing:
        """
        self.log.info('CALCULATING\t %i derivative with smoothing %i' % (diff, smoothing))
        xy = getattr(self, direction)
        out = general.differentiate(xy, diff=diff, smoothing=smoothing, norm=norm)
        return out

    def get_derivatives(self, diff=1, smoothing=2, norm=False):
        out = []
        for i in ['cooling', 'warming']:
            aux = self.get_derivative(direction=i, diff=diff, smoothing=smoothing, norm=norm)
            out.append(aux)
        return out


    def get_data(self, direction='cooling', norm=None, norm2=None):
        xy = getattr(self, direction)
        normalization = {'max': max(xy[:, 1]),
                         'cooling': self.c_max,
                         'warming': self.w_max,
                         'temp': None}

        if norm or norm2:
            if norm or type(norm) == float or type(norm) == int:
                if type(norm) == float or type(norm) == int:
                    norm_HT = self.get_M_T(temperature=norm, direction=direction)[1]
                    normalization['temp'] = norm
                    norm = 'temp'
                else:
                    norm_HT = normalization[norm]
            else:
                norm_HT = 1
            if norm2:
                normalization_LT = self.get_M_T(temperature=norm2, direction=direction)
                norm_LT = normalization_LT[1]
                xy[:, 1] -= norm_LT
                norm_HT -= norm_LT
            if norm in normalization or type(norm) == float or type(norm) == int:
                self.log.info('NORMALIZING to << %s %.3e >>' % (norm, normalization[norm]))
                xy[:, 1] /= norm_HT
            else:
                self.log.error('NORMALIZATION not recognized or non float')
        return xy

    def get_M_T(self, temperature, direction='cooling'):
        self.log.info('GETTING moment of T(%i)' % temperature)
        xy = getattr(self, direction)
        idx = np.argmin(abs(xy[:, 0] - temperature))
        out = xy[idx]
        self.log.info('FOUND moment of T(%.2f)|M=%.3e' % (out[0], out[1]))
        return out


class Pseudo_Thellier(Measurement):
    def __init__(self, sample_obj,
                 af_obj, parm_obj,
                 mtype=None, mfile=None, machine=None,
                 **options):
        log = 'RockPy.MEASUREMENT.pseudo-thellier'
        super(Pseudo_Thellier, self).__init__(sample_obj, mtype, mfile, machine, log=log, **options)

        parm_obj.subtract_af3()
        parm_data = parm_obj.data

        self.components = {'fields': 0, 'x': 1, 'y': 2, 'z': 3, 'm': 4}

        self.lab_field = np.mean(parm_obj.dc_field)

        self.trm = af_obj.data[0, :]
        self.nrm = af_obj.data[0, :]

        self.parm = np.column_stack((
            parm_obj.u_window_limit,
            np.array([np.sum(parm_data[:i, :], axis=0) for i in range(len(parm_data))])[:, 1:]))
        self.af = np.array([af_obj.data[i] for i in range(len(af_obj.data)) if af_obj.data[i][0] in self.parm[:, 0]])

        self.sum = np.c_[self.parm[:, 0], self.parm[:, 1:] + self.af[:, 1:]]

    @property
    def th(self):
        return self.af

    @property
    def ptrm(self):
        return self.parm

    def plot(self, norm='nrm', component='m'):
        RPplt.Pseudo_Thellier(self.sample_obj, norm=norm)

    def _get_af_arm_data(self, b_min=0, b_max=200, **options):
        """
        Returns af, arm data with [x,y,z] if option 'm'=True is given, it outputs [x,y,z,m]

        :param b_min: float
        :param b_max: float
        :return: ndarray
        """
        m = options.get('m')

        if m:
            ran = 5
        else:
            ran = 4
        xy_data = np.array([[i[1:ran], j[1:ran]] for i in self.af for j in self.parm
                            if i[0] == j[0]
                            if b_min <= i[0] <= b_max
                            if b_min <= j[0] <= b_max])
        x = xy_data[:, 0]
        y = xy_data[:, 1]
        return x, y


    def calculate_slope(self, b_min=0, t_max=200, **options):
        self.log.info('CALCULATING\t arai line fit << t_min=%.1f , t_max=%.1f >>' % (b_min, t_max))

        y, x = self._get_af_arm_data(t_min=b_min, t_max=t_max, m=True)

        x = np.fabs(x)
        y = np.fabs(y)

        ''' calculate averages '''

        x_mean = np.mean(x, axis=0)
        y_mean = np.mean(y, axis=0)

        ''' calculate differences '''
        x_diff = x - x_mean
        y_diff = y - y_mean

        ''' square differences '''
        x_diff_sq = x_diff ** 2
        y_diff_sq = y_diff ** 2

        ''' sum squared differences '''
        x_sum_diff_sq = np.sum(x_diff_sq, axis=0)
        y_sum_diff_sq = np.sum(y_diff_sq, axis=0)

        mixed = x_diff * y_diff
        mixed_sum = np.sum(x_diff * y_diff, axis=0)
        ''' calculate slopes '''
        N = len(x)

        slopes = np.sqrt(y_sum_diff_sq / x_sum_diff_sq) * np.sign(mixed_sum)

        sigmas = np.sqrt((2 * y_sum_diff_sq - 2 * slopes * mixed_sum) / ( (N - 2) * x_sum_diff_sq))

        y_intercept = y_mean + abs(slopes * x_mean)
        x_intercept = - y_intercept / slopes

        return slopes, sigmas, y_intercept, x_intercept


