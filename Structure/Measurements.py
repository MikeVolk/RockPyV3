# coding=utf-8
__author__ = 'Mike'
from RockPyV3.Functions import general
from RockPyV3.ReadIn import machines, helper
from RockPyV3.fit import fitting, distributions
import treatments
import logging
import numpy as np
import copy
import scipy
from scipy.interpolate import interp1d, splrep
from scipy.interpolate import UnivariateSpline
from scipy import stats, interpolate
import matplotlib.pyplot as plt
from pprint import pprint
from RockPyV3.Plots.general import Hys_Fabian2003
class Measurement(object):
    general.create_logger(__name__)

    def __init__(self, sample, mtype, mfile, machine, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.MEASUREMENT')
        implemented = {
            'af-demag': '',
            'parm-spectra': '',
            'hys': '',
            'palint': '',
            'zfc': '',
            'fc': '',
            'irm':'',
            'coe':''
        }

        if mtype.lower() in implemented:
            self.log.debug('FOUND\t measurement type: << %s >>' % mtype.lower())
            self.mtype = mtype.lower()
            self.machine = machine.lower()
            self.mfile = mfile
            self.sample = sample
            self.import_data()
        else:
            self.log.error('UNKNOWN\t measurement type: << %s >>' % mtype)
            return None

        self.treatment = None

    def import_data(self):
        implemented = {'sushibar': {'af-demag': machines.sushibar,
                                    'parm-spectra': machines.sushibar},
                       'vsm': {'hys': machines.vsm,
                               'irm': machines.vsm,
                               'coe': machines.vsm},
                       'cryo_nl': {'palint': machines.cryo_nl},
                       'mpms': {'zfc': machines.mpms,
                                'fc': machines.mpms, }}
        self.log.info(' IMPORTING << %s , %s >> data' % (self.machine, self.mtype))
        if self.machine in implemented:
            if self.mtype in implemented[self.machine]:
                self.raw_data = implemented[self.machine][self.mtype](self.mfile, self.sample.name)
                if self.raw_data == None:
                    self.log.error('IMPORTING\t did not transfer data - CHECK sample name and data file')
                    return
            else:
                self.log.error('IMPORTING UNKNOWN\t measurement type << %s >>' % (self.mtype))
        else:
            self.log.error('UNKNOWN\t machine << %s >>' % (self.machine))

    def add_treatment(self, ttype, options=None):
        self.ttype = ttype.lower()
        self.log.info('ADDING\t treatment to measurement << %s >>' % (self.ttype))

        implemented = {
            'pressure': treatments.Pressure,
        }

        if ttype.lower() in implemented:
            self.treatment = implemented[ttype.lower()](self.sample, self, options)
        else:
            self.log.error('UNKNOWN\t treatment type << %s >> is not know or not implemented' % ttype)

    def interpolate(self, what, x_new=None):
        self.log.info('INTERPOLATIN << %s >> using interp1d (Scipy)' % (what))
        xy = np.sort(getattr(self, what), axis=0)
        mtype = getattr(self, 'mtype')

        if mtype == 'coe':
            xy = np.array([[-i[0], i[1]] for i in xy])
            xy = xy[::-1]
        x = xy[:, 0]
        y = xy[:, 1]

        if x_new == None:
            x_new = np.linspace(min(x), max(x), 50000)

        # fc = interp1d(x, y, kind='cubic')
        # y_new = fc(x_new)

        fc = splrep(x, y, s=0)
        y_new = interpolate.splev(x_new, fc, der=0)

        out = np.array([[x_new[i], y_new[i]] for i in range(len(x_new))])

        if mtype == 'coe':
            out = np.array([[-i[0], i[1]] for i in out])
        self.log.debug('INTERPOLATION READY')
        return out

class Af_Demag(Measurement):
    general.create_logger('RockPy.MEASUREMENT.af-demag')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.af-demag')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)
        self.mag_method = mag_method
        try:
            self.raw_data.pop('sample')
            self.__dict__.update(self.raw_data)
            self.__dict__['fields'] = self.__dict__.pop('par1')
        except AttributeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

        except TypeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None


class pARM_spectra(Measurement):
    general.create_logger('RockPy.MEASUREMENT.PARM-SPECTRA')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.PARM-SPECTRA')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)
        try:
            self.raw_data.pop('sample')
            self.__dict__.update(self.raw_data)
            self.__dict__['ac_field'] = self.__dict__.pop('par1')
            self.__dict__['dc_field'] = self.__dict__.pop('par2')
            self.__dict__['u_window_limit'] = self.__dict__.pop('par3')
            self.__dict__['l_window_limit'] = self.__dict__.pop('par4')
            self.__dict__['window_size'] = self.__dict__.pop('par5')
            self.windows = np.array([[self.__dict__['l_window_limit'][i], self.__dict__['u_window_limit'][i]] for i in
                                     range(len(self.u_window_limit))])
        except AttributeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

        except TypeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

    def subtract_af3(self):
        self.log.debug('CREATING\t COPY of %s' % self)
        self.log.info('SUBTRACTING\t AF3 data of pARM measurement')
        self_copy = copy.copy(self)
        self_copy.x -= self.x[0]
        self_copy.y -= self.y[0]
        self_copy.z -= self.z[0]

        self_copy.m = np.array([np.sqrt(self.x[i] ** 2 + self.y[i] ** 2 + self.z[i] ** 2) for i in range(len(self.x))])
        self.log.info('RETURNING\t COPY with subtracted AF3')
        return self_copy


class Coe(Measurement):
    general.create_logger('RockPy.MEASUREMENT.COE')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.COE')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)
        self.info = self.raw_data[-1]
        self.measurement_settings = self.info['SCRIPT']


        if self.measurement_settings['Include IRM?'] == 'Yes':
            self.fields = self.raw_data[1][1][:,0]
            self.rem = self.raw_data[1][1][:,1]
            if self.measurement_settings['Include direct moment?']== 'Yes':
                self.dmom = self.raw_data[1][1][:,2]
                self.direct_moment = np.column_stack((self.fields,self.dmom))

        else:
            self.fields = self.raw_data[1][0][:,0]
            self.rem = self.raw_data[1][0][:,1]
            if self.measurement_settings['Include direct moment?']== 'Yes':
                self.dmom = self.raw_data[1][0][:,2]
                self.direct_moment = np.column_stack((self.fields,self.dmom))

        self.remanence = np.column_stack((self.fields,self.rem))



        self.remanence_interpolated = self.interpolate('remanence')

        if self.measurement_settings['Include direct moment?']== 'Yes':
            self.direct_moment_interpolated = self.interpolate('remanence')
        else:
            self.direct_moment_interpolated = []

        ''' Bcr calculation '''
        bcri_idx = np.argmin(abs(self.remanence_interpolated[:,1]))
        self.bcri = self.remanence_interpolated[bcri_idx][0] # in T
        # print self.remanence_interpolated

class Irm(Measurement):
    general.create_logger('RockPy.MEASUREMENT.IRM')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.IRM')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)

        self.fields = self.raw_data[1][0][:,0]
        self.rem = self.raw_data[1][0][:,1]
        self.dmom = self.raw_data[1][0][:,2]

        self.remanence = np.column_stack((self.fields,self.rem))
        self.direct_moment = np.column_stack((self.fields,self.dmom))


class Hysteresis(Measurement):
    general.create_logger('RockPy.MEASUREMENT.HYSTERESIS')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.HYSTERESIS')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)
        self.info = self.raw_data[2]
        self.calibration_factor = float(self.info['SETTINGS']['Calibration factor'])

        adj = self.info.get('adjusted', False)

        if adj:
            self.fields = np.array([[j[2] for j in i] for i in self.raw_data[1]])
            self.moments = np.array([[j[3] for j in i] for i in self.raw_data[1]])
        else:
            self.fields = np.array([[j[0] for j in i] for i in self.raw_data[1]])
            self.moments = np.array([[j[1] for j in i] for i in self.raw_data[1]])

        # ## initialize
        self.virgin = None
        self.msi = None
        self.up_field = []
        self.down_field = []

        for i in range(len(self.fields)):
            if self.fields[i][0] > self.fields[i][-1]:
                self.down_field = np.column_stack((self.fields[i], self.moments[i]))
                self.log.debug('FOUND\t << down_field branch >> saved as measurement.down_field [field,moment]')
            else:
                if abs(self.fields[i][0]) < self.fields[i][len(self.fields[i]) / 2]:
                    self.log.debug('FOUND\t << virgin branch >> saved as measurement.virgin [field,moment]')
                    self.virgin = np.column_stack((self.fields[i], self.moments[i]))
                else:
                    self.log.debug('FOUND\t << up_field branch >> saved as measurement.up_field [field,moment]')
                    self.up_field = np.column_stack((self.fields[i], self.moments[i]))[::-1]

        self.log.debug('CALCULATING\t << irreversible >> saved as measurement.irrev [field,moment]')

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
        ### check for msi branch ###
        ### msi is found, if the first virgin value is 0.9*mrs ###

        if self.virgin != None:

            self.mrs_msi = self.virgin[0,1] / self.mrs[0]

            if self.virgin[0,1] >= 0.9 * self.mrs[0]: #todo change back to 0.9
                self.log.debug('FOUND\t << Msi branch >> saved as measurement.msi [field,moment]')
                self.msi = self.virgin
            else:
                self.log.error('NO\t << Msi branch >> FOUND\t virgin(0)/MRS\t value: %.2f' %(self.mrs_msi))

            self.log.debug('MSI/MRS\t value: %.2f' %(self.mrs_msi))

    def E_delta_t(self):
        '''
        Method calculates the :math:`E^{\Delta}_t` value for the hysteresis.
        It uses scipy.integrate.simps for calculation of the area under the downfield branch for positive fields
         and later subtracts the area under the Msi curve.
        '''
        if self.msi != None:
            df_positive = np.array([i for i in self.down_field if i[0]>0])[::-1]
            df_energy =  scipy.integrate.simps(df_positive[:,1], x=df_positive[:,0])
            virgin_energy = scipy.integrate.simps(self.msi[:,1], x=self.msi[:,0])
            out = 2 * (df_energy-virgin_energy)
            return out

        else:
            self.log.error('UNABLE\t to find Msi branch')
            return 0.0

    def E_hys(self, check=False):
        '''
        '''
        # shifting hysteresis by Ms
        df = np.array([[i[0], i[1]-min(self.down_field[:,1])] for i in self.down_field])[::-1]
        uf = np.array([[i[0], i[1]-min(self.up_field[:,1])] for i in self.up_field])
        df_energy =  scipy.integrate.simps(df[:,1], x=df[:,0])
        uf_energy =  scipy.integrate.simps(uf[:,1], x=uf[:,0])

        out = df_energy - uf_energy

        if check:
            plt.plot(df[:,0], df[:,1])
            plt.plot(uf[:,0], uf[:,1])
            plt.plot([-1,1], [0,0])
            plt.show()

        return out

    def calculate_mrs(self, check = None):
        '''
        '''
        ''' interpoalted '''
        dfi_idx = np.argmin(self.df_interp[:,0]**2)
        ufi_idx = np.argmin(self.uf_interp[:,0]**2)
        mrs_dfi = self.df_interp[dfi_idx,1]
        mrs_ufi = self.uf_interp[ufi_idx,1]
        mrs_i = (abs(mrs_dfi)+abs(mrs_ufi) ) /2
        self.log.info('CALCULATING\t Mrs from interpolated data: df: %.1e, uf %.1e, mean:%.2e' %(mrs_dfi, mrs_ufi, mrs_i))

        ''' non interpolated '''
        df_idx = np.argmin(self.down_field[:,0]**2)
        uf_idx = np.argmin(self.up_field[:,0]**2)
        mrs_df = self.down_field[df_idx,1]
        mrs_uf = self.up_field[uf_idx,1]
        mrs = (abs(mrs_df)+abs(mrs_uf) ) /2
        self.log.info('CALCULATING\t Mrs from NON interpolated data: df: %.1e, uf %.1e, mean: %.2e' %(mrs_df, mrs_uf, mrs))


        out = [mrs_i, mrs_dfi, mrs_ufi, mrs, mrs_df, mrs_uf]
        self.log.info('RETURNING\t mrs_i, mrs_dfi, mrs_ufi, mrs, mrs_df, mrs_uf')

        if check:
            plt.axhline(0, color='black')
            plt.axvline(0, color='black')
            plt.plot(self.down_field[:,0], self.down_field[:,1], 'b.-')
            plt.plot(self.up_field[:,0], self.up_field[:,1], 'b.-')
            plt.plot(self.down_field[df_idx,0],self.down_field[df_idx,1], 'ro')
            plt.plot(self.up_field[uf_idx,0],self.up_field[uf_idx,1], 'ro')
            plt.plot(self.down_field_interpolated[:,0], self.down_field_interpolated[:,1], 'g--')
            plt.plot(self.up_field_interpolated[:,0], self.up_field_interpolated[:,1], 'g--')
            plt.plot(self.down_field_interpolated[dfi_idx,0],self.down_field_interpolated[dfi_idx,1], 'go')
            plt.plot(self.up_field_interpolated[ufi_idx,0],self.up_field_interpolated[ufi_idx,1], 'go')
            plt.show()

        return out

    def calculate_ms(self):
        df_max = max(self.down_field[:,1])
        uf_max = max(self.up_field[:,1])

        df_min = min(self.down_field[:,1])
        uf_min = min(self.up_field[:,1])

        ms = np.mean([abs(df_max),abs(uf_max),abs(df_min),abs(uf_min)])
        ms_sigm = np.std([abs(df_max),abs(uf_max),abs(df_min),abs(uf_min)])

        self.log.info('CALCULATING\t Ms: df:(%.2e,%.2e), uf:(%.2e,%.2e), mean: %.2e'%(df_max, df_min, uf_max, uf_min, ms))
        if df_max/df_min >= 0.1 or df_min/df_max >= 0.1:
            print 'test'
        out = [ms, ms_sigm]
        return out

    def approach_sat(self, branch='up_field', step=0.01, check=False, corrected=False):
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
        aux2 = np.array([i for i in out2 if i[0] > -1 if i[0] <= 0. if not np.isnan(i[1])])
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
        params = fitting.normal_skewed(self.irr[:, 0], self.irr[:, 1])
        self.field_fit = np.linspace(min(self.irr[:, 0]), max(self.irr[:, 0]), 500)
        self.irrev_fit = distributions.normal_skew(self.field_fit, parameters=params, check=False)
        self.irrev_fit /= max(self.irrev_fit)

    def correction_holder(self):
        # todo
        pass

    def plot_Fabian2003(self, norm = 'mass', out='show', folder=None):
        Hys_Fabian2003(self.sample, norm=norm, log=None, value=None, out=out, measurement = self).show(out=out, folder=folder)

    def plot(self, norm='mass', out='rtn', virg=False):

        factor = {'mass': self.sample.mass()}
        nrm = factor[norm]
        plot, = plt.plot(self.up_field[:, 0], self.up_field[:, 1] / nrm, label=self.sample.name)
        plt.plot(self.down_field[:, 0], self.down_field[:, 1] / nrm, color=plot.get_color())

        if self.virgin != None and virg:
            virg = plt.plot(self.virgin[:, 0], self.virgin[:, 1] / nrm, '.-', color=plot.get_color())
        if out == 'rtn':
            plt.xlim([-2, 2])
            return plot
        if out == 'show':
            plt.axhline(0, color='black')
            plt.axvline(0, color='black')
            plt.show()

    def Brh(self):
        MRS = self.mrs[0]
        idx = np.argmin(abs(self.up_field_interpolated[:,1] - MRS))
        out = self.up_field_interpolated[idx]
        return out


class PalInt(Measurement):
    '''

    '''
    general.create_logger('RockPy.MEASUREMENT.PALEOINTENSITY')

    def __init__(self, sample, mtype, mfile, machine, mag_method=None):
        self.log = logging.getLogger('RockPy.MEASUREMENT.PALEOINTENSITY')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)

        self.correction = ''

        self.holder = self.raw_data['acryl']

        self.nrm = helper.get_type(self.raw_data, 'NRM')
        self.trm = helper.get_type(self.raw_data, 'TRM')
        # check if trm if not replace with NRM for P0 steps
        if self.trm == None:
            self.trm = self.nrm

        self.th = helper.get_type(self.raw_data, 'TH')

        if self.nrm != None:
            self.th = np.vstack((self.nrm, self.th))

        self.pt = np.vstack((self.th[0], helper.get_type(self.raw_data, 'PT')))

        self.ck = helper.get_type(self.raw_data, 'CK')
        self.ac = helper.get_type(self.raw_data, 'AC')
        self.tr = helper.get_type(self.raw_data, 'TR')

        self.ptrm = self.calculate_ptrm()
        self.sum = self.calculate_sum()
        self.difference = self.calculate_difference()

    def get_data(self, step, type):
        if step not in self.__dict__.keys():
            self.log.error('CANT FIND\t << %s >> in data' % (step))
        types = {'temps': 0,
                 'x': 1, 'y': 2, 'z': 3,
                 'm': 4}
        out = getattr(self, step)[:, types[type]]
        return out


    def calculate_ptrm(self):
        self.log.info('CALCULATING\t PTRM')
        out = np.array([[i[0],
                         j[1] - i[1],  # X
                         j[2] - i[2],  # Y
                         j[3] - i[3],  # Z
                         np.linalg.norm([j[1] - i[1], j[2] - i[2], j[3] - i[3]])
                        ]
                        for i in self.th for j in self.pt
                        if i[0] == j[0]])
        return out

    def calculate_sum(self):

        self.log.info('CALCULATING\t SUM')
        out = np.array([[i[0],
                         j[1] + i[1],  # X
                         j[2] + i[2],  # Y
                         j[3] + i[3],  # Z
                         np.linalg.norm([j[1] + i[1], j[2] + i[2], j[3] + i[3]])
                        ]
                        for i in self.th for j in self.ptrm
                        if i[0] == j[0]])
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

        for i in ['th', 'pt', 'nrm', 'trm', 'ac', 'tr', 'ck']:
            if getattr(self, i) == None:
                self.log.debug('CANT FIND\t << %s >> in data structure' % (i))
                continue
            last_step = np.array([self.th[-1] * 0.95 for n in range(len(getattr(self, i)))])
            aux = getattr(self, i)
            d[i] = helper.claculate_difference(aux, last_step)
        self.__dict__.update(d)
        self.sum = self.calculate_sum()
        self.ptrm = self.calculate_ptrm()
        self.differnence = self.calculate_difference()


class Zfc_Fc(Measurement):
    general.create_logger('RockPy.MEASUREMENT.zfc-fc')

    def __init__(self, sample, mtype, mfile, machine, mag_method='IRM'):
        self.log = logging.getLogger('RockPy.MEASUREMENT.zfc-fc')
        super(Zfc_Fc, self).__init__(sample, mtype, mfile, machine, self.log)
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
            y_label = 'M(T)/M(T=%i)' % (norm)
            if type(norm2) == float or type(norm2) == int:
                y_label += '-M(T=%i)' % (norm2)
        else:
            if type(norm2) == float or type(norm2) == int:
                y_label = 'M(T)-M(T=%i)' % (norm2)
            else:
                if norm != None:
                    y_label = y_label[norm]
                else:
                    y_label = 'Am^2'
        if direction in implemented:
            self.log.info('PLOTTING\t %s curve' % (direction))
            xy = self.get_data(direction=direction, norm=norm, norm2=norm2)

            plot = plt.plot(xy[10:, 0], xy[10:, 1],
                            label=self.sample.name + ' ' + direction + '\t %.1fT' % (self.field),
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
        self.log.info('GETTING moment of T(%i)' % (temperature))
        xy = getattr(self, direction)
        idx = np.argmin(abs(xy[:, 0] - temperature))
        out = xy[idx]
        self.log.info('FOUND moment of T(%.2f)|M=%.3e' % (out[0], out[1]))
        return out