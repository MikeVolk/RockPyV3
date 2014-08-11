# coding=utf-8
__author__ = 'Mike'
import logging
import numpy as np
import scipy as sp
from RockPyV3.Functions import general, convert
from RockPyV3.Plots import general as RPplt
from RockPyV3.Plots.helper import get_colors
import measurements
import matplotlib.pyplot as plt
import csv


class SampleGroup():
    general.create_logger('RockPy.SAMPLEGROUP')

    def __init__(self, sample_dict=None, sample_list=None, name='', log=None):
        if not log:
            self.log = logging.getLogger('RockPy.SAMPLEGROUP')
        else:
            self.log = log

        self.log.info('CREATING\t new sample group')

        self.samples = []
        self.sample_names = []
        self.name = name

        if sample_dict:
            if type(sample_dict) == list:
                self.log.error('NO LIST found\t can`t add samples << %s >>' % sorted(sample_dict.keys()))
                return
            self.log.info('CREATING\t adding samples << %s >>' % sorted(sample_dict.keys()))
            self.add_sample_dict(sample_dict)

        if sample_list:
            if type(sample_list) == dict:
                self.log.error('NO DICT found\t can`t add samples << %s >>' % sorted(sample_dict.keys()))
                return
            self.log.info('CREATING\t adding samples << %s >>' % sorted(sample_dict.keys()))
            self.add_sample_list(sample_list)

    def add_sample(self, sample):
        self.samples.append(sample)
        self.sample_names.append(sample.name)
        self.sample_names = sorted(self.sample_names)

    def add_sample_list(self, sample_list):
        for sample in sample_list:
            self.add_sample(sample)

    def add_sample_dict(self, sample_dict):
        for sample in sample_dict:
            self.add_sample(sample_dict[sample])

    def pint_average_type(self, step='th'):
        test = []
        n_temps = []
        for sample in self.samples:
            for measurement in sample.measurements:
                t_aux = getattr(measurement, step)[:, 0]
                for i in getattr(measurement, step):
                    norm = getattr(measurement, 'trm')
                    aux = [sample.name, str(measurement.treatment), i[0], i[1] / norm[0, 1], i[2] / norm[0, 2],
                           i[3] / norm[0, 3], i[4] / norm[0, 4]]
                    test.append(aux)
                n_temps.append(set(t_aux))

        test = np.array(test)
        temps = sorted(list(set.intersection(*n_temps)))
        samples = sorted(list(set(test[:, 0])))
        measurements = sorted(list(set(test[:, 1])))

        ave, std = [], []
        for measurement in measurements:
            measurement_ave, measurement_std = [], []
            for t in temps:
                Taux = [float(i[2]) for i in test if float(i[2]) == t if i[1] == measurement]
                Xaux = [float(i[3]) for i in test if float(i[2]) == t if i[1] == measurement]
                Yaux = [float(i[4]) for i in test if float(i[2]) == t if i[1] == measurement]
                Zaux = [float(i[5]) for i in test if float(i[2]) == t if i[1] == measurement]
                Maux = [float(i[6]) for i in test if float(i[2]) == t if i[1] == measurement]
                t_ave = np.array([np.mean(Taux), np.mean(Xaux), np.mean(Yaux), np.mean(Zaux), np.mean(Maux)])
                t_std = np.array([np.std(Taux), np.std(Xaux), np.std(Yaux), np.std(Zaux), np.std(Maux)])
                measurement_ave.append(t_ave)
                measurement_std.append(t_std)
            ave.append(measurement_ave)
            std.append(measurement_std)
        return np.array(ave), np.array(std)

    def get_all_measurements(self, mtype):
        measurements = [measurement for sample in self.samples for measurement in sample.measurements if
                        measurement.mtype == mtype]
        treatments = [measurement.treatment for sample in self.samples for measurement in sample.measurements if
                      measurement.mtype == mtype]
        out = measurements
        self.log.info('FOUND\t %i measurements with << %s >>' % (len(out), mtype))
        return out

    def get_treatment_for_mtype(self, mtype):
        treatments = [measurement.treatment.label for sample in self.samples for measurement in sample.measurements if
                      measurement.mtype == mtype]
        out = sorted(list(set(treatments)))
        return out

    def plot(self, mtype, norm='mass', value=None, rtn='show'):
        if mtype == 'hys':
            OUT = RPplt.Hys(samples=self.samples, norm=norm, value=value, rtn=rtn).show()
            return OUT


class TTGroup(SampleGroup):
    general.create_logger('RockPy.TT-SAMPLEGROUP')

    def get_data_mean(self, measurements, step='th'):
        data = self.get_data_with_temp(measurements=measurements, step=step)
        mdat = np.ma.masked_array(data,np.isnan(data)) # masking nan entries
        mm = np.mean(mdat,axis=0)
        out = mm.filled(np.nan)
        return out

    def get_data_max(self, measurements, step='th'):
        data = self.get_data_with_temp(step=step)
        out = np.nanmax(np.fabs(data), axis=0)
        return out

    def get_data_stdev(self, measurements, step='th', *options):
        '''
        return stdev of all samples. First row is not usable (stdev of temperatures)
        '''
        data = self.get_data_with_temp(measurements=measurements, step=step)
        temps = self.get_temps(measurements=measurements, step=step)
        mdat = np.ma.masked_array(data,np.isnan(data))
        mm = np.std(mdat,axis=0)
        out = mm.filled(np.nan)
        out[:, 0] = temps
        return out

    def get_fill_nan(self, measurement, temps, step):
        aux = getattr(measurement, step.lower())
        out = []
        for i in temps:
            if i in aux[:, 0]:
                idx = np.where(aux[:, 0] == i)[0]
                out.append(aux[idx][0])
            else:
                out.append(np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]))

        out = np.array(out)
        return out

    def get_data_with_temp(self, measurements, step, norm='trm', **options):
        '''
        returns only the data (normalized), where a measurement was done at a common temperature.
        Replaces T-step with nan
        '''

        temps = self.get_temps(measurements=measurements, step=step)
        data = np.array([self.get_fill_nan(measurement, temps, step)[:, 0:6] for measurement in measurements])

        norm_data = np.array([np.c_[1,
                                    np.linalg.norm(getattr(measurement, norm)[:, 1:4]),
                                    np.linalg.norm(getattr(measurement, norm)[:, 1:4]),
                                    np.linalg.norm(getattr(measurement, norm)[:, 1:4]),
                                    np.linalg.norm(getattr(measurement, norm)[:, 1:4]),
                                    1
                              ]
                              for measurement in measurements])

        data /= norm_data

        return data

    def get_temps(self, measurements, step):
        '''
        gets the temperature steps common in all measurements
        '''
        temps = np.array([i[0] for measurement in measurements for i in getattr(measurement, step)])
        temps = sorted(list(set(temps.flatten())))  # temps in all measurements

        if step in ['pt', 'th']:
            temps = np.array([i for measurement in measurements for i in (getattr(measurement, 'temps'))])
            temps = sorted(list(set(temps)))  # temps in all measurements

        return temps

    def plot_arai(self, component='m', show_std=True, out='show', plt_opt={}):
        colors = get_colors()
        components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}

        if not 'color' in plt_opt:
            plt_opt.update({'color': colors[0]})

        idx = components[component]

        th_mean = self.get_data_mean(step='th')
        ptrm_mean = self.get_data_mean(step='ptrm')

        th_max = self.get_data_max(step='th')
        ptrm_max = self.get_data_max(step='ptrm')

        th_stdev = self.get_data_stdev(step='th')
        ptrm_stdev = self.get_data_stdev(step='ptrm')

        plt.plot(ptrm_mean[:, idx], th_mean[:, idx], '.-', **plt_opt)

        if show_std:
            plt.fill_between(ptrm_mean[:, idx], th_mean[:, idx] + th_stdev[:, idx], th_mean[:, idx] - th_stdev[:, idx],
                             color=colors[0], alpha=0.1, antialiased=True)
            # plt.fill_betweenx(th_mean[:,idx]-th_stdev[:,idx], ptrm_mean[:,idx], ptrm_mean[:,idx]-ptrm_stdev[:,idx],
            # where= th_stdev[:,idx]-th_stdev[:,idx] < th_stdev[:,idx],
            #                   color=colors[1], alpha = 0.1, antialiased=True)
            # plt.fill_betweenx(th_mean[:,idx]+th_stdev[:,idx], ptrm_mean[:,idx], ptrm_mean[:,idx]+ptrm_stdev[:,idx],
            #                   where= th_stdev[:,idx]+th_stdev[:,idx] > th_stdev[:,idx],
            #                   color=colors[1], alpha = 0.1, antialiased=True)

        if out == 'show':
            plt.xlabel(r'pTRM gained')
            plt.ylabel(r'NRM remaining')
            plt.show()
        if out == 'rtn':
            return


    def get_sample_obj(self):
        sample_obj = Sample(name=self.name)
        measurements = self.get_all_measurements(mtype='palint')
        treatments = self.get_treatment_for_mtype(mtype='palint')

        for treat in treatments:  # getting measurements with the same treatments
            measurements_treat = [i for i in measurements if i.treatment.label == treat]
            TT_obj = sample_obj.add_measurement(mtype='palint', mfile='', machine='')
            TT_obj.nrm = self.get_data_mean(measurements=measurements_treat, step='nrm')
            TT_obj.trm = self.get_data_mean(measurements=measurements_treat, step='trm')

            TT_obj.th = self.get_data_mean(measurements=measurements_treat, step='th')
            TT_obj.ptrm = self.get_data_mean(measurements=measurements_treat, step='ptrm')

            TT_obj.sum = np.c_[TT_obj.th[:, 0], (TT_obj.th[:, 1:] + TT_obj.ptrm[:, 1:])/TT_obj.trm[0,4]]

            TT_obj.pt = self.get_data_mean(measurements=measurements_treat, step='pt')

            TT_obj.difference = np.c_[TT_obj.th[:, 0], TT_obj.th[:, 1:] - TT_obj.ptrm[:, 1:]]

            TT_obj.th_stdev = self.get_data_stdev(measurements=measurements_treat, step='th')
            TT_obj.ptrm_stdev = self.get_data_stdev(measurements=measurements_treat, step='ptrm')
            TT_obj.sum_stdev = np.c_[TT_obj.th[:, 0], TT_obj.th_stdev[:, 1:] + TT_obj.ptrm_stdev[:, 1:]]
            TT_obj.difference_stdev = np.c_[TT_obj.th[:, 0], TT_obj.th_stdev[:, 1:] - TT_obj.ptrm_stdev[:, 1:]]

            TT_obj.ac = self.get_data_mean(measurements=measurements_treat, step='ac')
            TT_obj.ck = self.get_data_mean(measurements=measurements_treat, step='ck')
            TT_obj.tr = self.get_data_mean(measurements=measurements_treat, step='tr')

            TT_obj.temps = self.get_temps(measurements=measurements_treat, step='th')

            TT_obj.add_treatment(ttype='pressure', options={'p_max': measurements_treat[0].treatment.p_max,
                                                            'p_seen': measurements_treat[0].treatment.p_seen})

        return sample_obj

    def get_TT_obj(self):
        TT_obj = measurements.Thellier(sample=self.name, mtype=None, mfile=None, machine=None)

        TT_obj.nrm = self.get_data_mean(step='nrm')
        TT_obj.trm = self.get_data_mean(step='trm')

        TT_obj.th = self.get_data_mean(step='th')
        TT_obj.ptrm = self.get_data_mean(step='ptrm')

        TT_obj.sum = np.c_[TT_obj.th[:, 0], TT_obj.th[:, 1:] + TT_obj.ptrm[:, 1:]]
        TT_obj.difference = np.c_[TT_obj.th[:, 0], TT_obj.th[:, 1:] - TT_obj.ptrm[:, 1:]]

        TT_obj.ac = self.get_data_mean(step='ac')
        TT_obj.ck = self.get_data_mean(step='ck')
        TT_obj.tr = self.get_data_mean(step='tr')

        TT_obj.temps = self.get_temps('th')

        return TT_obj

    def plot_dunlop(self, component='m', show_std=True, out='show', plt_opt={}):
        colors = get_colors()
        components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
        idx = components[component]

        th_mean = self.get_data_mean(step='th')
        ptrm_mean = self.get_data_mean(step='ptrm')
        sum_mean = np.c_[th_mean[:, 0], th_mean[:, 1:] + ptrm_mean[:, 1:]]

        th_max = self.get_data_max(step='th')
        ptrm_max = self.get_data_max(step='ptrm')
        sum_max = self.get_data_max(step='sum')

        th_stdev = self.get_data_stdev(step='th')
        ptrm_stdev = self.get_data_stdev(step='ptrm')
        sum_stdev = np.c_[th_stdev[:, 0], th_stdev[:, 1:] + ptrm_stdev[:, 1:]]

        plt.plot(th_mean[:, 0], th_mean[:, idx], '.-', color=colors[1], **plt_opt)
        plt.plot(ptrm_mean[:, 0], ptrm_mean[:, idx], '.-', color=colors[2], **plt_opt)
        plt.plot(sum_mean[:, 0], sum_mean[:, idx], '.-', color=colors[0], **plt_opt)

        if show_std:
            plt.fill_between(th_mean[:, 0], th_mean[:, idx] - th_stdev[:, idx], th_mean[:, idx] + th_stdev[:, idx],
                             color=colors[1], alpha=0.1)
            plt.fill_between(ptrm_mean[:, 0], ptrm_mean[:, idx] - ptrm_stdev[:, idx],
                             ptrm_mean[:, idx] + ptrm_stdev[:, idx], color=colors[2], alpha=0.1)
            plt.fill_between(sum_mean[:, 0], sum_mean[:, idx] - sum_stdev[:, idx], sum_mean[:, idx] + sum_stdev[:, idx],
                             color=colors[0], alpha=0.1)

        if out == 'show':
            plt.xlabel(r'Temperature [$^\circ C$]')
            plt.ylabel(r'Normalized Magnetic Moment')
            plt.show()
        if out == 'rtn':
            return


class Sample():
    general.create_logger('RockPy.SAMPLE')

    def __init__(self, name, mass=1.0, mass_unit=None, height=None, diameter=None, length_unit=None):
        """

        :param name: str - name of the sample.
        :param mass: float - mass of the sample. If not kg then please specify the mass_unit. It is stored in kg.
        :param mass_unit: str - has to be specified in order to calculate the sample mass properly.
        :param height: float - sample height - stored in 'm'
        :param diameter: float - sample diameter - stored in 'm'
        :param length_unit: str - if not 'm' please specify
        """
        self.name = name
        self.log = logging.getLogger('RockPy.SAMPLE')
        self.log.info('CREATING\t new sample << %s >>' % self.name)

        self.measurements = []

        if mass_unit:
            mass_factor = convert.convert2(mass_unit, 'kg', 'mass')
        else:
            self.log.info(' MISSING\t length_unit: assuming << mm >>')
            mass_factor = convert.convert2('mg', 'kg', 'mass')

        if mass:
            self.mass_kg = mass * mass_factor
            self.log.debug(' ADDING\t << mass >> input: %.1f [%s] stored: %1f [kg]' % (mass, mass_unit, self.mass_kg))
        else:
            self.log.debug('MISSING\t << mass >>')
            self.mass_kg = None

        # get length _unit for conversion
        if length_unit:
            length_factor = convert.convert2(length_unit, 'm', 'length')
        else:
            self.log.info(' MISSING\t mass_unit: assuming << mg >>')
            length_factor = convert.convert2('mm', 'm', 'length')

        if height:
            self.height_m = float(height) * length_factor
            self.log.debug(
                ' ADDING\t << height >> input: %.1f [%s] stored: %1f [m]' % (height, length_unit, self.height_m))

        else:
            self.log.debug('MISSING\t << height >>')
            self.height_m = None

        if diameter:
            self.diameter_m = float(diameter) * length_factor
            self.log.debug(
                ' ADDING\t << diameter >> input: %.1f [%s] stored: %1f [m]' % (diameter, length_unit, self.diameter_m))

        else:
            self.log.debug('MISSING\t << diameter >>')
            self.diameter_m = None

    def __repr__(self):
        return '<< %s - Structure.sample.Sample >>' %self.name
    ''' ADD FUNCTIONS '''

    def add_mass(self, mass, mass_unit='mg'):
        '''
        Adds mass to the sample object

        :param mass: float
        :param mass_unit: str

        changes sample.mass_kg from its initial value to a new value.
        mass unit should be given, if other than 'mg'. See convert2 helper function
        #todo point to helper
        '''
        self.mass_kg = float(mass) * convert.convert2(mass_unit, 'kg', 'mass')
        self.log.debug(' ADDING\t << mass >> input: %.1f [%s] stored: %1f [kg]' % (mass, mass_unit, self.mass_kg))

    def add_height(self, height, length_unit='mm'):
        '''
        Adds height in m to sample

        :param height: float
        :param length_unit: str

        changes sample.height_m from its initial value to a new value.
        Length unit should be given, if other than 'mm'. See convert2 helper function
        #todo point to helper
        '''
        self.height_m = float(height) * convert.convert2(length_unit, 'm', 'length')
        self.log.debug(' ADDING\t << height >> input: %.1f [%s] stored: %1f [m]' % (height, length_unit, self.height_m))

    def add_diameter(self, diameter, length_unit='mm'):
        '''
        Adds diameter in m to sample

        :param diameter: float
        :param length_unit: str

        changes sample.diameter_m from its initial value to a new value.
        Length unit should be given, if other than 'mm'. See convert2 helper function
        #todo point to helper
        '''
        self.diameter_m = float(diameter) * convert.convert2(length_unit, 'm', 'length')
        self.log.debug(
            ' ADDING\t << diameter >> input: %.1f [%s] stored: %1f [m]' % (diameter, length_unit, self.diameter_m))

    def add_measurement(self, mtype, mfile, machine, mag_method=None):
        '''

        :param mtype: str - the type of measurement
        :param mfile: str -  the measurement file
        :param machine: str - the machine from which the file is output
        :param mag_method: str - only used for af-demag
        :return: RockPyV3.measurement object

        :mtypes:

        - af-demagnetization - 'af-demag'
        - hysteresis - 'hys'
        - Thellier-Thellier - 'thellier'
        - Zero Field Cooling-  'zfc'
        - IRM acquisition - 'irm'
        - Backfield - 'coe'
        - Viscosity - 'visc'

        '''

        implemented = {'af-demag': measurements.Af_Demag,
                       'hys': measurements.Hysteresis,
                       'palint': measurements.Thellier,
                       'thellier': measurements.Thellier,
                       'zfc': measurements.Zfc_Fc,
                       'irm': measurements.Irm,
                       'coe': measurements.Coe,
                       'visc': measurements.Viscosity,
        }

        if mtype.lower() in implemented:
            self.log.info(' ADDING\t << measurement >> %s' % mtype)
            measurement = implemented[mtype.lower()](self, mtype, mfile, machine, mag_method)
            # if measurement.raw_data:
            self.measurements.append(measurement)
            return measurement
        else:
            self.log.error(' << %s >> not implemented, yet' % (mtype))


    ''' RETURN FUNCTIONS '''

    def mass(self, mass_unit='mg'):
        '''
        Returns mass in specified mass unit
        :param mass_unit: str - unit to be output
        :return: float - mass in mass_unit
        '''
        OUT = self.mass_kg * convert.convert2('kg', mass_unit, 'mass')
        return OUT

    def height(self, length_unit='mm'):
        OUT = self.height_m * convert.convert2('m', length_unit, 'length')
        return OUT

    def diameter(self, length_unit='mm'):
        OUT = self.diameter_m * convert.convert2('m', length_unit, 'length')
        return OUT

    def volume(self, length_unit='mm'):
        out = np.pi * (self.diameter(length_unit=length_unit) / 2) ** 2 * self.height(length_unit=length_unit)
        return out

    def infos(self, mass_unit='mg', length_unit='mm', header=True):
        text = '%s\t ' % self.name
        hdr = 'Sample\t '
        if self.mass_kg:
            text += '%.2f\t ' % (self.mass(mass_unit=mass_unit))
            hdr += 'mass [%s]\t ' % mass_unit
        if self.height_m:
            text += '%.2f\t ' % self.height(length_unit=length_unit)
            hdr += 'length [%s]\t ' % length_unit
        if self.diameter_m:
            text += '%.2f\t ' % self.height(length_unit=length_unit)
            hdr += 'diameter [%s]\t ' % length_unit
        if header:
            print hdr
        print text


    ''' FIND FUNCTIONS '''

    def find_measurement(self, mtype):
        self.log.info('SEARCHING\t measurements with mtype << %s >>' % (mtype.lower()))
        out = [m for m in self.measurements if m.mtype == mtype.lower()]
        if len(out) != 0:
            self.log.info('FOUND\t sample << %s >> has %i measurements with mtype << %s >>' % (
                self.name, len(out), mtype.lower()))
        else:
            self.log.error('UNKNOWN\t mtype << %s >> or no measurement found' % mtype.lower())
        return out


    ''' ADDITIONAL '''

    def Experiment(self):
        experiment_matrix = {
            'af-demag': {},
        }

    def plot(self):
        # def plot(self, norm='mass', out='show', virgin=False, folder=None, name='output.pdf', figure=None):

        fig = plt.figure()

        mtypes = list(set([i.mtype for i in self.measurements]))

        print mtypes


def sample_import(sample_file, mass_unit='mg', length_unit='mm'):
    '''
    imports a csv list with mass, diameter and height data.
    has to be tab separated e.g.:

        Sample	Mass	Height	Diameter
        1a	320.0	5.17	5.84

    example
    -------

    samples = RockPyV3.sample_import(sample_file = the_data_file.txt, mass_unit = 'mg', length_unit='mm')

    :param sample_file: str
    :param mass_unit: str
    :param length_unit: str
    '''
    log = logging.getLogger('RockPy.READIN.get_data')
    reader_object = csv.reader(open(sample_file), delimiter='\t')
    list = [i for i in reader_object]
    header = list[0]

    dict = {i[0]: {header[j].lower(): float(i[j]) for j in range(1, len(i))} for i in list[1:]}

    out = {}
    for sample in dict:
        mass = dict[sample].get('mass', None)
        height = dict[sample].get('height', None)
        diameter = dict[sample].get('diameter', None)
        aux = {sample: Sample(sample, mass=mass, height=height, diameter=diameter, mass_unit=mass_unit,
                              length_unit=length_unit)}
        out.update(aux)

    return out


