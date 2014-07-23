# coding=utf-8
__author__ = 'Mike'
import logging
import numpy as np
import scipy as sp
from RockPyV3.Functions import general, convert
import measurements
from RockPyV3.Plots import general as RPplt
import matplotlib.pyplot as plt
import csv


class SampleGroup():
    general.create_logger('RockPy.SAMPLEGROUP')

    def __init__(self, sample_dict=None, sample_list=None, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.SAMPLEGROUP')
        else:
            self.log = log

        self.log.info('CREATING\t new sample group')

        self.samples = []
        self.sample_names = []

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


    def plot(self, mtype, norm='mass', value=None, rtn='show'):
        if mtype == 'hys':
            OUT = RPplt.Hys(samples=self.samples, norm=norm, value=value, rtn=rtn).show()
            return OUT


# class TT_Group(SampleGroup):
#
# def average_intensity(self):
# for sample in self.samples:
# p0 =

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


    ''' ADD FUNCTIONS '''

    def add_mass(self, mass, mass_unit='mg'):
        self.mass_kg = float(mass) * convert.convert2(mass_unit, 'kg', 'mass')
        self.log.debug(' ADDING\t << mass >> input: %.1f [%s] stored: %1f [kg]' % (mass, mass_unit, self.mass_kg))

    def add_height(self, height, length_unit='mm'):
        self.height_m = float(height) * convert.convert2(length_unit, 'm', 'length')
        self.log.debug(' ADDING\t << height >> input: %.1f [%s] stored: %1f [m]' % (height, length_unit, self.height_m))

    def add_diameter(self, diameter, length_unit='mm'):
        self.diameter_m = float(diameter) * convert.convert2(length_unit, 'm', 'length')
        self.log.debug(
            ' ADDING\t << diameter >> input: %.1f [%s] stored: %1f [m]' % (diameter, length_unit, self.diameter_m))

    def add_measurement(self, mtype, mfile, machine, mag_method=None):
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
            if measurement.raw_data:
                self.measurements.append(measurement)
            return measurement
        else:
            self.log.error(' << %s >> not implemented, yet' % (mtype))

    ''' RETURN FUNCTIONS '''

    def mass(self, mass_unit='mg'):
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


# ## tests

def run_6c_demag_test():
    test_sample = Sample('6c', mass=324., diameter=14)
    test_sample.add_height(33, 'mm')
    test_sample.mass()
    m = test_sample.add_measurement(mtype='AF-demag',
                                    mfile='/Users/Mike/PycharmProjects/RockPy V3/test_data/LF4c-6c_pdemag_P00.AF',
                                    machine='sushibar', mag_method='IRM')
    m = test_sample.add_measurement(mtype='AF-demag',
                                    mfile='/Users/Mike/PycharmProjects/RockPy V3/test_data/LF4c-6c_pdemag_P10.AF',
                                    machine='sushibar', mag_method='IRM')
    m.add_treatment(ttype='Pressure', options={'p_max': 0, 'p_seen': 0.2})
    m = test_sample.add_measurement(mtype='AF-demag',
                                    mfile='/Users/Mike/PycharmProjects/RockPy V3/test_data/LF4c-6c_pdemag_P11.AF',
                                    machine='sushibar', mag_method='IRMp')
    m.add_treatment(ttype='Pressure', options={'p_max': 0.2, 'p_seen': 0.2})
    RPplt.Af_Demag(test_sample).show()


def run_14a():
    sample1 = Sample('14a')
    sample2 = Sample('14b')
    m = sample1.add_measurement(mtype='AF-demag',
                                mfile='/Users/Mike/PycharmProjects/RockPy V3/test_data/LF6C-TRM(50muT).AF',
                                machine='sushibar',
                                mag_method='TRM@50µT')
    m = sample2.add_measurement(mtype='AF-demag',
                                mfile='/Users/Mike/PycharmProjects/RockPy V3/test_data/LF6C-TRM(50muT).AF',
                                machine='sushibar',
                                mag_method='TRM@50µT')
    RPplt.Af_Demag([sample1, sample2]).show()


# if __name__ == '__main__':
# run_14a()

