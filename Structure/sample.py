# coding=utf-8
__author__ = 'Mike'
import logging
import numpy as np
import scipy as sp
from Functions import General, Convert
import Measurements
from Plots import General as RPplt


class Sample():
    General.create_logger('RockPy.SAMPLE')

    def __init__(self, name, mass=1, mass_unit=None, height=None, diameter=None, length_unit=None):
        self.name = name
        self.log = logging.getLogger('RockPy.SAMPLE')
        self.log.info('CREATING\t new sample << %s >>' % self.name)

        self.measurements = []

        if mass_unit:
            mass_factor = Convert.convert2(mass_unit, 'kg', 'mass')
        else:
            self.log.info(' MISSING\t length_unit: assuming << mm >>')
            mass_factor = Convert.convert2('mg', 'kg', 'mass')

        if mass:
            self.mass_kg = mass * mass_factor
            self.log.debug(' ADDING\t << mass >> input: %.1f [%s] stored: %1f [kg]' % (mass, mass_unit, self.mass_kg))
        else:
            self.log.debug('MISSING\t << mass >>')
            self.mass_kg = None

        # get length _unit for conversion
        if length_unit:
            length_factor = Convert.convert2(length_unit, 'm', 'length')
        else:
            self.log.info(' MISSING\t mass_unit: assuming << mg >>')
            length_factor = Convert.convert2('mm', 'm', 'length')

        if height:
            self.height_m = float(height) * length_factor
            self.log.debug(
                ' ADDING\t << height >> input: %.1f [%s] stored: %1f [kg]' % (height, length_unit, self.height_m))

        else:
            self.log.debug('MISSING\t << height >>')
            self.height_m = None

        if diameter:
            self.diameter_m = float(diameter) * length_factor
            self.log.debug(
                ' ADDING\t << diameter >> input: %.1f [%s] stored: %1f [kg]' % (diameter, length_unit, self.diameter_m))

        else:
            self.log.debug('MISSING\t << diameter >>')
            self.diameter_m = None


    ''' ADD FUNCTIONS '''

    def add_mass(self, mass, mass_unit='mg'):
        self.mass_kg = float(mass) * Convert.convert2(mass_unit, 'kg', 'mass')
        self.log.debug(' ADDING\t << mass >> input: %.1f [%s] stored: %1f [kg]' % (mass, mass_unit, self.mass_kg))

    def add_height(self, height, length_unit='mm'):
        self.height_m = float(height) * Convert.convert2(length_unit, 'm', 'length')
        self.log.debug(' ADDING\t << height >> input: %.1f [%s] stored: %1f [m]' % (height, length_unit, self.height_m))

    def add_diameter(self, diameter, length_unit='mm'):
        self.diameter_m = float(diameter) * Convert.convert2(length_unit, 'm', 'length')
        self.log.debug(
            ' ADDING\t << diameter >> input: %.1f [%s] stored: %1f [m]' % (diameter, length_unit, self.diameter_m))

    def add_measurement(self, mtype, mfile, machine, mag_method=None):
        implemented = {'af-demag': Measurements.Af_Demag,
                       'hys':Measurements.Hysteresis}
        #todo hys
        #todo coe
        #todo irm
        #todo palint
        if mtype.lower() in implemented:
            measurement = implemented[mtype.lower()](self, mtype, mfile, machine, mag_method)
            if measurement.raw_data:
                self.measurements.append(measurement)
            return measurement

    ''' RETURN FUNCTIONS '''

    def mass(self, mass_unit='mg'):
        OUT = self.mass_kg * Convert.convert2('kg', mass_unit, 'mass')
        return OUT

    def height(self, mass_unit='mm'):
        OUT = self.height_m * Convert.convert2('m', mass_unit, 'length')
        return OUT

    def diamter(self, mass_unit='mm'):
        OUT = self.diameter_m * Convert.convert2('m', mass_unit, 'length')
        return OUT

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


if __name__ == '__main__':
    run_14a()
