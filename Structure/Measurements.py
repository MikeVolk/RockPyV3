# coding=utf-8
__author__ = 'Mike'
from Functions import General
import Treatments
import ReadIn
import logging
import numpy as np
import copy


class Measurement(object):
    General.create_logger(__name__)

    def __init__(self, sample, mtype, mfile, machine, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.MEASUREMENT')

        implemented = {
            'af-demag': '',
            'parm-spectra':'',
            'hysteresis' : '',
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
        implemented = {'sushibar': {'af-demag': ReadIn.Machines.sushibar,
                                    'parm-spectra': ReadIn.Machines.sushibar},
                       'vsm':{'hysteresis':ReadIn.Machines.vsm,
                              'irm': ReadIn.Machines.vsm}}

        self.log.info('IMPORTING << %s , %s >> data' % (self.machine, self.mtype))
        if self.machine in implemented:
            if self.mtype in implemented[self.machine]:
                self.raw_data = implemented[self.machine][self.mtype](self.mfile, self.sample.name)
                if self.raw_data == None:
                    self.log.error('IMPORTING\t did not transfer data - CHECK sample name and data file')
                    return
            else:
                self.log.error('UNKNOWN\t measurement type << %s >>' % (self.mtype))
        else:
            self.log.error('UNKNOWN\t machine << %s >>' % (self.machine))

    def add_treatment(self, ttype, options=None):
        self.ttype = ttype.lower()
        self.log.info('ADDING\t treatment to measurement << %s >>' % (self.ttype))

        implemented = {
            'pressure': Treatments.Pressure,
        }

        if ttype.lower() in implemented:
            self.treatment = implemented[ttype.lower()](self.sample, self, options)
        else:
            self.log.error('UNKNOWN\t treatment type << %s >> is not know or not implemented' % ttype)


class Af_Demag(Measurement):
    General.create_logger('RockPy.MEASUREMENT.AF-DEMAG')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.AF-DEMAG')
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
    General.create_logger('RockPy.MEASUREMENT.PARM-SPECTRA')

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
            self.windows = np.array([[self.__dict__['l_window_limit'][i], self.__dict__['u_window_limit'][i]] for i in range(len(self.u_window_limit))])
        except AttributeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

        except TypeError:
            self.log.error('SOMETHING IS NOT RIGHT - raw data not transfered')
            return None

    def subtract_af3(self):
        self.log.debug('CREATING\t COPY of %s' %self)
        self.log.info('SUBTRACTING\t AF3 data of pARM measurement')
        self_copy = copy.copy(self)
        self_copy.x -= self.x[0]
        self_copy.y -= self.y[0]
        self_copy.z -= self.z[0]

        self_copy.m = np.array([np.sqrt(self.x[i]**2+self.y[i]**2+self.z[i]**2) for i in range(len(self.x))])
        self.log.info('RETURNING\t COPY with subtracted AF3')
        return self_copy

class Hysteresis(Measurement):
    General.create_logger('RockPy.MEASUREMENT.HYSTERESIS')

    def __init__(self, sample, mtype, mfile, machine, mag_method):
        self.log = logging.getLogger('RockPy.MEASUREMENT.HYSTERESIS')
        Measurement.__init__(self, sample, mtype, mfile, machine, self.log)

        implemented = {'vsm':}