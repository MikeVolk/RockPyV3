# coding=utf-8
__author__ = 'Mike'
import logging
from RockPyV3.Functions import general


class Treatment(object):
    general.create_logger('RockPy.TREATMENT')

    def __init__(self, sample, measurement, ttype, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.TREATMENT')
            self.sample = sample
            self.measurement = measurement
            self.ttype = ttype


class Pressure(Treatment):
    def __init__(self, sample, measurement, options):
        super(Pressure, self).__init__(sample, measurement, ttype='pressure')
        self.log.info('CREATING\t new treatment << %s >> for sample << %s >> | measurement << %s >>' % (
            'pressure', self.sample.name, self.measurement.mtype))
        p_max = options.get('p_max', 0)
        p_seen = options.get('p_seen', 0)
        self.p_max = p_max
        self.p_seen = p_seen

    def get_label(self):
        label = 'P' + str(self.p_seen) + '|' + str(self.p_max)
        return label

    def __repr__(self):
        return 'P' + str(self.p_seen) + '|' + str(self.p_max)