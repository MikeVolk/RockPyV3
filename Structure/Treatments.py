# coding=utf-8
__author__ = 'Mike'
import logging
from RockPyV3.Functions import general


class Treatment(object):
    # general.create_logger('RockPy.TREATMENT')

    def __init__(self, sample, measurement, ttype, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.TREATMENT')
        self.sample = sample
        self.measurement = measurement
        self.ttype = ttype


class Pressure(Treatment):
    # general.create_logger('RockPy.TREATMENT.pressure')

    def __init__(self, sample, measurement, options):
        self.log = logging.getLogger('RockPy.TREATMENT.pressure')
        super(Pressure, self).__init__(sample, measurement, ttype='pressure')
        self.log.info('CREATING\t new treatment << %s >> for sample << %s >> | measurement << %s >>' % (
            'pressure', self.sample.name, self.measurement.mtype))
        self.log.info('CREATING\t                                                    with options: %s' % options)

        p_max = options.get('p_max', 0)
        p_seen = options.get('p_seen', 0)
        if p_seen < p_max:
            p_seen = p_max
            self.log.warning('P_SEEN < P_MAX --> impossible --> setting << P_SEEN = P_MAX >>')

        self.p_unit = options.get('p_unit', 'GPa')
        self.p_max = p_max
        self.p_seen = p_seen
        self.label = self.get_label()

    def get_label(self):
        if self.p_seen != self.p_max:
            label = 'P ' + str(self.p_seen) + '|' + str(self.p_max) + ' ' + self.p_unit
        else:
            label = 'P ' + str(self.p_seen) + ' ' + self.p_unit
        return label

    def __repr__(self):
        return 'P' + str(self.p_seen) + '|' + str(self.p_max)