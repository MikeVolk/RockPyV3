__author__ = 'Mike'
import logging
import matplotlib.pyplot as plt


class Plot(object):
    log = logging.getLogger('RockPy.PLOTTING')

    def __init__(self, samples, norm=None, log=None):
        if not log:
            self.log = logging.getLogger('RockPy.PLOTTING')
        if type(samples) is not list:
            self.log.debug('CONVERTING Sample Instance to Samples List')
            samples = [samples]

        self.samples = [i for i in samples]
        self.fig1 = plt.figure(figsize=(11.69, 8.27), dpi=100)
        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        self.plot_data = []


class Af_Demag(Plot):
    def show(self, dtype='m'):
        for sample in self.samples:
            for measurement in sample.find_measurement('af-demag'):
                label = sample.name + measurement.mag_method
                if measurement.treatment:
                    label += '\t' + measurement.treatment.get_label()
                self.ax.plot(measurement.fields, getattr(measurement, dtype), label=label)
                handles, labels = self.ax.get_legend_handles_labels()
                self.ax.legend(handles, labels, prop={'size': 8})
        plt.show()