__author__ = 'Mike'
# coding=utf-8
from RockPyV3.Functions import general
import logging
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc, lines
import numpy as np
import backfield, hysteresis

class Plot(object):
    log = general.create_logger('RockPy.PLOTTING')

    def __init__(self, samples, norm=None, log=None, value=None, out='show'):
        matplotlib.rcParams.update({'font.size': 10})

        if not log:
            self.log = logging.getLogger('RockPy.PLOTTING')

        if type(samples) is not list:
            self.log.debug('CONVERTING Sample Instance to Samples List')
            samples = [samples]
        self.out = out
        self.norm = norm
        self.samples = [i for i in samples]
        self.fig1 = plt.figure(figsize=(8, 8), dpi=100)
        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        self.plot_data = []
        self.ax.xaxis.major.formatter._useMathText = True
        self.ax.ticklabel_format(style='sci', scilimits=(1, 4), axis='both')

        # rc('text', usetex=True)
        # rc('font', family='serif')


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


class Hys(Plot):
    def __init__(self):
        super(Hys, self).__init__(self, samples, norm=None, log=None, value=None, out='show')
        self.ax.axhline(0, color='#555555')
        self.ax.axvline(0, color='#555555')

    def show(self):
        for sample in self.samples:
            for measurement in sample.find_measurement('hys'):
                label = sample.name
                if measurement.treatment:
                    label += '\t' + measurement.treatment.get_label()

                factor = {'mass': measurement.sample.mass_kg,
                          'max': max(measurement.up_field[:, 1]),
                          'calibration': measurement.calibration_factor,
                          None: 1}

                self.norm_factor = factor[self.norm]

                self.log.info('NORMALIZING\t by: %s %s' % (self.norm, str(norm_factor)))
                self.log.info('PLOT\t of: %s' % (sample.name))

                std, = self.ax.plot(measurement.up_field[:, 0], measurement.up_field[:, 1] / self.norm_factor, label=label)
                self.ax.plot(measurement.down_field[:, 0], measurement.down_field[:, 1] / self.norm_factor,
                             color=std.get_color())

                handles, labels = self.ax.get_legend_handles_labels()
                self.ax.legend(handles, labels, prop={'size': 8}, loc='best')
                plt.xlabel('Field [T]')
                if self.norm == None:
                    plt.ylabel('Magnetic Moment $Am^2$')
                if self.norm == 'mass':
                    plt.ylabel('Magnetic Moment $Am^2/kg$')
                else:
                    plt.ylabel('Magnetic Moment normalized')

                plt.xlim([-2, 2])
                # plt.ylim([-1.1,1.1])
        if self.rtn == 'show':
            plt.show()
        if self.rtn == 'rtn':
            return self.ax
        else:
            plt.savefig(self.rtn, dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype='a4', format='pdf',
                        transparent=False, bbox_inches='tight', pad_inches=0.1,
                        frameon=None)


class Hys_Fabian2003(Plot):
    def __init__(self, sample, measurement, norm='mass', log=None, value=None, out='show'):
        super(Hys_Fabian2003, self).__init__(samples=sample, norm=norm, log=log, value=value, out=out)

        self.hys = measurement
        try:
            self.coe = self.samples[0].find_measurement(mtype='coe')[0]
        except IndexError:
            self.log.error('NOT FOUND\t << coe >> measurement')
            self.coe = None

        self.factor = {'mass': measurement.sample.mass_kg,
                       'max': max(measurement.up_field[:, 1]),
                       'calibration': measurement.calibration_factor,
                       None: 1}
        self.norm_factor = self.factor[norm]
        self.ax.xaxis.major.formatter._useMathText = True
        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        # self.ax_low = plt.subplot2grid((3, 1), (2, 0), sharex=self.ax, colspan=1, rowspan=1)
        self.ax.axhline(0, color='#555555')
        self.ax.axvline(0, color='#555555')

        plt.suptitle(sample.name)



    def add_label(self, ax=None):
        if ax == None:
            ax = self.ax

        if self.norm == None:
            ax.set_ylabel('Magnetic Moment $Am^2$', fontsize=10)
        if self.norm == 'mass':
            ax.set_ylabel('Magnetic Moment $Am^2/kg$', fontsize=10)
        else:
            ax.set_ylabel('Magnetic Moment normalized', fontsize=10)

        ax.set_xlabel('Field [$T$]', fontsize=10)


    def show(self, out='show', folder=None, name='output.pdf'):
        norm_factor = self.factor[self.norm]  # # NORMALIZATION FACTOR

        self.ax.set_ylim(
            [min(self.hys.down_field[:, 1]) / norm_factor * 1.1, max(self.hys.down_field[:, 1]) / norm_factor * 1.1])
        self.ax.set_xlim([min(self.hys.down_field[:, 0]), max(self.hys.down_field[:, 0])])

        backfield.plot_coe(coe_obj=self.coe, ax=self.ax, norm_factor=self.norm_factor)

        hysteresis.plot_hys(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_virgin(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_rev(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_irrev(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')

        hysteresis.fill_hys(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor)
        hysteresis.fill_virgin(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor)

        hysteresis.add_virgin_info(self.hys, ax=self.ax, norm_factor=self.norm_factor)

        hysteresis.add_05ms_line(hys_obj=self.hys,ax=self.ax, norm_factor=self.norm_factor, text=True)
        hysteresis.add_ms_line(hys_obj=self.hys,ax=self.ax, norm_factor=self.norm_factor, text=True)

        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # handles, labels = self.ax.get_legend_handles_labels()
        # self.ax_low.legend(handles, labels, prop={'size': 8})
        # self.ax_low.legend(loc='best', frameon = False, prop={'size':8})

        if out == 'show':
            plt.show()

        if out == 'save':
            if folder != None:
                plt.savefig(folder + self.samples[0].name + '_' + name, dpi=300)

        if out == 'rtn':
            return self.fig1

class IRM(Plot):
    def show(self):
        for sample in self.samples:
            for measurement in sample.find_measurement('irm'):
                label = sample.name
                if measurement.treatment:
                    label += '\t' + measurement.treatment.get_label()
                factor = {'mass': measurement.sample.mass_kg,
                          'max': max(measurement.up_field[:, 1]),
                          None: 1}
                norm_factor = factor[self.norm]
                self.log.info('NORMALIZING\t by: %s %s' % (self.norm, str(norm_factor)))
                self.log.info('PLOT\t of: %s' % (sample.name))

                std, = self.ax.plot(measurement.up_field[:, 0], measurement.up_field[:, 1] / norm_factor, label=label)
                self.ax.plot(measurement.down_field[:, 0], measurement.down_field[:, 1] / norm_factor,
                             color=std.get_color())

                handles, labels = self.ax.get_legend_handles_labels()
                self.ax.legend(handles, labels, prop={'size': 8})
        plt.show()

