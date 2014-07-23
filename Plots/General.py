# coding=utf-8

__author__ = 'Mike'
import RockPyV3
from RockPyV3.Functions import general
from RockPyV3.Plots import paleointensity
# from RockPyV3.Structure.measurements import Thellier
import logging
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc, lines
import numpy as np
import backfield, hysteresis


class Plot(object):
    log = general.create_logger('RockPy.PLOTTING')

    def __init__(self, samples, norm=None, log=None, value=None, out='show', folder=None, name='output.pdf'):
        matplotlib.rcParams.update({'font.size': 10})
        params = {'backend': 'ps',
                  'text.latex.preamble': [r"\usepackage{upgreek}",
                                          r"\usepackage[nice]{units}"],
                  'axes.labelsize': 12,
                  'text.fontsize': 12,
                  'legend.fontsize': 8,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10,
                  # 'text.usetex': True,
                  'axes.unicode_minus': True}
        plt.rcParams.update(params)

        if not log:
            self.log = logging.getLogger('RockPy.PLOTTING')

        if type(samples) is not list:
            self.log.debug('CONVERTING Sample Instance to Samples List')
            samples = [samples]

        self.name = name

        if folder == None:
            from os.path import expanduser

            folder = expanduser("~") + '/Desktop/'
        self.folder = folder

        self.norm = norm
        self.samples = [i for i in samples]
        self.fig1 = plt.figure(figsize=(8, 8), dpi=100)
        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        self.plot_data = []
        self.ax.xaxis.major.formatter._useMathText = True
        self.ax.yaxis.major.formatter._useMathText = True
        self.ax.ticklabel_format(style='sci', scilimits=(1, 4), axis='both')
        plt.rc('axes', color_cycle=['k', 'm', 'y', 'c'])

    def out(self, out, *args):
        out_options = {'show': plt.show,
                       'rtn': self.get_fig,
                       'save': self.save_fig}

        if out in ['show', 'save']:
            if not 'nolable' in args:
                self.ax.set_xlabel(self.x_label)
                self.ax.set_ylabel(self.y_label)
        out_options[out]()

    def get_ax(self):
        return self.ax

    def get_fig(self):
        return self.fig1

    def save_fig(self):
        plt.savefig(self.folder + self.samples[0].name + '_' + self.name, dpi=300)


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


class Hysteresis(Plot):
    def __init__(self, sample, norm='mass', log=None, value=None, out='show'):
        super(Hysteresis, self).__init__(samples=sample, norm=norm, log=log, value=value, out=out)
        self.measurements = sample.find_measurement('hys')

    def show(self, out='show', folder=None, name='output.pdf'):
        for measurement in self.measurements:
            factor = {'mass': measurement.sample.mass_kg,
                      'max': max(measurement.up_field[:, 1]),
                      'calibration': measurement.calibration_factor,
                      None: 1}
            norm_factor = factor[self.norm]  # # NORMALIZATION FACTOR

            hysteresis.plot_hys(measurement, ax=self.ax, norm_factor=norm_factor, out='rtn')

        if out == 'show':
            plt.ylim([-1.1, 1.1])
            plt.show()

        if out == 'save':
            if folder != None:
                plt.savefig(folder + self.samples[0].name + '_' + name, dpi=300)

        if out == 'rtn':
            return self.fig1


class Hys(Plot):
    def __init__(self):
        super(Hys, self).__init__(self, samples, norm=None, log=None, value=None, out='show')


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

                std, = self.ax.plot(measurement.up_field[:, 0], measurement.up_field[:, 1] / self.norm_factor,
                                    label=label)
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
        backfield.add_bcr_line(coe_obj=self.coe, ax=self.ax, norm_factor=self.norm_factor, method='fit', text=True)

        hysteresis.plot_hys(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_virgin(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_rev(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')
        hysteresis.plot_irrev(self.hys, ax=self.ax, norm_factor=self.norm_factor, out='rtn')

        hysteresis.fill_hys(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor)
        hysteresis.fill_virgin(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor)

        hysteresis.add_virgin_info(self.hys, ax=self.ax, norm_factor=self.norm_factor)

        hysteresis.add_05ms_line(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor, text=True)
        hysteresis.add_ms_line(hys_obj=self.hys, ax=self.ax, norm_factor=self.norm_factor, text=True)

        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.add_label(ax=self.ax)
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


class Dunlop(Plot):
    def __init__(self, sample_obj, component='m', norm='mass', log=None, value=None, out='show', folder=None,
                 name='dunlop plot', **plt_opt):


        measurement_test = False

        if isinstance(sample_obj, RockPyV3.Structure.measurements.Thellier):
            measurement_test = True
            self.measurements = [sample_obj]

        if isinstance(sample_obj, RockPyV3.Structure.measurements.Thellier):
            sample_obj = sample_obj.sample

        super(Dunlop, self).__init__(samples=sample_obj, norm=norm, log=log, value=value, out=out, folder=folder,
                                     name=name)

        if not measurement_test:
            self.measurements = sample_obj.find_measurement('palint')

        for measurement in self.measurements:
            components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
            factors = {'mass': sample_obj.mass_kg,
                       'max': measurement.sum[:, components[component]],
                       'trm': measurement.trm[:, components[component]],
            }

            norm_factor = factors[norm]
            idx = self.measurements.index(
                measurement)  # index of measurement used to distinguish repeated measurements on the same sample

            lines = ['-', '--', ':']

            plt_opt.update({'linestyle': lines[idx]})
            paleointensity.dunlop(palint_object=measurement, ax=self.ax,
                                  plt_idx=idx,
                                  norm_factor=norm_factor, component=component,
                                  plt_opt=plt_opt)

            paleointensity.add_dunlop_labels(palint_object=measurement, ax=self.ax,
                                             norm=norm, norm_factor=norm_factor,
                                             text=True, plt_idx=idx)
            self.fig1.suptitle('%s Dunlop Plot' % sample_obj.name, fontsize=16)

        self.out(out, 'nolable')


class Arai(Plot):
    def __init__(self, sample_obj, component='m', norm='mass', log=None, value=None,
                 t_min=20, t_max=700, line=True, check=True,
                 out='show', folder=None,
                 name='arai plot', plt_opt={}):
        super(Arai, self).__init__(samples=sample_obj, norm=norm, log=log, value=value, out=out, folder=folder,
                                   name=name)
        self.measurements = sample_obj.find_measurement('palint')
        self.t_min = t_min
        self.t_max = t_max

        for measurement in self.measurements:
            components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
            factors = {'mass': [sample_obj.mass_kg, sample_obj.mass_kg],
                       'max': [max(measurement.sum[:, components[component]]),
                               max(measurement.sum[:, components[component]])],
                       'trm': [measurement.trm[:, components[component]], measurement.trm[:, components[component]]],
                       'th': [measurement.th[0, components[component]], measurement.th[0, components[component]]],
                       'dnorm': [measurement.ptrm[-1, components[component]], measurement.th[0, components[component]]],
            }

            norm_factor = factors[norm]

            idx = self.measurements.index(
                measurement)  # index of measurement used to distinguish repeated measurements on the same sample

            lines = ['-', '--', ':']

            if not 'ls' in plt_opt:
                if not 'linestyle' in plt_opt:
                    plt_opt.update({'linestyle': lines[idx]})

            paleointensity.arai(palint_object=measurement, ax=self.ax,
                                t_min=self.t_min, t_max=self.t_max,
                                line=line, check=check,
                                plt_idx=idx,
                                norm_factor=norm_factor, component=component,
                                plt_opt=plt_opt)

            self.fig1.suptitle('%s Arai Plot' % sample_obj.name, fontsize=16)

        self.x_label = 'pTRM gained'
        self.y_label = 'NRM remaining'

        self.out(out)