# coding=utf-8

__author__ = 'Mike'
import RockPyV3
from RockPyV3.Functions import general
from RockPyV3.Plots import paleointensity, helper
import logging
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc, lines
import numpy as np
import backfield, hysteresis


class Plot(object):
    log = general.create_logger('RockPy.PLOTTING')

    def __init__(self, samples_list, norm=None, log=None, value=None,
                 plot='show', folder=None, name='output.pdf',
                 plt_opt={},
                 **options):

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

        self.plot = plot

        if not log:
            self.log = logging.getLogger('RockPy.PLOTTING')
        else:
            self.log = logging.getLogger(log)

        if type(samples_list) is not list:
            self.log.debug('CONVERTING Sample Instance to Samples List')
            samples_list = [samples_list]

        self.name = name

        if folder == None:
            from os.path import expanduser

            folder = expanduser("~") + '/Desktop/'
        self.folder = folder
        self.colors = helper.get_colors()
        self.norm = norm
        self.samples = [i for i in samples_list]

        self.fig1 = options.get('fig', plt.figure(figsize=(8, 8), dpi=100))

        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
        self.plot_data = []
        self.ax.xaxis.major.formatter._useMathText = True
        self.ax.yaxis.major.formatter._useMathText = True
        self.ax.ticklabel_format(style='sci', scilimits=(1, 4), axis='both')
        plt.rc('axes', color_cycle=['k', 'm', 'y', 'c'])

    def out(self, *args):
        out_options = {'show': plt.show,
                       'rtn': self.get_fig,
                       'save': self.save_fig}
        if self.plot in ['show', 'save']:
            if not 'nolable' in args:
                self.ax.set_xlabel(self.x_label)
                self.ax.set_ylabel(self.y_label)
        out_options[self.plot]()

    def get_ax(self):
        return self.ax

    def get_fig(self):
        return self.fig1

    def save_fig(self):
        plt.savefig(self.folder + self.samples[0].name + '_' + self.name, dpi=300)


class Af_Demag(Plot):

    def __init__(self, samples_list, norm='mass', value=None,
                 plot='show', folder=None, name='hysteresis',
                 plt_opt={}, **options):
        log = 'RockPy.PLOT.af-demag'
        dtype = options.get('dtype', 'm')
        super(Af_Demag, self).__init__(samples_list=samples_list,
                                 norm=norm, log=log, value=value,
                                 plot=plot, folder=folder, name=name,
                                 **options)
        try:
            self.show(dtype=dtype)
            self.out()
        except AttributeError:
            self.log.error('can\'t plot')

    def show(self, dtype='m'):
        self.x_label = 'field'
        self.y_label = 'moment'
        for sample in self.samples:
            for measurement in sample.find_measurement('af-demag'):
                label = sample.name + measurement.mag_method
                if measurement.treatment:
                    label += '  ' + measurement.treatment.get_label()
                self.ax.plot(measurement.fields, getattr(measurement, dtype), label=label)
                handles, labels = self.ax.get_legend_handles_labels()

        self.ax.legend(handles, labels)
        return self.ax

class Hysteresis(Plot):
    log = general.create_logger('RockPy.PLOTTING.hysteresis')

    def __init__(self, samples_list, norm='mass', log=None, value=None,
                 plot='show', folder=None, name='hysteresis',
                 plt_opt={}, **options):

        if not log:
            log = logging.getLogger('RockPy.PLOTTING.hysteresis')

        super(Hysteresis, self).__init__(samples_list=samples_list,
                                         norm=norm, log=log, value=value,
                                         plot=plot, folder=folder, name=name,
                                         **options)

        self._get_labels()
        try:
            self.show()
            self.out()
        except AttributeError:
            pass

    def show(self):
        for sample in self.samples:

            self.measurements = sample.find_measurement('hys')
            for measurement in self.measurements:
                factor = {'mass': measurement.sample.mass_kg,
                          'max': max(measurement.up_field[:, 1]),
                          'calibration': measurement.calibration_factor,
                          None: 1}
                norm_factor = factor[self.norm]  # # NORMALIZATION FACTOR

                self.ax = hysteresis.plot_hys(measurement, ax=self.ax, norm_factor=norm_factor, out='rtn')
                self.ax = hysteresis.plot_virgin(measurement, ax=self.ax, norm_factor=norm_factor, out='rtn')

    def _get_labels(self):
        self.x_label = 'Field [T]'

        if self.norm == None:
            self.y_label = 'Magnetic Moment $Am^2$'
        if self.norm == 'mass':
            self.y_label = 'Magnetic Moment $Am^2/kg$'
        else:
            self.y_label = 'Magnetic Moment normalized'


class Hys_Fabian2003(Hysteresis):
    log = general.create_logger('RockPy.PLOTTING.fabian2003_hys')

    def __init__(self, samples_list, norm='mass', log=None, value=None,
                 plot='show', folder=None, name='Fabian2003_hysteresis.pdf',
                 plt_opt={}, **options):
        super(Hys_Fabian2003, self).__init__(samples_list=samples_list,
                                             norm=norm, log=log, value=value,
                                             plot=plot, folder=folder, name=name,
                                             **options)

        if len(self.samples) > 1:
            self.log.warning('MORE THAN one sample found, using first one for plot')

        self.hys = self.samples[0].find_measurement('hys')[0]
        self.samples[0].find_measurement('hys')[0]

        try:
            self.coe = self.samples[0].find_measurement(mtype='coe')[0]
        except IndexError:
            self.log.error('NOT FOUND\t << coe >> measurement')
            self.coe = None

        self.factor = {'mass': self.hys.sample.mass_kg,
                       'max': max(self.hys.up_field[:, 1]),
                       'calibration': self.hys.calibration_factor,
                       None: 1}

        self.norm_factor = self.factor[norm]

        self.ax.xaxis.major.formatter._useMathText = True
        self.ax = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)

        plt.suptitle(self.hys.sample.name)

        self._get_labels()
        self.show()
        self.out()

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


    def show(self):
        self.factor = {'mass': self.hys.sample.mass_kg,
                       'max': max(self.hys.up_field[:, 1]),
                       'calibration': self.hys.calibration_factor,
                       None: 1}

        norm_factor = self.factor[self.norm]  # # NORMALIZATION FACTOR
        self.ax.set_ylim(
            [min(self.hys.down_field[:, 1]) / norm_factor * 1.1, max(self.hys.down_field[:, 1]) / norm_factor * 1.1])
        self.ax.set_xlim([min(self.hys.down_field[:, 0]), max(self.hys.down_field[:, 0])])

        if self.coe:
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
    # self, samples_list, norm='mass', log=None, value=None,
    #              plot='show', folder=None, name='hysteresis',
    #              plt_opt={}, **options)
    def __init__(self, sample_obj, component='m', norm='mass', log=None, value=None, plot='show', folder=None,
                 name='dunlop plot', **plt_opt):


        measurement_test = False

        if isinstance(sample_obj, RockPyV3.Structure.measurements.Thellier):
            measurement_test = True
            self.measurements = [sample_obj]

        if isinstance(sample_obj, RockPyV3.Structure.measurements.Thellier):
            sample_obj = sample_obj.sample

        super(Dunlop, self).__init__(samples_list=sample_obj, norm=norm, log=log, value=value, plot=plot, folder=folder,
                                     name=name)
        for sample_obj in self.samples:
            if not measurement_test:
                self.measurements = sample_obj.find_measurement('palint')

            for measurement in self.measurements:
                components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
                factors = {'mass': sample_obj.mass_kg,
                           'max': measurement.sum[:, components[component]],
                           'trm': measurement.trm[:, components[component]],
                           'None': 1,
                           }

                norm_factor = factors[norm]

                idx = self.measurements.index(measurement)  # index of measurement used to distinguish repeated measurements on the same sample
                if len(self.samples) >1:
                    idx = self.samples.index(sample_obj)  # index of measurement used to distinguish samples


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

        self.out(plot, 'nolable')


class Arai(Plot):
    def __init__(self, sample_obj, component='m', norm='mass', log=None, value=None,
                 t_min=20, t_max=700, line=True, check=True,
                 plot='show', folder=None,
                 name='arai plot', plt_opt={},
                 **options):
        super(Arai, self).__init__(samples_list=sample_obj, norm=norm, log=log, value=value, plot=plot, folder=folder,
                                   name=name)

        for sample_obj in self.samples:
            self.measurements = sample_obj.find_measurement('palint')
            self.t_min = t_min
            self.t_max = t_max
            for measurement in self.measurements:
                components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
                factors = {'mass': [sample_obj.mass_kg, sample_obj.mass_kg],
                           'max': [max(measurement.sum[:, components[component]]),
                                   max(measurement.sum[:, components[component]])],
                           'trm': [measurement.trm[:, components[component]][0],
                                   measurement.trm[:, components[component]][0]],
                           'th': [measurement.th[0, components[component]], measurement.th[0, components[component]]],
                           'dnorm': [measurement.ptrm[-1, components[component]], measurement.th[0, components[component]]],
                           'none': [1, 1],
                           }

                norm_factor = factors[norm.lower()]

                idx = self.measurements.index(
                    measurement)  # index of measurement used to distinguish repeated measurements on the same sample
                if len(self.samples) >1:
                    idx = self.samples.index(sample_obj)  # index of measurement used to distinguish samples


                lines = ['-', '--', ':']

                if not 'ls' in plt_opt:
                    if not 'linestyle' in plt_opt:
                        plt_opt.update({'linestyle': lines[idx]})

                paleointensity.arai(palint_object=measurement, ax=self.ax,
                                    t_min=self.t_min, t_max=self.t_max,
                                    line=line, check=check,
                                    plt_idx=idx,
                                    norm_factor=norm_factor, component=component,
                                    plt_opt=plt_opt, **options)

                if measurement.th_stdev != None:
                    paleointensity.arai_stdev(palint_object=measurement, ax=self.ax,
                                    t_min=self.t_min, t_max=self.t_max,
                                    plt_idx=idx,
                                    norm_factor=norm_factor, component=component,
                                    plt_opt=plt_opt, **options)


                self.fig1.suptitle('%s Arai Plot' % sample_obj.name, fontsize=16)

        self.x_label = 'pTRM gained'
        self.y_label = 'NRM remaining'
        plt.gca().set_ylim(bottom=0)
        self.out(plot)