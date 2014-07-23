__author__ = 'mike'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, lines
from scipy import stats
from RockPyV3 import ReadIn
import helper


def dunlop(palint_object, ax, component='m', norm_factor=1, plt_idx=0, out='show', folder=None, name='output.pdf',
           plt_opt={}):
    components = {'x': 1, 'y': 2, 'z': 3, 'm': 4}
    marker = helper.get_marker()
    ax.plot(palint_object.th[:, 0], palint_object.th[:, components[component]] / norm_factor, '.-',
            marker=marker[plt_idx],
            color='r', **plt_opt)
    ax.plot(palint_object.sum[:, 0], palint_object.sum[:, components[component]] / norm_factor, '.-',
            marker=marker[plt_idx],
            color='b', **plt_opt)
    ax.plot(palint_object.ptrm[:, 0], palint_object.ptrm[:, components[component]] / norm_factor, '.-',
            marker=marker[plt_idx],
            color='g', **plt_opt)
    ax.plot(palint_object.tr[:, 0], palint_object.tr[:, components[component]] / norm_factor,
            color='k', marker=marker[plt_idx], alpha=0.8, ls='')
    return ax


def arai(palint_object, ax, component='m', norm=None, norm_factor=[1, 1], plt_idx=0, t_min=20, t_max=700,
         out='show', folder=None, name='output.pdf',
         line=True, check=True, plt_opt={}):
    idx = palint_object.components[component]
    markers = helper.get_marker()
    colors = helper.get_colors()
    components = palint_object.components
    factors = {'mass': [palint_object.sample.mass_kg, palint_object.sample.mass_kg],
               'max': [max(palint_object.sum[:, components[component]]),
                       max(palint_object.sum[:, components[component]])],
               'trm': [palint_object.trm[:, components[component]], palint_object.trm[:, components[component]]],
               'th': [palint_object.th[0, components[component]], palint_object.th[0, components[component]]],
               'dnorm': [palint_object.ptrm[-1, components[component]], palint_object.th[0, components[component]]],
               'none': [1, 1]
    }

    if norm == 'th':
        norm_factor_ptrm = palint_object.th[0, idx]
        norm_factor_th = palint_object.th[0, idx]
    if norm == None:
        norm_factor_ptrm = 1
        norm_factor_th = 1

    xy = np.array(
        [[i[idx], j[idx]] for i in palint_object.ptrm for j in palint_object.th if i[0] == j[0]]) / norm_factor

    th = xy[:, 1]
    ptrm = xy[:, 0]

    if not 'ls' in plt_opt or 'linestyle' in plt_opt:
        plt_opt.update({'ls': ''})

    ax.plot(np.fabs(ptrm), np.fabs(th),
            marker=markers[plt_idx], color=colors[plt_idx], **plt_opt)

    if line:
        add_arai_line(palint_object=palint_object, ax=ax,
                      t_min=t_min, t_max=t_max,
                      plt_idx=plt_idx,
                      norm_factor=norm_factor,
                      component=component,
                      plt_opt=plt_opt)

    if check:
        add_ck_check(palint_object=palint_object, ax=ax,
                     plt_idx=plt_idx,
                     norm_factor=norm_factor,
                     component=component,
                     plt_opt=plt_opt)
        add_ac_check(palint_object=palint_object, ax=ax,
                     plt_idx=plt_idx,
                     norm_factor=norm_factor,
                     component=component,
                     plt_opt=plt_opt)
    return ax


def add_arai_line(palint_object, ax, component='m', t_min=20, t_max=700, norm=None, norm_factor=[1, 1],
                  plt_idx=0,
                  plt_opt={}):
    colors = helper.get_colors()
    idx = palint_object.components[component]

    slopes, sigmas, y_intercept, x_intercept = palint_object.calculate_slope(t_min=t_min, t_max=t_max)

    y, x = palint_object._get_th_ptrm_data(t_min=t_min, t_max=t_max)
    x_m = [np.linalg.norm(i) for i in x]
    x = np.c_[x, x_m]

    x_fit = np.linspace(x[0, idx - 1], x[-1, idx - 1], 2)
    y_fit = (y_intercept[idx - 1] + slopes[idx - 1] * x_fit)

    ax.plot(x_fit / norm_factor[0], y_fit / norm_factor[1], color=colors[plt_idx])
    add_slope_text(palint_object=palint_object, ax=ax, slope=slopes[idx - 1], plt_idx=plt_idx)


def add_ac_check(palint_object, ax, component='m', norm=None,
                 norm_factor=[1, 1],
                 plt_idx=0, plt_opt={},
                 **options):
    colors = helper.get_colors()
    idx = palint_object.components[component]

    check_data = palint_object._get_ac_data()
    for i in check_data:
        ac_i = i[0][idx] / norm_factor[0]
        th_i = i[1][idx] / norm_factor[0]
        ptrm_j = i[2][idx] / norm_factor[1]
        th_j = i[3][idx] / norm_factor[1]

        vline = lines.Line2D([ptrm_j, ptrm_j], [th_j, th_i], color=colors[plt_idx], ls='--')
        hline = lines.Line2D([ac_i, ptrm_j], [th_i, th_i], color=colors[plt_idx], ls='--')
        ax.add_line(hline)
        ax.add_line(vline)
        ax.plot(ac_i, th_i, 's', markeredgecolor=colors[plt_idx], markerfacecolor='w')


def add_ck_check(palint_object, ax, component='m', norm=None,
                 norm_factor=[1, 1],
                 plt_idx=0, plt_opt={},
                 **options):
    colors = helper.get_colors()
    idx = palint_object.components[component]

    check_data = palint_object._get_ck_data()
    for i in check_data:
        ck_i = i[0][idx] / norm_factor[0]
        th_i = i[1][idx] / norm_factor[0]
        ptrm_j = i[2][idx] / norm_factor[1]
        th_j = i[3][idx] / norm_factor[1]

        hline = lines.Line2D([ck_i, ptrm_j], [th_j, th_j], color=colors[plt_idx])
        vline = lines.Line2D([ck_i, ck_i], [th_j, th_i], color=colors[plt_idx])
        ax.add_line(hline)
        ax.add_line(vline)
        ax.plot(ck_i, th_i, '^', markeredgecolor=colors[plt_idx], markerfacecolor='w')


def add_dunlop_labels(palint_object, ax, norm=None, norm_factor=1, plt_idx=0, text=False, **options):
    ax.set_xlabel('Temperature [$^\\circ C$]')

    ylabel = helper.get_moment_label(norm)
    ax.set_ylabel(ylabel)

    if text:
        if norm.lower() in ['max', 'trm']:
            add_norm_text(palint_object=palint_object, ax=ax, norm_factor=norm_factor, plt_idx=plt_idx)

    return ax


def add_norm_text(palint_object, ax, norm_factor=1, plt_idx=0, **options):
    ax.text(0.01, 0.99 - plt_idx * 0.03, '%i: 1= %.3e' % (plt_idx, norm_factor),
            horizontalalignment='left',
            verticalalignment='top',
            color='#808080',
            transform=ax.transAxes)
    return ax


def add_slope_text(palint_object, ax, slope=1, plt_idx=0, **options):
    ax.text(0.01, 0.01 + plt_idx * 0.04, '%i: $B_{\text{pal}}$ = %.2f' % (plt_idx, slope * -palint_object.lab_field),
            horizontalalignment='left',
            verticalalignment='bottom',
            color='#808080',
            transform=ax.transAxes)
    return ax