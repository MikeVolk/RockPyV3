__author__ = 'mike'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, lines
from scipy import stats
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


def arai(palint_object, ax, component='m', norm=None, norm_factor=1, plt_idx=0, out='show', folder=None,
         name='output.pdf', line=True, plt_opt={}):
    idx = palint_object.components[component]
    markers = helper.get_marker()
    xy = np.array([[i[idx], j[idx]] for i in palint_object.ptrm for j in palint_object.th if i[0] == j[0]])

    if norm == 'th':
        norm_factor_th = palint_object.th[0, idx]
        norm_factor_ptrm = palint_object.th[0, idx]
    if norm == None:
        norm_factor_th = 1
        norm_factor_ptrm = 1

    if not 'ls' in plt_opt or 'linestyle' in plt_opt:
        plt_opt.update({ls:''})
    ax.plot(xy[:, 0] / norm_factor_ptrm, xy[:, 1] / norm_factor_th, marker=markers[plt_idx], **plt_opt)

    if line:
        add_arai_line(palint_object=palint_object, ax=ax,
                      plt_idx=plt_idx,
                      norm_factor_th=norm_factor_th, norm_factor_ptrm=norm_factor_ptrm,
                      component=component,
                      plt_opt=plt_opt)
    return ax


def add_arai_line(palint_object, ax, component='m', norm=None, norm_factor_th=1, norm_factor_ptrm=1, plt_idx=0,
                  out='show', folder=None, name='output.pdf', plt_opt={}):
    idx = palint_object.components[component]
    xy = np.array([[i[idx], j[idx]] for i in palint_object.ptrm for j in palint_object.th if i[0] == j[0]])
    slope, intercept, r_value, p_value, std_err = stats.linregress(xy[:, 0], xy[:, 1])
    x_fit = np.linspace(0, max(palint_object.th[:, idx]) / norm_factor_th, 2)
    y_fit = (intercept + slope * x_fit) / norm_factor_ptrm

    ax.plot(x_fit, y_fit, color='#500000')
    add_slope_text(palint_object=palint_object, ax=ax, slope=slope, plt_idx=plt_idx)


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