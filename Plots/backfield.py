__author__ = 'mike'
from matplotlib import rc, lines
import numpy as np


def plot_coe(coe_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    # todo implement out
    if coe_obj is not None:
        ax.plot(coe_obj.remanence_interpolated[:, 0], coe_obj.remanence_interpolated[:, 1] / norm_factor, '-',
                color='k')
        ax.plot(coe_obj.remanence[:, 0], coe_obj.remanence[:, 1] / norm_factor,
                'x',
                markersize=3,
                color='k')


def add_bcr_line(coe_obj, ax, norm_factor=1, method='fit', text=False):
    implemented = {'fit': coe_obj.bcr,
                   'interpolated': coe_obj.bcri}

    if method in implemented:
        data = implemented[method]

    Bcr_line = lines.Line2D([-data, -data], [-1 / norm_factor, 1 / norm_factor], lw=1.5, ls='--', color='r', alpha=0.5)
    ax.add_line(Bcr_line)

    if text:
        add_bcr_text(coe_obj=coe_obj, ax=ax, norm_factor=norm_factor, method=method)

    return ax


def add_bcr_text(coe_obj, ax, norm_factor=1, method='fit'):
    '''
    Adds a line to an ax object, where Bcr is
    :param coe_obj:
    :param ax:
    :param norm_factor:
    :param method:
    '''
    implemented = {'fit': coe_obj.bcr,
                   'interpolated': coe_obj.bcri}

    if method in implemented:
        data = implemented[method]

    ax.text(-data, - 2 * max(coe_obj.remanence[:, 1]) / norm_factor,
            '$B_{cr}=%.2f mT$' % (data * 1000),
            horizontalalignment='left',
            verticalalignment='bottom',
            fontsize=8,
    )
    return ax


def plot_henkel(coe_obj, irm_obj, ax,
                norm_factor=[1, 1],
                plt_opt={}, **options):
    '''
    Simple function that plots a henkel_plot

    :param coe_obj: Measurement.coe object
    :param irm_obj: Measurement.irm object
    :param ax: ax
    :param norm_factor:
    :param plt_opt: plotting options e.g. linestyle, alpha
    :param options:
    :return:
    '''
    s = (irm_obj.rem[0] - irm_obj.rem[-1]) / (coe_obj.rem[0] - coe_obj.rem[-1])
    x = coe_obj.rem
    y = x * s + irm_obj.rem[-1] / 2

    ax.plot(coe_obj.rem / norm_factor[0], irm_obj.rem / norm_factor[1], '.-', color='#505050')
    ax.plot(x / norm_factor[0], y/ norm_factor[1], color='r')
    ax.fill_between(coe_obj.rem/ norm_factor[0], y/ norm_factor[1], irm_obj.rem/ norm_factor[1], color='#808080', alpha=0.2)

    # ax.set_xlim([min(coe_obj.rem)/ norm_factor[0], max(coe_obj.rem)/ norm_factor[0]])
    # ax.set_ylim([min(irm_obj.rem)/ norm_factor[1], max(irm_obj.rem)/ norm_factor[1]])
    return ax
