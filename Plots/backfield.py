__author__ = 'mike'
from matplotlib import rc, lines


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