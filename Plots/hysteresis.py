import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, lines

def plot_hys(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    std, = ax.plot(hys_obj.up_field[:, 0], hys_obj.up_field[:, 1] / norm_factor,
                   color='k')
    ax.plot(hys_obj.down_field[:, 0], hys_obj.down_field[:, 1] / norm_factor,
            color=std.get_color())

    if out == 'show':
        plt.show()
    if out == 'rtn':
        return ax


def plot_virgin(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    ''' VIRGIN '''
    # todo out options

    if hys_obj.virgin != None:
        ax.plot(hys_obj.virgin[:, 0], hys_obj.virgin[:, 1] / norm_factor, color='#808080')
    if out == 'show':
        plt.show()
    if out == 'rtn':
        return ax


def plot_mrs_shift(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    # todo out options
    MRS = hys_obj.mrs[0]
    M_B_Mrs = np.array([[i[0], (i[1] + MRS) / norm_factor] for i in hys_obj.up_field_interpolated if i[0] <= 0])
    ax.plot(M_B_Mrs[:, 0], M_B_Mrs[:, 1], '--', color='k')


def plot_rev(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    # todo out options
    ax.plot(hys_obj.rev[:, 0], hys_obj.rev[:, 1] / norm_factor, 'k--', label='irreversible')
    return ax


def plot_irrev(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    # todo out options
    ax.plot(hys_obj.irr[:, 0], hys_obj.irr[:, 1] / norm_factor, 'k--', label='irreversible')
    return ax


def plot_irrev_assymetry(hys_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    '''
    '''
    pos = np.array([[abs(i[0]), i[1]] for i in hys_obj.irr if i[0] > 0])
    neg = np.array([[abs(i[0]), i[1]] for i in hys_obj.irr if i[0] <= 0])
    # diff =
    ax.plot(pos[:,0], pos[:,1], '-', label='pos. fields')
    ax.plot(neg[:,0], neg[:,1], '-', label='neg. fields')

    if out == 'show':
        plt.legend()
        plt.show()



def add_virgin_info(hys_obj, ax, norm_factor=1):
    if hys_obj.virgin != None:
        mn = ax.get_xlim()[0]
        mn -= 0.2 * mn

        ax.axhline(y=hys_obj.virgin[0, 1] / norm_factor, xmax=0.15, color='#808080')
        ax.axhline(y=hys_obj.mrs[0] / norm_factor, xmax=0.15, color='#808080')

        ax.axhspan(hys_obj.virgin[0, 1] / norm_factor, hys_obj.mrs[0] / norm_factor,
                   facecolor='0', alpha=0.1, xmax=0.1)

        ax.text(mn, hys_obj.mrs[0] / norm_factor, '$M_{rs}=%.2e$' % (hys_obj.mrs[0] / norm_factor),
                horizontalalignment='left',
                verticalalignment='bottom',
                fontsize=8,
        )
        ax.text(mn, hys_obj.virgin[0, 1] / norm_factor, '$M_{si}(0)=%.2e$' % (hys_obj.virgin[0, 1] / norm_factor),
                horizontalalignment='left',
                verticalalignment='bottom',
                fontsize=8,
        )
        return ax


# ### FILLS


def fill_hys(hys_obj, ax=None, norm_factor=1):
    ax.fill_between(hys_obj.up_field[:, 0], hys_obj.up_field[:, 1] / norm_factor,
                    hys_obj.down_field[:, 1] / norm_factor,
                    color='#555555', alpha=0.1)
    return ax


def fill_virgin(hys_obj, ax=None, norm_factor=1):
    if hys_obj.virgin != None:
        vlen = len(hys_obj.virgin)

        ax.fill_between(hys_obj.virgin[:, 0][::-1], hys_obj.down_field[0:vlen, 1] / norm_factor,
                        hys_obj.virgin[:, 1][::-1] / norm_factor,
                        color='#555555', alpha=0.2)


####   LINES

def add_brh_line(hys_obj, ax, norm_factor=1):
    ''' Brh '''


    Brh = hys_obj.Brh()
    YMX = max(hys_obj.up_field_interpolated[:, 1])

    Brh_line = lines.Line2D([[-Brh[0], -Brh[0]]], [min(hys_obj.down_field[:, 1]) / norm_factor,
                                                   max(hys_obj.down_field[:, 1]) / norm_factor],
                            lw=1.5, ls='-', color='#808080', alpha=0.5)
    ax.add_line(Brh_line)

    ax.plot(-Brh[0], -Brh[1] / norm_factor, ls='.', color='r')
    return ax


def add_bcr_line(hys_obj,ax, norm_factor = 1, text=False):
    Bcr_line = lines.Line2D([[self.coe.bcri, self.coe.bcri]],
                            [-max(hys_obj.down_field[:, 1]) / norm_factor, 0],
                            lw=1., ls=':', color='k', alpha=0.5)
    ax.add_line(Bcr_line)
    return ax


def add_mdf_line(hys_obj,ax, norm_factor = 1):
    if ax == None:
        ax = self.ax
    MDF_line = lines.Line2D([[-2, 2]], [0.5 * max(hys_obj.irr[:, 1] / norm_factor),
                                        0.5 * max(hys_obj.irr[:, 1] / norm_factor)],
                            lw=1., ls='--', color='k', alpha=0.5)
    ax.add_line(MDF_line)
    return ax


def add_05ms_line(hys_obj,ax, norm_factor = 1, text = False):

    MDF_line = lines.Line2D([[-2, 2]], [0.5 * max(hys_obj.rev[:, 1] / norm_factor),
                                        0.5 * max(hys_obj.rev[:, 1] / norm_factor)],
                            lw=1., ls='--', color='k', alpha=0.5)
    ax.add_line(MDF_line)

    if text:
        add_05ms_text(hys_obj=hys_obj ,ax=ax , norm_factor = norm_factor)
    return ax


def add_ms_line(hys_obj,ax, norm_factor = 1, text = False):

    MS_line = lines.Line2D([[-2, 2]], [max(hys_obj.rev[:, 1] / norm_factor),
                                       max(hys_obj.rev[:, 1] / norm_factor)],
                           lw=1., ls='--', color='k', alpha=0.5)
    ax.add_line(MS_line)

    if text:
        add_ms_text(hys_obj=hys_obj, ax=ax, norm_factor=norm_factor)
    return ax

####   TEXT

def add_brh_text(hys_obj,ax):

    Brh = hys_obj.Brh()
    ax.text(-Brh[0], 0, '$B_{rh}=%.1f mT$' % (Brh[0] * 1000),
            horizontalalignment='right',
            verticalalignment='bottom',
            fontsize=8,
    )
    return ax


def add_bcr_text(hys_obj,ax):

    ax.text(self.coe.bcri, -max(hys_obj.down_field[:, 1]) / norm_factor,
            '$B_{cr}=%.2f mT$' % (self.coe.bcri * 1000),
            horizontalalignment='left',
            verticalalignment='bottom',
            fontsize=8,
    )
    return ax


def add_05ms_text(hys_obj,ax = None, norm_factor=1):
    ax.text(min(hys_obj.down_field[:, 0]), 0.5 * hys_obj.ms[0] / norm_factor,
            '$\\frac{1}{2}M_{s}=%.2e$' % (0.5 * hys_obj.ms[0] / norm_factor),
            horizontalalignment='left',
            verticalalignment='bottom',
            # fontsize=8,
    )
    return ax


def add_ms_text(hys_obj,ax = None, norm_factor=1):

    ax.text(min(hys_obj.down_field[:, 0]), hys_obj.ms[0] / norm_factor,
            '$M_{s}=%.2e$' % (hys_obj.ms[0] / norm_factor),
            horizontalalignment='left',
            verticalalignment='bottom',
            # fontsize=8,
    )
    return ax
