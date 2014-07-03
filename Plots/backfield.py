__author__ = 'mike'

def plot_coe(coe_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf'):
    #todo implement out
    if coe_obj != None:
        ax.plot(coe_obj.remanence_interpolated[:, 0], coe_obj.remanence_interpolated[:, 1] / norm_factor, '-',
                color='k')
        ax.plot(coe_obj.remanence[:, 0], coe_obj.remanence[:, 1] / norm_factor,
                'x',
                markersize = 3,
                color='k')