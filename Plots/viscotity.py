__author__ = 'mike'

def plot_visc(visc_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf', plt_opt={}):
    '''
    Main plotting of a hysteresis.
    :param visc_obj: hysteresis object

    :param ax: pyplot.ax object
    :param norm_factor: normalization of the y-data
    :param out: choice of:

       - 'show' : shows the plot
       - 'rtn' : returns the plot
       - 'save' : saves a pdf of the plot. a folder has to be specified. Name can be specified
    :param folder: str
       where the pdf should be saved
    :param name: str
       name of the pdf

    implemented plot_options:
        color
        linestyle
        markerstyle
        alpha
    '''

    # ax.axhline(0, color='#555555') # 0 line horizontally
    # ax.axvline(0, color='#555555') # 0 line vertically

    # color_cycle = ax._get_lines.color_cycle
    #
    # color = popt.get('color', next(color_cycle))
    # ls = popt.get('linestyle', '-')
    # markerstyle = popt.get('markerstyle', 'x')
    # alpha = popt.get('alpha', 1.0)

    std, = ax.plot(visc_obj.times, visc_obj.moments / norm_factor, label='t',**plt_opt)


    if out == 'show':
        plt.show()
    if out == 'rtn':
        return ax
    if out == 'save':
        if folder != None:
            plt.savefig(folder + self.samples[0].name + '_' + name, dpi=300)


def plot_log_visc(visc_obj, ax, norm_factor=1, out='show', folder=None, name='output.pdf', plt_opt={}):
    '''
    Main plotting of a hysteresis.
    :param visc_obj: hysteresis object

    :param ax: pyplot.ax object
    :param norm_factor: normalization of the y-data
    :param out: choice of:

       - 'show' : shows the plot
       - 'rtn' : returns the plot
       - 'save' : saves a pdf of the plot. a folder has to be specified. Name can be specified
    :param folder: str
       where the pdf should be saved
    :param name: str
       name of the pdf
    '''

    # ax.axhline(0, color='#555555') # 0 line horizontally
    # ax.axvline(0, color='#555555') # 0 line vertically

    std, = ax.plot(visc_obj.log_times, visc_obj.moments / norm_factor, label='log(t)', **plt_opt)


    if out == 'show':
        plt.show()
    if out == 'rtn':
        return ax
    if out == 'save':
        if folder != None:
            plt.savefig(folder + self.samples[0].name + '_' + name, dpi=300)


def add_label(visc_object, ax=None, log=True):
    if visc_object.norm == None:
        ax.set_ylabel('Magnetic Moment $Am^2$')
    if visc_object.norm == 'mass':
        ax.set_ylabel('Magnetic Moment $Am^2/kg$')
    else:
        ax.set_ylabel('Magnetic Moment normalized')

    if log:
        ax.set_xlabel('log(time [$s$])')
    else:
        ax.set_xlabel('time [$s$]')
    return ax


def add_formula_text(visc_obj,ax = None, norm_factor=1):
    text = '$M = %.3e \cdot x %.3e \; R^2 = %.3e \; \sigma = %.3e$' % (visc_obj.intercept, visc_obj.slope, visc_obj.r_value**2, visc_obj.std_err)
    ax.text(0, 0.01,
            text,
            horizontalalignment='left',
            verticalalignment='bottom',
            transform=ax.transAxes,
            # fontsize=8,
    )
    return ax