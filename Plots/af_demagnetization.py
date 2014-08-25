__author__ = 'mike'


def af_plot(af_demag_obj, ax,
            component='m', norm_factor=[1, 1],
            out='show', folder=None, name='output.pdf',
            plt_opt=dict(),
            **options):
    if af_demag_obj is not None:
        x = getattr(af_demag_obj.data, 'variable')
        y = getattr(af_demag_obj.data, component)
        ax.plot(x, y / norm_factor[1], '-',
                **plt_opt)
    return ax

def af_diff_plot(af_demag_obj, ax,
            component='m', norm_factor=[1, 1],
            out='show', folder=None, name='output.pdf',
            plt_opt=dict(),
            **options):
    strength = options.get('strength', 3)
    print norm_factor
    if af_demag_obj is not None:
        data_object = af_demag_obj.data.diff(strength=strength)
        x = getattr(data_object, 'variable')
        y = getattr(data_object, component)
        ax.plot(x, y / norm_factor[1], '-',
                **plt_opt)
    return ax