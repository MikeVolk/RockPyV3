__author__ = 'mike'

def get_moment_label(nom_factor):
    # if not norm_factor:
    #     norm_factor = 'none'

    texts = {'max': 'magnetic moment normalized',
             'mass': 'magnetic moment [$\\frac{Am^2}{kg}$]',
             None: 'magnetic moment [$Am^2$]',
             'trm':'magnetic moment normalized',
             }

    out = texts[nom_factor]
    return out

def get_marker():
    marker = ['o','v','^','<','>','1','2', '3','4','8','s','p','*','h','H']
    return marker