__author__ = 'mike'


def get_moment_label(nom_factor):
    # if not norm_factor:
    # norm_factor = 'none'

    texts = {'max': 'magnetic moment normalized',
             'mass': 'magnetic moment [$\\frac{Am^2}{kg}$]',
             None: 'magnetic moment [$Am^2$]',
             'trm': 'magnetic moment normalized',
    }

    out = texts[nom_factor]
    return out


def get_marker():
    marker = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H']
    return marker


def get_colors(version=2):
    if version == 1:
        out = ['#FF0000', '#E4E400', '#00FF00', '#00FFFF', '#B0B0FF', '#FF00FF', '#E4E4E4',
               '#B00000', '#BABA00', '#00B000', '#00B0B0', '#8484FF', '#B000B0', '#BABABA',
               '#870000', '#878700', '#008700', '#008787', '#4949FF', '#870087', '#878787',
               '#550000', '#545400', '#005500', '#005555', '#0000FF', '#550055', '#545454']
    else:
        out = ["#023FA5", "#7D87B9", "#BEC1D4", "#D6BCC0", "#BB7784", "#4A6FE3", "#8595E1", "#B5BBE3", "#E6AFB9",
               "#E07B91", "#D33F6A", "#11C638", "#8DD593", "#C6DEC7", "#EAD3C6", "#F0B98D", "#EF9708", "#0FCFC0",
               "#9CDED6", "#D5EAE7", "#F3E1EB", "#F6C4E1", "#F79CD4"]
    out *= 10

    return out