__author__ = 'mike'


def get_moment_label(norm_factor):
    # if not norm_factor:
    # norm_factor = 'none'

    texts = {'max': 'magnetic moment normalized',
             'mass': 'magnetic moment [$\\frac{Am^2}{kg}$]',
             'None': 'magnetic moment [$Am^2$]',
             'trm': 'magnetic moment normalized',
             'nrm': 'magnetic moment normalized',
    }

    out = texts[norm_factor]
    return out


def get_marker():
    marker = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', None]
    return marker


def get_colors(version=1):
    if version == 1:
        out = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
    else:
        out = [
            "#8DD593", "#C6DEC7", "#EAD3C6",
            "#023FA5", "#7D87B9", "#BEC1D4",
            "#8595E1", "#B5BBE3", "#E6AFB9",
            "#E07B91", "#D33F6A", "#11C638",
            "#F0B98D", "#EF9708", "#0FCFC0",
            "#9CDED6", "#D5EAE7", "#F3E1EB",
            "#D6BCC0", "#BB7784", "#4A6FE3",
            "#F6C4E1", "#F79CD4"
        ]
    v3 = ['#F0A3FF', '#0075DC', '#993F00', '#4C005C',
          '#191919', '#005C31', '#2BCE48', '#FFCC99',
          '#808080', '#94FFB5', '#8F7C00', '#9DCC00',
          '#C20088', '#003380', '#FFA405', '#FFA8BB',
          '#426600', '#FF0010', '#5EF1F2', '#00998F',
          '#E0FF66', '#740AFF', '#990000', '#FFFF80',
          '#FFFF00', '#FF5005']
    out *= 10

    return out