# coding=utf-8
__author__ = 'mike'
def __get_M(self, step='th'):
    implemented = {'th': self.th,
                   'pt': self.pt,
                   'ck': self.ck,
                   'ac': self.ac,
                   'tr': self.tr,
                   'ptrm': self.ptrm,
                   'sum': self.sum}
    OUT = [np.sqrt(i[1] ** 2 + i[2] ** 2 + i[3] ** 2) for i in implemented[step]]
    return np.array(OUT, dtype=np.ndarray)

def __get_D(self, step='th'):
    """
    :Parameter:
       step : str [default = 'th']
            The paleomagnetic step
    :Return:

    """
    from math import degrees

    if step not in ['th', 'pt', 'ptrm', 'ac', 'ck', 'tr', 'sum']:
        print 'No such step: %s' % step
        return

    implemented = {'th': self.th,
                   'pt': self.pt,
                   'ck': self.ck,
                   'ac': self.ac,
                   'tr': self.tr,
                   'ptrm': self.ptrm,
                   'sum': self.sum}

    aux = [np.arctan2(i[2], i[1]) for i in implemented[step]]
    D = map(degrees, aux)
    D = np.array(D)

    for i in range(len(D)):
        if D[i] < 0:
            D[i] += 360
        if D[i] > 360:
            D[i] -= 360
    return D

def __get_I(self, step='th'):
    """
    Calculates the Inclination from a given step.

    :Parameter:
       step : str [default = 'th']
            The paleomagnetic step
    :Return:
       I : inclination Data

    Inclination is calculated with,

    .. math::

       I = \\tan^{-1} \\left( \\sqrt{\\frac{z}{x^2 + y^2} } \\right)
    """
    from math import degrees

    implemented = {'th': self.th,
           'pt': self.pt,
           'ck': self.ck,
           'ac': self.ac,
           'tr': self.tr,
           'ptrm': self.ptrm,
           'sum': self.sum}

    aux = [np.arctan2(i[3], np.sqrt(i[2] ** 2 + i[1] ** 2)) for i in implemented[step]]
    I = map(degrees, aux)
    I = np.array(I)

    return I
