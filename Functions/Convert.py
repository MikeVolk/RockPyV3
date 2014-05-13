__author__ = 'Mike'
import logging


def convert2(in_unit, out_unit, unit_type):
    """

    :param in_unit: str the input unit
    :param out_unit: str the desired outpu unit
    :param unit_type: what kind of units they are.
    :return: conversion factor

    :example:

        assume you want to convert a samples mass from *mg* to *kg*

        ``sample_mass = 0.00345  # sample weight in kg``

        ``sample_mass *= convert2('kg', 'mg', 'mass')  # gives new mass in mg``

        ``>> sample_mass = 3450.0``

    :implemented units:

         * 'volume' [only metric units]
         * 'pressure'
         * 'length' [only metric units]
         * 'mass' [only metric units]
         * 'area' [only metric units]

    """
    conversion_table = {'mass':
                            {'T': 1E-9,
                             'kg': 1E-6,
                             'g': 1E-3,
                             'mg': 1,
                             'mug': 1E3,
                             'ng': 1E6},
                        'length':
                            {'km': 1E-6,
                             'm': 1E-3,
                             'dm': 1E-2,
                             'cm': 1E-1,
                             'mm': 1,
                             'micron': 1E3,
                             'mum': 1E6,
                             'nm': 1E9},
                        'area':
                            {'km2': 1E-12,
                             'm2': 1E-6,
                             'dm2': 1E-4,
                             'cm2': 1E-2,
                             'mm2': 1,
                             'micron2': 1E6,
                             'nm2': 1E12},
                        'volume':
                            {'km3': 1E-18,
                             'm3': 1E-9,
                             'dm3': 1E-6,
                             'cm3': 1E-3,
                             'mm3': 1,
                             'micron3': 1E9,
                             'nm3': 1E18},
                        'pressure':
                            {'GPa': 0.001,
                             'MPa': 1.0,
                             'N/cm2': 1000.0,
                             'N/m2': 1000000.0,
                             'N/mm2': 1.0,
                             'Pa': 1000000.0,
                             'TPa': 1e-06,
                             'at': 10.197162129779,
                             'atm': 9.86923266716,
                             'bar': 10.0,
                             'cmHg': 750.063755419211,
                             'dyn/cm2': 10000000.0,
                             'hPa': 10000.0,
                             'kN/m2': 1000.0,
                             'kPa': 1000.0,
                             'kgf/cm2': 10.197162129779,
                             'kgf/m2': 101971.621297793,
                             'kgf/mm2': 0.101971621298,
                             'lbf/ft2': 20885.434233120002,
                             'lbf/in2': 145.03773773,
                             'mbar': 10000.0,
                             'mmHg': 7500.63755419211,
                             'mubar': 10000000.0,
                             'psi': 145.03773773,
                             'torr': 7500.61682704},
    }
    log = logging.getLogger('RockPy.HELPER_FUNCTIONS')

    if unit_type not in conversion_table:
        log.error('UNKNOWN\t << unit_type >>')

    conversion = conversion_table[unit_type]
    if in_unit not in conversion:
        log.error('UNKNOWN\t << in_unit >>')
        return
    if out_unit not in conversion:
        log.error('UNKNOWN\t << out_unit >>')
        return

    factor = conversion[out_unit] / conversion[in_unit]

    log.debug('CONVERSION %s -> %s : FACTOR = %.2e' % (in_unit, out_unit, factor))

    return factor


def conversion(unit):
    '''
    This is a conversion table. Ask for the unit like:
    conversion['mass'] returns the appropriate conversion factor

    :param unit: str
    '''

    conversion_table = {'mass':
                            {'T': 1E-9,
                             'kg': 1E-6,
                             'g': 1E-3,
                             'mg': 1,
                             'mug': 1E3,
                             'ng': 1E6},
                        'length':
                            {'km': 1E-6,
                             'm': 1E-3,
                             'dm': 1E-2,
                             'cm': 1E-1,
                             'mm': 1,
                             'micron': 1E3,
                             'mum': 1E6,
                             'nm': 1E9},
                        'area':
                            {'km2': 1E-12,
                             'm2': 1E-6,
                             'dm2': 1E-4,
                             'cm2': 1E-2,
                             'mm2': 1,
                             'micron2': 1E6,
                             'nm2': 1E12},
                        'volume':
                            {'km3': 1E-18,
                             'm3': 1E-9,
                             'dm3': 1E-6,
                             'cm3': 1E-3,
                             'mm3': 1,
                             'micron3': 1E9,
                             'nm3': 1E18},
                        'pressure':
                            {'GPa': 0.001,
                             'MPa': 1.0,
                             'N/cm2': 1000.0,
                             'N/m2': 1000000.0,
                             'N/mm2': 1.0,
                             'Pa': 1000000.0,
                             'TPa': 1e-06,
                             'at': 10.197162129779,
                             'atm': 9.86923266716,
                             'bar': 10.0,
                             'cmHg': 750.063755419211,
                             'dyn/cm2': 10000000.0,
                             'hPa': 10000.0,
                             'kN/m2': 1000.0,
                             'kPa': 1000.0,
                             'kgf/cm2': 10.197162129779,
                             'kgf/m2': 101971.621297793,
                             'kgf/mm2': 0.101971621298,
                             'lbf/ft2': 20885.434233120002,
                             'lbf/in2': 145.03773773,
                             'mbar': 10000.0,
                             'mmHg': 7500.63755419211,
                             'mubar': 10000000.0,
                             'psi': 145.03773773,
                             'torr': 7500.61682704},
    }

    if unit == 'all':
        return conversion_table
    else:
        if not unit in conversion_table:
            print('unit not recognized')
            return
        else:
            return conversion_table[unit]
