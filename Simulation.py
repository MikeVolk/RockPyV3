from RockPyV3.Functions.general import rotate
from Structure.sample import Sample
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt


class Simulation(object):
    def __init__(self, parameters=None):
        if not parameters: parameters = {}
        self.parameters = parameters


class Sim_Thellier(Simulation):
    def __init__(self, blocking_params=None, aniso_params=None, errors=None, **options):

        if not blocking_params: blocking_params = {}
        if not aniso_params: aniso_params = {}
        if not errors: errors = {}
        super(Sim_Thellier, self).__init__()

        default_blocking_params = {'lab_field': 35.0,
                                   'paleo_field': 35.0,
                                   'mean_nrm': 510.0,
                                   'skewness_nrm': -10,
                                   'skewness_ptrm': -10,  # - skewed to left
                                   'w_nrm': 50.0,
                                   'w_ptrm': 50.0,
                                   'ptrm_t_diff': 2.0}

        default_aniso_params = {'atrm_matrix': np.array([[1, 0., 0.],
                                                         [0., 1., 0.],
                                                         [0., 0., 1.]]),
                                'lab_field_direction': np.array([0, 0, 1]),
                                'paleo_field_direction': np.array([0, 0, 1]),
        }

        default_errors = {'alignment_error': 0,
                          'field_error': 0,
                          'temperature_error': 0}
        default_steps = [20.0, 60.0, 100.0, 140.0, 185.0, 225.0, 265.0, 310.0, 350.0, 390.0, 430.0, 475.0, 515.0, 555.0,
                         600.0, 680.0]

        ''' OPTIONS '''
        for key in default_blocking_params:
            if key not in blocking_params:
                blocking_params[key] = default_blocking_params[key]

        for key in default_aniso_params:
            if key not in aniso_params:
                aniso_params[key] = default_aniso_params[key]

        for key in default_errors:
            if key not in errors:
                errors[key] = default_errors[key]

        # lists = {key: value for key, value in parameters.iteritems() if type(parameters[key]) == list}

        self.__dict__.update(blocking_params)
        self.__dict__.update(aniso_params)
        self.__dict__.update(errors)

        self.t_steps = options.get('t_steps', default_steps)


        # this are the put in values for plots oven etc. There is no error yet
        ''' zero field steps '''
        self.th_steps = self.t_steps

        self.ac_steps = [self.th_steps[i] for i in range(1, len(self.t_steps), 3)]
        self.tr_steps = [self.th_steps[i] for i in range(1, len(self.t_steps), 3)]

        ''' in field steps '''
        self.pt_steps = self.t_steps
        self.ck_steps = [self.th_steps[i] for i in range(1, len(self.t_steps), 2)]

        ''' adding temperature error '''

        if not self.temperature_error <= 0:
            self.th_temps = self.th_steps + np.random.normal(0, self.temperature_error, len(self.th_steps))
            self.ac_temps = self.ac_steps + np.random.normal(0, self.temperature_error, len(self.ac_steps))
            self.tr_temps = self.tr_steps + np.random.normal(0, self.temperature_error, len(self.tr_steps))

            self.pt_temps = self.pt_steps + np.random.normal(0, self.temperature_error, len(self.pt_steps))
            self.ck_temps = self.ck_steps + np.random.normal(0, self.temperature_error, len(self.ck_steps))

        else:  # error cannot be <=0 -> no arror added
            self.th_temps = self.th_steps
            self.ac_temps = self.ac_steps
            self.tr_temps = self.tr_steps

            self.pt_temps = self.pt_steps
            self.ck_temps = self.ck_steps

        ''' field errors '''

        if not self.field_error <= 0:
            self.ptrm_fields = np.random.normal(0, self.field_error, len(self.pt_steps)) + self.lab_field
            self.ck_fields = np.random.normal(0, self.field_error, len(self.ck_steps)) + self.lab_field
        else:
            self.ptrm_fields = np.ones(len(self.pt_steps)) * self.lab_field
            self.ck_fields = np.ones(len(self.ck_steps)) * self.lab_field

        ''' field errors '''

        if not self.alignment_error <= 0:
            self.ptrm_field_directions = [rotate(xyz=self.lab_field_direction, degree=i) for i in
                                          np.random.normal(0, self.alignment_error, len(self.pt_steps))]
            self.ck_field_directions = [rotate(xyz=self.lab_field_direction, degree=i) for i in
                                        np.random.normal(0, self.alignment_error, len(self.ck_steps))]
        else:
            self.ptrm_field_directions = np.array([self.lab_field_direction for i in range(len(self.pt_steps))])
            self.ck_field_directions = np.array([self.lab_field_direction for i in range(len(self.ck_steps))])

        ''' generate sample and measurements '''
        self.sample = Sample(name='Thellier Simulation')
        self.sample.add_measurement(mtype='palint', mfile='', machine='simulation')

        N = len(self.th_steps) + len(self.ck_steps) + len(self.pt_steps) + len(self.ac_steps) + len(self.tr_steps)

        ''' TH Steps '''
        self.nrm = np.dot(self.atrm_matrix, self.paleo_field_direction) * self.paleo_field

        m_nrm, m_ptrm = self.data()
        self.th = np.array([self.nrm * m_nrm[np.argmin(abs(i - m_nrm[:, 0])), 1] for i in self.th_temps])
        th_m = [np.linalg.norm(i) for i in self.th]
        self.th = np.c_[self.th_steps, self.th[:, 0], self.th[:, 1], self.th[:, 2], th_m]

        ''' PTRM '''
        m_ptrm_calc = np.array([m_ptrm[np.argmin(abs(i - m_ptrm[:, 0])), 1] for i in self.pt_temps])
        self.ptrm = np.array(
            [np.dot(self.atrm_matrix, self.ptrm_field_directions[i]) * self.ptrm_fields[i] * m_ptrm_calc[i]
             for i in range(len(self.ptrm_field_directions))])
        self.ptrm = np.c_[
            self.pt_steps, self.ptrm[:, 0], self.ptrm[:, 1], self.ptrm[:, 2], map(np.linalg.norm, self.ptrm)]

        m_ck_calc = np.array([m_ptrm[np.argmin(abs(i - m_ptrm[:, 0])), 1] for i in self.ck_temps])
        ck = np.array([np.dot(self.atrm_matrix, self.ck_field_directions[i]) * self.ck_fields[i] * m_ck_calc[i]
                       for i in range(len(self.ck_field_directions))])

        ck = [ck[i] + self.th[j, 1:4] for i in range(len(self.ck_steps)) for j in range(len(self.th_steps))
              if self.th_steps[j] == self.ck_steps[i]]

        ''' SUM '''

        self.sum = self.th[:, 1:4] + self.ptrm[:, 1:4]
        self.sum = np.c_[self.pt_steps, self.sum[:, 0], self.sum[:, 1], self.sum[:, 2], map(np.linalg.norm, self.sum)]

        # self.ptrm = np.array([])
        # plt.plot(self.th[:, 0], self.th[:, 4])
        # plt.plot(self.ptrm[:, 0], self.ptrm[:, 4])
        # plt.plot(self.sum[:, 0], self.sum[:, 4])
        # plt.plot([35, 0], [0, 35], '--')
        # plt.plot(self.ptrm[:, 4], self.th[:, 4])
        # plt.show()

    def skew(self, temps, mean=0, w=1, a=0):
        t = (temps - mean) / w
        return 2 * norm.pdf(t) * norm.cdf(a * t)


    def data(self, check=False):

        mean_ptrm = self.mean_nrm + self.ptrm_t_diff
        t = np.arange(0, 700, 0.1)

        p_nrm = self.skew(t, self.mean_nrm, self.w_nrm, self.skewness_nrm)
        p_ptrm = self.skew(t, mean_ptrm, self.w_ptrm, self.skewness_ptrm)

        p_nrm /= np.sum(p_nrm)
        p_ptrm /= np.sum(p_ptrm)

        m_nrm = np.array([1 - np.sum(p_nrm[:i]) for i in range(len(t))])
        m_ptrm = np.array([np.sum(p_ptrm[:i]) for i in range(len(t))])

        if check:
            plt.plot(m_nrm)
            plt.plot(m_ptrm)
            plt.show()
        m_nrm = np.c_[t, m_nrm]
        m_ptrm = np.c_[t, m_ptrm]
        return m_nrm, m_ptrm


# todo outputs a sample object

if __name__ == '__main__':
    test = Sim_Thellier(blocking_params={'ptrm_t_diff': 0, 'lab_field': 35.0, 'w_ptrm': 50},
                        errors={'field_error': 0.0, 'temperature_error': 0.0, 'alignment_error': 0.0},
                        t_steps=np.arange(0, 700, 50))