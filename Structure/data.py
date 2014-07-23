__author__ = 'mike'
import numpy as np


class tensor():
    '''

    '''

    def __init__(self, tensor):
        '''

        '''
        self.tensor = tensor
        self.get_eigensystem()


    def get_eigensystem(self, norm=False):
        '''
        :Parameters:
           norm : boolean
              if True :math:`\\sum_i \\tau_i =1`


        '''
        evals, evecs = np.linalg.eigh(self.tensor)
        idx = np.argsort(-evals)
        self.evals = evals[idx]
        evecs = evecs[:, idx]
        norm_evals = [i / sum(evals) for i in evals]
        self.evecs = np.array(evecs)
        self.evals_norm = norm_evals
