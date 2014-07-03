__author__ = 'mike'
import numpy as np


class tensor():
    '''

    '''

    def __init__(self, tensor):
        '''

        '''
        self.tensor = tensor

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
        self.norm_evals = norm_evals

        if norm:
            return self.norm_evals, self.evecs
        else:
            return self.evals, self.evecs

    def eigenvalues(self, norm=False):
        self.get_eigensystem()
        if norm:
            return self.norm_evals
        else:
            return self.evals

    def eigenvectors(self):
        self.get_eigensystem()