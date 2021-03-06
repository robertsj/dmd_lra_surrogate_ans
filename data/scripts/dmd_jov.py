"""
Derived module from dmdbase.py for classic dmd.
"""
import numpy as np
import scipy as sp
from pydmd import DMDBase


class DMD_jov(DMDBase):
    """
    Dynamic Mode Decomposition

    :param svd_rank: the rank for the truncation; If 0, the method computes the
        optimal rank and uses it for truncation; if positive interger, the
        method uses the argument for the truncation; if float between 0 and 1,
        the rank is the number of the biggest singular values that are needed
        to reach the 'energy' specified by `svd_rank`; if -1, the method does
        not compute truncation.
    :type svd_rank: int or float
    :param int tlsq_rank: rank truncation computing Total Least Square. Default
        is 0, that means TLSQ is not applied.
    :param bool exact: flag to compute either exact DMD or projected DMD.
        Default is False.
    :param bool opt: flag to compute optimized DMD. Default is False.
    """
    @staticmethod
    def _compute_amplitudes(modes,snapshots,eigs,opt, s, V,W):
        """
        Compute the amplitude coefficients. If `opt` is False the amplitudes
        are computed by minimizing the error between the modes and the first
        snapshot; if `opt` is True the amplitudes are computed by minimizing
        the error between the modes and all the snapshots, at the expense of
        bigger computational cost.
        :param numpy.ndarray modes: 2D matrix that contains the modes, stored
            by column.
        :param numpy.ndarray snapshots: 2D matrix that contains the original
            snapshots, stored by column.
        :param numpy.ndarray eigs: array that contains the eigenvalues of the
            linear operator.
        :param bool opt: flag for optimized dmd.
        :return: the amplitudes array
        :rtype: numpy.ndarray
        """
        if opt =='Jov':
            
            Vand = np.zeros((s.shape[0], V.shape[0]),dtype=complex); # Vandermonde matrix
            for k in range(V.shape[0]):
                Vand[:, k] = eigs**(k)
             
    
        # the next 5 lines follow Jovanovic et al, 2014 code:
            G = np.diag(s).dot( V.conj().T)
            P = (W.conj().T.dot(W))*(Vand.dot(Vand.conj().T)).conj()
            q = (np.diag(Vand.dot(G.conj().T).dot(W))).conj()
            Pl = sp.linalg.cholesky(P,lower=True)
            b = np.linalg.inv(Pl.conj().T).dot((np.linalg.inv(Pl)).dot(q)) # Optimal vector of amplitudes b
        elif opt == True:
            L = np.concatenate(
                [
                    modes.dot(np.diag(eigs**i))
                    for i in range(snapshots.shape[1])
                ],
                axis=0)
            a = np.reshape(snapshots, (-1, ), order='F')

            b = np.linalg.lstsq(L, a)[0]
        elif opt == False:
            b = np.linalg.lstsq(modes, snapshots.T[0])[0]
        else:
            print('opt must be True, False, or jov')
            return
        return b

    def fit(self, X):
        """
        Compute the Dynamic Modes Decomposition to the input data.

        :param X: the input snapshots.
        :type X: numpy.ndarray or iterable
        """
        self._snapshots, self._snapshots_shape = self._col_major_2darray(X)

        n_samples = self._snapshots.shape[1]
        X = self._snapshots[:, :-1]
        Y = self._snapshots[:, 1:]

        X, Y = self._compute_tlsq(X, Y, self.tlsq_rank)

        U, s, V = self._compute_svd(X, self.svd_rank)

        self._Atilde = self._build_lowrank_op(U, s, V, Y)

        self._eigs, self._modes = self._eig_from_lowrank_op(
            self._Atilde, Y, U, s, V, self.exact)

        _, lowrank_eigenvectors = np.linalg.eig(self._Atilde)
        
        self._b = self._compute_amplitudes(self._modes, self._snapshots,
                                           self._eigs, self.opt,s,V,lowrank_eigenvectors)

        # Default timesteps
        self.original_time = {'t0': 0, 'tend': n_samples - 1, 'dt': 1}
        self.dmd_time = {'t0': 0, 'tend': n_samples - 1, 'dt': 1}

        return self
