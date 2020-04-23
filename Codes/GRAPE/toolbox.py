import numpy as np
import scipy.linalg

def expm(A, method="eigen"):
        """
        Calculates the exponentials of a given matrix list
        :param A: List of matrices (2D ndarraays) [A1, A2, A3, ...]
        :return: List of matrices (2D ndarrays) [exp(A1), exp(A2), ...]
        """
        if method == "eigen":
            vals, vects = np.linalg.eig(A)
            return np.einsum('...ik, ...k, ...kj -> ...ij',
                             vects, np.exp(vals), np.linalg.inv(vects))
        if method == "direct":
            for i in range(len(A)):
                A[i] = scipy.linalg.expm(A[i])
            return A