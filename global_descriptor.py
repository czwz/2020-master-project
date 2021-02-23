import functools
import numpy as np

def sqeuclidean_distances(XA, XB):
    """
        Fast computation of squared euclidean distances
        between all samples in XA and XB (for the Gaussian kernel)
    """
    
    XA2 = np.sum(XA**2, axis=1).reshape((-1, 1))
    XB2 = np.sum(XB**2, axis=1).reshape((1, -1))
    D = XA2 + XB2 - 2*np.matmul(XA, XB.T)
    return D

def kernel_decorator(kernel_func):
    """
        A decorator for the kernel functions that
        does the structure averaging
    """
    
    @functools.wraps(kernel_func)
    def kernel_wrapper(XA, XB, **kwargs):
        """
            Wrapper for the kernel functions.
            If XA and/or XB is a list, it is assumed that each list
            entry contains the feature matrix 
            of local atomic environments for a single structure,
            and the kernel values are averaged over the structures.
            If only one of XA or XB is a list, then we perform
            the averaging so that the end result is a kernel
            between a set of structures and a set of local atom-based features.
            If both XA and XB are arrays, no summation is performed
            and the kernel is that between the local atom-based features
            instead of global features.
        """
        
        K = np.zeros((len(XA), len(XB)))
        
        # XA and XB are lists of structures
        if isinstance(XA, list) and isinstance(XB, list):
            for adx, a in enumerate(XA):
                for bdx, b in enumerate(XB):
                    K[adx, bdx] = np.mean(kernel_func(a, b, **kwargs))
                    
        # XA is a list of structures, XB is an array of atomic environments
        elif isinstance(XA, list):
            for adx, a in enumerate(XA):
                K[adx, :] = np.mean(kernel_func(a, XB, **kwargs), axis=0)
                
        # XA is an array of atomic environments, XB is a list of structures
        elif isinstance(XB, list):
            for bdx, b in enumerate(XB):
                K[:, bdx] = np.mean(kernel_func(XA, b, **kwargs), axis=1)
                
        # XA and XB are arrays of atomic environments
        else:
            K = kernel_func(XA, XB, **kwargs)                 
        return K
    return kernel_wrapper

@kernel_decorator
def gaussian_kernel(XA, XB, gamma=1.0):
    """
        Computes a Gaussian kernel between two matrices XA and XB
        with samples as rows and features as columns
    """
    D = sqeuclidean_distances(XA, XB)
    K = np.exp(-gamma*D)
    return K

@kernel_decorator
def linear_kernel(XA, XB, zeta=1):
    """
        Computes a linear kernel between two matrices XA and XB
        with samples as rows and features as columns
    """
    K = np.matmul(XA, XB.T)**zeta
    return K

def center_kernel(K, K_ref=None):
    """
        Centers a kernel matrix K relative to a reference kernel K_ref
    """
    if K_ref is None:
        K_ref = K
        
    if K.shape[1] != K_ref.shape[0] or K_ref.shape[0] != K_ref.shape[1]:
        print("Error: kernels must have compatible shapes " \
              "and the reference kernel must be square")
    else:
        col_mean = np.mean(K_ref, axis=0)
        row_mean = np.reshape(np.mean(K, axis=1), (-1, 1))
        k_mean = np.mean(K_ref)
        
        Kc = K - row_mean - col_mean + k_mean
        
        return Kc
