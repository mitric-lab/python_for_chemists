#!/usr/bin/env python

### ANCHOR: import
import numpy as np
### ANCHOR_END: import


### ANCHOR: qr_algorithm_function
def qr_algorithm(a, eps, maxiter=5000):
    n = a.shape[0]

    a_k = a.copy()
    eigvecs = np.eye(n)
    
    for i in range(0, maxiter):
        q_k, r_k = np.linalg.qr(a_k)
        a_k = np.dot(r_k, q_k)
        eigvecs = eigvecs @ q_k
        if np.max(np.abs(np.tril(a_k, k=-1))) < eps:
            break
    if i == maxiter - 1:
        print('QR algorithm did not converge.')
        print('Max residual:', np.max(np.abs(np.tril(a_k, k=-1))))
    
    eigvals = np.diag(a_k)
    
    sort_idx = np.argsort(eigvals)
    eigvals = eigvals[sort_idx]
    eigvecs = eigvecs[:, sort_idx]

    return eigvals, eigvecs
### ANCHOR_END: qr_algorithm_function


### ANCHOR: test_qr_algorithm
NDIM = 10

a_mat = np.random.rand(NDIM, NDIM)
a_mat = a_mat + a_mat.T

eigvals, eigvecs = qr_algorithm(a_mat, 1e-8)
### ANCHOR_END: test_qr_algorithm

### ANCHOR: verify_results
eigvals_ref, eigvecs_ref = np.linalg.eigh(a_mat)

print(eigvals)
print(eigvals_ref)
assert np.allclose(eigvals, eigvals_ref)
assert np.allclose(np.abs(eigvecs.T @ eigvecs_ref), np.eye(NDIM))
### ANCHOR_END: verify_results

