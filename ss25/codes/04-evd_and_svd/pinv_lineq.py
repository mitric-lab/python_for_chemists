#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
### ANCHOR_END: imports


### ANCHOR: pinv_function
def pinv(mat, rcond=1e-12):
    u, s, vh = np.linalg.svd(mat)
    nrows, ncols = mat.shape
    s_inv = np.zeros((ncols, nrows))
    for i in range(0, min(nrows, ncols)):
        if s[i] > rcond:
            s_inv[i, i] = 1 / s[i]
    return vh.T @ s_inv @ u.T
### ANCHOR_END: pinv_function

### ANCHOR: define_lineq
a_mat = np.array([
    [1, 1],
    [2, 1],
    [3, 1],
    [4, 1],
])
b_vec = np.array([4.0, 3.5, 4.5, 6.5])
### ANCHOR_END: define_lineq

### ANCHOR: solve_lineq
a_pinv = pinv(a_mat)
x0 = a_pinv @ b_vec
### ANCHOR_END: solve_lineq

### ANCHOR: verify_results
print(x0)

assert np.isclose(x0[0], 0.85)
assert np.isclose(x0[1], 2.50)
### ANCHOR_END: verify_results
