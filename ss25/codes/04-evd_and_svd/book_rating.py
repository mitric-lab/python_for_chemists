#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
### ANCHOR_END: imports

### ANCHOR: svd
a_mat = np.array([
    [1, 1, 1, 0, 0],
    [3, 3, 3, 0, 0],
    [4, 4, 4, 0, 0],
    [5, 5, 5, 0, 0],
    [0, 2, 0, 4, 4],
    [0, 0, 0, 5, 5],
    [0, 1, 0, 2, 2],
])
u, s, vh = np.linalg.svd(a_mat, full_matrices=True)

print(s)
print(u)
print(vh)
### ANCHOR_END: svd

### ANCHOR: exact_recon
rank = 3
a_recon = u[:, :rank] @ np.diag(s[:rank]) @ vh[:rank]
assert np.allclose(a_recon, a_mat, atol=1e-12)
### ANCHOR_END: exact_recon

### ANCHOR: approx_recon
rank = 2
a_approx_recon = u[:, :rank] @ np.diag(s[:rank]) @ vh[:rank]

mae = np.mean(np.abs(a_mat - a_approx_recon))
rmse = np.linalg.norm(a_mat - a_approx_recon) / np.sqrt(a_mat.size)
maxae = np.max(np.abs(a_mat - a_approx_recon))
assert np.close(mae, 0.12496405)
assert np.close(rmse, 0.22744110)
assert np.close(maxae, 0.73442940)
### ANCHOR_END: exact_recon
