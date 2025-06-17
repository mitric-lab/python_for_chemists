#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
### ANCHOR_END: imports

### ANCHOR: random_matrix
np.random.seed(42)
a_mat = np.random.rand(9, 6)
p = np.min(a_mat.shape)
### ANCHOR_END: random_matrix

### ANCHOR: double_evd
ata = np.dot(a_mat.T, a_mat)
aat = np.dot(a_mat, a_mat.T)

sigma_sq_v, v_mat = np.linalg.eigh(ata)
sigma_sq_u, u_mat = np.linalg.eigh(aat)

assert np.allclose(sigma_sq_v[-p:], sigma_sq_u[-p:])
s_devd = np.sqrt(sigma_sq_v[-p:][::-1])
u_devd = u_mat[:, -p:][:, ::-1]
v_devd = v_mat[:, -p:][:, ::-1]
### ANCHOR_END: double_evd

### ANCHOR: svd
u, s, vh = np.linalg.svd(a_mat, full_matrices=False)
print(u.shape)  # (9, 6)
print(s.shape)  # (6,)
print(vh.shape)  # (6, 6)
### ANCHOR_END: svd

### ANCHOR: comparison_decomposition
assert np.allclose(s_devd, s)
assert np.allclose(np.abs(u_devd.T @ u), np.eye(p))
assert np.allclose(np.abs(v_devd.T @ vh.T), np.eye(p))
### ANCHOR_END: comparison_decomposition

### ANCHOR: comparison_reconstruction
a_mat_recon_svd = u @ np.diag(s) @ vh
assert np.allclose(a_mat_recon_svd, a_mat)

a_mat_recon_devd = u_devd @ np.diag(s_devd) @ v_devd.T
print(np.allclose(a_mat_recon_devd, a_mat))  # usually false
### ANCHOR_END: comparison_reconstruction

