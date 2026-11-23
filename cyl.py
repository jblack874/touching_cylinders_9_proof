import json, numpy as np
from mpmath import mp, mpf, matrix, lu_solve
raw = json.load(open('c9.txt'))
P = np.array([c[0] for c in raw], float); U = np.array([c[1] for c in raw], float); R = np.array([c[2] for c in raw], float)
U /= np.linalg.norm(U, axis=1)[:, None]
P = P - (P * U).sum(1)[:, None] * U          # standardize p.u = 0
n = len(R)
def dist(i, j):
    nn = np.cross(U[i], U[j]); return abs((P[j] - P[i]) @ nn) / np.linalg.norm(nn)
res = [dist(i, j) - (R[i] + R[j]) for i in range(n) for j in range(i + 1, n)]
print('n', n, 'radii', np.round(R, 4), 'max |residual|', max(map(abs, res)))
# gauge transform: p1 -> 0, u1 -> e3, rotate about z so standardized p2 has y = 0
def rot_to_z(u):
    z = np.array([0, 0, 1.]); v = np.cross(u, z); s = np.linalg.norm(v); c = u @ z
    if s < 1e-14: return np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)
P = P - P[0]; Rm = rot_to_z(U[0]); P = P @ Rm.T; U = U @ Rm.T
P = P - (P * U).sum(1)[:, None] * U
a = np.arctan2(P[1][1], P[1][0]); Rz = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])
P = P @ Rz.T; U = U @ Rz.T
s = R[0]; P /= s; R = R / s                       # scale r1 = 1
res = [dist(i, j) - (R[i] + R[j]) for i in range(n) for j in range(i + 1, n)]
print('after gauge: p1', np.round(P[0], 12), 'u1', np.round(U[0], 12), 'p2', np.round(P[1], 6), 'max res', max(map(abs, res)))
np.save('gauge.npy', np.concatenate([P.ravel(), U.ravel(), R]))
