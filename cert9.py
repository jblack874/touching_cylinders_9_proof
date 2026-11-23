import numpy as np, json, sys, time
import flint
from flint import arb, arb_mat, ctx
from mpmath import mp, mpf
n = 9
g0 = np.load('gauge.npy'); P0 = g0[:27].reshape(9, 3); U0 = g0[27:54].reshape(9, 3); R0 = g0[54:]
# full parameter vector layout: for each cylinder i: p(3), u(3), r(1)
def full_index(i, k): return 7 * i + k          # k: 0-2 p, 3-5 u, 6 r
FIX0 = {full_index(0, k) for k in range(7)} | {full_index(1, 1)}   # cyl 0 fully fixed; p1_y = 0 (rotation about z)
def unpack(full, T):
    P = [[full[7*i+k] for k in range(3)] for i in range(n)]
    U = [[full[7*i+3+k] for k in range(3)] for i in range(n)]
    R = [full[7*i+6] for i in range(n)]
    return P, U, R
def sub(a, b): return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def F_and_J(full, zero, cols):
    """Residuals (52) and Jacobian w.r.t. full parameter vector (63 columns), restricted to cols."""
    P, U, R = unpack(full, None)
    rows, Jrows = [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = sub(P[j], P[i]); nn = cross(U[i], U[j]); g = dot(d, nn); s = R[i] + R[j]; N2 = dot(nn, nn)
            rows.append(g * g - s * s * N2)
            grad = {}
            for k in range(3):
                grad[7*j+k] = 2 * g * nn[k]; grad[7*i+k] = -2 * g * nn[k]
            gi = cross(U[j], d); gj = cross(d, U[i]); uu = dot(U[i], U[j]); ui2 = dot(U[i], U[i]); uj2 = dot(U[j], U[j])
            for k in range(3):
                grad[7*i+3+k] = 2 * g * gi[k] - s * s * (2 * uj2 * U[i][k] - 2 * uu * U[j][k])
                grad[7*j+3+k] = 2 * g * gj[k] - s * s * (2 * ui2 * U[j][k] - 2 * uu * U[i][k])
            grad[7*i+6] = -2 * s * N2; grad[7*j+6] = -2 * s * N2
            Jrows.append(grad)
    for i in range(1, n):
        rows.append(dot(U[i], U[i]) - 1); Jrows.append({7*i+3+k: 2 * U[i][k] for k in range(3)})
        rows.append(dot(P[i], U[i])); gr = {7*i+k: U[i][k] for k in range(3)}
        for k in range(3): gr[7*i+3+k] = P[i][k]
        Jrows.append(gr)
    J = [[Jr.get(c, zero) for c in cols] for Jr in Jrows]
    return rows, J
full0 = np.concatenate([np.concatenate([P0[i], U0[i], [R0[i]]]) for i in range(n)])
free = [c for c in range(63) if c not in FIX0]               # 55 candidates
_, Jf = F_and_J(list(full0), 0.0, free)
Jf = np.array(Jf, float)
import scipy.linalg as sl
q, r, piv = sl.qr(Jf, pivoting=True)
keep = sorted(free[p] for p in piv[:52]); fixextra = sorted(free[p] for p in piv[52:])
print('extra fixed params (3):', fixextra, 'cond of kept columns: %.3g' % np.linalg.cond(Jf[:, [free.index(c) for c in keep]]))
# dyadic rounding of every fixed coordinate (exact in arb and mpmath)
mp.prec = 400
def dy(v, b=40): return mpf(int(round(v * 2**b))) / 2**b
full = [mpf(0)] * 63
for c in range(63): full[c] = mpf(full0[c])
for c in fixextra: full[c] = dy(full0[c])
full[full_index(0,0)] = full[full_index(0,1)] = full[full_index(0,2)] = mpf(0)
full[full_index(0,3)] = full[full_index(0,4)] = mpf(0); full[full_index(0,5)] = mpf(1); full[full_index(0,6)] = mpf(1); full[full_index(1,1)] = mpf(0)
from mpmath import matrix, lu_solve, norm
for it in range(30):
    Fv, Jv = F_and_J(full, mpf(0), keep)
    nr = max(abs(x) for x in Fv)
    if it % 3 == 0 or nr < mpf(10)**-110: print('newton', it, 'max|F| = %s' % mp.nstr(nr, 5), flush=True)
    if nr < mpf(10)**-110: break
    dx = lu_solve(matrix(Jv), matrix([-x for x in Fv]))
    for a, c in enumerate(keep): full[c] += dx[a]
json.dump({'keep': keep, 'fixed_extra': fixextra, 'x': [mp.nstr(v, 120) for v in full]}, open('sol9.json', 'w'))
P, U, R = unpack([float(v) for v in full], None)
print('radii:', [round(r, 6) for r in R])
