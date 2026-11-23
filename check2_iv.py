"""Independent checker: mpmath.iv intervals, forward-mode automatic differentiation (no hand-derived gradients),
float64 preconditioner, written separately from cert9.py / kraw_arb.py. Reads only sol9.json."""
import json, numpy as np
from mpmath import iv, mp
iv.prec = 320; mp.prec = 320
d = json.load(open('sol9.json')); keep = d['keep']; xs = d['x']; m = len(keep); n = 9
idx = {c: a for a, c in enumerate(keep)}
class Dual:
    __slots__ = ('v', 'g')
    def __init__(s, v, g=None): s.v = v; s.g = g or {}
    def _c(o): return o if isinstance(o, Dual) else Dual(o)
    def __add__(s, o):
        o = Dual._c(o); g = dict(s.g)
        for k, w in o.g.items(): g[k] = g[k] + w if k in g else w
        return Dual(s.v + o.v, g)
    __radd__ = __add__
    def __neg__(s): return Dual(-s.v, {k: -w for k, w in s.g.items()})
    def __sub__(s, o): return s + (-Dual._c(o))
    def __rsub__(s, o): return Dual._c(o) - s
    def __mul__(s, o):
        o = Dual._c(o); g = {k: w * o.v for k, w in s.g.items()}
        for k, w in o.g.items(): g[k] = g[k] + w * s.v if k in g else w * s.v
        return Dual(s.v * o.v, g)
    __rmul__ = __mul__
def residuals(vals):
    V = [Dual(vals[c], {idx[c]: iv.mpf(1)}) if c in idx else Dual(vals[c]) for c in range(63)]
    p = lambda i: V[7*i:7*i+3]; u = lambda i: V[7*i+3:7*i+6]; r = lambda i: V[7*i+6]
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = u(i), u(j); w = [p(j)[k] - p(i)[k] for k in range(3)]
            # distance^2 between skew axes times |a x b|^2, via the 3x3 determinant det(w, a, b)
            det = w[0]*(a[1]*b[2]-a[2]*b[1]) - w[1]*(a[0]*b[2]-a[2]*b[0]) + w[2]*(a[0]*b[1]-a[1]*b[0])
            aa = a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; bb = b[0]*b[0]+b[1]*b[1]+b[2]*b[2]; ab = a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
            gram = aa*bb - ab*ab                                   # |a x b|^2 by Lagrange's identity
            s = r(i) + r(j)
            out.append(det*det - s*s*gram)
    for i in range(1, n):
        a = u(i); out.append(a[0]*a[0]+a[1]*a[1]+a[2]*a[2] - 1)
        q = p(i); out.append(q[0]*a[0]+q[1]*a[1]+q[2]*a[2])
    return out
mid = [iv.mpf(mp.mpf(x)) for x in xs]
mid = [iv.mpf(x.mid) for x in mid]
rho = iv.mpf(2) ** -150
box = list(mid)
for c in keep: box[c] = mid[c] + iv.mpf([-1, 1]) * rho
Fm = [e.v for e in residuals(mid)]
JX = [[e.g.get(a, iv.mpf(0)) for a in range(m)] for e in residuals(box)]
Jf = np.array([[float(e.g.get(a, iv.mpf(0)).mid) for a in range(m)] for e in residuals(mid)])
Cf = np.linalg.inv(Jf)
C = [[iv.mpf(float(Cf[i, j])) for j in range(m)] for i in range(m)]
K = []
for i in range(m):
    t = -sum((C[i][k] * Fm[k] for k in range(m)), iv.mpf(0))
    for j in range(m):
        cij = (iv.mpf(1) if i == j else iv.mpf(0)) - sum((C[i][k] * JX[k][j] for k in range(m)), iv.mpf(0))
        t += cij * (iv.mpf([-1, 1]) * rho)
    K.append(t)
ok = all(k.a > -rho and k.b < rho for k in K)
worst = max(max(abs(k.a), abs(k.b)) for k in K) / rho
print('checker 2 (mpmath.iv, %d bits, automatic differentiation): box radius 2^-150, max|K|/rho = %s, strictly inside: %s' % (iv.prec, mp.nstr(worst.b, 4), ok))
radii = [box[7*i+6] for i in range(n)]
print('radii lower bounds all > 0:', all(r.a > 0 for r in radii))
