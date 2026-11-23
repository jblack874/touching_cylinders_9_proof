import json, sys
from flint import arb, arb_mat, ctx
ctx.prec = 600
sys.argv = sys.argv[:1]
exec(open('cert9.py').read().split('full0 = np.concatenate')[0])   # reuse F_and_J (generic over number type)
d = json.load(open('sol9.json')); keep = d['keep']
xt = [arb(s) for s in d['x']]
xt = [arb(x.mid()) for x in xt]                      # exact midpoints
rho = arb(2) ** -330
X = list(xt)
for c in keep: X[c] = xt[c] + arb(0, rho)            # box: midpoint +/- rho on the 52 unknowns
Fm, _ = F_and_J(xt, arb(0), keep)
_, JX = F_and_J(X, arb(0), keep)
_, Jm = F_and_J(xt, arb(0), keep)
m = len(keep)
Jmid = arb_mat([[arb(e.mid()) for e in row] for row in Jm])
C = Jmid.inv(); C = arb_mat([[arb(C[i, j].mid()) for j in range(m)] for i in range(m)])   # exact preconditioner
I = arb_mat([[arb(1) if i == j else arb(0) for j in range(m)] for i in range(m)])
dX = arb_mat([[arb(0, rho)] for _ in range(m)])
K = -(C * arb_mat([[f] for f in Fm])) + (I - C * arb_mat(JX)) * dX
worst = max(float((abs(K[i, 0]).upper())) / float(rho.upper()) for i in range(m))
inside = all(abs(K[i, 0]) < rho for i in range(m))
print('Krawczyk (arb, %d bits): box radius 2^-330, max |K_i|/rho = %.3e, strictly inside: %s' % (ctx.prec, worst, inside))
P, U, R = unpack(X, None)
pos = all(r > 0 for r in R)
nonpar = all(dot(cross(U[i], U[j]), cross(U[i], U[j])) > 0 for i in range(n) for j in range(i + 1, n))
print('all radii > 0:', pos, '| no two axes parallel:', nonpar, '| radii:', [r.str(8, radius=False) for r in R])
