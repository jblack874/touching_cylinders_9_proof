# Nine mutually touching cylinders exist

**A computer-assisted proof that nine infinite circular cylinders, with unequal radii, can touch pairwise in three-dimensional space.**

Littlewood asked how many infinite cylinders can touch each other pairwise. For equal radii the answer is seven (Bozóki–Lee–Rónyai 2015 constructed seven; Matić–Radoičić–Spomer 2026 proved no more). For arbitrary radii the same 2026 paper proves at most ten, and asks whether ten is possible. Nine had only been found numerically (Pikhitsa–Choi 2014; Georgiev, Gómez-Serrano, Tao and Wagner, the AlphaEvolve paper, 2025). This repository turns the numerical nine into a proof.

## What it shows

1. **Mathematical objective:** prove existence of 9 cylinders whose axes are pairwise skew lines at distance exactly r_i + r_j (so every pair touches and no interiors overlap).
2. **Baseline:** 9 numerically (no rigorous proof found); 7 rigorously (equal radii); upper bound 10 for arbitrary radii (arXiv:2610.26940, v2, 2 Nov 2026).
3. **Candidate result:** a configuration of 9 touching cylinders exists within 2^-330 of the point in `sol9.json`. Radii (cylinder 0 scaled to 1): 1, 0.7178, 1.8436, 0.8511, 0.6894, 2.1131, 6.4365, 2.5414, 6.1258.
4. **Verification:** two independent Krawczyk interval proofs (arb balls at 600 bits; mpmath intervals at 320 bits with automatic differentiation). Both pass. No Lean proof.
5. **Record eligibility:** an existence result, not a numeric record. Raises the rigorous lower bound for arbitrary radii from 7 to 9.
6. **Novelty:** new mathematics (apparently new within the searched scope: 2026 arXiv papers on touching cylinders, the Fejes Tóth–Kuperberg survey, Pikhitsa's papers, teorth/optimizationproblems and its open PRs, web search, 23 Nov 2026). The numerical configuration is AlphaEvolve's.
7. **Search economics:** about 4 minutes from the published float configuration to both proofs; the whole pipeline re-runs in under 2 seconds.
8. **Decision:** shipped. The natural next question is ten cylinders (45 equations, 43 unknowns: overdetermined, so this method does not apply).

## Numbers

| Quantity | Value |
| --- | --- |
| Unknowns after removing rigid motions, scale and the 2-parameter solution family | 52 |
| Equations (36 tangency, 8 unit-norm, 8 orthogonality) | 52 |
| Condition number of the Jacobian at the solution | 4.9 × 10³ |
| Residual of AlphaEvolve's float configuration | 4.4 × 10⁻¹³ |
| Newton residual after 4 steps | 10⁻¹¹⁹ |
| Krawczyk contraction, checker 1 (box radius 2⁻³³⁰) | 6.1 × 10⁻²⁰ |
| Krawczyk contraction, checker 2 (box radius 2⁻¹⁵⁰) | 2.9 × 10⁻¹³ |

## Verification rules

- Cylinder i: axis {p_i + t·u_i}, radius r_i > 0. Two cylinders touch when their axes are skew at distance r_i + r_j; then their interiors are disjoint and their boundaries meet.
- Tangency equation: det(p_j − p_i, u_i, u_j)² = (r_i + r_j)²·|u_i × u_j|². With r_i + r_j > 0 and |u_i × u_j| > 0 this is exactly distance = r_i + r_j.
- Gauge: p_0 = 0, u_0 = (0,0,1), r_0 = 1, the y-coordinate of p_1 = 0, and three more coordinates fixed at exact dyadic values chosen by column-pivoted QR.
- A passing Krawczyk test (image strictly inside the box) proves a unique exact solution in the box. Both checkers also confirm every radius is positive; checker 1 confirms no two axes are parallel.

## Honest limits

- Computer-assisted, not peer reviewed, not formalized in Lean (the Lean toolchain here cannot fetch Mathlib).
- The two checkers share the equation set and the gauge choice; they differ in library, precision, derivative method (hand-derived vs automatic) and preconditioner.
- Proves one configuration (a point of a 2-parameter family). Says nothing about ten cylinders.

## Run it

```sh
pip install python-flint mpmath numpy scipy
python3 cyl.py        # load AlphaEvolve's configuration, fix the gauge -> gauge.npy
python3 cert9.py      # square system, Newton to 1e-119 -> sol9.json
python3 kraw_arb.py   # checker 1
python3 check2_iv.py  # checker 2
```

## How it works

| Piece | Where |
| --- | --- |
| AlphaEvolve's numerical 9 cylinders | `c9.txt` |
| Gauge fixing (translate, rotate, scale) | `cyl.py` |
| Square system, QR gauge choice, Newton | `cert9.py` |
| Refined point (120 digits) and gauge | `sol9.json` |
| Checker 1 (arb, hand-derived Jacobian) | `kraw_arb.py` |
| Checker 2 (mpmath.iv, automatic differentiation) | `check2_iv.py` |
| Logs | `results/` |
| Computation manifest | `computation-manifest.json` |

## Method

Start from the float configuration (pairwise distance errors about 4 × 10⁻¹³). Translate, rotate and scale so cylinder 0 is the z-axis with radius 1, and rotate about z so cylinder 1's closest point has y = 0. The remaining 55 coordinates satisfy 52 equations; column-pivoted QR picks the 52 best-conditioned unknowns and the other three are fixed at 40-bit dyadic roundings. Newton's method in 400-bit arithmetic converges quadratically. Krawczyk's operator K = x̃ − C·F(x̃) + (I − C·J(X))(X − x̃) lands strictly inside X, which proves an exact solution exists in X.

## Stack

Python 3.13, python-flint 0.9 (arb), mpmath 1.4, NumPy, SciPy.

## Sources and credit

- Numerical configuration: B. Georgiev, J. Gómez-Serrano, T. Tao, A. Zs. Wagner, *Mathematical exploration and discovery at scale*, arXiv:2511.02864; data from [google-deepmind/alphaevolve_repository_of_problems](https://github.com/google-deepmind/alphaevolve_repository_of_problems) (`experiments/kissing_cylinders`).
- I. Matić, R. Radoičić, A. Spomer, *Mutually touching cylinders: Littlewood's seven, and at most ten in general*, arXiv:2610.26940.
- S. Bozóki, T.-L. Lee, L. Rónyai, *Seven mutually touching infinite cylinders*, Comput. Geom. 48 (2015).
- P. V. Pikhitsa, M. Choi, *Seven, eight, and nine mutually touching infinitely long straight round cylinders*, arXiv:1312.6207.

Prepared with AI assistance (Claude, Anthropic).
