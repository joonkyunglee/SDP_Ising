# Understanding How `B` Is Parametrized in Our SDP (Beginner-Friendly)

Date: 2026-04-22  
Workspace: `/Users/joonkyunglee/SDP_Ising`

## 1. Why this note

You asked: "How is `B` parametrized in our SDP, and is there a better way?"

This file explains that in plain steps, assuming only basic familiarity with SDP.

---

## 2. Very short SDP primer

An SDP (semidefinite program) typically looks like:

$$
\text{find } y \text{ such that } F(y) \succeq 0,
$$

or optimization form:

$$
\max c^T y \quad \text{s.t. } F_0 + \sum_i y_i F_i \succeq 0.
$$

Here, "`\succeq 0`" means the matrix is positive semidefinite (PSD), i.e. all eigenvalues are nonnegative.

---

## 3. SOS as SDP (the key bridge)

For a polynomial `p(x)`, saying `p(x)` is sum-of-squares (SOS) means:

$$
p(x) = z(x)^T Q z(x), \quad Q \succeq 0,
$$

where `z(x)` is a chosen monomial basis (for example `[1, x, x^2]`).

Unknowns are entries of `Q`. Matching coefficients of both sides gives linear equations. PSD constraint on `Q` gives the SDP constraint.

So:
- polynomial nonnegativity claim -> SOS representation
- SOS representation -> SDP in Gram matrices `Q`

---

## 4. Constrained domain via localizing multipliers

Usually we need nonnegativity only on a region, not all space.

If region is:

$$
\Omega = \{x : g_1(x) \ge 0, g_2(x) \ge 0\},
$$

we use a Putinar-style certificate:

$$
f(x) - \gamma = \sigma_0(x) + \sigma_1(x) g_1(x) + \sigma_2(x) g_2(x),
$$

with each `sigma_i` SOS.

If this identity holds, then on `Omega` the RHS is nonnegative, so `f(x) >= gamma` on `Omega`.

---

## 5. What this project does for fixed `B`

In fixed-`B` script (`run_d2_sos_fixed_B.jl`), `B` is numeric. Unknowns are only in `(u,v)` variables.

The core identity is:

$$
\hat f(u,v) - \gamma
= \sigma_0 + \sigma_u\,u(1-u) + \sigma_v\,v(1-v) + \tau\,\hat g(u,v).
$$

Interpretation:
- `u(1-u)` and `v(1-v)` enforce domain `u,v in [0,1]`
- `tau * g` handles the equality constraint `g=0`

So this certifies one chosen `B` at a time.

---

## 6. What this project does for interval/global `B`

In piecewise script (`run_d2_sos_piecewise_B.jl`), `B` is promoted to a variable `b`.

For one interval `[L,U]`, we enforce `b in [L,U]` by localizer:

$$
g_I(b) = (b-L)(U-b) \ge 0.
$$

Then the certificate is:

$$
\hat f(b,u,v) - \gamma =
\sigma_0
+ \sigma_I\,g_I(b)
+ \sigma_u\,u(1-u)
+ \sigma_v\,v(1-v)
+ \tau\,\tilde g(b,u,v)
\;(+\;\sigma_B\,b(1-b)\text{ optionally}).
$$

So this is **not sampling `B` values**. It is one SDP proving the claim for all `b` in that whole interval.

---

## 7. Toy example 1: fixed parameter vs interval parameter

### Fixed parameter check (sampling style)

Suppose we want `p(x,t) >= 0` for all `x in [-1,1]` and all `t in [0,1]`.

Sampling style:
- pick `t = 0, 0.25, 0.5, 0.75, 1`
- prove each fixed-`t` problem

This only proves those points.

### Interval-parameter check (our piecewise style)

Treat `t` as variable and certify:

$$
p(x,t) = \sigma_0(x,t) + \sigma_x(x,t)(1-x^2) + \sigma_t(x,t)t(1-t).
$$

This certifies all `t in [0,1]` at once (subject to relaxation degree limits).

---

## 8. Toy example 2: why a localizer works

If `x in [0,1]`, then `x(1-x) >= 0`.

If we can write:

$$
f(x) = s_0(x) + s_1(x)\,x(1-x),
$$

with `s_0, s_1` SOS, then for `x in [0,1]`:
- `s_0(x) >= 0`
- `s_1(x) >= 0`
- `x(1-x) >= 0`

Hence `f(x) >= 0` on `[0,1]`.

Exactly the same logic is used for `b` with `(b-L)(U-b)`.

---

## 9. Why numerics get hard at small `B`

For your hard interval near `B=0.05`, common issues are:

1. Large coefficient scale spread in `fhat`, `gtilde` after clearing denominators.
2. Certificate optimum near boundary (`gamma` very close to `0`), which is numerically delicate.
3. Large Gram blocks (`500+` basis size) and large `tau` basis.
4. Interior-point stagnation (`SLOW_PROGRESS`) before reaching strict optimality.

This is why you can see feasible near-zero `gamma` but no formal `OPTIMAL` status.

---

## 10. Is there a better `B` parametrization?

Yes. The most practical improvements are:

### A) Affine remap of interval to `[-1,1]`

For each interval `[L,U]`, set

$$
b = c + h z, \quad c = (L+U)/2, \; h = (U-L)/2, \; z \in [-1,1].
$$

Then use localizer `1-z^2` instead of `(b-L)(U-b)`.

Why often better:
- normalized domain
- typically better-conditioned moments and Gram systems

### B) Two linear localizers instead of one quadratic

Instead of only `(b-L)(U-b)`, use

$$
\sigma_L (b-L) + \sigma_U (U-b)
$$

with SOS multipliers.

This can strengthen representation and improve stability in some instances.

### C) Orthogonal basis in the parameter variable

Use Chebyshev/Legendre basis in the transformed variable `z`, not raw monomials in `b`.

This usually reduces conditioning problems for high degree.

---

## 11. What we currently use (summary)

Current piecewise approach in this repo:

- `b` as polynomial variable (good)
- interval localizer `(b-L)(U-b)` (valid)
- optional `b(1-b)` multiplier (sometimes helps)
- scaling by max coefficient (`fscale`, `gscale`)
- adaptive splitting and feasibility mode added

So the framework is mathematically sound. Remaining difficulty is mostly numerical conditioning.

---

## 12. Practical recommendation for next step

If the target is hard low-`B` intervals, best next engineering move is:

1. implement affine remap `b = c + h z`,
2. keep adaptive splitting,
3. keep feasibility mode (`gamma = -eps`) for tight intervals.

That combination is usually stronger than just increasing time limits.

---

## 13. Glossary

- `SOS`: sum-of-squares polynomial.
- `Gram matrix`: PSD matrix `Q` in representation `z^T Q z`.
- `Localizer`: polynomial nonnegative on domain, multiplied by SOS.
- `tau multiplier`: unconstrained polynomial multiplier for equality constraint.
- `SLOW_PROGRESS`: solver termination due insufficient numerical progress (not proof of falsehood).

