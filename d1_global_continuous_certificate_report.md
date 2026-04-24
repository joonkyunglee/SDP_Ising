# d = 1 Global Continuous Certificate (SOS)

Generated at: `2026-04-24 17:04:51`

Workspace: `/Users/joonkyunglee/SDP_Ising`

## 1. Certificate target

We certify the `d=1` inequality for all `B in (0,1)` (continuous range, not a sampled grid).

We use variables `(b,u,t)` with

$$
b=B, \qquad u\in(0,1), \qquad t=\alpha B^2(1-u), \qquad u=\frac{s-1}{B^{-2}-1}.
$$

Set `v=1-u`, `A=u+b^2v`, and define

$$
\mathrm{den}_R=b^2v^2+1-v^2, \quad
\mathrm{num}_L=b^2v^2+2tv+t^2, \quad
\mathrm{den}_L=b^6v^3+3tb^2v^2+3t^2v+t^3.
$$

$$
N(b,u,t)=\mathrm{den}_L^2 - \mathrm{num}_L^3\,\mathrm{den}_R^3.
$$

This is exactly the cross-multiplied form of `rhs^6-lhs^6`,
so proving `N >= 0` on the feasible set proves the original inequality.

The coupling equality is

$$
q_t(b,u,t):=(A^2-1)t^2 + 2v(A^2-b^2)t + v^2b^2(A^2-b^4)=0.
$$


## 2. SOS model

We solve

$$
N - \gamma = \sigma_0 + \sigma_b\,b(1-b) + \sigma_u\,u(1-u) + \sigma_t\,t(1-t) + \tau\,q_t,
$$

with `sigma_*` SOS and `tau` free polynomial.

We additionally use box localizers `b(1-b)`, `u(1-u)`, `t(1-t)`.

Why `t(1-t)` is valid: for `v=1-u`, `r=1-b^2`,
- `q_t(0)=v^2b^2(A^2-b^4) >= 0`,
- `q_t(1)=r v^2 h(r,v)` where
  `h=-2-2v-2r+4rv+rv^2+r^2-r^2v^2 <= -1` on `[0,1]^2`,
so `q_t(1) <= 0` and the nonnegative root lies in `[0,1]`.

Model settings:
- `level = 0`
- `max_tau_terms = 0`
- `time_limit = 600.0 s`
- `threads = 8`

Support/basis sizes:
- `|EN| = 281`, `|Eq| = 28`
- `|S0| = 119`, `|Sb| = 128`, `|Su| = 110`, `|St| = 112`, `|Sq| = 301`

## 3. Solver result

- `status = OPTIMAL`
- `primal_status = FEASIBLE_POINT`
- `dual_status = FEASIBLE_POINT`
- `gamma_lower_bound_scaled = -0.0000000000000000e+00`
- `gamma_lower_bound_original = -0.0000000000000000e+00`
- `certificate_result = NONNEGATIVE (within tolerance)`
- `min_eig(Q0) = 1.420e-10`
- `min_eig(Qb) = 5.151e-10`
- `min_eig(Qu) = 4.970e-10`
- `min_eig(Qt) = 3.681e-10`

## 4. Interpretation

This is a continuous semialgebraic certificate over the whole parameter region.
It is not a discrete sampling argument.

## 5. Reproduction

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \
  /Users/joonkyunglee/SDP_Ising/scripts/run_d1_global_continuous_sos.jl \
  --time=600 --threads=8 --level=0 --max_tau=0 \
  --report=/Users/joonkyunglee/SDP_Ising/d1_global_continuous_certificate_report.md
```
