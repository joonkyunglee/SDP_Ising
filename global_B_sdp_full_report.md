# Global-$B$ SOS/SDP Report for $d=2$

Date: 2026-04-19  
Workspace: `/Users/joonkyunglee/SDP_Ising`  
Main script: `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_global_B.jl`

## 1) Goal

We want a **single SOS certificate valid for all $B \in (0,1)$** (not fixed $B$).

For $d=2$, define

$$
P_t(x)=\sum_{k=0}^{t}\binom{t}{k}B^{\binom{k}{2}+\binom{t-k}{2}}x^k,\qquad
\Phi_t(x)=P_t(x)^{1/t}.
$$

Using

$$
u=\frac{\lambda}{1+\lambda},\qquad
v=\frac{\alpha}{1+\alpha},\qquad
b=B,
$$

the target is to certify globally:

$$
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}\le
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}.
$$

Equivalent polynomial form:

$$
f(b,u,v)\ge 0 \quad \text{subject to} \quad g(b,u,v)=0,\;\; b,u,v\in(0,1).
$$

---

## 2) Polynomial model used in code

The script builds:

1. `fhat` = cleared polynomial equivalent of the inequality.
2. `gtilde` = cleared equality constraint (`gtilde = b^3 * ghat`) to remove $1/b$ denominators.
3. box multipliers:
   - $g_u = u(1-u)$,
   - $g_v = v(1-v)$,
   - $g_b = b(1-b)$.

Certificate searched:

$$
f_{\text{model}}-\gamma
=
\sigma_0
+\sigma_u\,g_u
+\sigma_v\,g_v
+\sigma_b\,g_b
+\tau\,g_{\text{model}},
$$

with $\sigma_\bullet$ SOS, $\tau$ polynomial.

`--no-bmult` disables $\sigma_b g_b$.

---

## 3) Efficiency techniques implemented

The global run is large, so the script uses multiple structural reductions.

### 3.1 Sparse support extraction (instead of dense monomial basis)

Functions:
- `support_tuples`
- `half_support`
- `shifted_half_support`
- `tau_support`

Key idea:
- build monomial exponent support sets in $(b,u,v)$,
- derive candidate Gram bases from half-support rules,
- avoid dense `monomials([b,u,v], 0:d)` construction.

Observed raw supports:
- $|E_f| = 2381$,
- $|E_g| = 456$.

### 3.2 Multiplier-specific bases

Separate sparse supports are built for:
- `S0` for $\sigma_0$,
- `Su` for $\sigma_u$,
- `Sv` for $\sigma_v$,
- `Sb` for $\sigma_b$,
- `St` for $\tau$.

Default level-0 sizes from run logs:
- $|S0|=524$,
- $|Su|=543$,
- $|Sv|=561$,
- $|Sb|=565$,
- $|St|=5306$ before truncation.

### 3.3 Tau-basis capping

Function:
- `cap_tau_support`

Option:
- `--max_tau=<N>`

This keeps only the first $N$ monomials (sorted by total degree then lexicographic style) to reduce linear variable count in $\tau$.

Tried values:
- full: `max_tau=0` (no cap, 5306),
- reduced: `max_tau=2500`,
- aggressive: `max_tau=1200`.

### 3.4 Controlled support expansion (optional)

Function:
- `expand_support`

Option:
- `--level=<k>`

Purpose:
- add neighboring monomials in a bounded way to increase flexibility without switching to dense bases.

Current runs used:
- `level=0`.

### 3.5 Coefficient scaling for conditioning

Option:
- default ON (disable via `--no-scale`)

Technique:
- $f_{\text{model}} = fhat / \max|c_f|$,
- $g_{\text{model}} = gtilde / \max|c_g|$.

Observed scales:
- `fscale = 1.048e+08`,
- `gscale = 5.520e+03`.

This improves numeric conditioning while preserving feasibility structure.

### 3.6 Solver controls

MOSEK settings from script:
- `MSK_DPAR_OPTIMIZER_MAX_TIME = --time`,
- `MSK_IPAR_NUM_THREADS = --threads` (if > 0),
- silent mode by default (`--verbose` to disable silence).

### 3.7 Diagnostic output

Script prints:
- support and basis sizes,
- termination status / primal status / dual status,
- `gamma_lower_bound_*` when optimal,
- candidate gamma for feasible-but-not-optimal stops,
- minimum eigenvalue of Gram matrices when optimal.

---

## 4) Runs performed and outcomes

All commands used absolute path to Julia:

`/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia`

| Run | Command highlights | Result |
|---|---|---|
| 1 | `--time=900 --threads=8` | `status = SLOW_PROGRESS`, no certificate |
| 2 | `--time=900 --threads=8 --max_tau=2500` | `status = SLOW_PROGRESS`, no certificate |
| 3 | `--time=900 --threads=8 --max_tau=2500 --no-bmult` | `status = SLOW_PROGRESS`, no certificate |
| 4 | `--time=600 --threads=8 --max_tau=1200` | `status = SLOW_PROGRESS`, no certificate |
| 5 | scaled model, `--time=900 --threads=8 --max_tau=2500` | `status = SLOW_PROGRESS`, `primal_status = FEASIBLE_POINT`, `dual_status = FEASIBLE_POINT`, no optimal certificate |
| 6 | scaled model, `--time=300 --threads=8 --max_tau=2500` | `status = TIME_LIMIT`, no certificate |

Summary:
- The global-$B$ model is implemented and runs successfully end-to-end.
- On current settings/hardware, MOSEK did not return `OPTIMAL` for a full global certificate.
- Best behavior observed: feasible iterate reached, but solver terminated before optimality.

---

## 5) How to run

### Default global run

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia /Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_global_B.jl --time=1800 --threads=8 --max_tau=2500
```

### Reduced model (for quicker diagnostics)

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia /Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_global_B.jl --time=900 --threads=8 --max_tau=1200 --no-bmult
```

### Verbose solver logs

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia /Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_global_B.jl --time=900 --threads=8 --max_tau=2500 --verbose
```

---

## 6) Practical interpretation

The current code is already engineered for sparsity and conditioning (support-driven SOS bases, tau capping, optional multiplier removal, coefficient normalization).  
The limiting factor is now solver convergence on a very large global PSD system, not model construction correctness.

---

## 7) Files produced in this effort

- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_global_B.jl`
- `/Users/joonkyunglee/SDP_Ising/global_B_sdp_full_report.md` (this file)
