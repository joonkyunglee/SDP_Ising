# Piecewise SOS Success Report

Date: 2026-04-19  
Workspace: `/Users/joonkyunglee/SDP_Ising`  
Script: `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl`

## Scope

This report records the **successful** interval certifications obtained during piecewise runs for the global-$B$ problem.

Each success below has:
- `status = OPTIMAL`,
- `gamma_original = -0.000000000000e+00`,
- nonnegative Gram-matrix eigenvalue checks (numerically PSD).

## Successful intervals

| Interval | Command highlights | Solver result | Gram matrix checks |
|---|---|---|---|
| `[0.500000, 0.950000]` | `--intervals=1 --bmin=0.5 --bmax=0.95 --time=600 --tau_schedule=1200` (global-b multiplier enabled) | `OPTIMAL`, certified | `min_eig(Q0)=1.267e-10`, `min_eig(QI)=1.527e-10`, `min_eig(QB)=3.162e-10`, `min_eig(Qu)=4.186e-10`, `min_eig(Qv)=2.320e-10` |
| `[0.275000, 0.500000]` | `--intervals=2 --bmin=0.05 --bmax=0.5 --time=300 --tau_schedule=1200 --no-bmult` | `OPTIMAL`, certified | `min_eig(Q0)=4.655e-11`, `min_eig(QI)=1.567e-10`, `min_eig(Qu)=1.142e-10`, `min_eig(Qv)=5.430e-11` |
| `[0.162500, 0.275000]` | `--intervals=2 --bmin=0.05 --bmax=0.275 --time=300 --tau_schedule=1200 --no-bmult` | `OPTIMAL`, certified | `min_eig(Q0)=7.326e-12`, `min_eig(QI)=2.312e-11`, `min_eig(Qu)=1.681e-11`, `min_eig(Qv)=9.616e-12` |

## Coverage achieved by certified piecewise intervals

Union of certified subintervals:

$$
[0.1625,\,0.95].
$$

So the currently certified piecewise region excludes only:

$$
[0.05,\,0.1625].
$$

## Notes

- The unresolved low-\(B\) segment `[0.05, 0.1625]` was run with `tau_schedule=1200,1800,2500` and reached feasible near-zero gamma repeatedly, but ended with `SLOW_PROGRESS` before `OPTIMAL`.
- Detailed run information for that unresolved segment is in:
  - `/Users/joonkyunglee/SDP_Ising/remaining_interval_0.05_0.1625_report.md`
