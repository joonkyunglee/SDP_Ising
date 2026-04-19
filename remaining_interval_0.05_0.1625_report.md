# Piecewise SOS Run Report: Remaining Interval `[0.05, 0.1625]`

Date: 2026-04-19  
Workspace: `/Users/joonkyunglee/SDP_Ising`  
Script: `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl`

## Command executed

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia /Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl --intervals=1 --time=900 --threads=8 --tau_schedule=1200,1800,2500 --bmin=0.05 --bmax=0.1625 --no-bmult
```

## Setup summary

- Interval: `[0.05, 0.1625]`
- Time limit per attempt: `900` seconds
- Threads: `8`
- Tau schedule: `1200, 1800, 2500`
- Global-b multiplier: disabled (`--no-bmult`)
- Scaling: enabled (default)
- Printed scales:
  - `fscale = 1.048e+08`
  - `gscale = 5.520e+03`

Shared support/basis sizes (for each attempt):

- `|Ef| = 2381`, `|Eg| = 456`
- `|S0| = 524`, `|SI| = 565`, `|SB| = 565`, `|Su| = 543`, `|Sv| = 561`
- `|Sτ|` depended on tau cap (`1200`, `1800`, `2500`)

## Results

### Attempt 1: `tau = 1200`

- status: `SLOW_PROGRESS`
- feasible gamma (original scale): `-0.000000000000e+00`
- certification status: **not certified** (no optimal termination)

### Attempt 2: `tau = 1800`

- status: `SLOW_PROGRESS`
- feasible gamma (original scale): `-0.000000000000e+00`
- certification status: **not certified** (no optimal termination)

### Attempt 3: `tau = 2500`

- status: `SLOW_PROGRESS`
- feasible gamma (original scale): `-0.000000000000e+00`
- certification status: **not certified** (no optimal termination)

## Final summary from run

- Interval `[0.050000, 0.162500]`
- Final reported status: `SLOW_PROGRESS`
- Certified: `NO`
- Overall run result: `PARTIAL/NO coverage`

## Interpretation

This remaining low-$B$ interval appears numerically borderline:

- solver repeatedly reaches feasible points with near-zero gamma,
- but does not reach `OPTIMAL` before progress stalls.

This indicates the obstruction is mainly solver convergence/conditioning, not obvious primal infeasibility.
