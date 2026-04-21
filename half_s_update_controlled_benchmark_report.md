# Half-s Update and Controlled Benchmark Report

Date: 2026-04-21  
Workspace: `/Users/joonkyunglee/SDP_Ising`  
Main script: `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl`

## 1. Goal of this update

This report documents:

1. why the added `--half-s` terms were hurting runtime and certification,
2. what was changed in the SOS formulation,
3. controlled benchmark results on the interval `[0.1625, 0.275]`,
4. the practical conclusion for future runs.

## 2. Root-cause analysis

The observed slowdown was not a paradox. The previous `--half-s` formulation had two issues:

1. It introduced extra SOS cones (`\sigma_{uh}`, `\sigma_{vh}`) on top of existing cones, which enlarged the SDP significantly.
2. The old half multipliers used
   `u(0.5-u)` and `v(0.5-v)`, which are not nonnegative on all of `[0,1]`, so they were not globally valid localizing factors for the full-domain certificate.

So the mode increased model size and worsened conditioning, without being a clean symmetry-reduction reformulation.

## 3. Code updates applied

### 3.1 Certificate-safe center multipliers (experimental mode)

Inside `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl`, the half-mode multipliers were changed to globally nonnegative forms on `[0,1]`:

$$
 g_{uh}(u) = u(1-u)(2u-1)^2, \qquad g_{vh}(v) = v(1-v)(2v-1)^2.
$$

These appear in the script around:
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:256`

### 3.2 Dedicated supports for the new multipliers

Separate supports were added for the extra multipliers:

- `Suh_base = shifted_half_support(Ef, 2, [3, 4])`
- `Svh_base = shifted_half_support(Ef, 3, [3, 4])`

and threaded into `run_interval_certificate(...)`.

Relevant locations:
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:200`
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:390`
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:431`

### 3.3 Safe default behavior

To avoid accidental heavy runs, half-mode is now gated:

- `--half-s` alone: request noted, but **disabled by default**.
- `--half-s --half-s-experimental`: actually activates the extra center-multiplier cones.

Relevant locations:
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:369`
- `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl:408`

## 4. Controlled benchmark setup

Interval tested:

$$
B \in [0.1625, 0.275].
$$

Fixed command template (all modes):

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \
  /Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl \
  --intervals=1 --time=300 --threads=8 --tau_schedule=1200 \
  --bmin=0.1625 --bmax=0.275 --no-bmult
```

Modes benchmarked (3 repetitions each):

1. `baseline` (no `--half-s`)
2. `half_flag_only` (`--half-s` only; now still `use_half_s=false`)
3. `half_experimental` (`--half-s --half-s-experimental`; `use_half_s=true`)

Artifacts:
- `/Users/joonkyunglee/SDP_Ising/benchmarks/controlled_benchmark_01625_0275_20260421_210612.csv`
- `/Users/joonkyunglee/SDP_Ising/benchmarks/controlled_benchmark_01625_0275_20260421_210612.log`

## 5. Controlled benchmark results

### 5.1 Raw run table

| mode | rep | real (s) | status | certified | use_half_s |
|---|---:|---:|---|---|---|
| baseline | 1 | 174.25 | OPTIMAL | YES | false |
| baseline | 2 | 169.46 | OPTIMAL | YES | false |
| baseline | 3 | 179.69 | OPTIMAL | YES | false |
| half_flag_only | 1 | 196.68 | OPTIMAL | YES | false |
| half_flag_only | 2 | 189.07 | OPTIMAL | YES | false |
| half_flag_only | 3 | 224.54 | OPTIMAL | YES | false |
| half_experimental | 1 | 319.68 | SLOW_PROGRESS | NO | true |
| half_experimental | 2 | 325.31 | TIME_LIMIT | NO | true |
| half_experimental | 3 | 326.50 | TIME_LIMIT | NO | true |

### 5.2 Aggregates

- `baseline`: median `174.25s`, mean `174.47s`, certification `3/3`
- `half_flag_only`: median `196.68s`, mean `203.43s`, certification `3/3`
- `half_experimental`: median `325.31s`, mean `323.83s`, certification `0/3`

Relative to baseline median:

- `half_flag_only`: `+22.43s` (`1.129x`)
- `half_experimental`: `+151.06s` (`1.867x`)

## 6. Interpretation

The extra half-mode terms are redundant with existing localizers in a way that increases cone overlap and numerical burden.

In finite-time interior-point solves, that means:

1. larger linear algebra blocks,
2. worse conditioning in the KKT/Schur system,
3. slower convergence and more `TIME_LIMIT`/`SLOW_PROGRESS` outcomes.

So the observed degradation is consistent with the SDP geometry of the implemented formulation.

## 7. Practical recommendation

For production certification runs, keep half-mode off:

```bash
# recommended default path
... run_d2_sos_piecewise_B.jl ...
```

Use the heavy path only for research experiments:

```bash
# explicitly experimental
... run_d2_sos_piecewise_B.jl ... --half-s --half-s-experimental
```

At this stage, `--half-s` is documented as experimental and should not be expected to improve runtime or certification on this interval.
