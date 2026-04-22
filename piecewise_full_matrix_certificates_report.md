# Piecewise Full Matrix Certificates (Successful Intervals)

Date: 2026-04-22  
Workspace: `/Users/joonkyunglee/SDP_Ising`  
Solver script: `/Users/joonkyunglee/SDP_Ising/scripts/run_d2_sos_piecewise_B.jl`

## What was added

The piecewise solver now supports full certificate export on successful intervals via:

- `--dump_full_cert`
- `--cert_out_dir=...`
- `--cert_tag=...`

For each certified interval, it writes:

- full Gram matrices (`Q0.tsv`, `QI.tsv`, `Qu.tsv`, `Qv.tsv`, optional `QB.tsv`),
- monomial bases (`basis_*.txt`),
- `tau_basis.txt` and `tau_coeffs.tsv`,
- `metadata.txt` (status, interval, scaling factors, residual, min eigenvalues).

## Certified intervals and artifact directories

| Interval | Status | Tau cap used | Artifact directory | Run log |
|---|---|---:|---|---|
| `[0.1625, 0.275]` | `OPTIMAL` | `1200` | `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_01625_0275_i01_L0p162500_U0p275000_tau1200` | `/Users/joonkyunglee/SDP_Ising/benchmarks/cert_full_01625_0275_20260422_221433.log` |
| `[0.275, 0.5]` | `OPTIMAL` | `1200` | `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_0275_05_i01_L0p275000_U0p500000_tau1200` | `/Users/joonkyunglee/SDP_Ising/benchmarks/cert_full_0275_05_20260422_221745.log` |
| `[0.5, 0.95]` | `OPTIMAL` | `1200` | `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_05_095_i01_L0p500000_U0p950000_tau1200` | `/Users/joonkyunglee/SDP_Ising/benchmarks/cert_full_05_095_20260422_222103.log` |

## Matrix sizes

### `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_01625_0275_i01_L0p162500_U0p275000_tau1200`

- `Q0.tsv`: `524x524`
- `QI.tsv`: `565x565`
- `Qu.tsv`: `543x543`
- `Qv.tsv`: `561x561`
- `tau_basis.txt`: `1200` lines
- `tau_coeffs.tsv`: `1201` lines (header + 1200 coefficients)

### `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_0275_05_i01_L0p275000_U0p500000_tau1200`

- `Q0.tsv`: `524x524`
- `QI.tsv`: `565x565`
- `Qu.tsv`: `543x543`
- `Qv.tsv`: `561x561`
- `tau_basis.txt`: `1200` lines
- `tau_coeffs.tsv`: `1201` lines

### `/Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_05_095_i01_L0p500000_U0p950000_tau1200`

- `Q0.tsv`: `524x524`
- `QI.tsv`: `565x565`
- `QB.tsv`: `565x565`
- `Qu.tsv`: `543x543`
- `Qv.tsv`: `561x561`
- `tau_basis.txt`: `1200` lines
- `tau_coeffs.tsv`: `1201` lines

## Metadata highlights

From each `metadata.txt`:

- all three have `status=OPTIMAL`
- all three have `gamma_original=-0.00000000000000000e+00`
- min eigenvalues are positive (numerically PSD)
- `residual_scaled_max_abs_coeff` is small (`~1e-8` scale)

## How to inspect quickly

Example for one interval:

```bash
cd /Users/joonkyunglee/SDP_Ising/certificates/piecewise_full/succ_01625_0275_i01_L0p162500_U0p275000_tau1200
head -n 30 metadata.txt
head -n 10 basis_Q0.txt
head -n 5 Q0.tsv
head -n 20 tau_coeffs.tsv
```

