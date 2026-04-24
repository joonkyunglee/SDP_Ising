# d = 1 Full-B-Range Management Report

Generated at: `2026-04-24 16:05:53`

Workspace: `/Users/joonkyunglee/SDP_Ising`

## 1. Why this report exists

The previous `d=1` note showed explicit formulas and sample `B` values,
but not a concrete way to manage the **entire range** `B in (0,1)`.

This report provides that management strategy and a dense numerical check.

## 2. Rectangle parameterization for the full domain

For `d=1`, define

$$
s = 1 + (B^{-2}-1)u, \qquad u\in(0,1).
$$

Then `(B,s)` over

$$
B\in(0,1), \quad s\in[1,B^{-2})
$$

is equivalent to the fixed rectangle

$$
(B,u)\in(0,1)\times(0,1).
$$

In this variable, `lambda` is explicit and cancellation-free:

$$
\lambda = \frac{u}{B(1-u)}.
$$

For `alpha`, use `s^2 = H_2(alpha)` which gives the quadratic

$$
(s^2B^4-1)\alpha^2 + 2(s^2B^2-1)\alpha + B^2(s^2-1)=0.
$$

The code uses a numerically stable quadratic-root formula and keeps the nonnegative root.

## 3. Full-range scan setup

- `B`-grid: log-clustered near `0` and `1`, linear in the middle.
- `u`-grid: same strategy (because gap is smallest near `u=0` and `u=1`).
- Precision: `BigFloat` (256-bit).

Parameters:
- `epsB = 1.000000e-04`, `nlogB = 90`, `nlinB = 241`
- `epsU = 1.000000e-06`, `nlogU = 120`, `nlinU = 321`
- `|Bgrid| = 421`, `|Ugrid| = 561`
- Total points = `236181`
- Scan time = `6.818 s`

## 4. Main result from the dense full-range scan

Global minimum observed gap:

- `min gap = 1.333062686603e-30`
- at `B = 9.999000000000e-01`, `u = 9.999990000000e-01`

Top worst-case `B` rows (smallest per-`B` minimum gaps):

| B | min_u(B) | min gap on Ugrid |
|---:|---:|---:|
| 9.9990000000e-01 | 9.9999900000e-01 | 1.333062686603e-30 |
| 9.9989192932e-01 | 9.9999900000e-01 | 1.682546817009e-30 |
| 9.9988320729e-01 | 9.9999900000e-01 | 2.123651146679e-30 |
| 9.9987378133e-01 | 9.9999900000e-01 | 2.680393672867e-30 |
| 9.9986359463e-01 | 9.9999900000e-01 | 3.383088281720e-30 |
| 9.9985258580e-01 | 9.9999900000e-01 | 4.269994969573e-30 |
| 9.9984068848e-01 | 9.9999900000e-01 | 5.389402556687e-30 |
| 9.9982783096e-01 | 9.9999900000e-01 | 6.802257223990e-30 |
| 9.9981393576e-01 | 9.9999900000e-01 | 8.585479894140e-30 |
| 9.9979891912e-01 | 9.9999900000e-01 | 1.083615294518e-29 |
| 9.9978269053e-01 | 9.9999900000e-01 | 1.367680402380e-29 |
| 9.9976515219e-01 | 9.9999900000e-01 | 1.726207438449e-29 |

## 5. Algebra-to-numeric consistency checks

Residuals for `s=H_1(lambda)` and `s=sqrt(H_2(alpha))` remain near machine-zero at 256-bit precision.

| B | u | gap | |H1(lambda)-s| | |sqrt(H2(alpha))-s| |
|---:|---:|---:|---:|---:|
| 1.000000e-01 | 2.000000e-01 | 4.5159801567e-02 | 2.764e-76 | 2.764e-76 |
| 1.000000e-01 | 8.000000e-01 | 3.7419397963e-04 | 1.105e-75 | 0.000e+00 |
| 4.000000e-01 | 5.000000e-01 | 5.4891742725e-03 | 0.000e+00 | 3.454e-77 |
| 8.000000e-01 | 2.000000e-01 | 8.9416662276e-05 | 0.000e+00 | 0.000e+00 |
| 8.000000e-01 | 8.000000e-01 | 3.5631292762e-05 | 0.000e+00 | 0.000e+00 |
| 1.000000e-04 | 5.000000e-01 | 2.1667867342e-04 | 0.000e+00 | 0.000e+00 |
| 9.999000e-01 | 5.000000e-01 | 2.0835937520e-14 | 0.000e+00 | 0.000e+00 |
| 5.000000e-01 | 1.000000e-06 | 6.3638548332e-18 | 1.727e-77 | 0.000e+00 |
| 5.000000e-01 | 9.999990e-01 | 4.9718387288e-20 | 0.000e+00 | 3.454e-77 |

## 6. Endpoint diagnostics (how we manage open boundaries)

For each listed `B`, define

$$
C_0(B) \approx \frac{\Delta(B,u_\varepsilon)}{u_\varepsilon^3}, \qquad
C_1(B) \approx \frac{\Delta(B,1-u_\varepsilon)}{u_\varepsilon^3}, \quad u_\varepsilon=10^{-6}.
$$

Both coefficients stay positive in the tested `B`-range:

| B | C0(B) ~ Delta/u^3 near u=0 | C1(B) ~ Delta/(1-u)^3 near u=1 |
|---:|---:|---:|
| 1.000000e-04 | 3.0253872064e+18 | 1.6666653667e-03 |
| 5.000000e-04 | 1.4862334916e+18 | 3.7267743723e-03 |
| 1.000000e-03 | 3.4259803224e+17 | 5.2704430027e-03 |
| 5.000000e-03 | 1.2558606262e+14 | 1.1784220320e-02 |
| 1.000000e-02 | 1.5895785054e+12 | 1.6661654668e-02 |
| 5.000000e-02 | 4.7256186195e+07 | 3.6988961436e-02 |
| 1.000000e-01 | 5.1115091800e+05 | 5.1139208398e-02 |
| 2.000000e-01 | 5.1513016055e+03 | 6.5944274444e-02 |
| 4.000000e-01 | 3.8131580341e+01 | 6.2476424225e-02 |
| 7.000000e-01 | 2.2460464520e-01 | 1.8497267397e-02 |
| 9.000000e-01 | 2.2674176150e-03 | 1.0845004589e-03 |
| 9.900000e-01 | 1.4020964047e-06 | 1.3068456574e-06 |

- `min C0(B) over table = 1.4020964047e-06`
- `min C1(B) over table = 1.3068456574e-06`

Additional boundary stress checks:
- Small-`B` rays (`B` from `1e-6` to `1e-2`, `u in {0.2,0.5,0.8,0.95}`):
  - `min gap = 2.0151582338e-08`
  - `min gap/sqrt(B) = 2.0145366501e-05`
- Near-`B=1` rays (`B=1-delta`, `delta` from `1e-6` to `1e-2`, same `u` probes):
  - `min gap = 1.4289559077e-22`
  - `min gap/delta^2 = 1.4289559077e-10`

## 7. Conclusion

This gives a practical full-range management recipe for `d=1`:
1. map `(B,s)` to `(B,u)` rectangle,
2. use explicit/stable formulas for `lambda` and `alpha`,
3. run boundary-aware dense scans and endpoint diagnostics.

All checks above found nonnegative gaps.

This is a strong numerical certificate, but still not a symbolic proof.

## 8. Reproduction

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \
  /Users/joonkyunglee/SDP_Ising/scripts/generate_d1_full_B_range_report.jl \
  --epsB=1e-4 --epsU=1e-6 \
  --nlogB=90 --nlinB=241 --nlogU=120 --nlinU=321 \
  --output=/Users/joonkyunglee/SDP_Ising/d1_full_B_range_management_report.md
```
