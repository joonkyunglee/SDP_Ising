# d = 1 Inequality Report (Explicit Derivation + Worked Examples)

Generated at: `2026-04-24 14:55:28`

Workspace: `/Users/joonkyunglee/SDP_Ising`

Scan points per B: `5000`

## 1. Target statement and domain

We start from

$$
P_t(x)=\sum_{k=0}^{t}\binom{t}{k}B^{\binom{k}{2}+\binom{t-k}{2}}x^k, \qquad \Phi_t(x)=P_t(x)^{1/t},
$$

$$
H_k(x)=\frac{P_k(x/B)}{P_k(xB)}, \qquad 0<B<1.
$$

For `d=1`, we set

$$
s = H_1(\lambda) = H_2(\alpha)^{1/2},
$$

and the inequality to verify is

$$
\frac{\Phi_2(\alpha B)}{\Phi_3(\alpha)} \le \frac{\Phi_1(\lambda B)}{\Phi_2(\lambda)}.
$$

Define the numerical gap

$$
\Delta(B,s) := \frac{\Phi_1(\lambda B)}{\Phi_2(\lambda)} - \frac{\Phi_2(\alpha B)}{\Phi_3(\alpha)}.
$$

The inequality is equivalent to `Delta(B,s) >= 0`.

## 2. Full d = 1 algebra (step by step)

### 2.1 Expand `P_1`, `P_2`, `P_3` explicitly

By direct expansion of the defining sum:

$$
P_1(x)=1+x,
$$

$$
P_2(x)=B+2x+Bx^2,
$$

$$
P_3(x)=B^3+3Bx+3x^2+B^3x^3.
$$

### 2.2 Derive `H_1` and solve for `lambda`

$$
H_1(x)=\frac{P_1(x/B)}{P_1(xB)}=\frac{1+x/B}{1+xB}.
$$

Set `s = H_1(lambda)` and solve:

$$
s(1+\lambda B)=1+\lambda/B \;\Longrightarrow\; \lambda=\frac{B(s-1)}{1-sB^2}.
$$

Because `lambda > 0`, we need `1 <= s < B^{-2}`.

### 2.3 Derive `H_2` and the quadratic for `alpha`

$$
H_2(x)=\frac{P_2(x/B)}{P_2(xB)}
     =\frac{B+2x/B+B(x/B)^2}{B+2xB+B(xB)^2}
     =\frac{B^2+2x+x^2}{B^2(1+2x+B^2x^2)}.
$$

Set `s^2 = H_2(alpha)`. After cross-multiplication, `alpha` satisfies:

$$
(s^2B^4-1)\alpha^2 + 2(s^2B^2-1)\alpha + B^2(s^2-1)=0.
$$

In code we compute `alpha` by monotone bisection on `H_2(alpha)=s^2`;
then we cross-check it against the quadratic roots.

## 3. Fully worked examples

All values below are computed with `BigFloat` precision = 256 bits.

### 3.1 Example with B = 0.4000, s = 2.300000

1. Compute lambda from closed form:

$$
\lambda=\frac{B(s-1)}{1-sB^2}=8.227848101266e-01.
$$

2. Build the quadratic `a*alpha^2 + b*alpha + c = 0`:

- `a = s^2*B^4 - 1 = -8.645760000000e-01`
- `b = 2*(s^2*B^2 - 1) = -3.072000000000e-01`
- `c = B^2*(s^2 - 1) = 6.864000000000e-01`
- `discriminant = b^2 - 4ac = 2.468151705600e+00`

3. Quadratic roots and feasibility check:
- `root1 = -1.086217541155e+00`, nonnegative = `false`, residual `|H2(root1)-s^2| = 2.073e-76`
- `root2 = 7.308988647171e-01`, nonnegative = `true`, residual `|H2(root2)-s^2| = 0.000e+00`
- bisection solution `alpha = 7.308988647171e-01`

4. Evaluate both sides:
- `LHS = Phi2(alpha*B)/Phi3(alpha) = 8.617564656454e-01`
- `RHS = Phi1(lambda*B)/Phi2(lambda) = 8.732914754892e-01`
- `gap = RHS-LHS = 1.153500984377e-02`

5. Consistency check:
- `|H1(lambda)-s| = 3.454e-77`
- `|sqrt(H2(alpha))-s| = 3.454e-77`

### 3.2 Example with B = 0.1000, s = 90.000000

1. Compute lambda from closed form:

$$
\lambda=\frac{B(s-1)}{1-sB^2}=8.900000000000e+01.
$$

2. Build the quadratic `a*alpha^2 + b*alpha + c = 0`:

- `a = s^2*B^4 - 1 = -1.900000000000e-01`
- `b = 2*(s^2*B^2 - 1) = 1.600000000000e+02`
- `c = B^2*(s^2 - 1) = 8.099000000000e+01`
- `discriminant = b^2 - 4ac = 2.566155240000e+04`

3. Quadratic roots and feasibility check:
- `root1 = -5.058835971211e-01`, nonnegative = `false`, residual `|H2(root1)-s^2| = 2.200e-69`
- `root2 = 8.426111467550e+02`, nonnegative = `true`, residual `|H2(root2)-s^2| = 7.075e-74`
- bisection solution `alpha = 8.426111467550e+02`

4. Evaluate both sides:
- `LHS = Phi2(alpha*B)/Phi3(alpha) = 3.177874758423e-01`
- `RHS = Phi1(lambda*B)/Phi2(lambda) = 3.178370780184e-01`
- `gap = RHS-LHS = 4.960217608243e-05`

5. Consistency check:
- `|H1(lambda)-s| = 1.105e-75`
- `|sqrt(H2(alpha))-s| = 0.000e+00`

## 4. Additional sample points (quick table)

| B | s | lambda | alpha | LHS | RHS | gap | |H1(lambda)-s| | |sqrt(H2(alpha))-s| |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.4000 | 2.300000 | 8.2278481013e-01 | 7.3089886472e-01 | 8.6175646565e-01 | 8.7329147549e-01 | 1.1535009844e-02 | 3.454e-77 | 3.454e-77 |
| 0.1000 | 10.000000 | 1.0000000000e+00 | 1.0000000000e+00 | 6.4975717516e-01 | 7.4161984871e-01 | 9.1862673545e-02 | 1.382e-76 | 0.000e+00 |
| 0.1000 | 90.000000 | 8.9000000000e+01 | 8.4261114676e+02 | 3.1778747584e-01 | 3.1783707802e-01 | 4.9602176082e-05 | 1.105e-75 | 0.000e+00 |
| 0.0500 | 200.000000 | 1.9900000000e+01 | 2.6450407541e+02 | 2.5348382055e-01 | 2.5830680987e-01 | 4.8229893237e-03 | 0.000e+00 | 2.211e-75 |
| 0.7000 | 1.500000 | 1.3207547170e+00 | 1.3984662532e+00 | 9.0020683823e-01 | 9.0098765280e-01 | 7.8081457062e-04 | 0.000e+00 | 0.000e+00 |
| 0.9000 | 1.100000 | 8.2568807339e-01 | 8.1700327884e-01 | 9.7978923463e-01 | 9.7981239936e-01 | 2.3164723233e-05 | 1.727e-77 | 0.000e+00 |

## 5. Uniform grid check on the full s-interval

For each fixed `B`, scan

$$
s_i = 1+10^{-20} + i\frac{(B^{-2}-10^{-20})-(1+10^{-20})}{n}, \qquad i=0,\ldots,n,
$$

with `n = 5000`, and record the minimum observed gap.

| B | min gap over scan | argmin s | gap near s=1+1e-12 | gap near s=B^{-2}-1e-12 |
|---:|---:|---:|---:|---:|
| 0.0500 | 5.8230942368e-70 | 4.0000000000e+02 | 7.4535599250e-37 | 5.8230936914e-46 |
| 0.1000 | 5.2704627740e-68 | 1.0000000000e+02 | 5.2704627669e-37 | 5.2704627669e-44 |
| 0.2000 | 4.7702783521e-66 | 2.5000000000e+01 | 3.7267799625e-37 | 4.7702783520e-42 |
| 0.4000 | 4.3175630987e-64 | 6.2500000000e+00 | 2.6352313835e-37 | 4.3175630987e-40 |
| 0.7000 | 1.6405369244e-62 | 2.0408163265e+00 | 1.9920476822e-37 | 1.6405369244e-38 |
| 0.9000 | 8.4028200100e-62 | 1.2345679012e+00 | 1.7568209223e-37 | 8.4028200099e-38 |

Observed pattern: all scanned minima are nonnegative and typically smallest near `s -> B^{-2}`.

## 6. Why this helps for the d = 2 project

- It isolates the `(B,s) -> (lambda, alpha)` parameter map in a case with closed-form equations.
- It validates the numeric pipeline (monotone inversion + high-precision evaluation).
- It gives a clean baseline where no SDP relaxation is needed, so failures in `d=2` are easier to diagnose.

## 7. Reproduction

```bash
/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \
  /Users/joonkyunglee/SDP_Ising/scripts/generate_d1_inequality_report.jl \
  --nscan=5000 \
  --output=/Users/joonkyunglee/SDP_Ising/d1_inequality_explicit_report.md
```
