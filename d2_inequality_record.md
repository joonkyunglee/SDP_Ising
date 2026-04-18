# Record: \(d=2\) inequality check for \(\Phi_t\)

Date: 2026-04-18

## Setup

Define

\[
P_t(x)=\sum_{k=0}^t \binom{t}{k} B^{\binom{k}{2}+\binom{t-k}{2}}x^k,\quad
\Phi_t(x)=P_t(x)^{1/t},
\]
with \(0<B<1\), and
\[
H_k(x)=\frac{P_k(x/B)}{P_k(xB)}.
\]
For \(s>0\), \(x_d(s)\) is the unique positive solution of
\[
H_d(x_d(s))^{1/d}=s.
\]
For \(d=2\), set
\[
\lambda=x_2(s),\quad \alpha=x_3(s),
\]
so
\[
H_2(\lambda)^{1/2}=s,\quad H_3(\alpha)^{1/3}=s.
\]

Goal:
\[
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}\le
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}.
\]

## Explicit low-degree polynomials

\[
P_2(x)=B+2x+Bx^2,
\]
\[
P_3(x)=B^3+3Bx+3Bx^2+B^3x^3,
\]
\[
P_4(x)=B^6+4B^3x+6B^2x^2+4B^3x^3+B^6x^4.
\]

The matching-\(s\) constraints become
\[
\frac{B^2+2\lambda+\lambda^2}{B^2(1+2\lambda+B^2\lambda^2)}=s^2,
\]
\[
\frac{B^4+3B\alpha+3\alpha^2+B\alpha^3}
{B^4+3B^3\alpha+3B^4\alpha^2+B^7\alpha^3}=s^3.
\]

## Polynomial formulation for SOS/SDP

Define
\[
F(B,\lambda,\alpha):=
P_2(B\lambda)^6P_4(\alpha)^3-P_3(B\alpha)^4P_3(\lambda)^4.
\]
Since all quantities are positive, target inequality is equivalent to \(F\ge 0\).

Eliminate \(s\) via the same-\(s\) condition:
\[
G(B,\lambda,\alpha):=
P_2(\lambda/B)^3P_3(B\alpha)^2-P_2(B\lambda)^3P_3(\alpha/B)^2=0.
\]

So the problem is:
\[
F\ge 0
\quad\text{on}\quad
\{0<B<1,\ \lambda>0,\ \alpha>0,\ G=0\}.
\]

An SOS/SDP-friendly compactification:
\[
u=\frac{\lambda}{1+\lambda},\quad
v=\frac{\alpha}{1+\alpha},
\]
then search for
\[
\widehat F=
\sigma_0+\sigma_1\,u(1-u)+\sigma_2\,v(1-v)+\sigma_3\,B(1-B)+\tau\,\widehat G,
\]
with \(\sigma_i\) SOS.

## Numerical evidence (computed in workspace)

- Dense grid: \(21{,}600\) points over \(B\in(0.01,0.99)\), \(s\in(1,1/B^2)\).
- Random stress: \(12{,}000\) random points.
- Worst observed value of
  \[
  \frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}-
  \frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}
  \]
  was \(-7.738254481637341\times 10^{-14}\) (floating-point noise scale).

This is consistent with the inequality holding for \(d=2\).

## User-provided example check

For \(B=0.4,\ s=2.3\):

- \(\lambda=x_2(s)=0.7308988647170938\),
- \(\alpha=x_3(s)=0.6725534189593911\),
- \(Bs\,x_2(s)=0.6724269555397262\), so \(x_3(s)>Bs\,x_2(s)\) holds.

And

- \(\Phi_3(\alpha B)/\Phi_4(\alpha)=0.8562817533943582\),
- \(\Phi_2(\lambda B)/\Phi_3(\lambda)=0.8617564656454331\),

so
\[
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}
\le
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}
\]
is satisfied in this case.

## Notes

- This environment did not have SDP/SOS libraries preinstalled, so only algebraic reduction and numeric validation were performed here.
- Next step for a formal computational certificate: run the above SOS ansatz in a toolchain like `SumOfSquares.jl`, `SOSTOOLS`, or another SDP-capable framework.
