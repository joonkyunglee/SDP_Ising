# Record: d = 2 inequality for Phi_t

Date: 2026-04-18

## 1) Problem statement

Given

$$
P_t(x)=\sum_{k=0}^t \binom{t}{k} B^{\binom{k}{2}+\binom{t-k}{2}}x^k,\qquad
\Phi_t(x)=P_t(x)^{1/t},
$$

with \(0 < B < 1\), define

$$
H_k(x)=\frac{P_k(x/B)}{P_k(xB)}.
$$

For \(s>0\), let \(x_d(s)\) be the unique positive solution of

$$
H_d(x_d(s))^{1/d}=s.
$$

For \(d=2\), set

$$
\lambda=x_2(s),\qquad \alpha=x_3(s),
$$

so

$$
H_2(\lambda)^{1/2}=s,\qquad H_3(\alpha)^{1/3}=s.
$$

Target inequality:

$$
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}
\le
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}.
$$

## 2) Explicit polynomials (d = 2 case)

$$
P_2(x)=B+2x+Bx^2,
$$

$$
P_3(x)=B^3+3Bx+3Bx^2+B^3x^3,
$$

$$
P_4(x)=B^6+4B^3x+6B^2x^2+4B^3x^3+B^6x^4.
$$

The equal-\(s\) constraints become:

$$
\frac{B^2+2\lambda+\lambda^2}{B^2(1+2\lambda+B^2\lambda^2)}=s^2,
$$

$$
\frac{B^4+3B\alpha+3\alpha^2+B\alpha^3}
{B^4+3B^3\alpha+3B^4\alpha^2+B^7\alpha^3}=s^3.
$$

## 3) SOS/SDP-ready polynomial form

Define

$$
F(B,\lambda,\alpha)
=P_2(B\lambda)^6P_4(\alpha)^3
-P_3(B\alpha)^4P_3(\lambda)^4.
$$

Since all terms are positive, the target inequality is equivalent to \(F\ge 0\).

Eliminate \(s\) using the equal-\(s\) condition:

$$
G(B,\lambda,\alpha)
=P_2(\lambda/B)^3P_3(B\alpha)^2
-P_2(B\lambda)^3P_3(\alpha/B)^2
=0.
$$

So we need

$$
F\ge 0
\quad\text{on}\quad
\{\,0<B<1,\ \lambda>0,\ \alpha>0,\ G=0\,\}.
$$

### Compactification (for numerical SOS)

Use

$$
u=\frac{\lambda}{1+\lambda},\qquad
v=\frac{\alpha}{1+\alpha},
$$

and seek a certificate of form

$$
\widehat F
=\sigma_0
+\sigma_1\,u(1-u)
+\sigma_2\,v(1-v)
+\sigma_3\,B(1-B)
+\tau\,\widehat G,
$$

where each \(\sigma_i\) is SOS.

## 4) Numerical evidence (workspace computation)

- Dense grid: 21,600 points over \(B\in(0.01,0.99)\), \(s\in(1,1/B^2)\).
- Random stress: 12,000 random points.
- Worst observed value of

$$
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}
-\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}
$$

was

$$
-7.738254481637341\times 10^{-14},
$$

which is at floating-point noise scale.

This is consistent with the inequality holding for \(d=2\).

## 5) Check for user example

For \(B=0.4,\ s=2.3\):

- \(\lambda=x_2(s)=0.7308988647170938\),
- \(\alpha=x_3(s)=0.6725534189593911\),
- \(Bs\,x_2(s)=0.6724269555397262\), so \(x_3(s)>Bs\,x_2(s)\).

And

$$
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}=0.8562817533943582,
$$

$$
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}=0.8617564656454331.
$$

Hence

$$
\frac{\Phi_3(\alpha B)}{\Phi_4(\alpha)}
\le
\frac{\Phi_2(\lambda B)}{\Phi_3(\lambda)}
$$

holds in this example.

## 6) Notes

- No SDP/SOS package was available in this environment, so this record includes algebraic reduction + numerical validation.
- For a formal computational certificate, run the SOS formulation in `SumOfSquares.jl`, `SOSTOOLS`, or another SDP-capable stack.
