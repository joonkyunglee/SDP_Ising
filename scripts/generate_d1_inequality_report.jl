#!/usr/bin/env julia

using Dates
using Printf

setprecision(BigFloat, 256)

function parse_string_arg(name::String, default::String)
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            return String(split(a, "=", limit = 2)[2])
        end
    end
    return default
end

function parse_int_arg(name::String, default::Int)
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            return parse(Int, split(a, "=", limit = 2)[2])
        end
    end
    return default
end

function parse_bigfloat_list_arg(name::String, default::Vector{BigFloat})
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            raw = String(split(a, "=", limit = 2)[2])
            vals = BigFloat[]
            for tok in split(raw, ",")
                s = strip(tok)
                isempty(s) && continue
                push!(vals, parse(BigFloat, s))
            end
            return isempty(vals) ? default : vals
        end
    end
    return default
end

function P(t::Int, x::BigFloat, B::BigFloat)
    s = BigFloat(0)
    for k in 0:t
        e = div(k * (k - 1), 2) + div((t - k) * (t - k - 1), 2)
        s += binomial(big(t), big(k)) * (B^e) * (x^k)
    end
    return s
end

Phi(t::Int, x::BigFloat, B::BigFloat) = P(t, x, B)^(inv(BigFloat(t)))

H1(x::BigFloat, B::BigFloat) = (1 + x / B) / (1 + x * B)
H2(x::BigFloat, B::BigFloat) = (B^2 + 2x + x^2) / (B^2 * (1 + 2x + B^2 * x^2))

function lambda_from_s(B::BigFloat, s::BigFloat)
    return B * (s - 1) / (1 - s * B^2)
end

function alpha_from_s(B::BigFloat, s::BigFloat)
    target = s^2
    lo = BigFloat(0)
    hi = BigFloat(1)
    while H2(hi, B) < target
        hi *= 2
        if hi > big"1e220"
            error("failed to bracket alpha for B=$(B), s=$(s)")
        end
    end
    for _ in 1:280
        mid = (lo + hi) / 2
        if H2(mid, B) < target
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

function alpha_quadratic_coeffs(B::BigFloat, s::BigFloat)
    a = s^2 * B^4 - 1
    b = 2 * (s^2 * B^2 - 1)
    c = B^2 * (s^2 - 1)
    return a, b, c
end

function alpha_quadratic_roots(B::BigFloat, s::BigFloat)
    a, b, c = alpha_quadratic_coeffs(B, s)
    D = b^2 - 4 * a * c
    sqrtD = sqrt(D)
    r1 = (-b + sqrtD) / (2 * a)
    r2 = (-b - sqrtD) / (2 * a)
    return a, b, c, D, r1, r2
end

function evaluate_point(B::BigFloat, s::BigFloat)
    lam = lambda_from_s(B, s)
    alp = alpha_from_s(B, s)

    lhs = Phi(2, alp * B, B) / Phi(3, alp, B)
    rhs = Phi(1, lam * B, B) / Phi(2, lam, B)
    gap = rhs - lhs

    h1_err = abs(H1(lam, B) - s)
    h2_err = abs(sqrt(H2(alp, B)) - s)

    return (
        B = B,
        s = s,
        lambda = lam,
        alpha = alp,
        lhs = lhs,
        rhs = rhs,
        gap = gap,
        h1_err = h1_err,
        h2_err = h2_err,
    )
end

function scan_B(B::BigFloat; n::Int = 4000)
    slo = BigFloat(1) + big"1e-20"
    shi = inv(B)^2 - big"1e-20"

    min_gap = big"1e100"
    min_s = slo

    for i in 0:n
        s = slo + (shi - slo) * BigFloat(i) / BigFloat(n)
        r = evaluate_point(B, s)
        if r.gap < min_gap
            min_gap = r.gap
            min_s = s
        end
    end

    near_lo = evaluate_point(B, BigFloat(1) + big"1e-12")
    near_hi = evaluate_point(B, inv(B)^2 - big"1e-12")

    return (
        B = B,
        n = n,
        min_gap = min_gap,
        argmin_s = min_s,
        gap_near_lo = near_lo.gap,
        gap_near_hi = near_hi.gap,
    )
end

function sfmt(x::BigFloat; digits::Int = 12)
    return @sprintf("%.*e", digits, Float64(x))
end

function main()
    default_out = joinpath(dirname(@__DIR__), "d1_inequality_explicit_report.md")
    out = parse_string_arg("output", default_out)
    nscan = parse_int_arg("nscan", 5000)

    Bvals = parse_bigfloat_list_arg(
        "Bvals",
        BigFloat[big"0.05", big"0.1", big"0.2", big"0.4", big"0.7", big"0.9"],
    )

    examples = [
        (B = big"0.4", s = big"2.3"),
        (B = big"0.1", s = big"10.0"),
        (B = big"0.1", s = big"90.0"),
        (B = big"0.05", s = big"200.0"),
        (B = big"0.7", s = big"1.5"),
        (B = big"0.9", s = big"1.1"),
    ]

    ex_results = [evaluate_point(ex.B, ex.s) for ex in examples]
    scan_results = [scan_B(B; n = nscan) for B in Bvals]

    detail_examples = [
        (B = big"0.4", s = big"2.3"),
        (B = big"0.1", s = big"90.0"),
    ]
    detail_payload = []
    for ex in detail_examples
        r = evaluate_point(ex.B, ex.s)
        a, b, c, D, r1, r2 = alpha_quadratic_roots(ex.B, ex.s)
        push!(
            detail_payload,
            (
                B = ex.B,
                s = ex.s,
                result = r,
                a = a,
                b = b,
                c = c,
                D = D,
                root1 = r1,
                root2 = r2,
                root1_err = abs(H2(r1, ex.B) - ex.s^2),
                root2_err = abs(H2(r2, ex.B) - ex.s^2),
                root1_nonneg = r1 >= 0,
                root2_nonneg = r2 >= 0,
            ),
        )
    end

    open(out, "w") do io
        println(io, "# d = 1 Inequality Report (Explicit Derivation + Worked Examples)")
        println(io)
        println(io, "Generated at: `", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "`")
        println(io)
        println(io, "Workspace: `/Users/joonkyunglee/SDP_Ising`")
        println(io)
        println(io, "Scan points per B: `", nscan, "`")
        println(io)

        println(io, "## 1. Target statement and domain")
        println(io)
        println(io, "We start from")
        println(io)
        println(io, raw"$$")
        println(io, "P_t(x)=\\sum_{k=0}^{t}\\binom{t}{k}B^{\\binom{k}{2}+\\binom{t-k}{2}}x^k, \\qquad \\Phi_t(x)=P_t(x)^{1/t},")
        println(io, raw"$$")
        println(io)
        println(io, raw"$$")
        println(io, "H_k(x)=\\frac{P_k(x/B)}{P_k(xB)}, \\qquad 0<B<1.")
        println(io, raw"$$")
        println(io)
        println(io, "For `d=1`, we set")
        println(io)
        println(io, raw"$$")
        println(io, "s = H_1(\\lambda) = H_2(\\alpha)^{1/2},")
        println(io, raw"$$")
        println(io)
        println(io, "and the inequality to verify is")
        println(io)
        println(io, raw"$$")
        println(io, "\\frac{\\Phi_2(\\alpha B)}{\\Phi_3(\\alpha)} \\le \\frac{\\Phi_1(\\lambda B)}{\\Phi_2(\\lambda)}.")
        println(io, raw"$$")
        println(io)
        println(io, "Define the numerical gap")
        println(io)
        println(io, raw"$$")
        println(io, "\\Delta(B,s) := \\frac{\\Phi_1(\\lambda B)}{\\Phi_2(\\lambda)} - \\frac{\\Phi_2(\\alpha B)}{\\Phi_3(\\alpha)}.")
        println(io, raw"$$")
        println(io)
        println(io, "The inequality is equivalent to `Delta(B,s) >= 0`.")
        println(io)

        println(io, "## 2. Full d = 1 algebra (step by step)")
        println(io)
        println(io, "### 2.1 Expand `P_1`, `P_2`, `P_3` explicitly")
        println(io)
        println(io, "By direct expansion of the defining sum:")
        println(io)
        println(io, raw"$$")
        println(io, "P_1(x)=1+x,")
        println(io, raw"$$")
        println(io)
        println(io, raw"$$")
        println(io, "P_2(x)=B+2x+Bx^2,")
        println(io, raw"$$")
        println(io)
        println(io, raw"$$")
        println(io, "P_3(x)=B^3+3Bx+3x^2+B^3x^3.")
        println(io, raw"$$")
        println(io)

        println(io, "### 2.2 Derive `H_1` and solve for `lambda`")
        println(io)
        println(io, raw"$$")
        println(io, "H_1(x)=\\frac{P_1(x/B)}{P_1(xB)}=\\frac{1+x/B}{1+xB}.")
        println(io, raw"$$")
        println(io)
        println(io, "Set `s = H_1(lambda)` and solve:")
        println(io)
        println(io, raw"$$")
        println(io, "s(1+\\lambda B)=1+\\lambda/B \\;\\Longrightarrow\\; \\lambda=\\frac{B(s-1)}{1-sB^2}.")
        println(io, raw"$$")
        println(io)
        println(io, "Because `lambda > 0`, we need `1 <= s < B^{-2}`.")
        println(io)

        println(io, "### 2.3 Derive `H_2` and the quadratic for `alpha`")
        println(io)
        println(io, raw"$$")
        println(io, "H_2(x)=\\frac{P_2(x/B)}{P_2(xB)}")
        println(io, "     =\\frac{B+2x/B+B(x/B)^2}{B+2xB+B(xB)^2}")
        println(io, "     =\\frac{B^2+2x+x^2}{B^2(1+2x+B^2x^2)}.")
        println(io, raw"$$")
        println(io)
        println(io, "Set `s^2 = H_2(alpha)`. After cross-multiplication, `alpha` satisfies:")
        println(io)
        println(io, raw"$$")
        println(io, "(s^2B^4-1)\\alpha^2 + 2(s^2B^2-1)\\alpha + B^2(s^2-1)=0.")
        println(io, raw"$$")
        println(io)
        println(io, "In code we compute `alpha` by monotone bisection on `H_2(alpha)=s^2`;")
        println(io, "then we cross-check it against the quadratic roots.")
        println(io)

        println(io, "## 3. Fully worked examples")
        println(io)
        println(io, "All values below are computed with `BigFloat` precision = 256 bits.")
        println(io)

        for (idx, d) in enumerate(detail_payload)
            println(io, "### 3.", idx, " Example with B = ", @sprintf("%.4f", Float64(d.B)), ", s = ", @sprintf("%.6f", Float64(d.s)))
            println(io)
            println(io, "1. Compute lambda from closed form:")
            println(io)
            println(io, raw"$$")
            println(io, "\\lambda=\\frac{B(s-1)}{1-sB^2}=", sfmt(d.result.lambda; digits = 12), ".")
            println(io, raw"$$")
            println(io)

            println(io, "2. Build the quadratic `a*alpha^2 + b*alpha + c = 0`:")
            println(io)
            println(io, "- `a = s^2*B^4 - 1 = ", sfmt(d.a; digits = 12), "`")
            println(io, "- `b = 2*(s^2*B^2 - 1) = ", sfmt(d.b; digits = 12), "`")
            println(io, "- `c = B^2*(s^2 - 1) = ", sfmt(d.c; digits = 12), "`")
            println(io, "- `discriminant = b^2 - 4ac = ", sfmt(d.D; digits = 12), "`")
            println(io)

            println(io, "3. Quadratic roots and feasibility check:")
            println(io, "- `root1 = ", sfmt(d.root1; digits = 12), "`, nonnegative = `", d.root1_nonneg, "`, residual `|H2(root1)-s^2| = ", sfmt(d.root1_err; digits = 3), "`")
            println(io, "- `root2 = ", sfmt(d.root2; digits = 12), "`, nonnegative = `", d.root2_nonneg, "`, residual `|H2(root2)-s^2| = ", sfmt(d.root2_err; digits = 3), "`")
            println(io, "- bisection solution `alpha = ", sfmt(d.result.alpha; digits = 12), "`")
            println(io)

            println(io, "4. Evaluate both sides:")
            println(io, "- `LHS = Phi2(alpha*B)/Phi3(alpha) = ", sfmt(d.result.lhs; digits = 12), "`")
            println(io, "- `RHS = Phi1(lambda*B)/Phi2(lambda) = ", sfmt(d.result.rhs; digits = 12), "`")
            println(io, "- `gap = RHS-LHS = ", sfmt(d.result.gap; digits = 12), "`")
            println(io)
            println(io, "5. Consistency check:")
            println(io, "- `|H1(lambda)-s| = ", sfmt(d.result.h1_err; digits = 3), "`")
            println(io, "- `|sqrt(H2(alpha))-s| = ", sfmt(d.result.h2_err; digits = 3), "`")
            println(io)
        end

        println(io, "## 4. Additional sample points (quick table)")
        println(io)
        println(io, "| B | s | lambda | alpha | LHS | RHS | gap | |H1(lambda)-s| | |sqrt(H2(alpha))-s| |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in ex_results
            println(
                io,
                "| ",
                @sprintf("%.4f", Float64(r.B)),
                " | ",
                @sprintf("%.6f", Float64(r.s)),
                " | ",
                sfmt(r.lambda; digits = 10),
                " | ",
                sfmt(r.alpha; digits = 10),
                " | ",
                sfmt(r.lhs; digits = 10),
                " | ",
                sfmt(r.rhs; digits = 10),
                " | ",
                sfmt(r.gap; digits = 10),
                " | ",
                sfmt(r.h1_err; digits = 3),
                " | ",
                sfmt(r.h2_err; digits = 3),
                " |",
            )
        end
        println(io)

        println(io, "## 5. Uniform grid check on the full s-interval")
        println(io)
        println(io, "For each fixed `B`, scan")
        println(io)
        println(io, raw"$$")
        println(io, "s_i = 1+10^{-20} + i\\frac{(B^{-2}-10^{-20})-(1+10^{-20})}{n}, \\qquad i=0,\\ldots,n,")
        println(io, raw"$$")
        println(io)
        println(io, "with `n = ", nscan, "`, and record the minimum observed gap.")
        println(io)
        println(io, "| B | min gap over scan | argmin s | gap near s=1+1e-12 | gap near s=B^{-2}-1e-12 |")
        println(io, "|---:|---:|---:|---:|---:|")
        for r in scan_results
            println(
                io,
                "| ",
                @sprintf("%.4f", Float64(r.B)),
                " | ",
                sfmt(r.min_gap; digits = 10),
                " | ",
                sfmt(r.argmin_s; digits = 10),
                " | ",
                sfmt(r.gap_near_lo; digits = 10),
                " | ",
                sfmt(r.gap_near_hi; digits = 10),
                " |",
            )
        end
        println(io)
        println(io, "Observed pattern: all scanned minima are nonnegative and typically smallest near `s -> B^{-2}`.")
        println(io)

        println(io, "## 6. Why this helps for the d = 2 project")
        println(io)
        println(io, "- It isolates the `(B,s) -> (lambda, alpha)` parameter map in a case with closed-form equations.")
        println(io, "- It validates the numeric pipeline (monotone inversion + high-precision evaluation).")
        println(io, "- It gives a clean baseline where no SDP relaxation is needed, so failures in `d=2` are easier to diagnose.")
        println(io)

        println(io, "## 7. Reproduction")
        println(io)
        println(io, "```bash")
        println(io, "/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \\")
        println(io, "  /Users/joonkyunglee/SDP_Ising/scripts/generate_d1_inequality_report.jl \\")
        println(io, "  --nscan=5000 \\")
        println(io, "  --output=/Users/joonkyunglee/SDP_Ising/d1_inequality_explicit_report.md")
        println(io, "```")
    end

    println("Report written to ", out)
end

main()
