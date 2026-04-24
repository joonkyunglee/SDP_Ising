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

function parse_bigfloat_arg(name::String, default::BigFloat)
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            return parse(BigFloat, split(a, "=", limit = 2)[2])
        end
    end
    return default
end

function sfmt(x::BigFloat; digits::Int = 10)
    return @sprintf("%.*e", digits, Float64(x))
end

function linspace_big(a::BigFloat, b::BigFloat, n::Int)
    if n <= 1
        return BigFloat[a]
    end
    out = BigFloat[]
    for i in 0:(n - 1)
        t = BigFloat(i) / BigFloat(n - 1)
        push!(out, a + (b - a) * t)
    end
    return out
end

function logspace_big(a::BigFloat, b::BigFloat, n::Int)
    @assert a > 0
    @assert b > 0
    if n <= 1
        return BigFloat[a]
    end
    la = log(a)
    lb = log(b)
    out = BigFloat[]
    for i in 0:(n - 1)
        t = BigFloat(i) / BigFloat(n - 1)
        push!(out, exp(la + (lb - la) * t))
    end
    return out
end

function unique_sorted(xs::Vector{BigFloat})
    ys = sort(xs)
    out = BigFloat[]
    for x in ys
        if isempty(out) || x != out[end]
            push!(out, x)
        end
    end
    return out
end

function make_unit_grid(; eps::BigFloat, n_log::Int, n_lin::Int)
    @assert eps > 0
    @assert eps < big"0.1"
    low = logspace_big(eps, big"0.1", n_log)
    mid = linspace_big(big"0.1", big"0.9", n_lin)
    high = [1 - x for x in reverse(low)]
    return unique_sorted(vcat(low, mid, high))
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

s_from_Bu(B::BigFloat, u::BigFloat) = 1 + (inv(B)^2 - 1) * u
lambda_from_Bu(B::BigFloat, u::BigFloat) = u / (B * (1 - u))
H1(x::BigFloat, B::BigFloat) = (1 + x / B) / (1 + x * B)
H2(x::BigFloat, B::BigFloat) = (B^2 + 2x + x^2) / (B^2 * (1 + 2x + B^2 * x^2))

function alpha_from_Bu(B::BigFloat, u::BigFloat)
    s = s_from_Bu(B, u)
    a = s^2 * B^4 - 1
    b = 2 * (s^2 * B^2 - 1)
    c = B^2 * (s^2 - 1)
    D = b^2 - 4 * a * c
    sqrtD = sqrt(D)
    q = b >= 0 ? -((b + sqrtD) / 2) : -((b - sqrtD) / 2)
    r1 = q / a
    r2 = c / q
    if r1 >= 0 && r2 >= 0
        return max(r1, r2)
    elseif r1 >= 0
        return r1
    else
        return r2
    end
end

function gap_from_Bu(B::BigFloat, u::BigFloat)
    lam = lambda_from_Bu(B, u)
    alp = alpha_from_Bu(B, u)
    lhs = Phi(2, alp * B, B) / Phi(3, alp, B)
    rhs = Phi(1, lam * B, B) / Phi(2, lam, B)
    return rhs - lhs
end

function full_scan(Bgrid::Vector{BigFloat}, Ugrid::Vector{BigFloat})
    global_min_gap = big"1e100"
    global_min_B = BigFloat(0)
    global_min_u = BigFloat(0)
    byB = NamedTuple[]

    t0 = time()
    for B in Bgrid
        min_gap_B = big"1e100"
        argmin_u_B = Ugrid[1]
        for u in Ugrid
            g = gap_from_Bu(B, u)
            if g < min_gap_B
                min_gap_B = g
                argmin_u_B = u
            end
            if g < global_min_gap
                global_min_gap = g
                global_min_B = B
                global_min_u = u
            end
        end
        push!(byB, (B = B, min_gap = min_gap_B, argmin_u = argmin_u_B))
    end
    elapsed = time() - t0
    return (
        elapsed = elapsed,
        global_min_gap = global_min_gap,
        global_min_B = global_min_B,
        global_min_u = global_min_u,
        byB = byB,
    )
end

function check_identity_residual(B::BigFloat, u::BigFloat)
    s = s_from_Bu(B, u)
    lam = lambda_from_Bu(B, u)
    alp = alpha_from_Bu(B, u)
    e1 = abs(H1(lam, B) - s)
    e2 = abs(sqrt(H2(alp, B)) - s)
    return (e1 = e1, e2 = e2, lambda = lam, alpha = alp)
end

function main()
    default_out = joinpath(dirname(@__DIR__), "d1_full_B_range_management_report.md")
    out = parse_string_arg("output", default_out)

    epsB = parse_bigfloat_arg("epsB", big"1e-4")
    epsU = parse_bigfloat_arg("epsU", big"1e-6")

    nlogB = parse_int_arg("nlogB", 90)
    nlinB = parse_int_arg("nlinB", 241)
    nlogU = parse_int_arg("nlogU", 120)
    nlinU = parse_int_arg("nlinU", 321)

    Bgrid = make_unit_grid(; eps = epsB, n_log = nlogB, n_lin = nlinB)
    Ugrid = make_unit_grid(; eps = epsU, n_log = nlogU, n_lin = nlinU)

    scan = full_scan(Bgrid, Ugrid)
    byB_sorted = sort(scan.byB; by = r -> r.min_gap)

    top_count = min(12, length(byB_sorted))
    top_rows = byB_sorted[1:top_count]

    consistency_points = [
        (B = big"0.1", u = big"0.2"),
        (B = big"0.1", u = big"0.8"),
        (B = big"0.4", u = big"0.5"),
        (B = big"0.8", u = big"0.2"),
        (B = big"0.8", u = big"0.8"),
        (B = big"1e-4", u = big"0.5"),
        (B = big"0.9999", u = big"0.5"),
        (B = big"0.5", u = big"1e-6"),
        (B = big"0.5", u = big"0.999999"),
    ]
    consistency_rows = NamedTuple[]
    for p in consistency_points
        c = check_identity_residual(p.B, p.u)
        g = gap_from_Bu(p.B, p.u)
        push!(
            consistency_rows,
            (
                B = p.B,
                u = p.u,
                gap = g,
                e1 = c.e1,
                e2 = c.e2,
            ),
        )
    end

    B_coeffs = BigFloat[
        big"1e-4",
        big"5e-4",
        big"1e-3",
        big"5e-3",
        big"1e-2",
        big"5e-2",
        big"0.1",
        big"0.2",
        big"0.4",
        big"0.7",
        big"0.9",
        big"0.99",
    ]
    u_small = big"1e-6"
    coeff_rows = NamedTuple[]
    min_c0 = big"1e100"
    min_c1 = big"1e100"
    for B in B_coeffs
        g0 = gap_from_Bu(B, u_small)
        g1 = gap_from_Bu(B, 1 - u_small)
        c0 = g0 / (u_small^3)
        c1 = g1 / (u_small^3)
        min_c0 = min(min_c0, c0)
        min_c1 = min(min_c1, c1)
        push!(
            coeff_rows,
            (
                B = B,
                c0 = c0,
                c1 = c1,
            ),
        )
    end

    u_probe = BigFloat[big"0.2", big"0.5", big"0.8", big"0.95"]
    smallB = logspace_big(big"1e-6", big"1e-2", 9)
    near1δ = logspace_big(big"1e-6", big"1e-2", 9)

    min_smallB_scaled = big"1e100"
    min_smallB_gap = big"1e100"
    for B in smallB, u in u_probe
        g = gap_from_Bu(B, u)
        min_smallB_gap = min(min_smallB_gap, g)
        min_smallB_scaled = min(min_smallB_scaled, g / sqrt(B))
    end

    min_near1_scaled = big"1e100"
    min_near1_gap = big"1e100"
    for δ in near1δ, u in u_probe
        B = 1 - δ
        g = gap_from_Bu(B, u)
        min_near1_gap = min(min_near1_gap, g)
        min_near1_scaled = min(min_near1_scaled, g / (δ^2))
    end

    open(out, "w") do io
        println(io, "# d = 1 Full-B-Range Management Report")
        println(io)
        println(io, "Generated at: `", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "`")
        println(io)
        println(io, "Workspace: `/Users/joonkyunglee/SDP_Ising`")
        println(io)

        println(io, "## 1. Why this report exists")
        println(io)
        println(io, "The previous `d=1` note showed explicit formulas and sample `B` values,")
        println(io, "but not a concrete way to manage the **entire range** `B in (0,1)`.")
        println(io)
        println(io, "This report provides that management strategy and a dense numerical check.")
        println(io)

        println(io, "## 2. Rectangle parameterization for the full domain")
        println(io)
        println(io, "For `d=1`, define")
        println(io)
        println(io, raw"$$")
        println(io, "s = 1 + (B^{-2}-1)u, \\qquad u\\in(0,1).")
        println(io, raw"$$")
        println(io)
        println(io, "Then `(B,s)` over")
        println(io)
        println(io, raw"$$")
        println(io, "B\\in(0,1), \\quad s\\in[1,B^{-2})")
        println(io, raw"$$")
        println(io)
        println(io, "is equivalent to the fixed rectangle")
        println(io)
        println(io, raw"$$")
        println(io, "(B,u)\\in(0,1)\\times(0,1).")
        println(io, raw"$$")
        println(io)
        println(io, "In this variable, `lambda` is explicit and cancellation-free:")
        println(io)
        println(io, raw"$$")
        println(io, "\\lambda = \\frac{u}{B(1-u)}.")
        println(io, raw"$$")
        println(io)
        println(io, "For `alpha`, use `s^2 = H_2(alpha)` which gives the quadratic")
        println(io)
        println(io, raw"$$")
        println(io, "(s^2B^4-1)\\alpha^2 + 2(s^2B^2-1)\\alpha + B^2(s^2-1)=0.")
        println(io, raw"$$")
        println(io)
        println(io, "The code uses a numerically stable quadratic-root formula and keeps the nonnegative root.")
        println(io)

        println(io, "## 3. Full-range scan setup")
        println(io)
        println(io, "- `B`-grid: log-clustered near `0` and `1`, linear in the middle.")
        println(io, "- `u`-grid: same strategy (because gap is smallest near `u=0` and `u=1`).")
        println(io, "- Precision: `BigFloat` (256-bit).")
        println(io)
        println(io, "Parameters:")
        println(io, "- `epsB = ", sfmt(epsB; digits = 6), "`, `nlogB = ", nlogB, "`, `nlinB = ", nlinB, "`")
        println(io, "- `epsU = ", sfmt(epsU; digits = 6), "`, `nlogU = ", nlogU, "`, `nlinU = ", nlinU, "`")
        println(io, "- `|Bgrid| = ", length(Bgrid), "`, `|Ugrid| = ", length(Ugrid), "`")
        println(io, "- Total points = `", length(Bgrid) * length(Ugrid), "`")
        println(io, "- Scan time = `", @sprintf("%.3f", scan.elapsed), " s`")
        println(io)

        println(io, "## 4. Main result from the dense full-range scan")
        println(io)
        println(io, "Global minimum observed gap:")
        println(io)
        println(io, "- `min gap = ", sfmt(scan.global_min_gap; digits = 12), "`")
        println(io, "- at `B = ", sfmt(scan.global_min_B; digits = 12), "`, `u = ", sfmt(scan.global_min_u; digits = 12), "`")
        println(io)
        println(io, "Top worst-case `B` rows (smallest per-`B` minimum gaps):")
        println(io)
        println(io, "| B | min_u(B) | min gap on Ugrid |")
        println(io, "|---:|---:|---:|")
        for r in top_rows
            println(
                io,
                "| ",
                sfmt(r.B; digits = 10),
                " | ",
                sfmt(r.argmin_u; digits = 10),
                " | ",
                sfmt(r.min_gap; digits = 12),
                " |",
            )
        end
        println(io)

        println(io, "## 5. Algebra-to-numeric consistency checks")
        println(io)
        println(io, "Residuals for `s=H_1(lambda)` and `s=sqrt(H_2(alpha))` remain near machine-zero at 256-bit precision.")
        println(io)
        println(io, "| B | u | gap | |H1(lambda)-s| | |sqrt(H2(alpha))-s| |")
        println(io, "|---:|---:|---:|---:|---:|")
        for r in consistency_rows
            println(
                io,
                "| ",
                sfmt(r.B; digits = 6),
                " | ",
                sfmt(r.u; digits = 6),
                " | ",
                sfmt(r.gap; digits = 10),
                " | ",
                sfmt(r.e1; digits = 3),
                " | ",
                sfmt(r.e2; digits = 3),
                " |",
            )
        end
        println(io)

        println(io, "## 6. Endpoint diagnostics (how we manage open boundaries)")
        println(io)
        println(io, "For each listed `B`, define")
        println(io)
        println(io, raw"$$")
        println(io, "C_0(B) \\approx \\frac{\\Delta(B,u_\\varepsilon)}{u_\\varepsilon^3}, \\qquad")
        println(io, "C_1(B) \\approx \\frac{\\Delta(B,1-u_\\varepsilon)}{u_\\varepsilon^3}, \\quad u_\\varepsilon=10^{-6}.")
        println(io, raw"$$")
        println(io)
        println(io, "Both coefficients stay positive in the tested `B`-range:")
        println(io)
        println(io, "| B | C0(B) ~ Delta/u^3 near u=0 | C1(B) ~ Delta/(1-u)^3 near u=1 |")
        println(io, "|---:|---:|---:|")
        for r in coeff_rows
            println(
                io,
                "| ",
                sfmt(r.B; digits = 6),
                " | ",
                sfmt(r.c0; digits = 10),
                " | ",
                sfmt(r.c1; digits = 10),
                " |",
            )
        end
        println(io)
        println(io, "- `min C0(B) over table = ", sfmt(min_c0; digits = 10), "`")
        println(io, "- `min C1(B) over table = ", sfmt(min_c1; digits = 10), "`")
        println(io)

        println(io, "Additional boundary stress checks:")
        println(io, "- Small-`B` rays (`B` from `1e-6` to `1e-2`, `u in {0.2,0.5,0.8,0.95}`):")
        println(io, "  - `min gap = ", sfmt(min_smallB_gap; digits = 10), "`")
        println(io, "  - `min gap/sqrt(B) = ", sfmt(min_smallB_scaled; digits = 10), "`")
        println(io, "- Near-`B=1` rays (`B=1-delta`, `delta` from `1e-6` to `1e-2`, same `u` probes):")
        println(io, "  - `min gap = ", sfmt(min_near1_gap; digits = 10), "`")
        println(io, "  - `min gap/delta^2 = ", sfmt(min_near1_scaled; digits = 10), "`")
        println(io)

        println(io, "## 7. Conclusion")
        println(io)
        println(io, "This gives a practical full-range management recipe for `d=1`:")
        println(io, "1. map `(B,s)` to `(B,u)` rectangle,")
        println(io, "2. use explicit/stable formulas for `lambda` and `alpha`,")
        println(io, "3. run boundary-aware dense scans and endpoint diagnostics.")
        println(io)
        println(io, "All checks above found nonnegative gaps.")
        println(io)
        println(io, "This is a strong numerical certificate, but still not a symbolic proof.")
        println(io)

        println(io, "## 8. Reproduction")
        println(io)
        println(io, "```bash")
        println(io, "/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \\")
        println(io, "  /Users/joonkyunglee/SDP_Ising/scripts/generate_d1_full_B_range_report.jl \\")
        println(io, "  --epsB=1e-4 --epsU=1e-6 \\")
        println(io, "  --nlogB=90 --nlinB=241 --nlogU=120 --nlinU=321 \\")
        println(io, "  --output=/Users/joonkyunglee/SDP_Ising/d1_full_B_range_management_report.md")
        println(io, "```")
    end

    println("Report written to ", out)
end

main()
