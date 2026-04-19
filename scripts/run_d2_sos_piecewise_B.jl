#!/usr/bin/env julia

using JuMP
using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using Printf

const Exp3 = NTuple{3, Int}

function parse_float_arg(name::String, default::Float64)
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            return parse(Float64, split(a, "=", limit = 2)[2])
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

function parse_string_arg(name::String, default::String)
    prefix = "--$(name)="
    for a in ARGS
        if startswith(a, prefix)
            return split(a, "=", limit = 2)[2]
        end
    end
    return default
end

function has_flag(name::String)
    any(==(name), ARGS)
end

function parse_int_list_arg(name::String, default::Vector{Int})
    raw = parse_string_arg(name, "")
    if isempty(raw)
        return default
    end
    vals = Int[]
    for tok in split(raw, ",")
        s = strip(tok)
        isempty(s) && continue
        push!(vals, parse(Int, s))
    end
    isempty(vals) ? default : vals
end

function build_global_polys(b, u, v)
    omu = 1 - u
    omv = 1 - v

    n2_Bλ = b * omu^2 + 2b * u * omu + b^3 * u^2
    n3_λ = b^3 * omu^3 + 3b * u * omu^2 + 3b * u^2 * omu + b^3 * u^3

    n4_α = b^6 * omv^4 +
           4b^3 * v * omv^3 +
           6b^2 * v^2 * omv^2 +
           4b^3 * v^3 * omv +
           b^6 * v^4

    n3_Bα = b^3 * omv^3 +
            3b^2 * v * omv^2 +
            3b^3 * v^2 * omv +
            b^6 * v^3

    # Denominator-cleared equality:
    # A = b * n2(lambda/b), C = b * n3(alpha/b)
    A = b^2 * omu^2 + 2u * omu + u^2
    C = b^4 * omv^3 + 3b * v * omv^2 + 3v^2 * omv + b * v^3

    fhat = n2_Bλ^6 * n4_α^3 - n3_Bα^4 * n3_λ^4
    gtilde = A^3 * n3_Bα^2 - b * n2_Bλ^3 * C^2

    return fhat, gtilde
end

function support_tuples(poly)::Set{Exp3}
    S = Set{Exp3}()
    for m in monomials(poly)
        e = exponents(m)
        push!(S, (Int(e[1]), Int(e[2]), Int(e[3])))
    end
    return S
end

function half_support(S::Set{Exp3})::Set{Exp3}
    H = Set{Exp3}()
    for e in S
        push!(H, (fld(e[1], 2), fld(e[2], 2), fld(e[3], 2)))
        push!(H, (cld(e[1], 2), cld(e[2], 2), cld(e[3], 2)))
    end
    return H
end

function shifted_half_support(S::Set{Exp3}, idx::Int, shifts::Vector{Int})::Set{Exp3}
    H = Set{Exp3}()
    for e in S, s in shifts
        if e[idx] >= s
            r = ntuple(i -> i == idx ? e[i] - s : e[i], 3)
            push!(H, (fld(r[1], 2), fld(r[2], 2), fld(r[3], 2)))
            push!(H, (cld(r[1], 2), cld(r[2], 2), cld(r[3], 2)))
        end
    end
    return H
end

function tau_support(Ef::Set{Exp3}, Eg::Set{Exp3})::Set{Exp3}
    T = Set{Exp3}()
    for e in Ef, g in Eg
        if g[1] <= e[1] && g[2] <= e[2] && g[3] <= e[3]
            push!(T, (e[1] - g[1], e[2] - g[2], e[3] - g[3]))
        end
    end
    return T
end

function expand_support(S::Set{Exp3}, caps::Exp3, level::Int)::Set{Exp3}
    T = copy(S)
    for _ in 1:level
        U = copy(T)
        for t in T
            if t[1] < caps[1]
                push!(U, (t[1] + 1, t[2], t[3]))
            end
            if t[2] < caps[2]
                push!(U, (t[1], t[2] + 1, t[3]))
            end
            if t[3] < caps[3]
                push!(U, (t[1], t[2], t[3] + 1))
            end
        end
        T = U
    end
    return T
end

function cap_tau_support(S::Set{Exp3}, max_terms::Int)::Set{Exp3}
    if max_terms <= 0 || length(S) <= max_terms
        return S
    end
    ordered = sort(collect(S); by = t -> (t[1] + t[2] + t[3], t[1], t[2], t[3]))
    return Set(ordered[1:max_terms])
end

function tuple_to_monomial(t::Exp3, vars)
    m = one(vars[1])
    if t[1] != 0
        m *= vars[1]^t[1]
    end
    if t[2] != 0
        m *= vars[2]^t[2]
    end
    if t[3] != 0
        m *= vars[3]^t[3]
    end
    return m
end

function tuples_to_monomials(S::Set{Exp3}, vars)
    ordered = sort(collect(S); by = t -> (t[1] + t[2] + t[3], t[1], t[2], t[3]))
    return [tuple_to_monomial(t, vars) for t in ordered], ordered
end

function eigmin_str(M)
    vals = eigvals(Symmetric(M))
    return @sprintf("%.3e", minimum(vals))
end

function build_uniform_intervals(bmin::Float64, bmax::Float64, n::Int)
    edges = collect(range(bmin, bmax; length = n + 1))
    return [(edges[i], edges[i + 1]) for i in 1:n]
end

function run_interval_certificate(
    L::Float64,
    U::Float64,
    vars,
    # Precomputed objects:
    fmodel,
    gmodel,
    fscale::Float64,
    Ef::Set{Exp3},
    Eg::Set{Exp3},
    S0_base::Set{Exp3},
    Sb_base::Set{Exp3},
    Su_base::Set{Exp3},
    Sv_base::Set{Exp3},
    Sτ_base::Set{Exp3};
    use_global_b_multiplier::Bool = true,
    level::Int = 0,
    max_tau_terms::Int = 2500,
    max_time::Float64 = 600.0,
    threads::Int = 0,
    silent::Bool = true,
)
    b, u, v = vars

    S0 = copy(S0_base)
    # Use shifted basis (like global b(1-b) multiplier) to avoid degree-overflow artifacts.
    SI = copy(Sb_base)
    SB = copy(Sb_base)
    Su = copy(Su_base)
    Sv = copy(Sv_base)
    Sτ = copy(Sτ_base)

    if level > 0
        S0 = expand_support(S0, (20, 8, 8), level)
        SI = expand_support(SI, (21, 8, 8), level)
        SB = expand_support(SB, (21, 8, 8), level)
        Su = expand_support(Su, (20, 7, 8), level)
        Sv = expand_support(Sv, (20, 8, 7), level)
        Sτ = expand_support(Sτ, (34, 14, 14), level)
    end

    Sτ = cap_tau_support(Sτ, max_tau_terms)

    mons0, ord0 = tuples_to_monomials(S0, vars)
    monsI, ordI = tuples_to_monomials(SI, vars)
    monsB, ordB = tuples_to_monomials(SB, vars)
    monsu, ordu = tuples_to_monomials(Su, vars)
    monsv, ordv = tuples_to_monomials(Sv, vars)
    monsτ, ordτ = tuples_to_monomials(Sτ, vars)

    println(@sprintf("Interval [%.6f, %.6f] setup", L, U))
    println("support sizes: |Ef|=$(length(Ef)), |Eg|=$(length(Eg))")
    println("basis sizes: |S0|=$(length(ord0)), |SI|=$(length(ordI)), |SB|=$(length(ordB)), |Su|=$(length(ordu)), |Sv|=$(length(ordv)), |Sτ|=$(length(ordτ))")

    gi = (b - L) * (U - b)
    gb = b * (1 - b)
    gu = u * (1 - u)
    gv = v * (1 - v)

    model = Model(Mosek.Optimizer)
    if silent
        set_silent(model)
    end
    set_optimizer_attribute(model, "MSK_DPAR_OPTIMIZER_MAX_TIME", max_time)
    if threads > 0
        set_optimizer_attribute(model, "MSK_IPAR_NUM_THREADS", threads)
    end

    @variable(model, σ0, SOSPoly(mons0))
    @variable(model, σI, SOSPoly(monsI))
    σB = nothing
    if use_global_b_multiplier
        @variable(model, σB_local, SOSPoly(monsB))
        σB = σB_local
    end
    @variable(model, σu, SOSPoly(monsu))
    @variable(model, σv, SOSPoly(monsv))
    @variable(model, τ, Poly(monsτ))
    @variable(model, γ)

    rhs = σ0 + σI * gi + σu * gu + σv * gv + τ * gmodel
    if use_global_b_multiplier
        rhs += σB * gb
    end
    @constraint(model, fmodel - γ == rhs)
    @objective(model, Max, γ)
    optimize!(model)

    status = termination_status(model)
    pstatus = primal_status(model)
    dstatus = dual_status(model)

    gamma_scaled = NaN
    gamma_original = NaN
    cert = false

    if status in (OPTIMAL, LOCALLY_SOLVED, ALMOST_OPTIMAL)
        gamma_scaled = value(γ)
        gamma_original = gamma_scaled * fscale
        cert = gamma_original >= -1e-7
        println("status = ", status, ", gamma_original = ", @sprintf("%.12e", gamma_original))
        if cert
            q0 = value.(Matrix(σ0.Q))
            qI = value.(Matrix(σI.Q))
            qu = value.(Matrix(σu.Q))
            qv = value.(Matrix(σv.Q))
            if use_global_b_multiplier
                qB = value.(Matrix(σB.Q))
                println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(QB) = ", eigmin_str(qB), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv))
            else
                println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv))
            end
        end
    elseif pstatus in (FEASIBLE_POINT, NEARLY_FEASIBLE_POINT)
        gamma_scaled = value(γ)
        gamma_original = gamma_scaled * fscale
        println("status = ", status, ", feasible gamma_original = ", @sprintf("%.12e", gamma_original), " (not certified)")
    else
        println("status = ", status, ", no feasible certificate data")
    end

    return (
        L = L,
        U = U,
        status = status,
        primal_status = pstatus,
        dual_status = dstatus,
        gamma_original = gamma_original,
        certified = cert,
        sizes = (S0 = length(ord0), SI = length(ordI), Su = length(ordu), Sv = length(ordv), Sτ = length(ordτ)),
    )
end

function main()
    bmin = parse_float_arg("bmin", 0.001)
    bmax = parse_float_arg("bmax", 0.999)
    nint = parse_int_arg("intervals", 8)
    level = parse_int_arg("level", 0)
    max_time = parse_float_arg("time", 600.0)
    threads = parse_int_arg("threads", 0)
    silent = !has_flag("--verbose")
    scale_polys = !has_flag("--no-scale")
    use_global_b_multiplier = !has_flag("--no-bmult")
    tau_schedule = parse_int_list_arg("tau_schedule", [1200, 1800, 2500])

    if !(0.0 < bmin < bmax < 1.0)
        error("Require 0 < bmin < bmax < 1.")
    end
    if nint < 1
        error("intervals must be >= 1")
    end

    @polyvar b u v
    fhat, gtilde = build_global_polys(b, u, v)
    Ef = support_tuples(fhat)
    Eg = support_tuples(gtilde)

    S0_base = half_support(Ef)
    Sb_base = shifted_half_support(Ef, 1, [1, 2])
    Su_base = shifted_half_support(Ef, 2, [1, 2])
    Sv_base = shifted_half_support(Ef, 3, [1, 2])
    Sτ_base = tau_support(Ef, Eg)

    fscale = 1.0
    gscale = 1.0
    fmodel = fhat
    gmodel = gtilde
    if scale_polys
        fscale = Float64(maximum(abs, coefficients(fhat)))
        gscale = Float64(maximum(abs, coefficients(gtilde)))
        fmodel = fhat / fscale
        gmodel = gtilde / gscale
    end

    println("Piecewise global-B run (interval multiplier on (b-L)(U-b))")
    println("range=[$bmin, $bmax], intervals=$nint, level=$level, tau_schedule=$(tau_schedule), time_per_try=$max_time")
    println("use_global_b_multiplier=$(use_global_b_multiplier)")
    println("global supports: |Ef|=$(length(Ef)), |Eg|=$(length(Eg))")
    println("base sizes: |S0|=$(length(S0_base)), |Sb|=$(length(Sb_base)), |Su|=$(length(Su_base)), |Sv|=$(length(Sv_base)), |Sτ|=$(length(Sτ_base))")
    if scale_polys
        println("scales: fscale=", @sprintf("%.3e", fscale), ", gscale=", @sprintf("%.3e", gscale))
    end

    intervals = build_uniform_intervals(bmin, bmax, nint)
    results = Any[]

    for (idx, (L, U)) in enumerate(intervals)
        println("\n=== Interval $(idx)/$(nint): [$(L), $(U)] ===")
        best = nothing
        for tau_cap in tau_schedule
            println("Trying tau cap = $(tau_cap)")
            r = run_interval_certificate(
                L, U,
                [b, u, v], fmodel, gmodel, fscale, Ef, Eg, S0_base, Sb_base, Su_base, Sv_base, Sτ_base;
                use_global_b_multiplier = use_global_b_multiplier,
                level = level,
                max_tau_terms = tau_cap,
                max_time = max_time,
                threads = threads,
                silent = silent,
            )
            best = r
            if r.certified
                println("Certified interval with tau cap = $(tau_cap)")
                break
            end
        end
        push!(results, best)
    end

    println("\n=== Summary ===")
    ncert = 0
    for (i, r) in enumerate(results)
        cmark = r.certified ? "YES" : "NO"
        if r.certified
            ncert += 1
        end
        gtxt = isnan(r.gamma_original) ? "NaN" : @sprintf("%.6e", r.gamma_original)
        println(@sprintf("[%02d] [%.6f, %.6f] status=%s certified=%s gamma=%s", i, r.L, r.U, string(r.status), cmark, gtxt))
    end
    println(@sprintf("Certified intervals: %d / %d", ncert, length(results)))
    if ncert == length(results)
        println("Overall result: FULL COVERAGE CERTIFIED on [$(bmin), $(bmax)].")
    else
        println("Overall result: PARTIAL/NO coverage. Increase time, level, or tau schedule.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
