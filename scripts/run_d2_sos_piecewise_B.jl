#!/usr/bin/env julia

using JuMP
using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using Printf
using DelimitedFiles
using Dates

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
            return String(split(a, "=", limit = 2)[2])
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

function eigmin_val(M)
    vals = eigvals(Symmetric(M))
    return Float64(minimum(vals))
end

function fmt_interval_token(x::Real)
    s = @sprintf("%.6f", Float64(x))
    s = replace(s, "-" => "m")
    s = replace(s, "." => "p")
    return s
end

function write_lines(path::String, lines)
    open(path, "w") do io
        for line in lines
            println(io, line)
        end
    end
end

function unique_dir(path::String)
    if !ispath(path)
        return path
    end
    i = 2
    while true
        cand = "$(path)_$(i)"
        if !ispath(cand)
            return cand
        end
        i += 1
    end
end

function max_abs_coeff(poly)
    cs = coefficients(poly)
    if isempty(cs)
        return 0.0
    end
    m = 0.0
    for c in cs
        v = abs(Float64(c))
        if v > m
            m = v
        end
    end
    return m
end

function dump_certificate_artifacts(
    cert_out_dir::String,
    cert_tag::String,
    L::Float64,
    U::Float64,
    tau_cap::Int,
    status,
    gamma_scaled::Float64,
    gamma_original::Float64,
    fscale::Float64,
    gscale::Float64,
    use_global_b_multiplier::Bool,
    use_half_s::Bool,
    q0,
    qI,
    qB,
    qu,
    qv,
    qUh,
    qVh,
    basis0::Vector{String},
    basisI::Vector{String},
    basisB::Vector{String},
    basisu::Vector{String},
    basisv::Vector{String},
    basisUh::Vector{String},
    basisVh::Vector{String},
    taubasis::Vector{String},
    taucoeffs::Vector{Float64},
    residual_scaled::Float64,
)
    tag = isempty(cert_tag) ? "cert" : cert_tag
    dname = "$(tag)_L$(fmt_interval_token(L))_U$(fmt_interval_token(U))_tau$(tau_cap)"
    outdir = unique_dir(joinpath(cert_out_dir, dname))
    mkpath(outdir)

    writedlm(joinpath(outdir, "Q0.tsv"), q0, '\t')
    writedlm(joinpath(outdir, "QI.tsv"), qI, '\t')
    writedlm(joinpath(outdir, "Qu.tsv"), qu, '\t')
    writedlm(joinpath(outdir, "Qv.tsv"), qv, '\t')
    write_lines(joinpath(outdir, "basis_Q0.txt"), basis0)
    write_lines(joinpath(outdir, "basis_QI.txt"), basisI)
    write_lines(joinpath(outdir, "basis_Qu.txt"), basisu)
    write_lines(joinpath(outdir, "basis_Qv.txt"), basisv)

    if use_global_b_multiplier && qB !== nothing
        writedlm(joinpath(outdir, "QB.tsv"), qB, '\t')
        write_lines(joinpath(outdir, "basis_QB.txt"), basisB)
    end
    if use_half_s && qUh !== nothing && qVh !== nothing
        writedlm(joinpath(outdir, "QUh.tsv"), qUh, '\t')
        writedlm(joinpath(outdir, "QVh.tsv"), qVh, '\t')
        write_lines(joinpath(outdir, "basis_QUh.txt"), basisUh)
        write_lines(joinpath(outdir, "basis_QVh.txt"), basisVh)
    end

    write_lines(joinpath(outdir, "tau_basis.txt"), taubasis)
    open(joinpath(outdir, "tau_coeffs.tsv"), "w") do io
        println(io, "index\tmonomial\tcoefficient")
        for i in eachindex(taubasis)
            println(io, i, '\t', taubasis[i], '\t', @sprintf("%.17e", taucoeffs[i]))
        end
    end

    open(joinpath(outdir, "metadata.txt"), "w") do io
        println(io, "generated_at=", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "status=", status)
        println(io, "interval_L=", @sprintf("%.12f", L))
        println(io, "interval_U=", @sprintf("%.12f", U))
        println(io, "tau_cap=", tau_cap)
        println(io, "gamma_scaled=", @sprintf("%.17e", gamma_scaled))
        println(io, "gamma_original=", @sprintf("%.17e", gamma_original))
        println(io, "fscale=", @sprintf("%.17e", fscale))
        println(io, "gscale=", @sprintf("%.17e", gscale))
        println(io, "residual_scaled_max_abs_coeff=", @sprintf("%.17e", residual_scaled))
        println(io, "use_global_b_multiplier=", use_global_b_multiplier)
        println(io, "use_half_s=", use_half_s)
        println(io, "min_eig_Q0=", @sprintf("%.17e", eigmin_val(q0)))
        println(io, "min_eig_QI=", @sprintf("%.17e", eigmin_val(qI)))
        println(io, "min_eig_Qu=", @sprintf("%.17e", eigmin_val(qu)))
        println(io, "min_eig_Qv=", @sprintf("%.17e", eigmin_val(qv)))
        if use_global_b_multiplier && qB !== nothing
            println(io, "min_eig_QB=", @sprintf("%.17e", eigmin_val(qB)))
        end
        if use_half_s && qUh !== nothing && qVh !== nothing
            println(io, "min_eig_QUh=", @sprintf("%.17e", eigmin_val(qUh)))
            println(io, "min_eig_QVh=", @sprintf("%.17e", eigmin_val(qVh)))
        end
    end

    return outdir
end

function build_uniform_intervals(bmin::Float64, bmax::Float64, n::Int)
    edges = collect(range(bmin, bmax; length = n + 1))
    return [(edges[i], edges[i + 1]) for i in 1:n]
end

function cert_status(status)
    status in (OPTIMAL, LOCALLY_SOLVED, ALMOST_OPTIMAL)
end

function run_interval_certificate(
    L::Float64,
    U::Float64,
    vars,
    # Precomputed objects:
    fmodel,
    gmodel,
    fscale::Float64,
    gscale::Float64,
    Ef::Set{Exp3},
    Eg::Set{Exp3},
    S0_base::Set{Exp3},
    Sb_base::Set{Exp3},
    Su_base::Set{Exp3},
    Sv_base::Set{Exp3},
    Suh_base::Set{Exp3},
    Svh_base::Set{Exp3},
    Sτ_base::Set{Exp3};
    use_global_b_multiplier::Bool = true,
    use_half_s::Bool = false,
    feas_eps_original::Float64 = 0.0,
    dump_full_cert::Bool = false,
    cert_out_dir::String = "",
    cert_tag::String = "",
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
    Suh = copy(Suh_base)
    Svh = copy(Svh_base)
    Sτ = copy(Sτ_base)

    if level > 0
        S0 = expand_support(S0, (20, 8, 8), level)
        SI = expand_support(SI, (21, 8, 8), level)
        SB = expand_support(SB, (21, 8, 8), level)
        Su = expand_support(Su, (20, 7, 8), level)
        Sv = expand_support(Sv, (20, 8, 7), level)
        Suh = expand_support(Suh, (20, 7, 8), level)
        Svh = expand_support(Svh, (20, 8, 7), level)
        Sτ = expand_support(Sτ, (34, 14, 14), level)
    end

    Sτ = cap_tau_support(Sτ, max_tau_terms)

    mons0, ord0 = tuples_to_monomials(S0, vars)
    monsI, ordI = tuples_to_monomials(SI, vars)
    monsB, ordB = tuples_to_monomials(SB, vars)
    monsu, ordu = tuples_to_monomials(Su, vars)
    monsv, ordv = tuples_to_monomials(Sv, vars)
    monsuh, orduh = tuples_to_monomials(Suh, vars)
    monsvh, ordvh = tuples_to_monomials(Svh, vars)
    monsτ, ordτ = tuples_to_monomials(Sτ, vars)
    basis0 = string.(mons0)
    basisI = string.(monsI)
    basisB = string.(monsB)
    basisu = string.(monsu)
    basisv = string.(monsv)
    basisUh = string.(monsuh)
    basisVh = string.(monsvh)
    taubasis = string.(monsτ)

    println(@sprintf("Interval [%.6f, %.6f] setup", L, U))
    println("support sizes: |Ef|=$(length(Ef)), |Eg|=$(length(Eg))")
    println("basis sizes: |S0|=$(length(ord0)), |SI|=$(length(ordI)), |SB|=$(length(ordB)), |Su|=$(length(ordu)), |Sv|=$(length(ordv)), |Sτ|=$(length(ordτ))")
    if use_half_s
        println("half-s bases: |Suh|=$(length(orduh)), |Svh|=$(length(ordvh))")
    end

    gi = (b - L) * (U - b)
    gb = b * (1 - b)
    gu = u * (1 - u)
    gv = v * (1 - v)
    # Nonnegative center multipliers on [0,1].
    guh = u * (1 - u) * (2u - 1)^2
    gvh = v * (1 - v) * (2v - 1)^2

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
    σuh = nothing
    σvh = nothing
    if use_half_s
        @variable(model, σuh_local, SOSPoly(monsuh))
        @variable(model, σvh_local, SOSPoly(monsvh))
        σuh = σuh_local
        σvh = σvh_local
    end
    @variable(model, τ, Poly(monsτ))
    @variable(model, γ)
    use_feas_mode = feas_eps_original > 0.0
    gamma_target_scaled = use_feas_mode ? (-feas_eps_original / fscale) : NaN

    rhs = σ0 + σI * gi + σu * gu + σv * gv + τ * gmodel
    if use_global_b_multiplier
        rhs += σB * gb
    end
    if use_half_s
        rhs += σuh * guh + σvh * gvh
    end
    @constraint(model, fmodel - γ == rhs)
    if use_feas_mode
        @constraint(model, γ == gamma_target_scaled)
        @objective(model, Max, 0.0)
    else
        @objective(model, Max, γ)
    end
    optimize!(model)

    status = termination_status(model)
    pstatus = primal_status(model)
    dstatus = dual_status(model)

    gamma_scaled = NaN
    gamma_original = NaN
    cert = false
    cert_dir = ""

    if cert_status(status)
        gamma_scaled = value(γ)
        gamma_original = gamma_scaled * fscale
        cert = use_feas_mode ? true : (gamma_original >= -1e-7)
        if use_feas_mode
            println(
                "status = ",
                status,
                ", feasibility target gamma_original = ",
                @sprintf("%.12e", -feas_eps_original),
                ", attained gamma_original = ",
                @sprintf("%.12e", gamma_original),
            )
        else
            println("status = ", status, ", gamma_original = ", @sprintf("%.12e", gamma_original))
        end
        if cert
            q0 = value.(Matrix(σ0.Q))
            qI = value.(Matrix(σI.Q))
            qu = value.(Matrix(σu.Q))
            qv = value.(Matrix(σv.Q))
            qB = nothing
            qUh = nothing
            qVh = nothing
            if use_half_s
                qUh = value.(Matrix(σuh.Q))
                qVh = value.(Matrix(σvh.Q))
            end
            if use_global_b_multiplier
                qB = value.(Matrix(σB.Q))
                if use_half_s
                    println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(QB) = ", eigmin_str(qB), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv), ", min_eig(QUh) = ", eigmin_str(qUh), ", min_eig(QVh) = ", eigmin_str(qVh))
                else
                    println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(QB) = ", eigmin_str(qB), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv))
                end
            else
                if use_half_s
                    println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv), ", min_eig(QUh) = ", eigmin_str(qUh), ", min_eig(QVh) = ", eigmin_str(qVh))
                else
                    println("min_eig(Q0) = ", eigmin_str(q0), ", min_eig(QI) = ", eigmin_str(qI), ", min_eig(Qu) = ", eigmin_str(qu), ", min_eig(Qv) = ", eigmin_str(qv))
                end
            end

            if dump_full_cert && !isempty(cert_out_dir)
                poly_rhs = value(σ0) + value(σI) * gi + value(σu) * gu + value(σv) * gv + value(τ) * gmodel
                if use_global_b_multiplier
                    poly_rhs += value(σB) * gb
                end
                if use_half_s
                    poly_rhs += value(σuh) * guh + value(σvh) * gvh
                end
                poly_lhs = fmodel - gamma_scaled
                residual_scaled = max_abs_coeff(poly_lhs - poly_rhs)
                taucoeffs = [Float64(value(MultivariatePolynomials.coefficient(τ, m))) for m in monsτ]

                cert_dir = dump_certificate_artifacts(
                    cert_out_dir,
                    cert_tag,
                    L,
                    U,
                    max_tau_terms,
                    status,
                    gamma_scaled,
                    gamma_original,
                    fscale,
                    gscale,
                    use_global_b_multiplier,
                    use_half_s,
                    q0,
                    qI,
                    qB,
                    qu,
                    qv,
                    qUh,
                    qVh,
                    basis0,
                    basisI,
                    basisB,
                    basisu,
                    basisv,
                    basisUh,
                    basisVh,
                    taubasis,
                    taucoeffs,
                    residual_scaled,
                )
                println("certificate_artifacts_dir = ", cert_dir)
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
        cert_dir = cert_dir,
        sizes = (S0 = length(ord0), SI = length(ordI), Su = length(ordu), Sv = length(ordv), Sτ = length(ordτ)),
    )
end

function run_interval_with_schedule(
    L::Float64,
    U::Float64,
    vars,
    fmodel,
    gmodel,
    fscale::Float64,
    gscale::Float64,
    Ef::Set{Exp3},
    Eg::Set{Exp3},
    S0_base::Set{Exp3},
    Sb_base::Set{Exp3},
    Su_base::Set{Exp3},
    Sv_base::Set{Exp3},
    Suh_base::Set{Exp3},
    Svh_base::Set{Exp3},
    Sτ_base::Set{Exp3},
    tau_schedule::Vector{Int};
    use_global_b_multiplier::Bool = true,
    use_half_s::Bool = false,
    feas_eps_original::Float64 = 0.0,
    dump_full_cert::Bool = false,
    cert_out_dir::String = "",
    cert_tag::String = "",
    level::Int = 0,
    max_time::Float64 = 600.0,
    threads::Int = 0,
    silent::Bool = true,
)
    best = nothing
    for tau_cap in tau_schedule
        println("Trying tau cap = $(tau_cap)")
        r = run_interval_certificate(
            L, U,
            vars, fmodel, gmodel, fscale, gscale, Ef, Eg, S0_base, Sb_base, Su_base, Sv_base, Suh_base, Svh_base, Sτ_base;
            use_global_b_multiplier = use_global_b_multiplier,
            use_half_s = use_half_s,
            feas_eps_original = feas_eps_original,
            dump_full_cert = dump_full_cert,
            cert_out_dir = cert_out_dir,
            cert_tag = cert_tag,
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
    return best
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
    half_s_requested = has_flag("--half-s")
    half_s_experimental = has_flag("--half-s-experimental")
    use_half_s = half_s_requested && half_s_experimental
    tau_schedule = parse_int_list_arg("tau_schedule", [1200, 1800, 2500])
    feas_eps_original = parse_float_arg("feas_eps", 0.0)
    adaptive_depth = parse_int_arg("adaptive_depth", 0)
    adaptive_min_width = parse_float_arg("adaptive_min_width", 0.0)
    cert_out_dir = parse_string_arg("cert_out_dir", "")
    cert_tag_base = parse_string_arg("cert_tag", "")
    dump_full_cert = has_flag("--dump_full_cert") || !isempty(cert_out_dir)
    if dump_full_cert && isempty(cert_out_dir)
        cert_out_dir = joinpath(dirname(@__DIR__), "certificates", "piecewise")
    end

    if !(0.0 < bmin < bmax < 1.0)
        error("Require 0 < bmin < bmax < 1.")
    end
    if nint < 1
        error("intervals must be >= 1")
    end
    if feas_eps_original < 0
        error("feas_eps must be >= 0")
    end
    if adaptive_depth < 0
        error("adaptive_depth must be >= 0")
    end
    if adaptive_min_width < 0
        error("adaptive_min_width must be >= 0")
    end

    @polyvar b u v
    fhat, gtilde = build_global_polys(b, u, v)
    Ef = support_tuples(fhat)
    Eg = support_tuples(gtilde)

    S0_base = half_support(Ef)
    Sb_base = shifted_half_support(Ef, 1, [1, 2])
    Su_base = shifted_half_support(Ef, 2, [1, 2])
    Sv_base = shifted_half_support(Ef, 3, [1, 2])
    Suh_base = shifted_half_support(Ef, 2, [3, 4])
    Svh_base = shifted_half_support(Ef, 3, [3, 4])
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
    if half_s_requested && !half_s_experimental
        println("NOTE: --half-s requested but disabled by default; use --half-s-experimental to enable the heavier center-multiplier mode.")
    end
    println("use_half_s=$(use_half_s)")
    println("feasibility_mode=$(feas_eps_original > 0.0), feas_eps_original=$(feas_eps_original)")
    println("adaptive_depth=$(adaptive_depth), adaptive_min_width=$(adaptive_min_width)")
    println("dump_full_cert=$(dump_full_cert)")
    if dump_full_cert
        mkpath(cert_out_dir)
        println("cert_out_dir=$(cert_out_dir)")
        println("cert_tag_base=$(isempty(cert_tag_base) ? "interval" : cert_tag_base)")
    end
    println("global supports: |Ef|=$(length(Ef)), |Eg|=$(length(Eg))")
    println("base sizes: |S0|=$(length(S0_base)), |Sb|=$(length(Sb_base)), |Su|=$(length(Su_base)), |Sv|=$(length(Sv_base)), |Sτ|=$(length(Sτ_base))")
    if use_half_s
        println("half-s extras: |Suh|=$(length(Suh_base)), |Svh|=$(length(Svh_base)) (via u(1-u)(2u-1)^2, v(1-v)(2v-1)^2)")
    end
    if scale_polys
        println("scales: fscale=", @sprintf("%.3e", fscale), ", gscale=", @sprintf("%.3e", gscale))
    end

    intervals = build_uniform_intervals(bmin, bmax, nint)
    results = Any[]

    if adaptive_depth == 0
        for (idx, (L, U)) in enumerate(intervals)
            println("\n=== Interval $(idx)/$(nint): [$(L), $(U)] ===")
            cert_tag = isempty(cert_tag_base) ? @sprintf("interval_i%02d", idx) : @sprintf("%s_i%02d", cert_tag_base, idx)
            r = run_interval_with_schedule(
                L, U,
                [b, u, v], fmodel, gmodel, fscale, gscale, Ef, Eg, S0_base, Sb_base, Su_base, Sv_base, Suh_base, Svh_base, Sτ_base, tau_schedule;
                use_global_b_multiplier = use_global_b_multiplier,
                use_half_s = use_half_s,
                feas_eps_original = feas_eps_original,
                dump_full_cert = dump_full_cert,
                cert_out_dir = cert_out_dir,
                cert_tag = cert_tag,
                level = level,
                max_time = max_time,
                threads = threads,
                silent = silent,
            )
            push!(results, merge(r, (depth = 0,)))
        end
    else
        queue = [(L = it[1], U = it[2], depth = 0) for it in intervals]
        leaf_idx = 0
        while !isempty(queue)
            task = popfirst!(queue)
            leaf_idx += 1
            println("\n=== Adaptive node $(leaf_idx): depth=$(task.depth), interval=[$(task.L), $(task.U)] ===")
            cert_tag = isempty(cert_tag_base) ? @sprintf("node_%03d_d%d", leaf_idx, task.depth) : @sprintf("%s_node_%03d_d%d", cert_tag_base, leaf_idx, task.depth)
            r = run_interval_with_schedule(
                task.L, task.U,
                [b, u, v], fmodel, gmodel, fscale, gscale, Ef, Eg, S0_base, Sb_base, Su_base, Sv_base, Suh_base, Svh_base, Sτ_base, tau_schedule;
                use_global_b_multiplier = use_global_b_multiplier,
                use_half_s = use_half_s,
                feas_eps_original = feas_eps_original,
                dump_full_cert = dump_full_cert,
                cert_out_dir = cert_out_dir,
                cert_tag = cert_tag,
                level = level,
                max_time = max_time,
                threads = threads,
                silent = silent,
            )
            width = task.U - task.L
            should_split = (!r.certified) && (task.depth < adaptive_depth) && (width > adaptive_min_width)
            if should_split
                mid = 0.5 * (task.L + task.U)
                println(@sprintf("Adaptive split triggered: [%.6f, %.6f] -> [%.6f, %.6f] and [%.6f, %.6f]", task.L, task.U, task.L, mid, mid, task.U))
                pushfirst!(queue, (L = mid, U = task.U, depth = task.depth + 1))
                pushfirst!(queue, (L = task.L, U = mid, depth = task.depth + 1))
            else
                push!(results, merge(r, (depth = task.depth,)))
            end
        end
        sort!(results; by = r -> r.L)
    end

    println("\n=== Summary ===")
    ncert = 0
    for (i, r) in enumerate(results)
        cmark = r.certified ? "YES" : "NO"
        if r.certified
            ncert += 1
        end
        gtxt = isnan(r.gamma_original) ? "NaN" : @sprintf("%.6e", r.gamma_original)
        println(@sprintf("[%02d] [%.6f, %.6f] depth=%d status=%s certified=%s gamma=%s", i, r.L, r.U, r.depth, string(r.status), cmark, gtxt))
        if !isempty(r.cert_dir)
            println("     cert_dir=$(r.cert_dir)")
        end
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
