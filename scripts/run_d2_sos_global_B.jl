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

function has_flag(name::String)
    any(==(name), ARGS)
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

    # Clear 1/b denominators from the equality constraint:
    # gtilde = b^3 * ghat
    A = b^2 * omu^2 + 2u * omu + u^2
    C = b^4 * omv^3 + 3b * v * omv^2 + 3v^2 * omv + b * v^3

    fhat = n2_Bλ^6 * n4_α^3 - n3_Bα^4 * n3_λ^4
    gtilde = A^3 * n3_Bα^2 - b * n2_Bλ^3 * C^2

    gu = u * (1 - u)
    gv = v * (1 - v)
    gb = b * (1 - b)

    return fhat, gtilde, gu, gv, gb
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

function max_exp(S::Set{Exp3})::Exp3
    mb = maximum(t -> t[1], S)
    mu = maximum(t -> t[2], S)
    mv = maximum(t -> t[3], S)
    return (mb, mu, mv)
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

function cap_tau_support(S::Set{Exp3}, max_terms::Int)::Set{Exp3}
    if max_terms <= 0 || length(S) <= max_terms
        return S
    end
    ordered = sort(collect(S); by = t -> (t[1] + t[2] + t[3], t[1], t[2], t[3]))
    return Set(ordered[1:max_terms])
end

function eigmin_str(M)
    vals = eigvals(Symmetric(M))
    return @sprintf("%.3e", minimum(vals))
end

function run_global_certificate(;
    level::Int = 0,
    use_b_multiplier::Bool = true,
    max_tau_terms::Int = 0,
    scale_polys::Bool = true,
    max_time::Float64 = 900.0,
    threads::Int = 0,
    silent::Bool = true,
)
    @polyvar b u v
    vars = [b, u, v]
    fhat, gtilde, gu, gv, gb = build_global_polys(b, u, v)

    Ef = support_tuples(fhat)
    Eg = support_tuples(gtilde)

    S0 = half_support(Ef)
    Su = shifted_half_support(Ef, 2, [1, 2])
    Sv = shifted_half_support(Ef, 3, [1, 2])
    Sb = shifted_half_support(Ef, 1, [1, 2])
    St = tau_support(Ef, Eg)

    # Controlled expansion: increases flexibility while keeping basis sparse.
    if level > 0
        S0 = expand_support(S0, (20, 8, 8), level)
        Su = expand_support(Su, (20, 7, 8), level)
        Sv = expand_support(Sv, (20, 8, 7), level)
        Sb = expand_support(Sb, (21, 8, 8), level)
        St = expand_support(St, (34, 14, 14), level)
    end

    St = cap_tau_support(St, max_tau_terms)

    mons0, ord0 = tuples_to_monomials(S0, vars)
    monsu, ordu = tuples_to_monomials(Su, vars)
    monsv, ordv = tuples_to_monomials(Sv, vars)
    monsb, ordb = tuples_to_monomials(Sb, vars)
    monst, ordt = tuples_to_monomials(St, vars)

    fscale = 1.0
    gscale = 1.0
    fmodel = fhat
    gmodel = gtilde
    if scale_polys
        fscale = maximum(abs, coefficients(fhat))
        gscale = maximum(abs, coefficients(gtilde))
        fmodel = fhat / fscale
        gmodel = gtilde / gscale
    end

    println("Global-B SOS setup")
    println("support sizes: |Ef|=$(length(Ef)), |Eg|=$(length(Eg))")
    println("basis sizes: |S0|=$(length(ord0)), |Su|=$(length(ordu)), |Sv|=$(length(ordv)), |Sb|=$(length(ordb)), |St|=$(length(ordt))")
    println("level=$(level), use_b_multiplier=$(use_b_multiplier), max_tau_terms=$(max_tau_terms), scale_polys=$(scale_polys)")
    if scale_polys
        println("scales: fscale=", @sprintf("%.3e", fscale), ", gscale=", @sprintf("%.3e", gscale))
    end

    model = Model(Mosek.Optimizer)
    if silent
        set_silent(model)
    end
    set_optimizer_attribute(model, "MSK_DPAR_OPTIMIZER_MAX_TIME", max_time)
    if threads > 0
        set_optimizer_attribute(model, "MSK_IPAR_NUM_THREADS", threads)
    end

    @variable(model, σ0, SOSPoly(mons0))
    @variable(model, σu, SOSPoly(monsu))
    @variable(model, σv, SOSPoly(monsv))
    σb = nothing
    if use_b_multiplier
        @variable(model, σb_local, SOSPoly(monsb))
        σb = σb_local
    end
    @variable(model, τ, Poly(monst))
    @variable(model, γ)

    rhs = σ0 + σu * gu + σv * gv + τ * gmodel
    if use_b_multiplier
        rhs += σb * gb
    end

    @constraint(model, fmodel - γ == rhs)
    @objective(model, Max, γ)

    optimize!(model)

    status = termination_status(model)
    pstatus = primal_status(model)
    dstatus = dual_status(model)
    println("status = ", status)
    println("primal_status = ", pstatus)
    println("dual_status = ", dstatus)

    if status in (OPTIMAL, LOCALLY_SOLVED, ALMOST_OPTIMAL)
        gval_scaled = value(γ)
        gval = gval_scaled * fscale
        println("gamma_lower_bound_scaled = ", @sprintf("%.12e", gval_scaled))
        println("gamma_lower_bound_original = ", @sprintf("%.12e", gval))
        if gval >= -1e-7
            println("certificate_result = NONNEGATIVE (within tolerance)")
        else
            println("certificate_result = NEGATIVE LOWER BOUND")
        end
        q0 = value.(Matrix(σ0.Q))
        qu = value.(Matrix(σu.Q))
        qv = value.(Matrix(σv.Q))
        println("min_eig(Q0) = ", eigmin_str(q0))
        println("min_eig(Qu) = ", eigmin_str(qu))
        println("min_eig(Qv) = ", eigmin_str(qv))
        if use_b_multiplier
            qb = value.(Matrix(σb.Q))
            println("min_eig(Qb) = ", eigmin_str(qb))
        end
    elseif pstatus in (FEASIBLE_POINT, NEARLY_FEASIBLE_POINT)
        gval_scaled = value(γ)
        gval = gval_scaled * fscale
        println("candidate_gamma_scaled = ", @sprintf("%.12e", gval_scaled))
        println("candidate_gamma_original = ", @sprintf("%.12e", gval))
        println("Solver stopped before optimality; this is not a valid certificate yet.")
    else
        println("No certificate found at this relaxation.")
    end
end

function main()
    level = parse_int_arg("level", 0)
    max_tau_terms = parse_int_arg("max_tau", 0)
    max_time = parse_float_arg("time", 900.0)
    threads = parse_int_arg("threads", 0)
    silent = !has_flag("--verbose")
    use_b_multiplier = !has_flag("--no-bmult")
    scale_polys = !has_flag("--no-scale")

    run_global_certificate(
        level = level,
        use_b_multiplier = use_b_multiplier,
        max_tau_terms = max_tau_terms,
        scale_polys = scale_polys,
        max_time = max_time,
        threads = threads,
        silent = silent,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
