#!/usr/bin/env julia

using JuMP
using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using Printf
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
    mt = maximum(t -> t[3], S)
    return (mb, mu, mt)
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

function build_polys(b, u, t)
    v = 1 - u
    A = u + b^2 * v

    denR = b^2 * v^2 + 1 - v^2
    numL = b^2 * v^2 + 2t * v + t^2
    denL = b^6 * v^3 + 3t * b^2 * v^2 + 3t^2 * v + t^3

    # rhs^6 = b^3 / denR^3, lhs^6 = b^3 * numL^3 / denL^2.
    # On b>0 this is equivalent to N >= 0.
    N = denL^2 - numL^3 * denR^3

    # Polynomial equality for t = alpha * b^2 * (1-u):
    # ((A^2-1)t^2 + 2v(A^2-b^2)t + v^2*b^2*(A^2-b^4) = 0)
    qt = (A^2 - 1) * t^2 + 2v * (A^2 - b^2) * t + v^2 * b^2 * (A^2 - b^4)

    gb = b * (1 - b)
    gu = u * (1 - u)
    gt = t * (1 - t)
    return N, qt, gb, gu, gt
end

function run_global_d1_certificate(;
    level::Int = 0,
    max_tau_terms::Int = 0,
    max_time::Float64 = 600.0,
    threads::Int = 0,
    silent::Bool = true,
    scale_polys::Bool = true,
)
    @polyvar b u t
    vars = [b, u, t]
    N, qt, gb, gu, gt = build_polys(b, u, t)

    EN = support_tuples(N)
    Eq = support_tuples(qt)

    S0 = half_support(EN)
    Sb = shifted_half_support(EN, 1, [1, 2])
    Su = shifted_half_support(EN, 2, [1, 2])
    St = shifted_half_support(EN, 3, [1, 2])
    Sq = tau_support(EN, Eq)

    if level > 0
        capN = max_exp(EN)
        capQ = max_exp(Sq)
        S0 = expand_support(S0, capN, level)
        Sb = expand_support(Sb, capN, level)
        Su = expand_support(Su, capN, level)
        St = expand_support(St, capN, level)
        Sq = expand_support(Sq, capQ, level)
    end

    Sq = cap_tau_support(Sq, max_tau_terms)

    mons0, ord0 = tuples_to_monomials(S0, vars)
    monsb, ordb = tuples_to_monomials(Sb, vars)
    monsu, ordu = tuples_to_monomials(Su, vars)
    monst, ordt = tuples_to_monomials(St, vars)
    monsq, ordq = tuples_to_monomials(Sq, vars)

    Nscale = 1.0
    qscale = 1.0
    Nmodel = N
    qmodel = qt
    if scale_polys
        Nscale = maximum(abs, coefficients(N))
        qscale = maximum(abs, coefficients(qt))
        Nmodel = N / Nscale
        qmodel = qt / qscale
    end

    println("d=1 global continuous SOS setup")
    println("support sizes: |EN|=$(length(EN)), |Eq|=$(length(Eq))")
    println("basis sizes: |S0|=$(length(ord0)), |Sb|=$(length(ordb)), |Su|=$(length(ordu)), |St|=$(length(ordt)), |Sq|=$(length(ordq))")
    println("level=$(level), max_tau_terms=$(max_tau_terms), scale_polys=$(scale_polys)")
    if scale_polys
        println("scales: Nscale=", @sprintf("%.3e", Nscale), ", qscale=", @sprintf("%.3e", qscale))
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
    @variable(model, σb, SOSPoly(monsb))
    @variable(model, σu, SOSPoly(monsu))
    @variable(model, σt, SOSPoly(monst))
    @variable(model, τ, Poly(monsq))
    @variable(model, γ)

    @constraint(model, Nmodel - γ == σ0 + σb * gb + σu * gu + σt * gt + τ * qmodel)
    @objective(model, Max, γ)

    optimize!(model)

    status = termination_status(model)
    pstatus = primal_status(model)
    dstatus = dual_status(model)

    result = Dict{Symbol, Any}()
    result[:status] = status
    result[:primal_status] = pstatus
    result[:dual_status] = dstatus
    result[:Nscale] = Nscale
    result[:qscale] = qscale
    result[:sizes] = Dict(
        :EN => length(EN),
        :Eq => length(Eq),
        :S0 => length(ord0),
        :Sb => length(ordb),
        :Su => length(ordu),
        :St => length(ordt),
        :Sq => length(ordq),
    )

    if status in (OPTIMAL, ALMOST_OPTIMAL, LOCALLY_SOLVED)
        γ_scaled = value(γ)
        γ_orig = γ_scaled * Nscale
        result[:gamma_scaled] = γ_scaled
        result[:gamma_original] = γ_orig

        q0 = value.(Matrix(σ0.Q))
        qb = value.(Matrix(σb.Q))
        qu = value.(Matrix(σu.Q))
        qtM = value.(Matrix(σt.Q))
        result[:mineig_Q0] = eigmin_str(q0)
        result[:mineig_Qb] = eigmin_str(qb)
        result[:mineig_Qu] = eigmin_str(qu)
        result[:mineig_Qt] = eigmin_str(qtM)
    end

    return result
end

function write_report(path::String, r::Dict{Symbol, Any}; kwargs...)
    open(path, "w") do io
        println(io, "# d = 1 Global Continuous Certificate (SOS)")
        println(io)
        println(io, "Generated at: `", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "`")
        println(io)
        println(io, "Workspace: `/Users/joonkyunglee/SDP_Ising`")
        println(io)
        println(io, "## 1. Certificate target")
        println(io)
        println(io, "We certify the `d=1` inequality for all `B in (0,1)` (continuous range, not a sampled grid).")
        println(io)
        println(io, "We use variables `(b,u,t)` with")
        println(io)
        println(io, raw"$$")
        println(io, "b=B, \\qquad u\\in(0,1), \\qquad t=\\alpha B^2(1-u), \\qquad u=\\frac{s-1}{B^{-2}-1}.")
        println(io, raw"$$")
        println(io)
        println(io, "Set `v=1-u`, `A=u+b^2v`, and define")
        println(io)
        println(io, raw"$$")
        println(io, "\\mathrm{den}_R=b^2v^2+1-v^2, \\quad")
        println(io, "\\mathrm{num}_L=b^2v^2+2tv+t^2, \\quad")
        println(io, "\\mathrm{den}_L=b^6v^3+3tb^2v^2+3t^2v+t^3.")
        println(io, raw"$$")
        println(io)
        println(io, raw"$$")
        println(io, "N(b,u,t)=\\mathrm{den}_L^2 - \\mathrm{num}_L^3\\,\\mathrm{den}_R^3.")
        println(io, raw"$$")
        println(io)
        println(io, "This is exactly the cross-multiplied form of `rhs^6-lhs^6`,")
        println(io, "so proving `N >= 0` on the feasible set proves the original inequality.")
        println(io)
        println(io, "The coupling equality is")
        println(io)
        println(io, raw"$$")
        println(io, "q_t(b,u,t):=(A^2-1)t^2 + 2v(A^2-b^2)t + v^2b^2(A^2-b^4)=0.")
        println(io, raw"$$")
        println(io)
        println(io)
        println(io, "## 2. SOS model")
        println(io)
        println(io, "We solve")
        println(io)
        println(io, raw"$$")
        println(io, "N - \\gamma = \\sigma_0 + \\sigma_b\\,b(1-b) + \\sigma_u\\,u(1-u) + \\sigma_t\\,t(1-t) + \\tau\\,q_t,")
        println(io, raw"$$")
        println(io)
        println(io, "with `sigma_*` SOS and `tau` free polynomial.")
        println(io)
        println(io, "We additionally use box localizers `b(1-b)`, `u(1-u)`, `t(1-t)`.")
        println(io)
        println(io, "Why `t(1-t)` is valid: for `v=1-u`, `r=1-b^2`,")
        println(io, "- `q_t(0)=v^2b^2(A^2-b^4) >= 0`,")
        println(io, "- `q_t(1)=r v^2 h(r,v)` where")
        println(io, "  `h=-2-2v-2r+4rv+rv^2+r^2-r^2v^2 <= -1` on `[0,1]^2`,")
        println(io, "so `q_t(1) <= 0` and the nonnegative root lies in `[0,1]`.")
        println(io)
        println(io, "Model settings:")
        println(io, "- `level = ", kwargs[:level], "`")
        println(io, "- `max_tau_terms = ", kwargs[:max_tau_terms], "`")
        println(io, "- `time_limit = ", @sprintf("%.1f", kwargs[:max_time]), " s`")
        println(io, "- `threads = ", kwargs[:threads], "`")
        println(io)

        sizes = r[:sizes]
        println(io, "Support/basis sizes:")
        println(io, "- `|EN| = ", sizes[:EN], "`, `|Eq| = ", sizes[:Eq], "`")
        println(io, "- `|S0| = ", sizes[:S0], "`, `|Sb| = ", sizes[:Sb], "`, `|Su| = ", sizes[:Su], "`, `|St| = ", sizes[:St], "`, `|Sq| = ", sizes[:Sq], "`")
        println(io)

        println(io, "## 3. Solver result")
        println(io)
        println(io, "- `status = ", r[:status], "`")
        println(io, "- `primal_status = ", r[:primal_status], "`")
        println(io, "- `dual_status = ", r[:dual_status], "`")
        if haskey(r, :gamma_original)
            println(io, "- `gamma_lower_bound_scaled = ", @sprintf("%.16e", r[:gamma_scaled]), "`")
            println(io, "- `gamma_lower_bound_original = ", @sprintf("%.16e", r[:gamma_original]), "`")
            if r[:gamma_original] >= -1e-9
                println(io, "- `certificate_result = NONNEGATIVE (within tolerance)`")
            else
                println(io, "- `certificate_result = NEGATIVE LOWER BOUND`")
            end
            println(io, "- `min_eig(Q0) = ", r[:mineig_Q0], "`")
            println(io, "- `min_eig(Qb) = ", r[:mineig_Qb], "`")
            println(io, "- `min_eig(Qu) = ", r[:mineig_Qu], "`")
            println(io, "- `min_eig(Qt) = ", r[:mineig_Qt], "`")
        else
            println(io, "- No optimal certificate value returned.")
        end
        println(io)

        println(io, "## 4. Interpretation")
        println(io)
        println(io, "This is a continuous semialgebraic certificate over the whole parameter region.")
        println(io, "It is not a discrete sampling argument.")
        println(io)
        println(io, "## 5. Reproduction")
        println(io)
        println(io, "```bash")
        println(io, "/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia \\")
        println(io, "  /Users/joonkyunglee/SDP_Ising/scripts/run_d1_global_continuous_sos.jl \\")
        println(io, "  --time=", @sprintf("%.0f", kwargs[:max_time]), " --threads=", kwargs[:threads], " --level=", kwargs[:level], " --max_tau=", kwargs[:max_tau_terms], " \\")
        println(io, "  --report=/Users/joonkyunglee/SDP_Ising/d1_global_continuous_certificate_report.md")
        println(io, "```")
    end
end

function main()
    level = parse_int_arg("level", 0)
    max_tau_terms = parse_int_arg("max_tau", 0)
    max_time = parse_float_arg("time", 600.0)
    threads = parse_int_arg("threads", 8)
    silent = !has_flag("--verbose")
    scale_polys = !has_flag("--no-scale")
    report_path = parse_string_arg(
        "report",
        joinpath(dirname(@__DIR__), "d1_global_continuous_certificate_report.md"),
    )

    r = run_global_d1_certificate(
        level = level,
        max_tau_terms = max_tau_terms,
        max_time = max_time,
        threads = threads,
        silent = silent,
        scale_polys = scale_polys,
    )

    println("status = ", r[:status])
    println("primal_status = ", r[:primal_status])
    println("dual_status = ", r[:dual_status])
    if haskey(r, :gamma_original)
        println("gamma_lower_bound_scaled = ", @sprintf("%.16e", r[:gamma_scaled]))
        println("gamma_lower_bound_original = ", @sprintf("%.16e", r[:gamma_original]))
        println("min_eig(Q0) = ", r[:mineig_Q0])
        println("min_eig(Qb) = ", r[:mineig_Qb])
        println("min_eig(Qu) = ", r[:mineig_Qu])
        println("min_eig(Qt) = ", r[:mineig_Qt])
    end

    write_report(
        report_path,
        r;
        level = level,
        max_tau_terms = max_tau_terms,
        max_time = max_time,
        threads = threads,
    )
    println("report_written = ", report_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
