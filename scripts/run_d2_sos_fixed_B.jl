#!/usr/bin/env julia

using JuMP
using SumOfSquares
using DynamicPolynomials
using MosekTools

# d=2 setup:
# Fhat(u,v) >= 0 on [0,1]^2 and Ghat(u,v)=0
# where u = λ/(1+λ), v = α/(1+α).

function build_cleared_polys(B::Float64, u, v)
    omu = 1 - u
    omv = 1 - v

    # Cleared numerators after substitution λ=u/(1-u), α=v/(1-v).
    n2_Bλ = B * omu^2 + 2B * u * omu + B^3 * u^2
    n3_λ = B^3 * omu^3 + 3B * u * omu^2 + 3B * u^2 * omu + B^3 * u^3

    n4_α = B^6 * omv^4 +
           4B^3 * v * omv^3 +
           6B^2 * v^2 * omv^2 +
           4B^3 * v^3 * omv +
           B^6 * v^4

    n3_Bα = B^3 * omv^3 +
            3B^2 * v * omv^2 +
            3B^3 * v^2 * omv +
            B^6 * v^3

    n2_λ_div_B = B * omu^2 + (2 / B) * u * omu + (1 / B) * u^2
    n3_α_div_B = B^3 * omv^3 +
                 3 * v * omv^2 +
                 (3 / B) * v^2 * omv +
                 v^3

    fhat = n2_Bλ^6 * n4_α^3 - n3_Bα^4 * n3_λ^4
    ghat = n2_λ_div_B^3 * n3_Bα^2 - n2_Bλ^3 * n3_α_div_B^2
    return fhat, ghat
end

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

function run_fixed_B_certificate(
    B::Float64;
    d0::Int = 8,
    d1::Int = 7,
    dt::Int = 8,
    silent::Bool = true,
)
    @polyvar u v
    fhat, ghat = build_cleared_polys(B, u, v)

    gu = u * (1 - u)
    gv = v * (1 - v)

    model = Model(Mosek.Optimizer)
    if silent
        set_silent(model)
    end

    mons0 = monomials([u, v], 0:d0)
    mons1 = monomials([u, v], 0:d1)
    monst = monomials([u, v], 0:dt)

    @variable(model, σ0, SOSPoly(mons0))
    @variable(model, σu, SOSPoly(mons1))
    @variable(model, σv, SOSPoly(mons1))
    @variable(model, τ, Poly(monst))
    @variable(model, γ)

    @constraint(model, fhat - γ == σ0 + σu * gu + σv * gv + τ * ghat)
    @objective(model, Max, γ)
    optimize!(model)

    status = termination_status(model)
    pstatus = primal_status(model)
    dual = dual_status(model)

    println("status = ", status)
    println("primal_status = ", pstatus)
    println("dual_status = ", dual)

    gamma_val = NaN
    result_label = "NO_CERTIFICATE"
    q0 = nothing
    qu = nothing
    qv = nothing

    if status in (OPTIMAL, LOCALLY_SOLVED, ALMOST_OPTIMAL)
        gv = value(γ)
        gamma_val = gv
        println("gamma_lower_bound = ", gv)
        if gv >= -1e-7
            result_label = "NONNEGATIVE (within tolerance)"
            println("certificate_result = ", result_label)
        else
            result_label = "NEGATIVE LOWER BOUND (increase degrees)"
            println("certificate_result = ", result_label)
        end
        q0 = value.(Matrix(σ0.Q))
        qu = value.(Matrix(σu.Q))
        qv = value.(Matrix(σv.Q))
    else
        println("No optimal solution. Try larger degree or inspect Mosek logs.")
    end

    return (
        B = B,
        d0 = d0,
        d1 = d1,
        dt = dt,
        status = status,
        primal_status = pstatus,
        dual_status = dual,
        gamma = gamma_val,
        certificate_result = result_label,
        basis0 = string.(σ0.basis.monomials),
        basis1 = string.(σu.basis.monomials),
        Q0 = q0,
        Qu = qu,
        Qv = qv,
    )
end

function main()
    B = parse_float_arg("B", 0.4)
    d0 = parse_int_arg("d0", 8)
    d1 = parse_int_arg("d1", 7)
    dt = parse_int_arg("dt", 8)
    silent = !has_flag("--verbose")

    if !(0.0 < B < 1.0)
        error("B must satisfy 0 < B < 1.")
    end

    println("Running fixed-B SOS certificate")
    println("B = ", B, ", d0 = ", d0, ", d1 = ", d1, ", dt = ", dt)
    run_fixed_B_certificate(B; d0 = d0, d1 = d1, dt = dt, silent = silent)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
