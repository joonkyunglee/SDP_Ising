#!/usr/bin/env julia

using Dates
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "run_d2_sos_fixed_B.jl"))

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

function fmt_sci(x::Real)
    @sprintf("%.6e", x)
end

function matrix_block_text(M::AbstractMatrix{<:Real}; maxn::Int = 12)
    n = size(M, 1)
    m = min(n, maxn)
    io = IOBuffer()
    for i in 1:m
        for j in 1:m
            @printf(io, "% .4e", M[i, j])
            if j < m
                print(io, "  ")
            end
        end
        print(io, "\n")
    end
    if m < n
        print(io, "... (truncated to top-left $(m)x$(m) block out of $(n)x$(n))\n")
    end
    return String(take!(io))
end

function basis_preview(basis::Vector{String}; max_terms::Int = 20)
    n = length(basis)
    if n <= max_terms
        return join(basis, ", ")
    end
    return join(vcat(basis[1:max_terms], ["..."]), ", ")
end

function summarize_matrix(M::AbstractMatrix{<:Real})
    eigmin = minimum(eigvals(Symmetric(M)))
    (
        n = size(M, 1),
        min_eig = eigmin,
        max_abs = maximum(abs, M),
        fro = sqrt(sum(abs2, M)),
    )
end

function append_summary_table(io::IO, results)
    println(io, "| B | status | gamma_lower_bound | certificate_result |")
    println(io, "|---:|:-------|------------------:|:-------------------|")
    for r in results
        gamma_text = isnan(r.gamma) ? "NaN" : fmt_sci(r.gamma)
        println(
            io,
            "| ",
            r.B,
            " | ",
            r.status,
            " | ",
            gamma_text,
            " | ",
            r.certificate_result,
            " |",
        )
    end
end

function append_detail_section(io::IO, r)
    println(io, "## Detailed Matrices for B = ", r.B)
    println(io)
    println(io, "- status: `", r.status, "`")
    println(io, "- primal_status: `", r.primal_status, "`")
    println(io, "- dual_status: `", r.dual_status, "`")
    println(io, "- gamma_lower_bound: `", fmt_sci(r.gamma), "`")
    println(io, "- certificate_result: `", r.certificate_result, "`")
    println(io)

    for (name, basis, M) in (("Q0 (sigma0)", r.basis0, r.Q0), ("Qu (sigmau)", r.basis1, r.Qu), ("Qv (sigmav)", r.basis1, r.Qv))
        s = summarize_matrix(M)
        println(io, "### ", name)
        println(io)
        println(io, "- size: `", s.n, "x", s.n, "`")
        println(io, "- min_eigenvalue: `", fmt_sci(s.min_eig), "`")
        println(io, "- max_abs_entry: `", fmt_sci(s.max_abs), "`")
        println(io, "- frobenius_norm: `", fmt_sci(s.fro), "`")
        println(io, "- basis_preview: `", basis_preview(basis), "`")
        println(io)
        println(io, "Top-left block:")
        println(io)
        println(io, "```text")
        print(io, matrix_block_text(M; maxn = 12))
        println(io, "```")
        println(io)
    end
end

function main()
    outfile = parse_string_arg("output", joinpath(dirname(@__DIR__), "d2_sos_matrices_report.md"))
    detailB = parse_float_arg("detailB", 0.4)
    d0 = parse_int_arg("d0", 12)
    d1 = parse_int_arg("d1", 11)
    dt = parse_int_arg("dt", 12)

    Bs = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    results = []
    for B in Bs
        println("Running B = ", B)
        r = run_fixed_B_certificate(B; d0 = d0, d1 = d1, dt = dt, silent = true)
        push!(results, r)
    end

    detail = nothing
    for r in results
        if isapprox(r.B, detailB; atol = 1e-12)
            detail = r
            break
        end
    end
    if isnothing(detail)
        error("detailB=$(detailB) not found in sweep values.")
    end

    open(outfile, "w") do io
        println(io, "# d=2 SOS Matrix Report")
        println(io)
        println(io, "- Generated at: `", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "`")
        println(io, "- Solver: `Mosek`")
        println(io, "- Relaxation degrees: `d0=", d0, ", d1=", d1, ", dt=", dt, "`")
        println(io, "- Sweep B values: `", join(Bs, ", "), "`")
        println(io)
        println(io, "## Summary")
        println(io)
        append_summary_table(io, results)
        println(io)
        append_detail_section(io, detail)
    end

    println("Report written to ", outfile)
end

main()
