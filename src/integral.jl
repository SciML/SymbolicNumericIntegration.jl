using LinearAlgebra
using Statistics: mean, std

Base.signbit(z::Complex{T}) where T<:Number = signbit(real(z))

"""
    integrate is the main entry point

    input:
    ------
    eq: a Symbolics expression to integrate
    abstol: the desired tolerance
    num_steps: the number of different steps with expanding basis to be tried
    num_trials: the number of trials in each step (no changes to the basis)
    lo and hi: the range used to generate random values of x (the independent variable)
    show_basis: if true, the basis is printed

    output:
    -------
    solved, unsolved

    a pair of expressions, solved is the solved integral and unsolved is the residual unsolved
    portion of the input
"""
function integrate(eq, x=nothing; abstol=1e-6, num_steps=2, num_trials=5, radius=1.0,
                   show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false,
                   symbolic=true, bypart=true, max_basis=110,
                   verbose=false, complex_plane=true, prune_basis=false)
    eq = expand(eq)
    eq = apply_div_rule(eq)

    if x == nothing
        x = var(eq)
        if x == nothing
            @syms 𝑥
            x = 𝑥
        end
    end

    l = Logger(verbose)

    # eq is a constant
    if !isdependent(eq, x)
        return x * eq, 0, 0
    end

    # check if eq is a rational function
    # if so, we perform a partial-fraction decomposition first (the first part of the Hermite's method)

    # q = to_rational(eq, x)
    # if q != nothing
    #     eq = q
    # end

    s₁, u₁, ϵ = integrate_sum(eq, x, l; bypass, abstol, num_trials, num_steps,
                              radius, show_basis, opt, symbolic,
                              max_basis, verbose, complex_plane, prune_basis)

    if isequal(u₁, 0) || !bypart
        return s₁, u₁, ϵ
    else
        s₂, u₂, ϵ = try_integration_by_parts(u₁, x, l; abstol, num_trials, num_steps,
                                             radius, show_basis, opt, symbolic,
                                             max_basis, verbose, complex_plane,
                                             prune_basis)
        return s₁ + s₂, u₂, ϵ
    end
end

"""
    ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
"""
function integrate_sum(eq, x, l; bypass=false, kwargs...)
    solved = 0
    unsolved = 0
    ϵ₀ = 0
    ts = bypass ? [eq] : terms(eq)

    if length(ts) > 1
        inform(l, "Integrating sum", ts)
    end

    for p in ts
        s, u, ϵ = integrate_term(p, x, l; kwargs...)
        solved += s
        unsolved += u
        ϵ₀ = max(ϵ₀, ϵ)
    end

    if !isequal(unsolved, 0)
        eq = apply_q_rules(apply_integration_rules(unsolved))

        if !isequal(eq, unsolved)
            # eq = expand(eq)
            unsolved = 0
            ϵ₀ = 0
            ts = bypass ? [eq] : terms(eq)

            if length(ts) > 1
                inform(l, "Integrating transformed sum", ts)
            else
                inform(l, "Transforming the expression", ts[1])
            end

            for p in ts
                s, u, ϵ = integrate_term(p, x, l; kwargs...)
                solved += s
                unsolved += u
                ϵ₀ = max(ϵ₀, ϵ)
            end
        end
    end

    return solved, unsolved, ϵ₀
end

function integrate_sum_fast(eq, x, l; kwargs...)
    solved = 0
    ϵ₀ = 0

    for p in terms(eq)
        s, u, ϵ = integrate_term(p, x, l; kwargs...)
        if !isequal(u, 0) return 0, eq, ϵ end
        solved += s
        ϵ₀ = max(ϵ₀, ϵ)
    end
    return solved, 0, ϵ₀
end

function test_point(complex_plane, radius)
    if complex_plane
        return radius * sqrt(rand()) * cis(2π*rand())
    else
        return Complex(radius * (2*rand() - 1))
    end
end

function accept_solution(eq, x, sol, radius; abstol=1e-6)
    try
        x₀ = test_point(true, radius)
        Δ = substitute(expand_derivatives(Differential(x)(sol)-eq), Dict(x => x₀))
        return abs(Δ) < abstol
    catch e
        #
    end
    return false
end

function integrate_term(eq, x, l; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis,
        radius, prune_basis = args[:abstol], args[:num_steps], args[:num_trials],
        args[:show_basis], args[:symbolic], args[:verbose], args[:max_basis],
        args[:radius], args[:prune_basis]

    attempt(l, "Integrating term", eq)

    if is_number(eq)
        y = eq * x
        result(l, "Successful", y)
        return y, 0, 0
    end

    # note that the order of the operations is important!
    # first, collecing hints, then applying transformation rules, and finally finding the basis.
    # basis = generate_basis(eq, x)
    basis = generate_by_parts(eq, x)

    if show_basis
        inform(l, "Generating basis (|β| = $(length(basis)))", basis)
    end

    if prune_basis || length(basis) > max_basis
        basis, ok = prune(basis, eq, x)
        if ok && show_basis
            inform(l, "Prunning the basis (|β| = $(length(basis)))", basis)
        end
    end

    if length(basis) > max_basis
        result(l, "|β| = $(length(basis)) is too large")
        return 0, eq, Inf
    end

    D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    for i = 1:num_steps
        Δbasis = [expand_derivatives(D(f)) for f in basis]

        # if show_basis println(basis) end

        if symbolic
            y, ϵ = try_symbolic(Float64, eq, x, basis, Δbasis; kwargs...)

            if !isequal(y, 0) && accept_solution(eq, x, y, radius; abstol)
                # if verbose printstyled("$i, symbolic\n"; color=:yellow) end
                result(l, "Successful symbolic", y)
                return y, 0, 0
            else
                inform(l, "Failed symbolic")
            end
        end

        for j = 1:num_trials
            r = radius*sqrt(2)^j
            y, ϵ = try_integrate(Float64, eq, x, basis, Δbasis, r; kwargs...)

            if ϵ < abstol && accept_solution(eq, x, y, r; abstol)
                # if verbose printstyled("$i, $j\n"; color=:yellow) end
                result(l, "Successful numeric (attempt $j out of $num_trials)", y)
                return y, 0, ϵ
            else
                ϵ₀ = min(ϵ, ϵ₀)
                y₀ = y
            end
        end

        inform(l, "Failed numeric")

        if i < num_steps
            basis = unique([basis; basis*x])
            if show_basis
                inform(l, "Expanding the basis (|β| = $(length(basis)))", basis)
            end
        end
    end

    if accept_solution(eq, x, y₀, radius; abstol=abstol*10)
        # if verbose printstyled("rescue\n"; color=:yellow) end
        result(l, "Accepting numeric (rescued)", y₀)
        return y₀, 0, ϵ₀
    else
        result(l, "Unsucessful", eq)
        return 0, eq, ϵ₀
    end
end

rms(x) = sqrt(sum(x.^2) / length(x))

"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol=1e-3)
    Q, R = qr(A)
    abs.(diag(R)) .> abstol
end

"""
    converts float to int or small rational numbers
"""
function nice_parameters(p; abstol=1e-3)
    c = lcm(collect(1:10)...)
    n = length(p)
    q = Array{Any}(undef, n)
    for i = 1:n
        den = 1
        while den < 10
            if abs(round(p[i]*den) - p[i]*den) < abstol
                a = round(Int, p[i]*den) // den
                q[i] = (denominator(a) == 1 ? numerator(a) : a)
                den = 10
            else
                q[i] = Float64(p[i])
            end
            den += 1
        end
    end
    q
end

function nice_parameter(u::T; abstol=1e-3, M=10) where T<:Real
    c = lcm(collect(1:M)...)
    for den = 1:M
        try
            if abs(round(u*den) - u*den) < abstol
                a = round(Int, u*den) // den
                return (denominator(a) == 1 ? numerator(a) : a)
            end
        catch e
        end
    end
    return u
end

function nice_parameter(u::Complex{T}; abstol=1e-3, M=10) where T<:Real
    α = nice_parameter(real(u))
    β = nice_parameter(imag(u))
    return β ≈ 0 ? α : Complex(α, β)
end

###############################################################################

"""
    the core of the randomized parameter-fitting algorithm

    `try_integrate` tries to find a linear combination of the basis, whose
    derivative is equal to eq

    output
    -------
    integral, error
"""
function try_integrate(T, eq, x, basis, Δbasis, radius; kwargs...)
    args = Dict(kwargs)
    abstol, opt, complex_plane, verbose =
        args[:abstol], args[:opt], args[:complex_plane], args[:verbose]

    basis = basis[2:end]    # remove 1 from the beginning
    Δbasis = Δbasis[2:end]
    n = length(basis)

    # A is an nxn matrix holding the values of the fragments at n random points
    A = zeros(Complex{T}, (n, n))
    X = zeros(Complex{T}, n)

    init_basis_matrix!(T, A, X, x, eq, Δbasis, radius, complex_plane; abstol)

    y₁, ϵ₁ = sparse_fit(T, A, x, basis, Δbasis, opt; abstol)
    if ϵ₁ < abstol
        return y₁, ϵ₁
    end

    y₂, ϵ₂ = find_singlet(T, A, basis; abstol)
    if ϵ₂ < abstol
        return y₂, ϵ₂
    end

    if n < 8    # 8 is arbitrary here and signifies a small basis
        y₃, ϵ₃ = find_dense(T, A, basis; abstol)
        if ϵ₃ < abstol
            return y₃, ϵ₃
        end
    end

    ∂eq = expand_derivatives(Differential(x)(eq))
    modify_basis_matrix!(T, A, X, x, eq, ∂eq, Δbasis, radius; abstol)
    y₄, ϵ₄ = sparse_fit(T, A, x, basis, Δbasis, opt; abstol)

    if ϵ₄ < abstol || ϵ₄ < ϵ₁
        # if verbose printstyled("improvement after moving toward poles\n"; color=:blue) end
        return y₄, ϵ₄
    else
        return y₁, ϵ₁
    end
end

function init_basis_matrix!(T, A, X, x, eq, Δbasis, radius, complex_plane; abstol=1e-6)
    n = size(A, 1)
    # X = zeros(Complex{T}, n)
    k = 1
    i = 1

    while k <= n
        try
            x₀ = test_point(complex_plane, radius)
            X[k] = x₀
            d = Dict(x => x₀)

            b₀ = Complex{T}(substitute(eq, d))
            if is_proper(b₀)
                for j = 1:n
                    A[k,j] = Complex{T}(substitute(Δbasis[j], d)) / b₀
                end
                if all(is_proper, A[k,:])
                    k += 1
                end
            end
        catch e
            println(e)
        end
    end
end

function modify_basis_matrix!(T, A, X, x, eq, ∂eq, Δbasis, radius; abstol=1e-6)
    n = size(A, 1)
    k = 1
    for k = 1:n
        d = Dict(x => X[k])
        # One Newton iteration toward the poles
        x₀ = X[k] + Complex{T}(substitute(eq, d)) / Complex{T}(substitute(∂eq, d))
        X[k] = x₀
        d = Dict(x => x₀)
        b₀ = Complex{T}(substitute(eq, d))
        for j = 1:n
            A[k,j] = Complex{T}(substitute(Δbasis[j], d)) / b₀
        end
    end
end

function sparse_fit(T, A, x, basis, Δbasis, opt; abstol=1e-6)
    n = length(basis)
    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, basis, Δbasis, n = A[l,l], basis[l], Δbasis[l], sum(l)

    try
        b = ones(n)
        # q₀ = A \ b
        q₀ = Optimize.init(opt, A, b)
        @views Optimize.sparse_regression!(q₀, A, permutedims(b)', opt, maxiter = 1000)
        ϵ = rms(A * q₀ - b)
        q = nice_parameter.(q₀)
        return sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x)), abs(ϵ)
    catch e
        println(e)
        return nothing, Inf
    end
end

function find_singlet(T, A, basis; abstol)
    σ = vec(std(A; dims=1))
    μ = vec(mean(A; dims=1))
    l = (σ .< abstol) .* (abs.(μ) .> abstol)
    if sum(l) == 1
        k = findfirst(l)
        return nice_parameter(1/μ[k]) * basis[k], σ[k]
    else
        return nothing, Inf
    end
end

function find_dense(T, A, basis; abstol=1e-6)
    n = size(A, 1)
    b = ones(T, n)

    try
        q = A \ b
        if minimum(abs.(q)) > abstol
            ϵ = maximum(abs.(A*q .- b))
            if ϵ < abstol
                y = sum(nice_parameter.(q) .* basis)
                return y, ϵ
            end
        end
    catch e
        #
    end
    return nothing, Inf
end
