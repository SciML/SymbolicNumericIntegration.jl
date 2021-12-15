using LinearAlgebra
using Statistics: mean, std

Base.signbit(z::Complex{T}) where T<:Number = signbit(real(z))

"""
    integrate is the main entry point

    input:
    ------
    eq: a Symbolics expression to integrate
    x: the independent variable (optional)

    abstol: the desired tolerance
    num_steps: the number of different steps with expanding basis to be tried
    num_trials: the number of trials in each step (no changes to the basis)
    radius: the radius of the disk in the complex plane to generate random test points
    show_basis: if true, the basis (list of candidate terms) is printed
    opt: the sparse regression optimizer
    bypass: if true, do not integrate terms separately but consider all at once
    symbolic: try symbolic integration first
    max_basis: the maximum number of candidate terms to consider
    verbose: print a detailed report
    complex_plane: generate random test points on the complex plane (if false, the points will be on real axis)
    homotomy: use the homotopy algorithm to generat the basis

    output:
    -------
    solved, unsolved, err

    solved is the solved integral and unsolved is the residual unsolved portion of the input
    err is the numerical error in reaching the solution
"""
function integrate(eq, x=nothing; abstol=1e-6, num_steps=2, num_trials=5, radius=1.0,
                   show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false,
                   symbolic=true, max_basis=100, verbose=false, complex_plane=true,
                   homotopy=true)
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

    # homotopy = homotopy && !bypass

    return integrate_sum(eq, x, l; bypass, abstol, num_trials, num_steps,
                              radius, show_basis, opt, symbolic,
                              max_basis, verbose, complex_plane, homotopy)
end

"""
    integrate_sum applies the integral summation rule ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx

    eq: the integrand
    x: the indepedent variable
    l: a logger

    output is the same as integrate
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
        eq = factor_rational(simplify_trigs(unsolved))

        if !isequal(eq, unsolved)
            eq = expand(eq)
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

                if !isequal(u, 0)   # premature termination on the first failure
                    return 0, eq, ϵ₀
                end

            end
        end
    end

    return expand(solved), unsolved, ϵ₀
end

"""
    integrate_term is the central part of the code that tries different
    methods to integrate eq

    eq: the integrand
    x: the indepedent variable
    l: a logger

    output is the same as integrate
"""
function integrate_term(eq, x, l; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis,
        radius, homotopy = args[:abstol], args[:num_steps],
        args[:num_trials], args[:show_basis], args[:symbolic], args[:verbose],
        args[:max_basis], args[:radius], args[:homotopy]

    attempt(l, "Integrating term", eq)

    if is_number(eq)
        y = eq * x
        result(l, "Successful", y)
        return y, 0, 0
    end

    basis = generate_basis(eq, x; homotopy)

    if show_basis
        inform(l, "Generating basis (|β| = $(length(basis)))", basis)
    end

    if length(basis) > max_basis
        result(l, "|β| = $(length(basis)) is too large")
        return 0, eq, Inf
    end

    D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    # rescue
    ϵᵣ = Inf
    yᵣ = 0

    for i = 1:num_steps
        if length(basis) > max_basis break end

        Δbasis = [expand_derivatives(D(f)) for f in basis]

        if symbolic
            y, ϵ = try_symbolic(Float64, eq, x, basis, Δbasis; kwargs...)

            if !isequal(y, 0) && accept_solution(eq, x, y, radius) < abstol
                result(l, "Successful symbolic", y)
                return y, 0, 0
            else
                inform(l, "Failed symbolic")
            end
        end

        for j = 1:num_trials
            r = radius #*sqrt(2)^j
            y, ϵ = try_integrate(Float64, eq, x, basis, Δbasis, r; kwargs...)

            ϵ = accept_solution(eq, x, y, r)
            if ϵ < abstol
                result(l, "Successful numeric (attempt $j out of $num_trials)", y)
                return y, 0, ϵ
            elseif ϵ < ϵᵣ
                ϵᵣ = ϵ
                yᵣ = y
            end
        end

        inform(l, "Failed numeric")

        if i < num_steps
            basis = expand_basis(basis, x)
            if show_basis
                inform(l, "Expanding the basis (|β| = $(length(basis)))", basis)
            elseif verbose
                inform(l, "Expanding the basis (|β| = $(length(basis)))")
            end
        end
    end

    if ϵᵣ < abstol * 10
        result(l, "Accepting numeric (rescued)", yᵣ)
        return yᵣ, 0, ϵᵣ
    else
        result(l, "Unsucessful", eq)
        return 0, eq, ϵ₀
    end
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

    if n < 8    # here, 8 is arbitrary and signifies a small basis
        y₃, ϵ₃ = find_dense(T, A, basis; abstol)
        if ϵ₃ < abstol
            return y₃, ϵ₃
        end
    end

    # moving toward the poles
    ∂eq = expand_derivatives(Differential(x)(eq))
    modify_basis_matrix!(T, A, X, x, eq, ∂eq, Δbasis, radius; abstol)
    y₄, ϵ₄ = sparse_fit(T, A, x, basis, Δbasis, opt; abstol)

    if ϵ₄ < abstol || ϵ₄ < ϵ₁
        return y₄, ϵ₄
    else
        return y₁, ϵ₁
    end
end

function init_basis_matrix!(T, A, X, x, eq, Δbasis, radius, complex_plane; abstol=1e-6)
    n = size(A, 1)
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
            println("Error from init_basis_matrix!: ", e)
        end
    end
end

function modify_basis_matrix!(T, A, X, x, eq, ∂eq, Δbasis, radius; abstol=1e-6)
    n = size(A, 1)
    k = 1
    for k = 1:n
        d = Dict(x => X[k])
        # One Newton iteration toward the poles
        # note the + sign instead of the usual - in Newton-Raphson's method. This is
        # because we are moving toward the poles and not zeros.
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
        if sum(iscomplex.(q)) > 2 return nothing, Inf end   # eliminating complex coefficients
        return sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x)), abs(ϵ)
    catch e
        println("Error from sparse_fit", e)
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
