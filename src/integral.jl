using LinearAlgebra
using Statistics: mean, std

Base.signbit(z::Complex{T}) where T<:Number = signbit(real(z))

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x, h=[])
    Δeq = expand_derivatives(Differential(x)(eq))
    kers = expand(eq + Δeq)
    return [one(x); candidates(kers, x); h]
end

"""
    candidates returns a list of candidate expressions to form the integration
    basis
"""
candidates(eq, x) = isdependent(eq,x) ? [eq] : []
candidates(eq::Num, x) = candidates(value(eq), x)

# the candidates of an Add is the union of the candidates of the terms
# ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
candidates(eq::SymbolicUtils.Add, x) = unique(∪([candidates(t,x) for t in arguments(eq)]...))

# the candidates of a Mul is the outer product of the candidates of the terms
# d(uv)/dx = u dv/dx + v + du/dx
function candidates(eq::SymbolicUtils.Mul, x)
    terms = [candidates(q,x) for q in arguments(eq)]
    n = length(terms)

    l = Any[one(x)]

    for j = 1:n
        m = length(l)
        for t in terms[j]
            for k = 1:m
                push!(l, l[k]*t)
            end
        end
    end

    unique(l[2:end])    # removing the initial 1
end

# the candidates of a Pow encode different integration rules
function candidates(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) return [one(x)] end

    p = arguments(eq)[1]    # eq = p ^ k
    k = arguments(eq)[2]

    if k < 0 && k ≈ round(k)
        return candidate_pow_minus(p, k, x)
    elseif k ≈ 0.5 || k ≈ -0.5
        if check_poly(p,x) == :real_poly && leading(p,x) < 0 p = -p end
        # ∫ √f(x) dx = ... + c * log(df/dx + √f) if deg(f) == 2
        Δ = expand_derivatives(Differential(x)(p))
        return [[p^k, p^(k+1)]; log(0.5*Δ + sqrt(p))]
    end

    # ∫ p^k dp = c * p^(k+1)
    return [p^k, p^(k+1)]
end

nice_abs2(u) = abs2(u)     # nice_parameter(abs2(u))

function candidate_pow_minus(p, k, x)
    check = check_poly(p, x)
    if check == :not_poly || check == :complex_poly
        return [p^k, p^(k+1), log(p)]
        # Δp = expand_derivatives(Differential(x)(p))
        # return [p^k, p^(k+1), log(p); Δp]
    end

    r, s = find_roots(p, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    # ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
    q = [[log(x - u) for u in r];
          [atan((x - real(u))/imag(u)) for u in s];
          [log(x^2 - 2*real(u)*x + nice_abs2(u)) for u in s]
         ]

    # return [[p^k, p^(k+1)]; candidates(q₁, x)]
    if k ≈ -1
        return [[p^k]; q]
    else
        return [[p^k, p^(k+1)]; q]
    end
end

###############################################################################

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
function integrate(eq, x=nothing; abstol=1e-6, num_steps=3, num_trials=2, radius=1.0,
                   show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false,
                   attempt_ratio=5, symbolic=true, bypart=true, max_basis=110,
                   verbose=false, margin=1.0, complex_plane=true)
    eq = expand(eq)

    if x == nothing
        x = var(eq)
        if x == nothing
            @syms 𝑥
            x = 𝑥
        end
    end

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

    s₁, u₁, ϵ = integrate_sum(eq, x; bypass, abstol, num_trials, num_steps,
                              radius, show_basis, opt, attempt_ratio, symbolic,
                              max_basis, verbose, margin, complex_plane)

    if isequal(u₁, 0) || !bypart
        return s₁, u₁, ϵ
    else
        s₂, u₂, ϵ = try_integration_by_parts(u₁, x; abstol, num_trials, num_steps,
                                             radius, show_basis, opt, attempt_ratio,
                                             symbolic, max_basis, verbose, margin,
                                             complex_plane)
        return s₁ + s₂, u₂, ϵ
    end
end

"""
    ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
"""
function integrate_sum(eq::SymbolicUtils.Add, x; bypass=false, kwargs...)
    if bypass
        return integrate_term(eq, x; kwargs...)
    else
        solved = 0
        unsolved = 0
        ϵ₀ = 0

        for p in arguments(eq)
            s, u, ϵ = integrate_term(p, x; kwargs...)
            solved += s
            unsolved += u
            ϵ₀ = max(ϵ₀, ϵ)
        end

        return solved, unsolved, ϵ₀
    end
end

function integrate_sum(eq, x; kwargs...)
    integrate_term(eq, x; kwargs...)
end

function accept_solution(eq, x, sol, radius, margin; abstol=1e-6)
    try
        Δ = substitute(expand_derivatives(Differential(x)(sol)-eq), Dict(x => test_point(radius, margin)))
        return abs(Δ) < abstol
    catch e
        #
    end
    return false
end

function integrate_term(eq, x; kwargs...)
    args = Dict(kwargs)
    abstol = args[:abstol]
    h = collect_hints(eq, x)
    s, u, ϵ = integrate_term_main(eq, x, h; kwargs...)
    if isequal(u, 0)
        return s, u, ϵ
    else
        eq = apply_integration_rules(eq)
        return integrate_term_main(eq, x, h; kwargs...)
    end
end

function integrate_term_main(eq, x, h; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis, radius, margin =
        args[:abstol], args[:num_steps], args[:num_trials], args[:show_basis],
        args[:symbolic], args[:verbose], args[:max_basis], args[:radius], args[:margin]

    # note that the order of the operations is important!
    # first, collecing hints, then applying transformation rules, and finally finding the basis.
    basis = generate_basis(eq, x, h)

    # basis = filter(u -> !(deg(u,x)>0), basis)

    if verbose printstyled("|β| = ", length(basis), '\n'; color=:yellow) end
    if length(basis) > max_basis return 0, eq, Inf end

    D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    for i = 1:num_steps
        basis = unique([basis; basis*x])
        Δbasis = [expand_derivatives(D(f)) for f in basis]

        if show_basis println(basis) end

        if symbolic
            y, ϵ = try_symbolic(Float64, eq, x, basis, Δbasis; kwargs...)

            if !isequal(y, 0) && accept_solution(eq, x, y, radius, margin; abstol)
                if verbose printstyled("$i, symbolic\n"; color=:yellow) end
                return y, 0, 0
            end
        end

        for j = 1:num_trials
            r = radius*sqrt(2)^j
            y, ϵ = try_integrate2(Float64, eq, x, basis, Δbasis, r, margin; kwargs...)

            if ϵ < abstol && accept_solution(eq, x, y, r, margin; abstol)
                if verbose printstyled("$i, $j\n"; color=:yellow) end
                return y, 0, ϵ
            else
                ϵ₀ = min(ϵ, ϵ₀)
                y₀ = y
            end
        end
    end

    if accept_solution(eq, x, y₀, radius, margin; abstol=abstol*10)
        if verbose printstyled("rescue\n"; color=:yellow) end
        return y₀, 0, ϵ₀
    else
        return 0, eq, ϵ₀
    end
end

rms(x) = sqrt(sum(x.^2) / length(x))

function test_point(radius, margin)
    # select a quadrant
    s1 = rand([-1,1])
    s2 = rand([-1,1])

    x = s1*(margin + rand()*radius)
    y = s2*(margin + rand()*radius)

    return Complex(x, y)
end

"""
    the core of the randomized parameter-fitting algorithm

    `try_integrate` tries to find a linear combination of the basis, whose
    derivative is equal to eq

    output
    -------
    integral, error
"""
function try_integrate(T, eq, x, basis, Δbasis, radius, margin=1.0; kwargs...)
    args = Dict(kwargs)
    abstol, opt, attempt_ratio = args[:abstol], args[:opt], args[:attempt_ratio]

    n = length(basis)
    # A is an nxn matrix holding the values of the fragments at n random points
    A = zeros(Complex{T}, (n, n))

    i = 1
    k = 1

    while i <= n
        x₀ = test_point(radius, margin)
        d = Dict(x => x₀)
        try
            b₀ = Complex{T}(substitute(eq, d))
            for j = 1:n
                A[i, j] = Complex{T}(substitute(Δbasis[j], d)) / b₀
            end
            i += 1
        catch e
            println("basis matrix error: ", e)
        end
        if k > attempt_ratio*n return nothing, 1e6 end
        k += 1
    end

    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, basis, Δbasis, n = A[l,l], basis[l], Δbasis[l], sum(l)

    if det(A) ≈ 0 return nothing, 1e6 end

    coefs = ones(Complex{T}, n)
    for j = 1:n
        coefs[j] = coef(Δbasis[j], x)
        A[:,j] /= coefs[j]
    end

    # q₀ = A \ b
    try
        b = ones(n)
        q₀ = Optimize.init(opt, A, b)
        @views Optimize.sparse_regression!(q₀, A, permutedims(b)', opt, maxiter = 1000)
        ϵ = rms(A * q₀ - b)
        q = nice_parameter.(q₀ ./ coefs)
        return sum(q[i]*basis[i] for i = 1:length(basis) if q[i] != 0; init=zero(x)), abs(ϵ)
    catch e
        return nothing, 1e6
    end
end

"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol=1e-3)
    Q, R = qr(A)
    abs.(diag(R)) .> abstol
end

function find_independent_subset2(A; abstol=1e-3)
    n = size(A, 1)
    l = BitVector(undef, n)
    for i = 1:n
        l[i] = 1
        if abs(det(A[l,l])) < abstol
            l[i] = 0
        end
    end
    l
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

function fibonacci_spiral(n, radius=1.0)
    z = zeros(Complex{Float64}, n)
    ϕ = (1 + sqrt(5)) / 2
    for i = 1:n
        θ, r = mod((i-1)/ϕ, 1), (i-1)/n
        z[i] = radius*sqrt(r)*cis(2π*θ)
    end
    return z
end

is_proper(x) = !isnan(x) && !isinf(x)

function try_integrate2(T, eq, x, basis, Δbasis, radius, margin=1.0; kwargs...)
    args = Dict(kwargs)
    abstol, opt, attempt_ratio, complex_plane, verbose =
        args[:abstol], args[:opt], args[:attempt_ratio], args[:complex_plane], args[:verbose]

    basis = basis[2:end]    # remove 1 from the beginning
    Δbasis = Δbasis[2:end]
    n = length(basis)

    # A is an nxn matrix holding the values of the fragments at n random points
    A = zeros(Complex{T}, (n, n))
    X = zeros(Complex{T}, n)

    init_basis_matrix!(T, A, X, x, eq, Δbasis, radius, complex_plane; abstol)

    y₀, ϵ₀ = find_singlet(T, A, basis; abstol)
    if ϵ₀ < abstol
        return y₀, ϵ₀
    end

    y₁, ϵ₁ = sparse_fit(T, A, x, basis, Δbasis, opt; abstol)
    if ϵ₁ < abstol
        return y₁, ϵ₁
    end

    ∂eq = expand_derivatives(Differential(x)(eq))
    modify_basis_matrix!(T, A, X, x, eq, ∂eq, Δbasis, radius; abstol)
    y₂, ϵ₂ = sparse_fit(T, A, x, basis, Δbasis, opt; abstol)
    if ϵ₂ < abstol || ϵ₂ < ϵ₁
        if verbose printstyled("improvement after moving toward poles\n"; color=:blue) end
        return y₂, ϵ₂
    else
        return y₁, ϵ₁
    end
end

function init_basis_matrix!(T, A, X, x, eq, Δbasis, radius, complex_plane; abstol=1e-6)
    n = size(A, 1)
    X = zeros(Complex{T}, n)
    k = 1
    while k <= n
        try
            if complex_plane
                x₀ = radius * sqrt(rand()) * cis(2π*rand())
            else
                x₀ = Complex(radius * (2*rand() - 1))
            end

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
            # println(e)
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

########################## Convert to a Polynomial? #############################

# is_multiple_x(p::SymbolicUtils.Add, x) = all(is_multiple_x(t,x) for t in arguments(p))
# is_multiple_x(p::SymbolicUtils.Mul, x) = any(is_multiple_x(t,x) for t in arguments(p))
# is_multiple_x(p::SymbolicUtils.Pow, x) = isequal(arguments(p)[1], x) && arguments(p)[2] >= 1
# is_multiple_x(p::SymbolicUtils.Sym, x) = isequal(p, x)
# is_multiple_x(p, x) = false


coef(eq::SymbolicUtils.Mul, x) = prod(t for t in arguments(eq) if !isdependent(t,x); init=1)
coef(eq::SymbolicUtils.Add, x) = minimum(abs(coef(t,x)) for t in arguments(eq))
coef(eq, x) = 1

##################### Special Functions ######################################

# using SpecialFunctions
#
# """
#     logarithmic integral
# """
# function li(x; n=10)
#     z = log(abs(x))
#     s = sum(z^k / (factorial(k) * k) for k = 1:n)
#     return SpecialFunctions.γ + log(z) + s
# end
