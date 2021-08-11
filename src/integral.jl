using LinearAlgebra
# using SpecialFunctions

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
function integrate(eq, x=nothing; abstol=1e-6, num_steps=2, num_trials=10, radius=1.0,
                   show_basis=false, opt = STLSQ(exp.(-10:1:0)), bypass=false,
                   attempt_ratio=5, symbolic=true, bypart=true, max_basis=110,
                   verbose=false)
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
                              max_basis, verbose)

    if isequal(u₁, 0) || !bypart
        return s₁, u₁, ϵ
    else
        s₂, u₂, ϵ = try_integration_by_parts(u₁, x; abstol, num_trials, num_steps,
                                             radius, show_basis, opt, attempt_ratio,
                                             symbolic, max_basis, verbose)
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

function accept_solution(eq, x, sol; abstol=1e-6)
    try
        Δ = substitute(expand_derivatives(Differential(x)(sol)-eq), Dict(x => Complex(rand())))
        return abs(Δ) < abstol
    catch e
        #
    end
    return false
end

function integrate_term(eq, x; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis, radius =
        args[:abstol], args[:num_steps], args[:num_trials], args[:show_basis],
        args[:symbolic], args[:verbose], args[:max_basis], args[:radius]

    # note that the order of the operations is important!
    # first, collecing hints, then applying transformation rules, and finally finding the basis.
    h = collect_hints(eq, x)
    eq = apply_integration_rules(eq)
    basis = generate_basis(eq, x, h)

    # basis = filter(u -> !(deg(u,x)>0), basis)

    if verbose printstyled("|β| = ", length(basis), ". "; color=:yellow) end
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
            if !isequal(y, 0) && accept_solution(eq, x, y; abstol)
                if verbose printstyled("$i, symbolic\n"; color=:yellow) end
                return y, 0, 0
            end
        end

        for j = 1:num_trials
            y, ϵ = try_integrate(Float64, eq, x, basis, Δbasis, radius*sqrt(2)^j, 1.0; kwargs...)
            if ϵ < abstol && accept_solution(eq, x, y; abstol)
                if verbose printstyled("$i, $j\n"; color=:yellow) end
                return y, 0, ϵ
            else
                ϵ₀ = min(ϵ, ϵ₀)
                y₀ = y
            end
        end
    end

    if accept_solution(eq, x, y₀; abstol)
        if verbose printstyled("rescue\n"; color=:yellow) end
        return y₀, 0, ϵ₀
    else
        return 0, eq, ϵ₀
    end
end

rms(x) = sqrt(sum(x.^2) / length(x))

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
    # b hold the value of the input function at those points
    A = zeros(Complex{T}, (n, n))
    b = zeros(Complex{T}, n)

    i = 1
    k = 1

    # select quadrant
    s1 = rand([-1,1])
    s2 = rand([-1,1])

    while i <= n
        x₀ = Complex{T}(s1*(margin + rand()*radius), s2*(margin + rand()*radius))

        d = Dict(x => x₀)
        try
            for j = 1:n
                A[i, j] = Complex{T}(substitute(Δbasis[j], d))
            end
            b[i] = Complex{T}(substitute(eq, d))
            A[i, :] .= A[i, :] ./ b[i]
            b[i] = one(T)
            i += 1
        catch e
            println("basis matrix error: ", e)
        end
        if k > attempt_ratio*n return nothing, 1e6 end
        k += 1
    end

    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)    
    A, b, basis, Δbasis, n = A[l,l], b[l], basis[l], Δbasis[l], sum(l)

    if det(A) ≈ 0 return nothing, 1e6 end

    coefs = ones(Complex{T}, n)
    for j = 1:n
        coefs[j] = coef(Δbasis[j], x)
        A[:,j] /= coefs[j]
    end

    # q₀ = A \ b
    try
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

########################## Transformation Rules ###############################

trig_rule1 = @rule tan(~x) => sin(~x) / cos(~x)
trig_rule2 = @rule sec(~x) => one(~x) / cos(~x)
# trig_rule2 = @rule sec(~x) => (tan(~x)/cos(~x) + 1/cos(~x)^2) / (tan(~x) + 1/cos(~x))
trig_rule3 = @rule csc(~x) => one(~x) / sin(~x)
# trig_rule3 = @rule csc(~x) => (cot(~x)/sin(~x) + 1/sin(~x)^2) / (cot(~x) + 1/sin(~x))
trig_rule4 = @rule cot(~x) => cos(~x) / sin(~x)

trig_rules = [trig_rule1, trig_rule2, trig_rule3, trig_rule4]

hyper_rule1 = @rule tanh(~x) => sinh(~x) / cosh(~x)
hyper_rule2 = @rule sech(~x) => one(~x) / cosh(~x)
hyper_rule3 = @rule csch(~x) => one(~x) / sinh(~x)
hyper_rule4 = @rule coth(~x) => cosh(~x) / sinh(~x)

hyper_rules = [hyper_rule1, hyper_rule2, hyper_rule3, hyper_rule4]

misc_rule1 = @rule sqrt(~x) => ^(~x, 0.5)
misc_rule2 = @acrule exp(~x) * exp(~y) => exp(~x + ~y)

misc_rules = [misc_rule1]

int_rules = [trig_rules; hyper_rules; misc_rules]
# int_rules = misc_rules

apply_integration_rules(eq) = Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq))

########################## Expansion Rules ####################################

x_rule_g1 = @rule sin(~x) => (exp(im * ~x) - exp(-im * ~x)) / 2im
x_rule_g2 = @rule cos(~x) => (exp(im * ~x) + exp(-im * ~x)) / 2
x_rule_g3 = @rule tan(~x) => -im * (exp(2im * ~x) - 1) / (exp(2im * ~x) + 1)
x_rule_g4 = @rule csc(~x) => -2im / (exp(im * ~x) - exp(-im * ~x))
x_rule_g5 = @rule sec(~x) => 2 / (exp(im * ~x) + exp(im * ~x))
x_rule_g6 = @rule cot(~x) => im * (exp(2im * ~x) + 1) / (exp(2im * ~x) - 1)

x_rule_h1 = @rule sinh(~x) => (exp(~x) - exp(-~x)) / 2
x_rule_h2 = @rule cosh(~x) => (exp(~x) + exp(-~x)) / 2
x_rule_h3 = @rule tanh(~x) => (exp(2 * ~x) - 1) / (exp(2 * ~x) + 1)
x_rule_h4 = @rule csch(~x) => 2 / (exp(~x) - exp(-~x))
x_rule_h5 = @rule sech(~x) => 2 / (exp(~x) + exp(-~x))
x_rule_h6 = @rule coth(~x) => (exp(2 * ~x) + 1) / (exp(2 * ~x) - 1)

x_rule_i1 = @rule asin(~x) => -im * log(sqrt(1-x^2) + im*x)
x_rule_i2 = @rule acos(~x) => -im * log(x + im*sqrt(1-x^2))
x_rule_i3 = @rule atan(~x) => im/2 * log((x + im) / (-x + im))
x_rule_i4 = @rule acsc(~x) => -im * log(sqrt(1-1.0/x^2) + im/x)
x_rule_i5 = @rule asec(~x) => -im * log(1.0/x + im*sqrt(1-1.0/x^2))
x_rule_i6 = @rule acot(~x) => im/2 * log((1.0/x + im) / (-1.0/x + im))

x_rule_j1 = @rule asin(~x) => log(sqrt(1+x^2) + x)
x_rule_j2 = @rule acos(~x) => log(x + sqrt(1-x^2))
x_rule_j3 = @rule atan(~x) => 1/2 * log((x + 1) / (-x + 1))
x_rule_j4 = @rule acsc(~x) => log(sqrt(1+1.0/x^2) + 1.0/x)
x_rule_j5 = @rule asec(~x) => log(1.0/x + sqrt(1-1.0/x^2))
x_rule_j6 = @rule acot(~x) => 1/2 * log((1.0/x + 1) / (-1.0/x + 1))

expansion_rules = [
    x_rule_g1,
    x_rule_g2,
    x_rule_g3,
    x_rule_g4,
    x_rule_g5,
    x_rule_g6,

    x_rule_h1,
    x_rule_h2,
    x_rule_h3,
    x_rule_h4,
    x_rule_h5,
    x_rule_h6,

    x_rule_i1,
    x_rule_i2,
    x_rule_i3,
    x_rule_i4,
    x_rule_i5,
    x_rule_i6,
]


apply_expansion(eq) = Fixpoint(Prewalk(PassThrough(Chain(expansion_rules))))(value(eq))


###############################################################################


function U(u...)
    u = map(x -> x isa AbstractArray ? x : [], u)
    return union(u...)
end

hints(eq::SymbolicUtils.Add, x, h) = map(t->hints(t,x,h), arguments(eq))
hints(eq::SymbolicUtils.Mul, x, h) = map(t->hints(t,x,h), arguments(eq))
hints(eq::SymbolicUtils.Pow, x, h) = hints(arguments(eq)[1],x,h)

function hints(eq::SymbolicUtils.Term, x, h)
    s = Symbol(operation(eq))
    u = arguments(eq)[1]

    if s == :sec
        push!(h, log(1/cos(u) + sin(u)/cos(u)))
    elseif s == :csc
        push!(h, log(1/sin(u) - cos(u)/sin(u)))
    elseif s == :tan
        push!(h, log(cos(u)))
    elseif s == :cot
        push!(h, log(sin(u)))
    elseif s == :tanh
        push!(h, log(cosh(u)))
    elseif s == :log
        if check_poly(u, x) == :real_poly && deg(u, x) == 2
            r, s = find_roots(u, x)
            if !isempty(r)
                push!(h, log(x - r[1]))
                push!(h, log(x - r[2]))
            else
                push!(h, log(u))
            end
        end
    end
end

function hints(eq, x, h)
end

function collect_hints(eq, x)
    h = []
    hints(eq, x, h)
    h
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

# """
#     logarithmic integral
# """
# function li(x; n=10)
#     z = log(abs(x))
#     s = sum(z^k / (factorial(k) * k) for k = 1:n)
#     return SpecialFunctions.γ + log(z) + s
# end
