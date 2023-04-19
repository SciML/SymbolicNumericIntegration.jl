using LinearAlgebra
using Statistics: mean, std

Base.signbit(z::Complex{T}) where {T <: Number} = signbit(real(z))
Base.signbit(x::SymbolicUtils.Sym{Number, Nothing}) = false

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
    homotopy: use the homotopy algorithm to generate the basis

    output:
    -------
    solved, unsolved, err

    solved is the solved integral and unsolved is the residual unsolved portion of the input
    err is the numerical error in reaching the solution
"""
function integrate(eq, x = nothing; abstol = 1e-6, num_steps = 2, num_trials = 5,
                   radius = 1.0,
                   show_basis = false, opt = STLSQ(exp.(-10:1:0)), bypass = false,
                   symbolic = true, max_basis = 100, verbose = false, complex_plane = true,
                   homotopy = true, use_optim=false)
    eq = expand(eq)
    eq = apply_div_rule(eq)

    if x == nothing
        x = var(eq)
        if x == nothing
            @syms 𝑥
            x = 𝑥
        end
    else
        x = value(x)    # needed for the transition from @syms to @variables
    end

    l = Logger(verbose)

    # eq is a constant
    if !isdependent(eq, x)
        return x * eq, 0, 0
    end

    # homotopy = homotopy && !bypass

    s, u, ϵ = integrate_sum(eq, x, l; bypass, abstol, num_trials, num_steps,
                            radius, show_basis, opt, symbolic,
                            max_basis, verbose, complex_plane, homotopy, use_optim)
    return simplify(s), u, ϵ
end

"""
    integrate_sum applies the integral summation rule ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx

    eq: the integrand
    x: the indepedent variable
    l: a logger

    output is the same as integrate
"""
function integrate_sum(eq, x, l; bypass = false, kwargs...)
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
                       args[:num_trials], args[:show_basis], args[:symbolic],
                       args[:verbose],
                       args[:max_basis], args[:radius], args[:homotopy]

    attempt(l, "Integrating term", eq)

    if is_number(eq)
        y = eq * x
        result(l, "Successful", y)
        return y, 0, 0
    end

    eq = cache(eq)
    basis = generate_basis(eq, x; homotopy)

    if show_basis
        inform(l, "Generating basis (|β| = $(length(basis)))", basis)
    end

    if length(basis) > max_basis
        result(l, "|β| = $(length(basis)) is too large")
        return 0, expr(eq), Inf
    end

    # D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    # rescue
    ϵᵣ = Inf
    yᵣ = 0

    for i in 1:num_steps
        if length(basis) > max_basis
            break
        end

        if symbolic
            y, ϵ = try_symbolic(Float64, expr(eq), x, expr.(basis), deriv!.(basis, x);
                                kwargs...)

            if !isequal(y, 0) && accept_solution(eq, x, y, radius) < abstol
                result(l, "Successful symbolic", y)
                return y, 0, 0
            else
                inform(l, "Failed symbolic")
            end
        end

        for j in 1:num_trials
            r = radius #*sqrt(2)^j
            y, ϵ = try_integrate(Float64, eq, x, basis, r; kwargs...)

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
        return 0, expr(eq), ϵ₀
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
function try_integrate(T, eq, x, basis, radius; kwargs...)
    args = Dict(kwargs)
	use_optim = args[:use_optim]
    basis = basis[2:end]    # remove 1 from the beginning
	
	if use_optim
		return solve_optim(T, eq, x, basis, radius; kwargs...)
	else
		return solve_sparse(T, eq, x, basis, radius; kwargs...)
	end
end

