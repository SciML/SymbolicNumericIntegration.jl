using LinearAlgebra
using Statistics: mean, std
using Symbolics

Base.signbit(z::Complex{T}) where {T <: Number} = signbit(real(z))
Base.signbit(x::SymbolicUtils.Sym{Number}) = false

"""
    integrate(eq, x; kwargs...)
   
is the main entry point to integrate a univariate expression `eq` with respect to `x' (optional). 

```julia
integrate(x * sin(2x))

# output

((1//4)*sin(2x) - (1//2)*x*cos(2x), 0, 0)
```

Arguments:
----------
- `eq`: a univariate expression
- `x`: the independent variable (optional)
- `domain`: the domain upon which to evaluate a definite integral (optional)

Keyword Arguments:
------------------
- `abstol` (default: `1e-6`): the desired tolerance
- `num_steps` (default: `2`): the number of different steps with expanding basis to be tried
- `num_trials` (default: `10`): the number of trials in each step (no changes to the basis)
- `show_basis` (default: `false`): if true, the basis (list of candidate terms) is printed
- `bypass` (default: `false`): if true do not integrate terms separately but consider all at once
- `symbolic` (default: `false`): try symbolic integration first
- `max_basis` (default: `100`): the maximum number of candidate terms to consider
- `verbose` (default: `false`): print a detailed report
- `complex_plane` (default: `true`): generate random test points on the complex plane (if false, the points will be on real axis)
- `radius` (default: `1.0`): the radius of the disk in the complex plane to generate random test points
- `opt` (default: `STLSQ(exp.(-10:1:0))`): the sparse regression optimizer (from DataDrivenSparse)
- `homotopy` (default: `true`): use the homotopy algorithm to generate the basis (*deprecated*, will be removed in a future version)
- `use_optim` (default: `false`): use Optim.jl `minimize` function instead of the STLSQ algorithm (*experimental*)

Output:
-------
- `solved`: the solved integral 
- `unsolved`: the residual unsolved portion of the input
- `err`: the numerical error in reaching the solution
"""
function integrate(eq, x = nothing, domain::Vector{<:Number} = nothing; abstol = 1e-6, num_steps = 2, num_trials = 10,
                   radius = 1.0,
                   show_basis = false, opt = STLSQ(exp.(-10:1:0)), bypass = false,
                   symbolic = true, max_basis = 100, verbose = false, complex_plane = true,
                   homotopy = true, use_optim = false)
    eq = expand(eq)

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

    s, u, ϵ = integrate_sum(eq, x, l; bypass, abstol, num_trials, num_steps,
                            radius, show_basis, opt, symbolic,
                            max_basis, verbose, complex_plane, use_optim)
    if domain != nothing
        s = substitute( s, Dict( [ x=>domain[1] ] ) ) - substitute( s, Dict( [ x=>domain[2] ] ) )
    end
    # return simplify(s), u, ϵ
    return s, u, ϵ
end

"""
    integrate_sum(eq, x, l; kwargs...) 

applies the integral summation rule ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx

inputs:
------
- eq: the integrand
- x: the indepedent variable
- l: a logger

The output is the same as `integrate`
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
    integrate_term(eq, x, l; kwargs...) 
    
is the central part of the code that tries different methods to integrate `eq`,
which is assume to be a single term.

inputs:
-------
- eq: the integrand
- x: the indepedent variable
- l: a logger

The output is the same as `integrate`
"""
function integrate_term(eq, x, l; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis,
    radius = args[:abstol], args[:num_steps],
             args[:num_trials], args[:show_basis], args[:symbolic],
             args[:verbose],
             args[:max_basis], args[:radius]

    attempt(l, "Integrating term", eq)

    if is_number(eq)
        y = eq * x
        result(l, "Successful", y)
        return y, 0, 0
    end

    eq = cache(eq)
    basis1 = generate_basis(eq, x, true)
    basis2 = generate_basis(eq, x, false)

    if show_basis
        inform(l, "Generating basis (|β| = $(length(basis1)))", basis1)
    end

    if length(basis1) > max_basis
        result(l, "|β| = $(length(basis1)) is too large")
        return 0, expr(eq), Inf
    end

    # D = Differential(x)
    ϵ₀ = Inf
    y₀ = 0

    # rescue
    ϵᵣ = Inf
    yᵣ = 0

    for i in 1:num_steps
        if length(basis1) > max_basis
            break
        end

        if symbolic
            y, ϵ = try_symbolic(Float64, expr(eq), x, expr.(basis1), deriv!.(basis1, x);
                                kwargs...)

            if !isequal(y, 0) && accept_solution(eq, x, y, radius) < abstol
                result(l, "Successful symbolic", y)
                return y, 0, 0
            else
                inform(l, "Failed symbolic")
            end
        end

        for j in 1:num_trials
            basis = isodd(j) ? basis1 : basis2
            r = radius #*sqrt(2)^j
            y, ϵ = try_integrate(eq, x, basis, r; kwargs...)

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
            basis1, ok1 = expand_basis(prune_basis(eq, x, basis1, radius; kwargs...), x)
            basis2, ok2 = expand_basis(prune_basis(eq, x, basis2, radius; kwargs...), x)

            if !ok1 && ~ok2
                break
            end

            if show_basis
                inform(l, "Expanding the basis (|β| = $(length(basis)))", basis1)
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
	try_integrate(eq, x, basis, radius; kwargs...) 
	
is the main dispatch point to call different sparse solvers. It tries to 
find a linear combination of the basis, whose derivative is equal to eq

output:
-------
- solved: the solved integration problem or 0 otherwise
- err: the numerical error in reaching the solution
"""
function try_integrate(eq, x, basis, radius; kwargs...)
    args = Dict(kwargs)
    use_optim = args[:use_optim]

    if isempty(basis)
        return 0, Inf
    end

    if use_optim
        return solve_optim(eq, x, basis, radius; kwargs...)
    else
        return solve_sparse(eq, x, basis, radius; kwargs...)
    end
end

#################################################################################

"""
	integrate_basis(eq, x; kwargs...)
	
is used for debugging and should not be called in the course of normal execution
"""
function integrate_basis(eq, x = var(eq); abstol = 1e-6, radius = 1.0, complex_plane = true)
    eq = cache(expand(eq))
    basis = generate_basis(eq, x, false)
    n = length(basis)
    A, X = init_basis_matrix(eq, x, basis, radius, complex_plane; abstol)
    return basis, A, X
end
