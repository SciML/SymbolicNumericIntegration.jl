using LinearAlgebra
using Statistics: mean, std

Base.signbit(z::Complex{T}) where {T <: Number} = signbit(real(z))
Base.signbit(x::SymbolicUtils.Sym{Number}) = false

"""
    integrate(eq, x; kwargs...)

is the main entry point to integrate a univariate expression `eq` with respect to `x' (optional).

```julia
julia> using Symbolics, SymbolNumericIntegration

julia> @variables x a

julia> integrate(x * sin(2x))
((1//4)*sin(2x) - (1//2)*x*cos(2x), 0, 0)

julia> integrate(x * sin(a * x), x; symbolic = true, detailed = false)
(sin(a*x) - a*x*cos(a*x)) / (a^2)

julia> integrate(x * sin(a * x), (x, 0, 1); symbolic = true, detailed = false)
(sin(a) - a*cos(a)) / (a^2)
```

## Arguments:

  - `eq`: a univariate expression
  - `x`: independent variable (optional if `eq` is univariate) or a tuple
    of (independent variable, lower bound, upper bound) for definite integration.

## Keyword Arguments:

  - `abstol` (default: `1e-6`): the desired tolerance
  - `num_steps` (default: `2`): the number of different steps with expanding basis to be tried
  - `num_trials` (default: `10`): the number of trials in each step (no changes to the basis)
  - `show_basis` (default: `false`): if true, the basis (list of candidate terms) is printed
  - `bypass` (default: `false`): if true do not integrate terms separately but consider all at once
  - `symbolic` (default: `false`): try symbolic integration first (will be forced if `eq` has symbolic constants)
  - `max_basis` (default: `100`): the maximum number of candidate terms to consider
  - `verbose` (default: `false`): print a detailed report
  - `complex_plane` (default: `true`): generate random test points on the complex plane (if false, the points will be on real axis)
  - `radius` (default: `1.0`): the radius of the disk in the complex plane to generate random test points
  - `opt` (default: `STLSQ(exp.(-10:1:0))`): the sparse regression optimizer (from DataDrivenSparse)
  - `homotopy` (default: `true`): *deprecated*, will be removed in a future version
  - `use_optim` (default: `false`): *deprecated*, will be removed in a future version
  - `detailed` (default: `true`): `(solved, unsolved, err)` output format. If `detailed=false`, only the final integral is returned.

Returns a tuple of (solved, unsolved, err) if `detailed == true`, where

    solved: the solved integral 
    unsolved: the residual unsolved portion of the input
    err: the numerical error in reaching the solution

Returns the resulting integral or nothing if `detailed == false`
"""
function integrate(eq, x = nothing;
        abstol = 1e-6,
        num_steps = 2,
        num_trials = 10,
        radius = 5.0,
        show_basis = false,
        opt = STLSQ(exp.(-10:1:0)),
        bypass = false,
        symbolic = false,
        max_basis = 100,
        verbose = false,
        complex_plane = true,
        homotopy = true,
        use_optim = false,
        detailed = true)
    deprecation_warnings(; homotopy, use_optim)

    eq = expand(eq)

    if x == nothing
        vars = get_variables(eq)
        if length(vars) > 1
            error("Multiple symbolic variables detect. Please pass the independent variable to `integrate`")
        elseif length(vars) == 1
            x = vars[1]
        else
            @syms 𝑥
            x = 𝑥
        end
    else
        x = value(x)    # needed for the transition from @syms to @variables
    end

    # eq is a constant
    if !isdependent(eq, x)
        if detailed
            return x * eq, 0, 0
        else
            return x * eq
        end
    end

    plan = NumericalPlan(abstol, radius, complex_plane, opt)

    s, u, ε = integrate_sum(eq, x; plan, bypass, num_trials, num_steps,
        show_basis, symbolic, max_basis, verbose, use_optim)

    s = beautify(s)

    if detailed
        return s, u, ε
    else
        if !isequal(s, 0) && !isequal(u, 0)
            @info("Integration is partially successful. Pass `detailed = true` to `integrate` for details")
        end

        return isequal(s, 0) || !isequal(u, 0) ? nothing : s
    end
end

# Definite integral
function integrate(eq, xx::Tuple; kwargs...)
    x, lo, hi = xx
    sol = integrate(eq, x; kwargs...)

    if sol isa Tuple
        if first(sol) != 0 && sol[2] == 0
            return substitute(first(sol), Dict(x => hi)) -
                   substitute(first(sol), Dict(x => lo))
        else
            return nothing
        end
    elseif sol != nothing
        return substitute(sol, Dict(x => hi)) - substitute(sol, Dict(x => lo))
    end

    return nothing
end

function get_solved(p, sol)
    if sol isa Tuple
        s = sol[1]
        return s == nothing ? 0 : s
    else
        return sol == nothing ? 0 : sol
    end
end

function get_unsolved(p, sol)
    if sol isa Tuple
        u = sol[2]
        return u == nothing ? 0 : u
    else
        return sol == 0 || sol == nothing ? p : 0
    end
end

function get_err(p, sol)
    if sol isa Tuple
        return sol[3]
    else
        return sol == 0 || sol == nothing ? Inf : 0
    end
end

# integrate_sum applies the integral summation rule ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
function integrate_sum(eq, x; bypass = false, kwargs...)
    solved = 0
    unsolved = 0
    ε₀ = 0
    ts = bypass ? [eq] : terms(eq)

    for p in ts
        sol = integrate_term(p, x; kwargs...)
        solved += get_solved(p, sol)
        unsolved += get_unsolved(p, sol)
        ε₀ = max(ε₀, get_err(p, sol))
    end

    if !isequal(unsolved, 0) && isempty(sym_consts(unsolved, x))
        eq = factor_rational(simplify_trigs(unsolved))

        if !isequal(eq, unsolved)
            eq = expand(eq)
            unsolved = 0
            ε₀ = 0
            ts = bypass ? [eq] : terms(eq)

            for p in ts
                sol = integrate_term(p, x; kwargs...)
                solved += get_solved(p, sol)
                unsolved += get_unsolved(p, sol)
                ε₀ = max(ε₀, get_err(p, sol))

                if !isequal(u, 0)   # premature termination on the first failure
                    return 0, eq, ε₀
                end
            end
        end
    end

    return expand(solved), unsolved, ε₀
end

# integrate_term is the central part of the univariate integration code that 
# tries different methods to integrate `eq`.
#
# note: this function will be replaced with solver(prob::Problem) in symbolic.jl
# in a future version
function integrate_term(eq, x; kwargs...)
    args = Dict(kwargs)
    plan, num_steps, num_trials, show_basis, symbolic, verbose, max_basis = args[:plan],
    args[:num_steps], args[:num_trials], args[:show_basis],
    args[:symbolic], args[:verbose], args[:max_basis]

    abstol = plan.abstol

    if is_number(eq)
        y = eq * x
        return y, 0, 0
    end

    params = sym_consts(eq, x)
    has_sym_consts = !isempty(params)

    if has_sym_consts && !symbolic
        @info("The input expression has constant parameters: [$(join(params, ", "))], forcing `symbolic = true`")
        symbolic = true
    end

    if symbolic
        return try_symbolic(eq, x, has_sym_consts, params)
    end

    eq = cache(eq)
    basis1 = generate_basis(eq, x, false)

    if has_sym_consts
        # kernel-based ansatz generator does not work correctly with sym consts
        basis2 = basis1
    else
        basis2 = generate_basis(eq, x, true)
    end

    if show_basis
        @info("Generating basis (|β| = $(length(basis1))): $basis1")
    end

    if length(basis1) > max_basis
        return 0, expr(eq), Inf
    end

    ε₀ = Inf
    y₀ = 0

    # rescue
    εᵣ = Inf
    yᵣ = 0

    for i in 1:num_steps
        if length(basis1) > max_basis
            break
        end

        for j in 1:num_trials
            basis = isodd(j) ? basis1 : basis2
            y, ε = try_integrate(eq, x, basis; plan)
            ε = accept_solution(eq, x, y; plan)

            if ε < abstol
                return y, 0, ε
            elseif ε < εᵣ
                εᵣ = ε
                yᵣ = y
            end
        end

        if i < num_steps
            basis1, ok1 = expand_basis(prune_basis(eq, x, basis1; plan), x)
            basis2, ok2 = expand_basis(prune_basis(eq, x, basis2; plan), x)

            if !ok1 && ~ok2
                break
            end
        end
    end

    if εᵣ < abstol * 10
        return yᵣ, 0, εᵣ
    else
        return 0, expr(eq), ε₀
    end
end

# try_integrate is the main dispatch point to call different sparse solvers. 
# It tries to find a linear combination of the basis, whose derivative is 
# equal to eq
function try_integrate(eq, x, basis; plan = default_plan())
    if isempty(basis)
        return 0, Inf
    end

    # return solve_optim(eq, x, basis; plan)
    return solve_sparse(eq, x, basis; plan)
end

function try_symbolic(eq, x, has_sym_consts = false, params = []; plan = default_plan())
    y = integrate_symbolic(eq, x; plan)

    if y == nothing
        if has_sym_consts && !isempty(params)
            @info("Symbolic integration failed. Try changing constant parameters ([$(join(params, ", "))]) to numerical values.")
        end

        return 0, eq, Inf
    else
        return y, 0, 0
    end
end

function deprecation_warnings(; use_optim = false, homotopy = true)
    if use_optim
        @warn("use_optim is deprecated and will be removed in a future version")
    end

    if !homotopy
        @warn("homotopy is deprecated and will be removed in a future version")
    end
end
