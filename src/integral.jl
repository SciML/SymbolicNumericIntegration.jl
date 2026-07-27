"""
    integrate(eq, x; kwargs...)

Compute an indefinite or definite integral of a scalar symbolic expression.

`integrate` combines symbolic candidate generation with numerical verification. It
accepts expressions constructed by Symbolics.jl and returns a detailed result by
default so callers can distinguish a complete solution from a partial one.

```julia
julia> using Symbolics, SymbolicNumericIntegration

julia> @variables x a

julia> integrate(x * sin(2x))
((1//4)*sin(2x) - (1//2)*x*cos(2x), 0, 0)

julia> integrate(x * sin(a * x), x; symbolic = true, detailed = false)
(sin(a*x) - a*x*cos(a*x)) / (a^2)

julia> integrate(x * sin(a * x), (x, 0, 1); symbolic = true, detailed = false)
(sin(a) - a*cos(a)) / (a^2)
```

## Arguments

  - `eq`: Scalar symbolic integrand. Arrays are rejected; use broadcasting for
    elementwise integration.
  - `x`: Independent symbolic variable. Omit it only for an expression with exactly
    one variable. Pass `(x, lower, upper)` to request a definite integral.

## Keyword Arguments

  - `abstol` (default: `1e-6`): the desired tolerance
  - `num_steps` (default: `2`): the number of different steps with expanding basis to be tried
  - `num_trials` (default: `10`): the number of trials in each step (no changes to the basis)
  - `show_basis` (default: `false`): Print the generated candidate basis.
  - `bypass` (default: `false`): Solve the complete expression rather than splitting
    it into terms.
  - `symbolic` (default: `false`): Attempt the symbolic solver first. This is forced
    when the integrand has symbolic constants.
  - `max_basis` (default: `100`): the maximum number of candidate terms to consider
  - `verbose` (default: `false`): Print solver diagnostics.
  - `complex_plane` (default: `true`): Sample verification points in the complex
    plane. Set to `false` to sample on the real axis.
  - `radius` (default: `5.0`): Radius of the verification-point disk.
  - `opt` (default: `STLSQ(exp.(-10:1:0))`): the sparse regression optimizer (from DataDrivenSparse)
  - `homotopy` (default: `true`): *deprecated*, will be removed in a future version
  - `use_optim` (default: `false`): *deprecated*, will be removed in a future version
  - `detailed` (default: `true`): Return `(solved, unsolved, err)`. When `false`,
    return only a complete solution or `nothing`.

## Returns

When `detailed = true`, return `(solved, unsolved, err)`, where

    solved: the solved integral
    unsolved: the residual unsolved portion of the input
    err: the numerical error in reaching the solution

When `detailed = false`, return the complete integral or `nothing`.
"""
function integrate(
        eq, x = nothing;
        abstol = 1.0e-6,
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
        detailed = true
    )
    deprecation_warnings(; homotopy, use_optim)

    # Check if eq is a vector/array expression and throw an error
    if eq isa AbstractArray || eq isa AbstractVector
        error("Vector expressions are not supported. Please use element-wise integration with `integrate.([expr1, expr2, ...], x)` instead.")
    end

    eq = expand(eq)

    if x === nothing
        vars = get_variables(eq)
        if length(vars) > 1
            error("Multiple symbolic variables detect. Please pass the independent variable to `integrate`")
        elseif length(vars) == 1
            x = first(vars)
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

    s, u,
        ε = integrate_sum(
        eq, x; plan, bypass, num_trials, num_steps,
        show_basis, symbolic, max_basis, verbose, use_optim
    )

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

# Helper to extract numeric value from symbolic expression
function extract_numeric_value(expr)
    # Already a number
    if expr isa Number
        return expr
    end

    # Unwrap Num type
    if expr isa Num
        expr = value(expr)
        if expr isa Number
            return expr
        end
    end

    # Try Float64 conversion (works for BasicSymbolic numeric literals)
    try
        result = Float64(expr)
        return result
    catch
    end

    # Try converting to Julia expression and evaluating
    # This handles cases like sin(0) that need to be evaluated
    try
        julia_expr = toexpr(expr)
        result = Base.invokelatest(eval, julia_expr)
        if result isa Number
            return result
        end
    catch
    end

    # Return as-is if conversion fails
    return expr
end

# Evaluate an expression at a bound, using limit for infinite bounds
function eval_at_bound(expr, x, bound)
    # Unwrap both expression and variable for consistent handling
    expr_unwrapped = value(expr)
    x_unwrapped = value(x)

    if isinf(bound)
        # Use symbolic limits for infinite bounds
        try
            result = limit(expr_unwrapped, x_unwrapped, bound)
            # limit returns a tuple (value, assumptions), extract the value
            if result isa Tuple
                result = first(result)
            end
            return extract_numeric_value(result)
        catch e
            # If limit computation fails, fall back to direct substitution
            # This may result in NaN for indeterminate forms
            return substitute(expr_unwrapped, Dict(x_unwrapped => bound))
        end
    else
        result = substitute(expr_unwrapped, Dict(x_unwrapped => bound))
        return extract_numeric_value(result)
    end
end

# Definite integral
function integrate(eq, xx::Tuple; kwargs...)
    x, lo, hi = xx
    sol = integrate(eq, x; kwargs...)

    if sol isa Tuple
        if !isequal(first(sol), 0) && isequal(sol[2], 0)
            hi_val = eval_at_bound(first(sol), x, hi)
            lo_val = eval_at_bound(first(sol), x, lo)
            result = extract_numeric_value(hi_val - lo_val)
            # Check if the result is valid (not NaN or undefined)
            if result isa Number && (isnan(result) || isinf(result))
                return nothing
            end
            return result
        else
            return nothing
        end
    elseif sol !== nothing
        hi_val = eval_at_bound(sol, x, hi)
        lo_val = eval_at_bound(sol, x, lo)
        result = extract_numeric_value(hi_val - lo_val)
        # Check if the result is valid (not NaN or undefined)
        if result isa Number && (isnan(result) || isinf(result))
            return nothing
        end
        return result
    end

    return nothing
end

function get_solved(p, sol)
    if sol isa Tuple
        s = sol[1]
        return s === nothing ? 0 : s
    else
        return sol === nothing ? 0 : sol
    end
end

function get_unsolved(p, sol)
    if sol isa Tuple
        u = sol[2]
        return u === nothing ? 0 : u
    else
        return isequal(sol, 0) || sol === nothing ? p : 0
    end
end

function get_err(p, sol)
    if sol isa Tuple
        return sol[3]
    else
        return isequal(sol, 0) || sol === nothing ? Inf : 0
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
    plan, num_steps,
        num_trials,
        show_basis,
        symbolic,
        verbose,
        max_basis = args[:plan],
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

    if y === nothing
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

    return if !homotopy
        @warn("homotopy is deprecated and will be removed in a future version")
    end
end
