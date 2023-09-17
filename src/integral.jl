using LinearAlgebra
using Statistics: mean, std

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

Keyword Arguments:
------------------
- `abstol` (default: `1e-6`): the desired tolerance
- `num_steps` (default: `2`): the number of different steps with expanding basis to be tried
- `num_trials` (default: `10`): the number of trials in each step (no changes to the basis)
- `show_basis` (default: `false`): if true, the basis (list of candidate terms) is printed
- `bypass` (default: `false`): if true do not integrate terms separately but consider all at once
- `symbolic` (default: `false`): try symbolic integration first (will be forced if `eq` has constant parameters)
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
function integrate(eq, x = nothing; abstol = 1e-6, num_steps = 2, num_trials = 10,
                   radius = 1.0,
                   show_basis = false, opt = STLSQ(exp.(-10:1:0)), bypass = false,
                   symbolic = false, max_basis = 100, verbose = false, complex_plane = true,
                   homotopy = true, use_optim = false, detailed = true)
    eq = expand(eq)

    if x == nothing
        if !isempty(sym_consts(eq, x))
            error("Multiple symbolic variables detect. Please pass the independent variable to `integrate`")            
        end
        x = var(eq)
        if x == nothing
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

    s, u, ε = integrate_sum(eq, x; bypass, abstol, num_trials, num_steps,
                            radius, show_basis, opt, symbolic,
                            max_basis, verbose, complex_plane, use_optim)
    
    if detailed
        return s, u, ε
    else
        return isequal(s, 0) ? nothing : s
    end
end

"""
    integrate_sum(eq, x; kwargs...) 

applies the integral summation rule ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx

inputs:
------
- eq: the integrand
- x: the indepedent variable

The output is the same as `integrate`
"""
function integrate_sum(eq, x; bypass = false, kwargs...)
    solved = 0
    unsolved = 0
    ε₀ = 0
    ts = bypass ? [eq] : terms(eq)

    for p in ts
        s, u, ε = integrate_term(p, x; kwargs...)
        solved += s
        unsolved += u
        ε₀ = max(ε₀, ε)
    end

    if !isequal(unsolved, 0) && isempty(const_params(unsolved, x))
        eq = factor_rational(simplify_trigs(unsolved))

        if !isequal(eq, unsolved)
            eq = expand(eq)
            unsolved = 0
            ε₀ = 0
            ts = bypass ? [eq] : terms(eq)

            for p in ts
                s, u, ε = integrate_term(p, x; kwargs...)
                solved += s
                unsolved += u
                ε₀ = max(ε₀, ε)

                if !isequal(u, 0)   # premature termination on the first failure
                    return 0, eq, ε₀
                end
            end
        end
    end

    return expand(solved), unsolved, ε₀
end

"""
    integrate_term(eq, x; kwargs...) 
    
is the central part of the code that tries different methods to integrate `eq`,
which is assume to be a single term.

inputs:
-------
- eq: the integrand
- x: the indepedent variable

The output is the same as `integrate`
"""
function integrate_term(eq, x; kwargs...)
    args = Dict(kwargs)
    abstol, num_steps, num_trials, show_basis, symbolic, verbose, max_basis,
    radius = args[:abstol], args[:num_steps],
             args[:num_trials], args[:show_basis], args[:symbolic],
             args[:verbose],
             args[:max_basis], args[:radius]

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
        y = integrate_symbolic(eq, x; abstol, radius)
        if y == nothing
            if has_sym_consts
                @info("Symbolic integration failed. Try changing constant parameters ([$(join(params, ", "))]) to numerical values.")
                return 0, eq, Inf
            end
        else
            return y, 0, 0
        end
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

    # D = Differential(x)
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
            r = radius 
            y, ε = try_integrate(eq, x, basis, r; kwargs...)

            ε = accept_solution(eq, x, y, r)
            if ε < abstol
                return y, 0, ε
            elseif ε < εᵣ
                εᵣ = ε
                yᵣ = y
            end
        end     

        if i < num_steps
            basis1, ok1 = expand_basis(prune_basis(eq, x, basis1, radius; kwargs...), x)
            basis2, ok2 = expand_basis(prune_basis(eq, x, basis2, radius; kwargs...), x)

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
