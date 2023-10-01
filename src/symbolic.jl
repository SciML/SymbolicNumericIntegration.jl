"""
    Convers floats to integers/rational numbers with small denominators 
    if possible
"""
function beautify(eq)
    if is_add(eq)
        return sum(beautify(t) for t in args(eq))
    elseif is_mul(eq)
        return prod(beautify(t) for t in args(eq))
    elseif is_pow(eq)
        p = args(eq)[1]
        k = args(eq)[2]
        return beautify(p)^k
    elseif is_div(eq)
        return beautify(numer(eq)) / beautify(denom(eq))
    elseif is_number(eq)
        return nice_parameter(eq)
    else
        return eq
    end
end

###################################################################

"""
    Removes numerical coefficients from terms
"""
function equiv(y, x)
    y = expand(value(y))

    if is_add(y)
        return sum(equiv(u, x) for u in args(y))
    elseif is_mul(y)
        return prod(isdependent(u, x) ? u : 1 for u in args(y))        
    elseif is_div(y)
        return expand_fraction(y, x)
    elseif is_number(y)
        return 1
    else
        return y
    end
end

"""
    Splits an expression into a list of terms
"""
function split_terms(eq, x)
    eq = value(eq)
    if is_add(eq)
        return [equiv(u, x) for u in args(eq)]
    else
        return [equiv(eq, x)]
    end
end

"""
    Checks whether y is a holonomic function, i.e., is closed 
    under differentiation w.r.t. x
    
    For our purpose, we define a holonomic function as one composed
    of the positive powers of x, sin, cos, exp, sinh, and cosh
    
    Args:
        y: the expression to check for holonomy
        x: independent variable
        
    Returns:
        true if y is holonomic
"""
function is_holonomic(y, x)
    y = value(y)

    # technically, polynomials are not holonomic, but practically we can include them
    if is_polynomial(y, x)
        return true
    end

    if is_number(y)
        return true
    end

    if is_term(y)
        return operation(y) in [sin, cos, sinh, cosh, exp] &&
               is_polynomial(args(y)[1], x)
    end

    if is_pow(y)
        p = numer(y)
        k = denom(y)
        return is_holonomic(p, x) && k isa Integer && k > 0
    end

    if is_add(y) || is_mul(y)
        return all(is_holonomic(t, x) for t in args(y))
    end

    return false
end

"""
    Generates a list of ansatzes based on repetative differentiation. 
    It works for holonomic functions, which are closed under differentiation.
"""
function blender(y, x; n = 3)
    basis = value(y)
    for i in 1:n
        new_basis = basis
        for t in split_terms(basis, x)
            new_basis += equiv(diff(t, x), x)
        end
        if isequal(new_basis, basis)
            break
        else
            basis = new_basis
        end
    end

    basis = expand(x + basis + x * basis)
    return split_terms(basis, x)
end

############################## Utils ###############################

"""
    Returns a dictionary of the symbolic constants in the expression
    and random real value assignments.

    Args:
        eq: the integrands
        x: independent variable

    Returns:
        a dict of sym : value pairs
"""
function subs_symbols(eq, x; include_x = false, radius = 5.0, as_complex=true)
    S = Dict()
    for v in get_variables(value(eq))
        if !isequal(v, x)
            S[v] = as_complex ? Complex(randn()) : randn()
        end
    end

    if include_x
        S[x] = as_complex ? Complex(randn()) : randn()
    end

    return S
end


"""
    Splits terms into a part dependent on x and a part constant w.r.t. x    
    For example, `atomize(a*sin(b*x), x)` is `(a, sin(b*x))`)
"""
function atomize(eq, xs...)
    eq = value(eq)

    if is_mul(eq)
        coef = 1
        atom = 1
        for t in args(eq)
            c, a = atomize(t, xs...)
            coef *= c
            atom *= a
        end
        return (coef, atom)
    elseif is_div(eq)
        c1, a1 = atomize(numer(eq), xs...)
        c2, a2 = atomize(denom(eq), xs...)
        return (c1 / c2, a1 / a2)
    elseif any(isdependent(eq, x) for x in xs)
        return (1, eq)
    else
        return (eq, 1)
    end
end


"""
    Generates the final integral based on the list of coefficients and
    a list of ansatzes (ker).
"""
function apply_coefs(q, ker, x)
    s = 0
    for (coef, expr) in zip(q, ker)
        s += beautify(coef) * expr
    end

    try
        s = simplify(s)
    catch e
        # println(e)
    end

    c, a = atomize(s, get_variables(s)...)
    return beautify(c) * a
end


function expand_basis_symbolic(basis, x)
    b = sum(basis)
    Δb = sum(diff(y, x) for y in basis)
    basis = split_terms(expand((1+x)*(b + Δb)), x)    
    
    return basis
end


complex_from_num(x) = Complex(value(real(x)), value(imag(x)))


################################ Refactoring ##########################################

"""
    The main entry point for symbolic integration.
    
    Argu:
        eq: the expression to integrate
        x: independent variable
        
    Returns:
        the integral or nothing if no solution
"""
function integrate_symbolic(eq, x; plan = default_plan(), num_steps=1)
    prob = Problem(eq, x; plan, num_steps)

    if prob != nothing
        return solver(prob)
    end
    
    return nothing
end


Expression = Union{Num, BasicSymbolic, Number}

struct Problem
    eq::Expression          # integrand
    x::Expression           # independent variable
    coef::Expression        # coefficient of the integrand
    ker::Array{Expression}  # the pruned list of the basis expressions (kernel)
    plan::NumericalPlan     # the numerial plan, containing various parameters
    num_steps::Int          # the number of expansion attempts
end


function Problem(eq, x; plan = default_plan(), num_steps=1)
    eq = expand(eq)
    coef, eq = atomize(eq, x)

    if is_holonomic(eq, x)
        basis = blender(eq, x)
    else
        basis = generate_homotopy(eq, x)
    end
    
    ker = best_hints(eq, x, basis; plan)

    if ker == nothing
        basis = expand_basis_symbolic(basis, x)
        ker = best_hints(eq, x, basis; plan)
        if ker == nothing
            return nothing
        end
    end

    ker = [atomize(y, x)[2] for y in ker]
    
    return Problem(
        eq,
        x,
        coef,
        ker,
        plan,
        num_steps
    )
end


# returns true if sol solves the integration problem represented by prob
function accept(prob, sol)
    if sol != nothing
        ε = accept_solution(prob.eq, prob.x, sol; prob.plan)
        return ε < prob.plan.abstol
    end
    return false
end


# `solve` would be a better name but is confused with the Symbolic solve function.
function solver(prob::Problem)
    # First, we attempt fully symbolic integration
    alg = SymbolicIntegrator(prob)
    sol = solver(prob, alg)

    if accept(prob, sol)
        return sol * prob.coef
    end
    
    # Next, dense numeric integration
    alg = make_square(prob, alg)
    if alg != nothing
        sol = solver(prob, alg)
    
        if accept(prob, sol)
            return sol * prob.coef
        end
    end
    
    # Finally, sparse numeric integration
    alg = NumericIntegrator(prob)
    sol = solver(prob, alg)

    if accept(prob, sol)
        return sol * prob.coef
    end
    
    return nothing
end


subs_symbols(prob::Problem) = subs_symbols(prob.eq, prob.x)


abstract type IntegrationAlgorithm end


struct SymbolicIntegrator <: IntegrationAlgorithm
   eqs::Vector{Equation}        # list of equations of the form Σ θ[i]*frag[i] ~ 0  
   vars::Vector{Expression}     # list of dummy variables (θ[1], θ[2], ...)
   frags::Dict                  # Dict of fragments of form frag => expression in θs
end


@variables θ[1:30]

function SymbolicIntegrator(prob::Problem)
    frags = Dict()
    x = prob.x

    for term in terms(prob.eq)
        c, a = atomize(term, x)        
        frags[a] = c
    end

    for (i, b) in enumerate(prob.ker)
        db = expand_fraction(diff(b, x), x)

        for term in terms(db)
            c, a = atomize(term, x)                        
            frags[a] = get(frags, a, 0) + c * θ[i]
        end        
    end    
    
    return SymbolicIntegrator(
        [y ~ 0 for y in values(frags)], 
        [θ[i] for i=1:length(prob.ker)], 
        frags
    )
end


function solver(prob::Problem, alg::SymbolicIntegrator)       
    plan = prob.plan
    
    try
        sys = linearize(alg)
        q = solver(sys)
        sol = apply_coefs(q, prob.ker, prob.x)
        return sol        
    catch e
        # println(e) 
    end

    return nothing
end


# If the linear system generated by the SymbolicIntegrator solver is
# under-determined, make_square uses semi-numerical methods to 
# generate a full-rank linear system.
# 
# The output is another SymbolicIntegrator.
function make_square(prob::Problem, alg::SymbolicIntegrator)
    n = length(alg.vars)
    eqs = []
    S = subs_symbols(prob)

    for i in 1:n
        S[prob.x] = Complex(randn())
        q = sum(substitute(k, S) * v for (k, v) in alg.frags)
        Q = (q ~ 0)
        if Q isa Array
            # Q returns a complex array
            # a different pathway is needed here!
            return nothing
        else
            push!(eqs, Q)
        end
    end
    
    SymbolicIntegrator(
        eqs, 
        alg.vars, 
        alg.frags
    )
end


struct NumericIntegrator
    A::AbstractMatrix           # training dataset as a normalized linear system 
    X::AbstractVector           # vector of test points
    V::AbstractMatrix           # verification dataset
    basis::Vector{Expression}   # expand basis 
end


function NumericIntegrator(prob::Problem; nv=1)
    eq, x, ker, plan = value(prob.eq), prob.x, prob.ker, prob.plan    
    params = sym_consts(eq, x)
    
    basis = []
    for p in [[1]; params]
        for y in ker            
            push!(basis, y*p)
        end
    end
    basis = unique(basis)    

    Δbasis = [diff(y, x) for y in basis]    
    n = length(basis)
    
    A = zeros(Complex{Float64}, (n+nv, n))    
    X = zeros(Complex{Float64}, n+nv)
    
    for i = 1:n+nv
        S = subs_symbols(eq, x; include_x = true, plan.radius)        
        X[i] = S[x]        
        b₀ = complex_from_num(substitute(eq, S))

        for j in 1:n            
            A[i, j] = complex_from_num(substitute(Δbasis[j], S)) / b₀
        end        
    end
    
    return NumericIntegrator(
        A[1:n,:], 
        X[1:n], 
        A[n+1:end,:], 
        basis
    )
end

function solver(prob::Problem, alg::NumericIntegrator)
    sol, _ = solve_sparse(prob.eq, prob.x, alg.basis; prob.plan, AX=(alg.A, alg.X, alg.V))
    return sol
end


struct LinearSystem
    # linear system is Ax = b
    A::AbstractMatrix   
    b::AbstractVector
end


function linearize(alg::SymbolicIntegrator)
    eqs, vars = alg.eqs, alg.vars

    n = length(eqs)
    m = length(vars)
    
    L = [eq.lhs for eq in eqs]
    S = Dict(v => 0 for v in vars)

    b = Array{Num, 1}(undef, n)
    for i in 1:n
        b[i] = substitute(L[i], S)
    end

    A = Array{Num, 2}(undef, (n, m))

    for j in 1:m
        S = Dict(v => (i == j ? 1 : 0) for (i, v) in enumerate(vars))
        for i in 1:n
            A[i, j] = substitute(L[i], S) - b[i]
        end
    end

    return LinearSystem(
        A,
        b
    )
end


function solver(sys::LinearSystem)
    q = sys.A \ sys.b
    return value.(q)    
end

