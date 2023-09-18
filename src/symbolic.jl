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
        return beautify(p) ^ k
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
    
    if is_polynomial(y, x)
        return false
    end
    
    if is_number(y) 
        return true
    end

    if is_term(y)
        return operation(y) in [sin, cos, sinh, cosh, exp]
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
function blender(y, x; n=3)
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

    return split_terms(basis, x)
end

#################################################################################

"""
    Returns a dictionary of the symbolic constants in the expression
    and random real value assignments.

    Args:
        eq: the integrands
        x: independent variable

    Returns:
        a dict of sym : value pairs
"""
function subs_symbols(eq, x; include_x=false, radius=5.0)
    S = Dict()
    for v in get_variables(value(eq))
        if !isequal(v, x)
            S[v] = randn()
        end
    end
    
    if include_x
        S[x] = randn()
    end

    return S
end


"""
    Splits terms into a part dependent on x and a part constant w.r.t. x    
    For example, `atomize(a*sin(b*x), x)` is `(a, sin(b*x))`)
"""
function atomize(eq, x)
    eq = value(eq)
    
    if is_mul(eq)
        coef = 1
        atom = 1
        for t in args(eq)
            c, a = atomize(t, x)
            coef *= c
            atom *= a
        end
        return (coef, atom)            
    elseif is_div(eq)
        c1, a1 = atomize(numer(eq), x)
        c2, a2 = atomize(denom(eq), x)
        return (c1/c2, a1/a2)
    elseif isdependent(eq, x)
        return (1, eq)
    else
        return (eq, 1)
    end   
end


function atomize(eq)
    eq = value(eq)
    
    if is_mul(eq)
        coef = 1
        atom = 1
        for t in args(eq)
            c, a = atomize(t)
            coef *= c
            atom *= a
        end
        return (coef, atom)            
    elseif is_div(eq)
        c1, a1 = atomize(numer(eq))
        c2, a2 = atomize(denom(eq))
        return (c1/c2, a1/a2)
    elseif is_number(eq)
        return (eq, 1)
    else
        return (1, eq)
    end   
end



############################## Numerical Utils ##############################

@variables θ[1:30]


"""
    Given an integral problem ∫eq dx = Σ θ[i]*basis[i], where 
    θs are coefficients and basis is a list of ansatzes, make_eqs
    generates a list of (hopefully linear!) equations in θs.
"""
function make_eqs(eq, x, basis)
    frags = Dict()
    
    for d in terms(eq)
        c, a = atomize(d, x)
        frags[a] = c
    end
    
    vars = []
    
    for (i, b) in enumerate(basis)
        db = expand_fraction(diff(b, x), x)

        for d in terms(db)
            c, a = atomize(d, x)
            frags[a] = get(frags, a, 0) + c * θ[i]
        end
        
        push!(vars, θ[i])
    end
  
    return [y ~ 0 for y in values(frags)], vars, frags
end


"""
    make_Ab transforms the output of make_eqs to a linear
    system.
    
    Args:
    --------
        eqs: a list of n linear equations in θs
        vars: the list of variables, i.e., θ[1], θ[2],...
        
    Returns:
    --------
        A: an n-by-n symbolic matrix
        b: a symbolic venctor of length n
"""
function make_Ab(eqs, vars)
    n = length(vars)
    L = [eq.lhs for eq in eqs]
    S = Dict(v => 0 for v in vars)
    
    b = Array{Num,1}(undef, n)
    for i = 1:n
        b[i] = substitute(L[i], S)
    end

    A = Array{Num,2}(undef, (n, n))

    for j = 1:n
        S = Dict(v => (i==j ? 1 : 0) for (i, v) in enumerate(vars))
        for i = 1:n
            A[i, j] = substitute(L[i], S) - b[i]
        end
    end

    return A, b
end


"""
    If the output of make_eqs is under-determined, make_square employs
    numerical methods to generate a solvable system. 
"""
function make_square(eq, x, vars, frags)
    n = length(vars)
    eqs = []
    
    for i = 1:n
        S = subs_symbols(eq, x; include_x=true)
        q = sum(substitute(k, S) * v for (k,v) in frags)
        push!(eqs, q ~ 0)
    end
    
    return eqs 
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
    s = simplify(s)
    c, a = atomize(s)
    return beautify(c) * a
end


"""
    The main entry point for symbolic integration.
    
    Argu:
        eq: the expression to integrate
        x: independent variable
        
    Returns:
        the integral or nothing if no solution
"""
function integrate_symbolic(eq, x; abstol=1e-6, radius=1.0)    
    eq = expand(eq)
    
    if is_holonomic(eq, x)
        basis = blender(eq, x)
    else
        basis = generate_homotopy(eq, x)        
    end

    ker = best_hints(eq, x, basis)
    
    if ker == nothing
        return nothing
    end
    
    ker = [atomize(y, x)[2] for y in ker]
    
    eqs, vars, frags = make_eqs(eq, x, ker)
    
    try
        if length(eqs) != length(vars)
            eqs = make_square(eq, x, vars, frags)
        end
        A, b = make_Ab(eqs, vars)        
        q = A \ b
        q = value.(q)        
        sol = apply_coefs(q, ker, x)

        # test if sol is the correct integral of eq
        S = subs_symbols(eq, x; include_x=true, radius)
        err = substitute(diff(sol, x) - eq, S)
        
        if abs(err) < abstol
            return sol
        else
            return nothing
        end
    catch e
        # println(e)
    end    

    return nothing
end





