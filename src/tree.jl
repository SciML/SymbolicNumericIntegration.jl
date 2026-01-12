# utility functions inspired from sympy

function is_add(eq)
    y = value(eq)
    return iscall(y) && operation(y) === (+)
end

function is_mul(eq)
    y = value(eq)
    return iscall(y) && operation(y) === (*)
end

function is_pow(eq)
    y = value(eq)
    return iscall(y) && operation(y) === (^)
end

function is_div(eq)
    y = value(eq)
    return iscall(y) && operation(y) === (/)
end

is_term(eq) = iscall(value(eq)) && !issym(value(eq))
# is_sym(eq) = SymbolicUtils.issym(value(eq))

function args(eq)
    eq = value(eq)
    return iscall(eq) ? arguments(eq) : []
end

diff(eq, x) = expand_derivatives(Differential(x)(eq))

numer(eq) = args(eq)[1]
denom(eq) = args(eq)[2]

function is_polynomial(p, x)
    if is_add(p) || is_mul(p)
        return all(is_polynomial(t, x) for t in args(p))
    elseif is_pow(p)
        q = args(p)[1]
        k = args(p)[2]
        return is_polynomial(q, x) && k isa Integer && k > 0
    else
        return !isdependent(p, x) || isequal(p, x)
    end
end

# Checks whether p is a univariate polynomial
function is_univar_poly(p)
    p = value(p)
    vars = get_variables(p)
    return length(vars) == 1 && is_polynomial(p, first(vars))
end

# Expands (Σ Ai) / B to Σ(Ai / B)
function expand_fraction(eq, x)
    if is_add(eq)
        return sum(expand_fraction(t, x) for t in args(eq))
    elseif is_div(eq)
        n = numer(eq)
        d = denom(eq)

        if is_add(n) && !isdependent(d, x)
            return sum(simplify(equiv(u, x) / d) for u in args(n))
        else
            return eq
        end
    else
        return eq
    end
end

# Returns the list of symbolic constants
function sym_consts(eq, x)
    return [v for v in get_variables(eq) if !isequal(v, x)]
end
