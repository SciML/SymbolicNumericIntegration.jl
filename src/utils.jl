
abstract type Op end

struct Add <: Op end
struct Mul <: Op end
struct Div <: Op end
struct Pow <: Op end
struct Sym <: Op end
struct Term <: Op end

ops(eq) = nothing, eq
ops(eq::Num) = ops(eq.val)

function ops(eq::BasicSymbolic)
    op = exprtype(eq)

    if op == SymbolicUtils.ADD
        return Add(), eq
    elseif op == SymbolicUtils.MUL
        return Mul(), eq
    elseif op == SymbolicUtils.DIV
        return Div(), eq
    elseif op == SymbolicUtils.POW
        return Pow(), eq
    elseif op == SymbolicUtils.SYM
        return Sym(), eq
    elseif op == SymbolicUtils.TERM
        return Term(), eq
    else
        return nothing, eq
    end
end

############################################################################

"""
    isdependent returns true if eq is dependent on x
"""
#isdependent(eq, x) = !isequal(expand_derivatives(Differential(x)(eq)), 0)

function isdependent(eq, x)
    return any(isequal(x, v) for v in get_variables(eq))
    # vars = get_variables(eq)
    # return length(vars) == 1 && isequal(x, vars[1])
end

"""
    is_number(x) returns true if x is a concrete numerical type
"""
is_number(x::T) where {T <: Integer} = true
# is_number(x::T) where {T <: Real} = true
is_number(x::T) where {T <: Float32} = true
is_number(x::T) where {T <: Float64} = true
is_number(x::T) where {T <: Complex} = true
is_number(x::T) where {T <: Rational} = true
is_number(x) = false

is_proper(x) = is_number(x) && !isnan(x) && !isinf(x)

"""
    check_poly checks if a Symbolic expression is a polynomial
"""
# function check_poly(eq, x)
#     if deg(eq, x) == 0 return :not_poly end
#     coefs = values(collect_powers(value(eq), x))
#     if any(z -> isdependent(restore_x(z,x),x), coefs) return :not_poly end
#     if any(z -> z isa Complex, coefs) return :complex_poly end
#     return :real_poly
# end

# is_linear_poly(eq, x) = check_poly(eq, x) != :not_poly && deg(eq, x) == 1

function leading(eq, x)
    l = 0
    k₀ = -1

    for (k, coef) in collect_powers(value(eq), x)
        if !is_number(coef)
            return 0
        end
        if k > k₀
            k₀ = k
            l = coef
        end
    end
    l
end

"""
    deg(p) returns the degree of p if p is a polynomial
"""
deg(::Add, p, x) = maximum(deg(t, x) for t in arguments(p))
deg(::Mul, p, x) = sum(deg(t, x) for t in arguments(p); init = 0)

function deg(::Pow, p, x)
    if !isequal(arguments(p)[1], x)
        return 0
    end
    return arguments(p)[2]
end

deg(::Sym, p, x) = isequal(p, x) ? 1 : 0
deg(::Term, p, x) = 0
deg(::Any, p, x) = 0

deg(p, x) = deg(ops(p)..., x)

"""
    var(p) returns the unique variable of an expression (is exists)
"""
function var(p)
    vars = get_variables(p)
    if length(vars) == 1
        return vars[1]
    elseif length(vars) == 0
        return nothing
    else
        error("expects only one variable")
    end
end

###############################################################################

# pox (power-of-x) is a symbolic function to keep track of the powers of x
# pox(k,n) means k*x^n
@syms pox(k, n)

is_pox(x) = istree(x) && operation(x) == pox
is_not_pox(x) = !is_pox(x)

get_coef(p) = is_pox(p) ? arguments(p)[1] : p
get_power(p) = is_pox(p) ? arguments(p)[2] : 0

replace_x(eq, x) = substitute(eq, Dict(x => pox(1, 1)))

function restore_x(eq, x)
    r = @rule pox(~k, ~n) => ~k * x^~n
    Prewalk(PassThrough(r))(eq)
end

iscomplex(x) = x isa Complex

count_rule1 = @rule ^(pox(~k, ~n1), ~n2) => isequal(~k, 1) ? pox(1, ~n1 * ~n2) :
                                            pox(^(~k, ~n2), ~n1 * ~n2)
count_rule2 = @rule pox(~k1, ~n1) * pox(~k2, ~n2) => pox(~k1 * ~k2, ~n1 + ~n2)
count_rule3 = @acrule pox(~k, ~n) * ~u::is_not_pox => pox(~k * ~u, ~n)

"""
    collect_powers separates the powers of x in eq (a polynomial) and returns
    a dictionary of power => term
"""
function collect_powers(eq, x)
    eq = expand(expand_derivatives(eq))
    eq = replace_x(eq, x)
    #eq = Prewalk(PassThrough(count_rule1))(eq)
    eq = Fixpoint(Prewalk(PassThrough(Chain([count_rule1, count_rule2, count_rule3]))))(eq)

    if !istree(eq)
        return Dict{Any, Any}(0 => eq)
    elseif is_pox(eq)
        return Dict{Any, Any}(get_power(eq) => get_coef(eq))
    else
        eqs = Dict{Any, Any}()
        for term in arguments(eq)
            n = get_power(term)
            if haskey(eqs, n)
                eqs[n] = eqs[n] + get_coef(term)
            else
                eqs[n] = get_coef(term)
            end
        end

        return eqs
    end
end
