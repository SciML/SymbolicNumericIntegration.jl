@syms 𝑥
@syms u[20]

mutable struct Transform
    k::Int
    sub::Dict
end

function next_variable!(f, eq)
    μ = u[f.k]
    f.k += 1
    f.sub[μ] = eq
    return μ
end

transformer(eq::SymbolicUtils.Add, f) = sum(transformer(t,f) for t in arguments(eq); init=0)
transformer(eq::SymbolicUtils.Mul, f) = prod(transformer(t,f) for t in arguments(eq); init=1)
transformer(eq::SymbolicUtils.Div, f) = transformer(arguments(eq)[1],f) * transformer(inv(arguments(eq)[2]),f)

function transformer(eq::SymbolicUtils.Pow, f)
    y, k = arguments(eq)

    if is_pos_int(k)
        μ = next_variable!(f, y)
        return μ ^ k
    elseif is_neg_int(k)
        μ = next_variable!(f, inv(y))
        return μ ^ -k
    else
        return next_variable!(f, y^k)
    end
end

function transformer(eq, f)
    if isdependent(eq, 𝑥)
        return next_variable!(f, eq)
    else
        return 1
    end
end

function transform(eq, x)
    eq = substitute(eq, Dict(x => 𝑥))
    f = Transform(1, Dict())
    q = transformer(eq, f)
    if !any(is_poly, values(f.sub))
        q *= next_variable!(f, 1)
    end
    return q, f.sub
end

##############################################################################

function substitute_x(eq, x, sub)
    eq = substitute(eq, sub)
    substitute(eq, Dict(𝑥 => x))
end

function generate_homotopy(eq, x)
    q, sub = transform(eq, x)    
    S = 0

    for i = 1:length(sub)
        μ = u[i]
        h₁, ∂h₁ = apply_partial_int_rules(sub[μ])
        h₂ = expand_derivatives(Differential(μ)(q))

        h₁ = substitute_x(h₁, x, sub)
        h₂ = substitute_x(h₂ * ∂h₁^-1, x, sub)

        S += expand((1 + h₁) * (1 + h₂))
    end

    unique([one(x); [equivalent(t,x) for t in terms(S)]])
end

##############################################################################

function ∂(x)
    d = expand_derivatives(Differential(𝑥)(x))
    return isequal(d, 0) ? 1 : d
end

partial_int_rules = [
    @rule 𝛷(^(sin(~x), ~k::is_neg)) => 𝛷(^(csc(~x), -~k))
    @rule 𝛷(^(cos(~x), ~k::is_neg)) => 𝛷(^(sec(~x), -~k))
    @rule 𝛷(^(tan(~x), ~k::is_neg)) => 𝛷(^(cot(~x), -~k))
    @rule 𝛷(^(csc(~x), ~k::is_neg)) => 𝛷(^(sin(~x), -~k))
    @rule 𝛷(^(sec(~x), ~k::is_neg)) => 𝛷(^(cos(~x), -~k))
    @rule 𝛷(^(cot(~x), ~k::is_neg)) => 𝛷(^(tan(~x), -~k))

    @rule 𝛷(^(sinh(~x), ~k::is_neg)) => 𝛷(^(csch(~x), -~k))
    @rule 𝛷(^(cosh(~x), ~k::is_neg)) => 𝛷(^(sech(~x), -~k))
    @rule 𝛷(^(tanh(~x), ~k::is_neg)) => 𝛷(^(coth(~x), -~k))
    @rule 𝛷(^(csch(~x), ~k::is_neg)) => 𝛷(^(sinh(~x), -~k))
    @rule 𝛷(^(sech(~x), ~k::is_neg)) => 𝛷(^(cosh(~x), -~k))
    @rule 𝛷(^(coth(~x), ~k::is_neg)) => 𝛷(^(tanh(~x), -~k))

    @rule 𝛷(sin(~x)) => (cos(~x), ∂(~x))
    @rule 𝛷(cos(~x)) => (sin(~x), ∂(~x))
    @rule 𝛷(tan(~x)) => (log(cos(~x)), ∂(~x))
    @rule 𝛷(csc(~x)) => (log(sin(~x)^-1 - cos(~x)*sin(~x)^-1), ∂(~x))
    @rule 𝛷(sec(~x)) => (log(cos(~x)^-1 + sin(~x)*cos(~x)^-1), ∂(~x))
    @rule 𝛷(cot(~x)) => (log(sin(~x)), ∂(~x))

    @rule 𝛷(sinh(~x)) => (cosh(~x), ∂(~x))
    @rule 𝛷(cosh(~x)) => (sinh(~x), ∂(~x))
    @rule 𝛷(tanh(~x)) => (log(cosh(~x)), ∂(~x))
    @rule 𝛷(csch(~x)) => (log(sinh(~x)^-1 + cosh(~x)*sinh(~x)^-1), ∂(~x))
    @rule 𝛷(sech(~x)) => (log(cosh(~x)^-1 + sinh(~x)*cosh(~x)^-1), ∂(~x))
    @rule 𝛷(coth(~x)) => (log(sinh(~x)), ∂(~x))

    @rule 𝛷(asin(~x)) => (~x*asin(~x) + sqrt(1 - ~x*~x), ∂(~x))
    @rule 𝛷(acos(~x)) => (~x*acos(~x) + sqrt(1 - ~x*~x), ∂(~x))
    @rule 𝛷(atan(~x)) => (~x*atan(~x) + log(~x*~x + 1), ∂(~x))
    @rule 𝛷(acsc(~x)) => (~x*acsc(~x) + acosh(~x), ∂(~x))     # needs an abs inside acosh
    @rule 𝛷(asec(~x)) => (~x*asec(~x) + acosh(~x), ∂(~x))     # needs an abs inside acosh
    @rule 𝛷(acot(~x)) => (~x*acot(~x) + log(~x*~x + 1), ∂(~x))

    @rule 𝛷(asinh(~x)) => (~x*asinh(~x) + sqrt(~x*~x + 1), ∂(~x))
    @rule 𝛷(acosh(~x)) => (~x*acosh(~x) + sqrt(~x*~x - 1), ∂(~x))
    @rule 𝛷(atanh(~x)) => (~x*atanh(~x) + log(~x + 1), ∂(~x))
    @rule 𝛷(acsch(~x)) => (acsch(~x), ∂(~x))
    @rule 𝛷(asech(~x)) => (asech(~x), ∂(~x))
    @rule 𝛷(acoth(~x)) => (~x*acot(~x) + log(~x + 1), ∂(~x))

    @rule 𝛷(log(~x)) => (~x + ~x * log(~x), ∂(~x))

    @rule 𝛷(^(~x, ~k::is_abs_half)) => (sum(candidate_sqrt(~x,~k); init=one(~x)), 1)
    @rule 𝛷(^(~x::is_poly, ~k::is_neg)) => (sum(candidate_pow_minus(~x, ~k); init=one(~x)), 1)
    @rule 𝛷(sqrt(~x)) => (sum(candidate_sqrt(~x,0.5); init=one(~x)), 1)
    @rule 𝛷(^(sqrt(~x),-1)) => 𝛷(^(~x,-0.5))

    @rule 𝛷(^(~x, -1)) => (log(~x), ∂(~x))
    @rule 𝛷(1 / ~x) => 𝛷(^(~x, -1))
    @rule 𝛷(^(~x, ~k)) => (^(~x, ~k+1), ∂(~x))

    @rule 𝛷(exp(~x)) => (exp(~x), ∂(~x))
    @rule 𝛷(1) => (𝑥, 1)
    @rule 𝛷(~x) => ((~x + ^(~x,2)), ∂(~x))
]

apply_partial_int_rules(eq) = expand(Fixpoint(Prewalk(Chain(partial_int_rules))))(𝛷(value(eq)))
