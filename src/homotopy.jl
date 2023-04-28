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

function transformer(eq::SymbolicUtils.Add, f)
    return sum(transformer(t, f) for t in arguments(eq); init = 0)
end
function transformer(eq::SymbolicUtils.Mul, f)
    return prod(transformer(t, f) for t in arguments(eq); init = 1)
end
function transformer(eq::SymbolicUtils.Div, f)
    return transformer(arguments(eq)[1], f) * transformer(arguments(eq)[2]^-1, f)
end

function transformer(eq::SymbolicUtils.Pow, f)
    y, k = arguments(eq)

    if is_pos_int(k)
        μ = next_variable!(f, y)
        return μ^k
    elseif is_neg_int(k)
        μ = next_variable!(f, inv(y))
        return μ^-k
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

Symbolics.@register_symbolic Ei(z)
Symbolics.@register_symbolic Si(z)
Symbolics.@register_symbolic Ci(z)
Symbolics.@register_symbolic Li(z)

Symbolics.derivative(::typeof(Ei), args::NTuple{1, Any}, ::Val{1}) = exp(args[1]) / args[1]
Symbolics.derivative(::typeof(Si), args::NTuple{1, Any}, ::Val{1}) = sin(args[1]) / args[1]
Symbolics.derivative(::typeof(Ci), args::NTuple{1, Any}, ::Val{1}) = cos(args[1]) / args[1]
Symbolics.derivative(::typeof(Li), args::NTuple{1, Any}, ::Val{1}) = 1 / log(args[1])

@syms si(𝑥) ci(𝑥) ei(𝑥) li(𝑥)

##############################################################################

function substitute_x(eq, x, sub)
    eq = substitute(eq, sub)
    substitute(eq, Dict(𝑥 => x))
end

function generate_homotopy(eq, x)
    eq = eq isa Num ? eq.val : eq
    x = x isa Num ? x.val : x

    q, sub = transform(eq, x)
    S = 0

    for i in 1:length(sub)
        μ = u[i]
        h₁, ∂h₁ = apply_partial_int_rules(sub[μ])
        h₁ = substitute(h₁, Dict(si => Si, ci => Ci, ei => Ei, li => Li))
        h₂ = expand_derivatives(Differential(μ)(q))

        h₁ = substitute_x(h₁, x, sub)
        h₂ = substitute_x(h₂ * ∂h₁^-1, x, sub)

        S += expand((1 + h₁) * (1 + h₂))
    end

    unique([one(x); [equivalent(t, x) for t in terms(S)]])
end

##############################################################################

function ∂(x)
    d = expand_derivatives(Differential(𝑥)(x))
    return isequal(d, 0) ? 1 : d
end

partial_int_rules = [
                     # trigonometric functions
                     @rule 𝛷(sin(~x)) => (cos(~x) + si(~x), ∂(~x))
                     @rule 𝛷(cos(~x)) => (sin(~x) + ci(~x), ∂(~x))
                     @rule 𝛷(tan(~x)) => (log(cos(~x)), ∂(~x))
                     @rule 𝛷(csc(~x)) => (log(csc(~x) + cot(~x)), ∂(~x))
                     @rule 𝛷(sec(~x)) => (log(sec(~x) + tan(~x)), ∂(~x))
                     @rule 𝛷(cot(~x)) => (log(sin(~x)), ∂(~x))
                     # hyperbolic functions
                     @rule 𝛷(sinh(~x)) => (cosh(~x), ∂(~x))
                     @rule 𝛷(cosh(~x)) => (sinh(~x), ∂(~x))
                     @rule 𝛷(tanh(~x)) => (log(cosh(~x)), ∂(~x))
                     @rule 𝛷(csch(~x)) => (log(tanh(~x / 2)), ∂(~x))
                     @rule 𝛷(sech(~x)) => (atan(sinh(~x)), ∂(~x))
                     @rule 𝛷(coth(~x)) => (log(sinh(~x)), ∂(~x))
                     # 1/trigonometric functions
                     @rule 𝛷(^(sin(~x), -1)) => (log(csc(~x) + cot(~x)), ∂(~x))
                     @rule 𝛷(^(cos(~x), -1)) => (log(sec(~x) + tan(~x)), ∂(~x))
                     @rule 𝛷(^(tan(~x), -1)) => (log(sin(~x)), ∂(~x))
                     @rule 𝛷(^(csc(~x), -1)) => (cos(~x), ∂(~x))
                     @rule 𝛷(^(sec(~x), -1)) => (sin(~x), ∂(~x))
                     @rule 𝛷(^(cot(~x), -1)) => (log(cos(~x)), ∂(~x))
                     # 1/hyperbolic functions
                     @rule 𝛷(^(sinh(~x), -1)) => (log(tanh(~x / 2)), ∂(~x))
                     @rule 𝛷(^(cosh(~x), -1)) => (atan(sinh(~x)), ∂(~x))
                     @rule 𝛷(^(tanh(~x), -1)) => (log(sinh(~x)), ∂(~x))
                     @rule 𝛷(^(csch(~x), -1)) => (cosh(~x), ∂(~x))
                     @rule 𝛷(^(sech(~x), -1)) => (sinh(~x), ∂(~x))
                     @rule 𝛷(^(coth(~x), -1)) => (log(cosh(~x)), ∂(~x))
                     # inverse trigonometric functions
                     @rule 𝛷(asin(~x)) => (~x * asin(~x) + sqrt(1 - ~x * ~x), ∂(~x))
                     @rule 𝛷(acos(~x)) => (~x * acos(~x) + sqrt(1 - ~x * ~x), ∂(~x))
                     @rule 𝛷(atan(~x)) => (~x * atan(~x) + log(~x * ~x + 1), ∂(~x))
                     @rule 𝛷(acsc(~x)) => (~x * acsc(~x) + atanh(1 - ^(~x, -2)), ∂(~x))
                     @rule 𝛷(asec(~x)) => (~x * asec(~x) + acosh(~x), ∂(~x))
                     @rule 𝛷(acot(~x)) => (~x * acot(~x) + log(~x * ~x + 1), ∂(~x))
                     # inverse hyperbolic functions
                     @rule 𝛷(asinh(~x)) => (~x * asinh(~x) + sqrt(~x * ~x + 1), ∂(~x))
                     @rule 𝛷(acosh(~x)) => (~x * acosh(~x) + sqrt(~x * ~x - 1), ∂(~x))
                     @rule 𝛷(atanh(~x)) => (~x * atanh(~x) + log(~x + 1), ∂(~x))
                     @rule 𝛷(acsch(~x)) => (acsch(~x), ∂(~x))
                     @rule 𝛷(asech(~x)) => (asech(~x), ∂(~x))
                     @rule 𝛷(acoth(~x)) => (~x * acot(~x) + log(~x + 1), ∂(~x))
                     # logarithmic and exponential functions
                     @rule 𝛷(log(~x)) => (~x + ~x * log(~x) +
                                          sum(candidate_pow_minus(~x, -1); init = one(~x)),
                                          ∂(~x))
                     @rule 𝛷(^(log(~x), -1)) => (log(log(~x)) + li(~x), ∂(~x))
                     @rule 𝛷(exp(~x)) => (exp(~x) + ei(~x), ∂(~x))
                     @rule 𝛷(^(exp(~x), ~k::is_neg)) => (^(exp(-~x), -~k), ∂(~x))
                     # square-root functions
                     @rule 𝛷(^(~x, ~k::is_abs_half)) => (sum(candidate_sqrt(~x, ~k);
                                                             init = one(~x)), 1);
                     @rule 𝛷(sqrt(~x)) => (sum(candidate_sqrt(~x, 0.5); init = one(~x)), 1);
                     @rule 𝛷(^(sqrt(~x), -1)) => 𝛷(^(~x, -0.5))
                     # rational functions                                                              
                     @rule 𝛷(^(~x::is_poly, ~k::is_neg)) => (sum(candidate_pow_minus(~x,
                                                                                     ~k);
                                                                 init = one(~x)), 1)
                     @rule 𝛷(^(~x, -1)) => (log(~x), ∂(~x))
                     @rule 𝛷(^(~x, ~k::is_neg_int)) => (sum(^(~x, i) for i in (~k + 1):-1),
                                                        ∂(~x))
                     @rule 𝛷(1 / ~x) => 𝛷(^(~x, -1))
                     @rule 𝛷(^(~x, ~k)) => (^(~x, ~k + 1), ∂(~x))
                     @rule 𝛷(1) => (𝑥, 1)
                     @rule 𝛷(~x) => ((~x + ^(~x, 2)), ∂(~x))]

function apply_partial_int_rules(eq)
    expand(Fixpoint(Prewalk(Chain(partial_int_rules))))(𝛷(value(eq)))
end
