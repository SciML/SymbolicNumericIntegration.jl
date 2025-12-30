transformer(eq) = transformer(ops(eq)...)

function transformer(::Mul, eq)
    return vcat([transformer(t) for t in arguments(eq)]...)
end

function transformer(::Div, eq)
    a = transformer(arguments(eq)[1])
    b = transformer(arguments(eq)[2])
    b = [(1 / q, k) for (q, k) in b]
    return [a; b]
end

function transformer(::Pow, eq)
    y, k = arguments(eq)
    if is_number(k)
        r = nice_parameter(k)
        if r isa Integer || r isa Rational
            if denominator(r) == 1
                return [(y, k)]
            else
                return [(y^(1 / denominator(r)), numerator(r))]
            end
        end
    end
    return [(eq, 1)]
end

function transformer(::Any, eq)
    return [(eq, 1)]
end

function transform(eq, x)
    p = transformer(eq)
    p = p[isdependent.(first.(p), x)]
    return p
end

const u = let
    @variables _u[1:20]
    Symbolics.scalarize(_u)
end

function rename_factors(p, ab = ())
    n = length(p)
    q = 1
    ks = Int[]
    sub = Dict()

    for (a, b) in ab
        sub[a] = b
    end

    for (i, (y, k)) in enumerate(p)
        μ = u[i]
        q *= μ^k
        sub[μ] = y
        push!(ks, k)
    end

    return q, sub, ks
end

##############################################################################

Symbolics.@register_symbolic Ei(z)
Symbolics.@register_symbolic Si(z)
Symbolics.@register_symbolic Ci(z)
Symbolics.@register_symbolic Li(z)
Symbolics.@register_symbolic Erfi(z)

Symbolics.derivative(::typeof(Ei), args::NTuple{1, Any}, ::Val{1}) = exp(args[1]) / args[1]
Symbolics.derivative(::typeof(Si), args::NTuple{1, Any}, ::Val{1}) = sin(args[1]) / args[1]
Symbolics.derivative(::typeof(Ci), args::NTuple{1, Any}, ::Val{1}) = cos(args[1]) / args[1]
Symbolics.derivative(::typeof(Li), args::NTuple{1, Any}, ::Val{1}) = 1 / log(args[1])
function Symbolics.derivative(::typeof(Erfi), args::NTuple{1, Any}, ::Val{1})
    2 / sqrt(2) * exp(args[1]^2)
end

@syms 𝑥 si(𝑥) ci(𝑥) ei(𝑥) li(𝑥) erfi_(𝑥)

##############################################################################

guard_zero(x) = isequal(x, 0) ? one(x) : x

# The core of ansatz generation.
# The name is not accurate and should be changed (not really homotopy, 
# just inspired from!)
function generate_homotopy(eq, x)
    eq = value(eq)
    x = value(x)

    if is_add(eq)
        return unique(union([generate_homotopy(t, x) for t in args(eq)]...))
    end

    p = transform(eq, x)
    q, sub, ks = rename_factors(p, (si => Si, ci => Ci, ei => Ei, li => Li, erfi_ => Erfi))
    S = 0

    for i in 1:length(ks)
        μ = u[i]
        y, dy = apply_partial_int_rules(sub[μ], x)

        y = substitute(y, sub)
        ∂y = guard_zero(diff(dy, x))

        for j in 1:ks[i]
            h = substitute((q / μ^j) / ∂y, sub)
            S += expand((ω + y) * (ω + h))
        end
    end

    S = substitute(S, Dict(ω => 1))
    unique([x; [equivalent(t, x) for t in terms(S)]])
end

########################## Main Integration Rules ##################################

@syms 𝛷(x, w)

partial_int_rules = [
                     # trigonometric functions
                     @rule 𝛷(~x, sin(~u)) => (cos(~u) + si(~u), ~u)
                     @rule 𝛷(~x, cos(~u)) => (sin(~u) + ci(~u), ~u)
                     @rule 𝛷(~x, tan(~u)) => (log(cos(~u)), ~u)
                     @rule 𝛷(~x, csc(~u)) => (log(csc(~u) + cot(~u)) + log(sin(~u)), ~u)
                     @rule 𝛷(~x, sec(~u)) => (log(sec(~u) + tan(~u)) + log(cos(~u)), ~u)
                     @rule 𝛷(~x, cot(~u)) => (log(sin(~u)), ~u)
                     # hyperbolic functions
                     @rule 𝛷(~x, sinh(~u)) => (cosh(~u), ~u)
                     @rule 𝛷(~x, cosh(~u)) => (sinh(~u), ~u)
                     @rule 𝛷(~x, tanh(~u)) => (log(cosh(~u)), ~u)
                     @rule 𝛷(~x, csch(~u)) => (log(tanh(~u / 2)), ~u)
                     @rule 𝛷(~x, sech(~u)) => (atan(sinh(~u)), ~u)
                     @rule 𝛷(~x, coth(~u)) => (log(sinh(~u)), ~u)
                     # 1/trigonometric functions
                     @rule 𝛷(~x, 1 /
                                 sin(~u)) => (log(csc(~u) + cot(~u)) + log(sin(~u)), ~u)
                     @rule 𝛷(~x, 1 /
                                 cos(~u)) => (log(sec(~u) + tan(~u)) + log(cos(~u)), ~u)
                     @rule 𝛷(~x, 1 / tan(~u)) => (log(sin(~u)) + log(tan(~u)), ~u)
                     @rule 𝛷(~x, 1 / csc(~u)) => (cos(~u) + log(csc(~u)), ~u)
                     @rule 𝛷(~x, 1 / sec(~u)) => (sin(~u) + log(sec(~u)), ~u)
                     @rule 𝛷(~x, 1 / cot(~u)) => (log(cos(~u)) + log(cot(~u)), ~u)
                     # 1/hyperbolic functions
                     @rule 𝛷(~x, 1 / sinh(~u)) => (log(tanh(~u / 2)) + log(sinh(~u)), ~u)
                     @rule 𝛷(~x, 1 / cosh(~u)) => (atan(sinh(~u)) + log(cosh(~u)), ~u)
                     @rule 𝛷(~x, 1 / tanh(~u)) => (log(sinh(~u)) + log(tanh(~u)), ~u)
                     @rule 𝛷(~x, 1 / csch(~u)) => (cosh(~u) + log(csch(~u)), ~u)
                     @rule 𝛷(~x, 1 / sech(~u)) => (sinh(~u) + log(sech(~u)), ~u)
                     @rule 𝛷(~x, 1 / coth(~u)) => (log(cosh(~u)) + log(coth(~u)), ~u)
                     # inverse trigonometric functions
                     @rule 𝛷(~x, asin(~u)) => (~u * asin(~u) + sqrt(1 - ~u * ~u), ~u)
                     @rule 𝛷(~x, acos(~u)) => (~u * acos(~u) + sqrt(1 - ~u * ~u), ~u)
                     @rule 𝛷(~x, atan(~u)) => (~u * atan(~u) + log(~u * ~u + 1), ~u)
                     @rule 𝛷(~x, acsc(~u)) => (~u * acsc(~u) + atanh(1 - ^(~u, -2)), ~u)
                     @rule 𝛷(~x, asec(~u)) => (~u * asec(~u) + acosh(~u), ~u)
                     @rule 𝛷(~x, acot(~u)) => (~u * acot(~u) + log(~u * ~u + 1), ~u)
                     # inverse hyperbolic functions
                     @rule 𝛷(~x, asinh(~u)) => (~u * asinh(~u) + sqrt(~u * ~u + 1), ~u)
                     @rule 𝛷(~x, acosh(~u)) => (~u * acosh(~u) + sqrt(~u * ~u - 1), ~u)
                     @rule 𝛷(~x, atanh(~u)) => (~u * atanh(~u) + log(~u + 1), ~u)
                     @rule 𝛷(~x, acsch(~u)) => (acsch(~u), ~u)
                     @rule 𝛷(~x, asech(~u)) => (asech(~u), ~u)
                     @rule 𝛷(~x, acoth(~u)) => (~u * acot(~u) + log(~u + 1), ~u)
                     # logarithmic and exponential functions
                     @rule 𝛷(~x,
                         log(~u)) => (
                         ~u + ~u * log(~u) +
                         sum(pow_minus_rule(~u, ~x, -1); init = one(~u)),
                         ~u)
                     @rule 𝛷(~x, 1 / log(~u)) => (log(log(~u)) + li(~u), ~u)
                     @rule 𝛷(~x, exp(~u)) => (exp(~u) + ei(~u) + erfi_(~x), ~u)
                     @rule 𝛷(~x, ^(exp(~u), ~k::is_neg)) => (^(exp(-~u), -~k), ~u)
                     # square-root functions
                     @rule 𝛷(~x, ^(~u, ~k::is_abs_half)) => (
                         sum(sqrt_rule(~u, ~x, ~k); init = one(~u)), ~u)
                     @rule 𝛷(~x, sqrt(~u)) => (
                         sum(sqrt_rule(~u, ~x, 0.5); init = one(~u)), ~u)
                     @rule 𝛷(~x, 1 /
                                 sqrt(~u)) => (
                         sum(sqrt_rule(~u, ~x, -0.5); init = one(~u)), ~u)
                     # rational functions                                                              
                     @rule 𝛷(~x,
                         1 / ^(~u::is_univar_poly, ~k::is_pos_int)) => (
                         sum(pow_minus_rule(~u,
                                 ~x,
                                 -~k);
                             init = one(~u)),
                         ~u)
                     @rule 𝛷(
                         ~x, 1 / ~u::is_univar_poly) => (
                         sum(pow_minus_rule(~u, ~x, -1); init = one(~u)),
                         ~u)
                     @rule 𝛷(~x, ^(~u, -1)) => (log(~u) + ~u * log(~u), ~u)
                     @rule 𝛷(~x, ^(~u, ~k::is_neg_int)) => (
                         sum(^(~u, i) for i in (~k + 1):-1), ~u)
                     @rule 𝛷(~x, 1 / ~u) => (log(~u), ~u)
                     @rule 𝛷(~x, ^(~u, ~k::is_pos_int)) => (
                         sum(^(~u, i + 1) for i in 1:(~k + 1)), ~u)
                     @rule 𝛷(~x, 1) => (𝑥, 1)
                     @rule 𝛷(~x, ~u) => ((~u + ^(~u, 2)), ~u)]

function apply_partial_int_rules(eq, x)
    y, dy = Chain(partial_int_rules)(𝛷(x, value(eq)))
    return y, dy
end

################################################################

function pow_minus_rule(p, x, k; abstol = 1e-8)
    if !is_univar_poly(p)
        return [p^k, p^(k + 1), log(p), p * log(p)]
    end

    # x = var(p)
    d = poly_deg(p)

    for j in 1:10  # will try 10 times to find the roots
        r, s = find_roots(p, x)
        if length(r) + length(s) >= d
            break
        end
    end
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    # applying ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
    q = Any[log(x - u) for u in r]
    for i in eachindex(s)
        β = s[i]
        if abs(imag(β)) > abstol
            push!(q, atan((x - real(β)) / imag(β)))
            push!(q, (1 + x) * log(x^2 - 2 * real(β) * x + abs2(β)))
        else
            push!(q, log(x - real(β)))
        end
    end
    q = unique(q)

    if k ≈ -1
        return [[p^k]; q]
    else
        return [[p^k, p^(k + 1)]; q]
    end
end

function sqrt_rule(p, x, k)
    h = Any[p ^ k, p ^ (k + 1)]

    Δ = diff(p, x)
    push!(h, log(Δ / 2 + sqrt(p)))

    if !is_univar_poly(p)
        return h
    end

    # x = var(p)

    if poly_deg(p) == 2
        r, s = find_roots(p, x)
        l = leading(p, x)
        if length(r) == 2
            if sum(r) ≈ 0
                r₁ = abs(r[1])
                if l > 0
                    push!(h, acosh(x / r₁))
                else
                    push!(h, asin(x / r₁))
                end
            end
        elseif real(s[1]) ≈ 0
            push!(h, asinh(x / imag.(s[1])))
        end
    end

    return h
end
