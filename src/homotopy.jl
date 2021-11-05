@syms 𝑥
@syms u[20]

mutable struct Transform
    k::Int
    sub::Dict
    deg::Int
    hasx::Bool
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

    # if is_poly(y)
    #     return next_variable!(f, y)^k
    # end

    # r = nice_parameter(k)
    # if r isa Rational || isinteger(r)
    if isinteger(k)
        a, b = k, 1
        # a, b = numerator(r), denominator(r)
        if k < 0
            y = inv(y)
        end
        f.deg = max(f.deg, abs(a))
        μ = next_variable!(f, b == 1 ?  y : y ^(1/b))
        return μ ^ abs(a)
    else
        return next_variable!(f, y^k)
    end
end

function transformer(eq, f)
    if isdependent(eq, 𝑥)
        f.hasx |= is_linear_poly(eq)
        return next_variable!(f, eq)
    else
        return 1
    end
end

function transform(eq, x)
    eq = substitute(eq, Dict(x => 𝑥))
    f = Transform(1, Dict(), 1, false)
    p = transformer(eq, f)
    if !any(is_poly, values(f.sub))
        p *= next_variable!(f, 1)
    end
    return p, f.sub, f.deg
end

function homotopy_integrand(eq, x)
    eq, sub, deg = transform(eq, x)
    I = (1 + x) * eq
    n = length(sub)

    for i = 1:n
        μ = u[i]
        H = apply_H_rules(sub[μ])
        I += H * expand_derivatives(Differential(μ)(I))
    end

    I = substitute(I, sub)
    I = substitute(I, Dict(𝑥 => x))
    return expand(I), deg
end

function expand_integrand(I, x, deg)
    # E = sum((Differential(x)^i)(I) for i=1:deg-1; init=I) #* (1+x)
    # S = Set{Any}()
    # enqueue_expr_ex!(S, expand(expand_derivatives(E)), x)
    # return [one(x); [s for s in S]]

    S = Set{Any}()
    # T = Set{Any}()
    Q₁ = Queue{Any}()

    enqueue_expr_ex!(S, Q₁, expand(I), x)

    D = Differential(x)

    for i = 1:deg
        Q₂ = Queue{Any}()
        while !isempty(Q₁) # && length(S) < max_terms
            y = dequeue!(Q₁)
            E = expand(expand_derivatives(D(y)))
            enqueue_expr_ex!(S, Q₂, E, x)
        end
        Q₁ = Q₂
    end

    return [one(x); [s for s in S]]
end

function expand_integrand(I, x, deg)
    E = sum((Differential(x)^i)(I) for i=1:deg-1; init=I) #* (1+x)
    S = Set{Any}()
    enqueue_expr_ex!(S, expand(expand_derivatives(E)), x)
    return [one(x); [s for s in S]]
end


function generate_homotopy(eq, x=var(eq))
    I, deg = homotopy_integrand(eq, x)
    expand_integrand(I, x, deg)
end

function substitute_x(eq, x, sub)
    eq = substitute(eq, sub)
    substitute(eq, Dict(𝑥 => x))
end

function generate_homotopy2(eq, x)
    q, sub = transform(eq, x)
    d = degree(q)
    n = length(sub)

    S = Set{Any}()

    for i = 1:n
        μ = u[i]
        h₁ = apply_H_rules(sub[μ])
        h₂ = expand_derivatives(Differential(μ)(q))

        h₁ = substitute_x(h₁, x, sub)
        h₂ = substitute_x(h₂, x, sub)

        H = sum((Differential(x)^i)(h₂) for i=1:d-1; init=(1 + h₂))
        I = expand(expand_derivatives((1 + h₁) * H))
        enqueue_expr_ex!(S, I, x)
    end

    # H = sum((Differential(x)^i)(eq) for i=1:deg-1; init=eq)
    # I = expand((1+x) * expand_derivatives(H))
    # enqueue_expr_ex!(S, I, x)

    # return [one(x); [s for s in S]]
    return [one(x); [s for s in S]]
end

##############################################################################

∂(x) = expand_derivatives(Differential(𝑥)(x))

H_rules = [
    # @rule 𝛷(+(~~xs)) => sum(map(𝛷, ~~xs))
    # @rule 𝛷(*(~~xs)) => prod(map(𝛷, ~~xs))

    @rule 𝛷(sin(~x)) => cos(~x) * ∂(~x)^-1
    @rule 𝛷(cos(~x)) => sin(~x) * ∂(~x)^-1
    @rule 𝛷(tan(~x)) => log(cos(~x)) * ∂(~x)^-1
    @rule 𝛷(csc(~x)) => log(sin(~x)^-1 - cos(~x)*sin(~x)^-1) * ∂(~x)^-1
    @rule 𝛷(sec(~x)) => log(cos(~x)^-1 + sin(~x)*cos(~x)^-1) * ∂(~x)^-1
    @rule 𝛷(cot(~x)) => log(sin(~x)) * ∂(~x)^-1

    @rule 𝛷(sinh(~x)) => cosh(~x) * ∂(~x)^-1
    @rule 𝛷(cosh(~x)) => sinh(~x) * ∂(~x)^-1
    @rule 𝛷(tanh(~x)) => log(cosh(~x)) * ∂(~x)^-1
    @rule 𝛷(csch(~x)) => log(sinh(~x)^-1 - cosh(~x)*sinh(~x)^-1) * ∂(~x)^-1
    @rule 𝛷(sech(~x)) => log(cosh(~x)^-1 + sinh(~x)*cosh(~x)^-1) * ∂(~x)^-1
    @rule 𝛷(coth(~x)) => log(sinh(~x)) * ∂(~x)^-1

    @rule 𝛷(asin(~x)) => (~x*asin(~x) + sqrt(1 - ~x*~x)) * ∂(~x)^-1
    @rule 𝛷(acos(~x)) => (~x*acos(~x) + sqrt(1 - ~x*~x)) * ∂(~x)^-1
    @rule 𝛷(atan(~x)) => (~x*atan(~x) + log(~x*~x + 1)) * ∂(~x)^-1
    @rule 𝛷(acsc(~x)) => acsc(~x) * ∂(~x)^-1
    @rule 𝛷(asec(~x)) => asec(~x) * ∂(~x)^-1
    @rule 𝛷(acot(~x)) => (~x*acot(~x) + log(~x*~x + 1)) * ∂(~x)^-1

    @rule 𝛷(asinh(~x)) => (~x*asinh(~x) + sqrt(~x*~x + 1)) * ∂(~x)^-1
    @rule 𝛷(acosh(~x)) => (~x*acosh(~x) + sqrt(~x*~x - 1)) * ∂(~x)^-1
    @rule 𝛷(atanh(~x)) => (~x*atanh(~x) + log(~x + 1)) * ∂(~x)^-1
    @rule 𝛷(acsch(~x)) => acsch(~x) * ∂(~x)^-1
    @rule 𝛷(asech(~x)) => asech(~x) * ∂(~x)^-1
    @rule 𝛷(acoth(~x)) => (~x*acot(~x) + log(~x + 1)) * ∂(~x)^-1

    @rule 𝛷(log(~x)) => (~x + ~x * log(~x)) * ∂(~x)^-1
    # @rule 𝛷(sqrt(~x)) => ~x * sqrt(~x) * ∂(~x)^-1

    # @rule 𝛷(^(-1 + ~x::is_sqr, -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + -(~x::is_sqr), -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + ~x::is_sqr, -1)) => atan(sqrt_of(~x))

    @rule 𝛷(^(~x, ~k::is_abs_half)) => sum(candidate_sqrt(~x,~k); init=one(~x))
    @rule 𝛷(^(~x::is_poly, ~k::is_neg)) => sum(candidate_pow_minus(~x, ~k); init=one(~x))
    @rule 𝛷(sqrt(~x)) => sum(candidate_sqrt(~x,0.5); init=one(~x))
    @rule 𝛷(^(sqrt(~x),-1)) => 𝛷(^(~x,-0.5))


    @rule 𝛷(^(~x, -1)) => log(~x) * ∂(~x)^-1
    @rule 𝛷(1 / ~x) => 𝛷(^(~x, -1))
    # @rule 𝛷(^(~x, ~k::is_pos_int)) => sum(^(~x, i) for i=1:~k+1) * ∂(~x)^-1
    # @rule 𝛷(^(~x, ~k::is_neg_int)) => sum(^(~x, i) for i=~k:-1) * ∂(~x)^-1
    @rule 𝛷(^(~x, ~k)) => ^(~x, ~k+1) * ∂(~x)^-1

    @rule 𝛷(exp(~x)) => exp(~x) * ∂(~x)^-1
    @rule 𝛷(1) => 𝑥
    @rule 𝛷(~x) => (~x + ^(~x,2)) * ∂(~x)^-1
]

apply_H_rules(eq) = expand(Fixpoint(Prewalk(Chain(H_rules))))(𝛷(value(eq)))
