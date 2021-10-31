@syms 𝜂

struct Part
    eq
    sub::Union{Dict,Nothing}
end

find_parts(eq::SymbolicUtils.Add, x, η) = Part(η, Dict(η => eq))

function find_parts(eq::SymbolicUtils.Mul, x, η)
    l = Part[]
    for t in arguments(eq)
        p = find_parts(t, x, η)
        if p != nothing
            push!(l, p)
        end
    end
    l
end

function find_parts(eq::SymbolicUtils.Div, x, η)
    a, b = arguments(eq)
    return find_parts(a*b^-1, x, η)
end

function find_parts(eq::SymbolicUtils.Pow, x, η)
    y, k = arguments(eq)

    # if k > 0
    #     return find_parts(y, x, η)
    # elseif k ≈ -1
    #     return Part(eq, Dict(η => x))
    # else
    #     # p = normal_form(y, x, η)
    #     p = find_parts(y, x, η)
    #     return Part(inv(p.eq), p.sub)
    # end

    if k ≈ -1
        return Part(eq, Dict(η => x))
    elseif k < 0
        p = find_parts(y, x, η)
        return Part(inv(p.eq), p.sub)
    else
        p = find_parts(y, x, η)
        return p
    end
end

function find_parts(eq::SymbolicUtils.Term, x, η)
    f = operation(eq)
    y = arguments(eq)[1]
    return Part(f(η), Dict(η => y))
end

function find_parts(eq, x, η)
    if isdependent(eq, x)
        return Part(η, Dict(η => eq))
    else
        return nothing
    end
end

function normal_form(eq, x, η)
    if !is_poly(eq) return Part(η, Dict(η => eq)) end

    d = poly_deg(eq)

    if d == 1
        return Part(η, Dict(η => eq))
    elseif d == 2
        r, s = find_roots(eq, x)

        if !isempty(r)
            if r[1] ≈ r[2]
                return Part(η^2, Dict(η => x - r[1]))
            else
                D = Dict(η => 2*(x-(r[1]+r[2])/2)/abs(r[1]-r[2]))
                l = leading(eq, x)
                return l > 0 ? Part(η^2 - 1, D) : Part(1 - η^2, D)
            end
        else
            D = Dict(η => (x-abs(real(s[1])))/abs(imag(s[1])))
            return Part(η^2 + 1, D)
        end
    end

    return Part(η, Dict(η => eq))
end

##############################################################################

is_sqr_rules = [
    @rule ~x * ~x => true
    @rule ^(~x, 2) => true
]

is_sqr(eq) = isequal(Chain(is_sqr_rules)(value(eq)), true)

sqrt_rules = [
    @rule ~x * ~x => ~x
    @rule ^(~x, 2) => ~x
]

sqrt_of(eq) = is_sqr(eq) ? Chain(sqrt_rules)(value(eq)) : sqrt(eq)

i_rules = [
    @rule 𝛷(+(~~xs)) => sum(map(𝛷, ~~xs))
    @rule 𝛷(*(~~xs)) => prod(map(𝛷, ~~xs))

    @rule 𝛷(sin(~x)) => cos(~x)
    @rule 𝛷(cos(~x)) => sin(~x)
    @rule 𝛷(tan(~x)) => log(cos(~x))
    @rule 𝛷(csc(~x)) => log(sin(~x)^-1 - cos(~x)*sin(~x)^-1)
    @rule 𝛷(sec(~x)) => log(cos(~x)^-1 + sin(~x)*cos(~x)^-1)
    @rule 𝛷(cot(~x)) => log(sin(~x))

    @rule 𝛷(sinh(~x)) => cosh(~x)
    @rule 𝛷(cosh(~x)) => sinh(~x)
    @rule 𝛷(tanh(~x)) => log(cosh(~x))
    @rule 𝛷(csch(~x)) => log(sinh(~x)^-1 - cosh(~x)*sinh(~x)^-1)
    @rule 𝛷(sech(~x)) => log(cosh(~x)^-1 + sinh(~x)*cosh(~x)^-1)
    @rule 𝛷(coth(~x)) => log(sinh(~x))

    @rule 𝛷(asin(~x)) => ~x*asin(~x) + sqrt(1 - ~x*~x)
    @rule 𝛷(acos(~x)) => ~x*acos(~x) + sqrt(1 - ~x*~x)
    @rule 𝛷(atan(~x)) => ~x*atan(~x) + log(~x*~x + 1)
    @rule 𝛷(acsc(~x)) => acsc(~x)
    @rule 𝛷(asec(~x)) => asec(~x)
    @rule 𝛷(acot(~x)) => ~x*acot(~x) + log(~x*~x + 1)

    @rule 𝛷(asinh(~x)) => ~x*asinh(~x) + sqrt(~x*~x + 1)
    @rule 𝛷(acosh(~x)) => ~x*acosh(~x) + sqrt(~x*~x - 1)
    @rule 𝛷(atanh(~x)) => ~x*atanh(~x) + log(~x + 1)
    @rule 𝛷(acsch(~x)) => acsch(~x)
    @rule 𝛷(asech(~x)) => asech(~x)
    @rule 𝛷(acoth(~x)) => ~x*acot(~x) + log(~x + 1)

    @rule 𝛷(log(~x)) => ~x + ~x * log(~x)
    @rule 𝛷(sqrt(~x)) => ~x * sqrt(~x)

    # @rule 𝛷(^(-1 + ~x::is_sqr, -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + -(~x::is_sqr), -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + ~x::is_sqr, -1)) => atan(sqrt_of(~x))

    # @rule 𝛷(^(~x, -1)) => log(~x)
    @rule 𝛷(^(~x, -1)) => sum(candidate_pow_minus(~x, -1))
    @rule 𝛷(1 / ~x) => log(~x)
    @rule 𝛷(^(~x, ~k)) => ^(~x, ~k+1)

    @rule 𝛷(exp(~x)) => exp(~x)
    @rule 𝛷(~x) => ~x * var(~x)
]

apply_i_rules(eq) = expand(Fixpoint(Prewalk(Chain(i_rules))))(𝛷(value(eq)))

##############################################################################

is_multiple_x(eq::SymbolicUtils.Mul, x) = any(z -> is_multiple_x(z,x), arguments(eq))
is_multiple_x(eq::SymbolicUtils.Pow, x) = is_multiple_x(arguments(eq)[1], x)
is_multiple_x(eq, x) = is_poly(eq)

sum_power(eq::SymbolicUtils.Add, x) = maximum(sum_power(t,x) for t in arguments(eq))
sum_power(eq::SymbolicUtils.Mul, x) = sum(sum_power(t,x) for t in arguments(eq))
sum_power(eq::SymbolicUtils.Pow, x) = arguments(eq)[2]
sum_power(eq, x) = 0

function generate_by_parts(eq, x=var(eq); max_terms=20)
    if !isdependent(eq, x) return [one(x)] end
    D = Differential(x)
    S = Set{Any}()
    T = Set{Any}()
    Q₁ = Queue{Any}()
    eq = eq / coef(eq, x)

    for y in terms(eq)
        w = x * y
        push!(S, w)
        if w ∉ T
            enqueue!(Q₁, w)
            push!(T, w)
        end

        ps = find_parts(y, x, 𝜂)

        if ps == nothing continue end
        if !(ps isa AbstractArray) ps = [ps] end

        for p in ps
            u = p.eq
            printstyled("integrating ", u, " with ", p.sub, '\n'; color=:blue)
            U = apply_i_rules(u)
            u = substitute(u, p.sub)
            U = substitute(U, p.sub)
            u′ = expand_derivatives(D(substitute(𝜂, p.sub)))
            v = simplify_fractions(y * (u * u′)^-1)
            for t in terms(expand(U*v))
                w = simplify_fractions(t / coef(t, x))
                push!(S, w)
            end
            v′ = expand_derivatives(D(v))
            for t in terms(expand(U*v′))
                w = simplify_fractions(t / coef(t, x))
                if w ∉ T
                    enqueue!(Q₁, w)
                    push!(T, w)
                end
            end
        end
    end

    for i = 1:1
        Q₂ = Queue{Any}()
        while !isempty(Q₁) # && length(S) < max_terms
            y = dequeue!(Q₁)
            enqueue_expr_ex!(S, Q₂, expand_derivatives(D(y)), x)
        end
        Q₁ = Q₂
    end

    unique([1; [s for s in S if isdependent(s,x)]])
end

##############################################################################

@syms 𝑥
@syms u[20]

mutable struct Transform
    k::Int
    sub::Dict
end

transformer(eq::SymbolicUtils.Add, f) = sum(transformer(t,f) for t in arguments(eq); init=0)
transformer(eq::SymbolicUtils.Mul, f) = prod(transformer(t,f) for t in arguments(eq); init=1)
transformer(eq::SymbolicUtils.Div, f) = transformer(arguments(eq)[1],f) * transformer(inv(arguments(eq)[2]),f)

function transformer(eq::SymbolicUtils.Pow, f)
    y, k = arguments(eq)
    if k > 0
        return transformer(y, f)^k
    else
        μ = u[f.k]
        f.k += 1
        f.sub[μ] = inv(y)
        return μ ^ -k
    end
end

function transformer(eq, f)
    if isdependent(eq, 𝑥)
        μ = u[f.k]
        f.k += 1
        f.sub[μ] = eq
        return μ
    else
        return 1
    end
end

function transform(eq, x)
    eq = substitute(eq, Dict(x => 𝑥))
    f = Transform(1, Dict())
    return transformer(eq, f), f.sub
end

function homotopy_integrand(eq, x)
    eq, sub = transform(eq, x)
    I = 0
    n = length(sub)

    for i = 1:n
        μ = u[i]
        H = apply_H_rules(sub[μ])
        I += H * expand_derivatives(Differential(μ)(eq))
    end

    I = substitute(I, sub)
    I = substitute(I, Dict(𝑥 => x))
    return expand(I)
end

function expand_integrand(I, x)
    S = Set{Any}()
    T = Set{Any}()
    Q₁ = Queue{Any}()

    enqueue_expr_ex!(S, Q₁, expand(I + x*I), x)

    println(S)
    D = Differential(x)

    for i = 1:2
        Q₂ = Queue{Any}()
        while !isempty(Q₁) # && length(S) < max_terms
            y = dequeue!(Q₁)
            enqueue_expr_ex!(S, Q₂, expand_derivatives(D(y)), x)
        end
        Q₁ = Q₂
    end

    return [one(x); [s for s in S]]
end

function generate_homotopy(eq, x=var(eq))
    I = homotopy_integrand(eq, x)
    expand_integrand(I, x)
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
    @rule 𝛷(sqrt(~x)) => ~x * sqrt(~x) * ∂(~x)^-1

    # @rule 𝛷(^(-1 + ~x::is_sqr, -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + -(~x::is_sqr), -1)) => log(sqrt_of(~x) - 1) + log(sqrt_of(~x) + 1)
    # @rule 𝛷(^(1 + ~x::is_sqr, -1)) => atan(sqrt_of(~x))

    @rule 𝛷(^(~x::is_poly, ~k::is_neg)) => sum(candidate_pow_minus(~x, ~k); init=one(~x))
    @rule 𝛷(sqrt(~x)) => sum(candidate_sqrt(~x,0.5); init=one(~x))
    @rule 𝛷(^(sqrt(~x),-1)) => 𝛷(^(~x,-0.5))
    @rule 𝛷(^(~x, ~k::is_abs_half)) => sum(candidate_sqrt(~x,~k); init=one(~x))    

    @rule 𝛷(^(~x, -1)) => log(~x) * ∂(~x)^-1
    @rule 𝛷(1 / ~x) => 𝛷(^(~x, -1))
    @rule 𝛷(^(~x, ~k)) => ^(~x, ~k+1) * ∂(~x)^-1

    @rule 𝛷(exp(~x)) => exp(~x) * ∂(~x)^-1
    @rule 𝛷(~x) => ^(~x,2) * ∂(~x)^-1
]

apply_H_rules(eq) = expand(Fixpoint(Prewalk(Chain(H_rules))))(𝛷(value(eq)))
