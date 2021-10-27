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

function find_parts(eq::SymbolicUtils.Pow, x, η)
    y, k = arguments(eq)

    if is_poly(y)
        return Part(η ^ k, Dict(η => y))
    end

    p = find_parts(y, x, η)
    if k > 0
        return p
    else
        return Part(inv(p.eq), p.sub)
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

##############################################################################

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

    @rule 𝛷(^(~x, -1)) => log(~x)
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

@syms 𝜂

function generate_by_parts(eq, x=var(eq); max_terms=50)
    if !isdependent(eq, x) return [one(x)] end
    D = Differential(x)
    S = Set{Any}()
    T = Set{Any}()
    Q = Queue{Any}()
    eq = eq / coef(eq, x)
    push!(T, eq)
    enqueue!(Q, eq)

    while !isempty(Q) && length(S) < max_terms
        y = dequeue!(Q)

        if !is_multiple_x(y, x)
            w = x * y
            push!(S, w)
            if w ∉ T
                enqueue!(Q, w)
                push!(T, w)
            end
        end

        ps = find_parts(y, x, 𝜂)

        if ps == nothing continue end
        if !(ps isa AbstractArray) ps = [ps] end

        for p in ps
            u = p.eq
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
                    enqueue!(Q, w)
                    push!(T, w)
                end
            end
        end
    end
    unique([one(x); [s for s in S if isdependent(s,x)]])
end
