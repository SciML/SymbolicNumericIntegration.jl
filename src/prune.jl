

function normal(eq, x)
    return eq / coef(eq, x)
end

function generate_fragments(C, x)
    D = Differential(x)
    S = Set{Any}[]

    for c in C
        s = Set{Any}()
        push!(S, s)

        for t in terms(expand_derivatives(D(c)))
            push!(s, normal(t, x))
        end
    end

    return S
end

function overlap(U, V; abstol=1e-6)
    for u in U
        for v in V
            if abs(u - v) < abstol
                return true
            end
        end
    end
    return false
end

function prune(C, eq, x=var(eq); abstol=1e-6)
    eq = normal(eq, x)
    n = length(C)

    S = generate_fragments(C, x)

    x₀ = test_point(true, 1.0)
    d = Dict(x => x₀)
    u = Complex(substitute(eq, d))

    V = Any[]
    for i = 1:n
        push!(V, [Complex(substitute(s, d)) for s in S[i]])
        push!(V[i], sum(V[i]))
    end

    W = []
    U = []
    selected = falses(n)

    for i = 1:n
        if overlap([u], V[i])
            W = W ∪ V[i]
            selected[i] = true
        end
    end

    U = copy(W)
    done = false

    while !done
        done = true
        for i = 1:n
            if !selected[i] && overlap(U, V[i])
                W = W ∪ V[i]
                selected[i] = true
                done = false
            end
        end
    end

    T = unique([one(x); [C[i] for i=1:n if selected[i]]])

    for u in W
        if !overlap([u], U) return C, false end
    end

    return T, true
end

function contains_value(l, val; abstol=1e-6)
    return minimum(abs.(l .- val)) < abstol
end

##############################################################################

@rule tan(~x) => sin(~x) * cos(~x)^-1
@rule csc(~x) => sin(~x)^-1
@rule sec(~x) => cos(~x)^-1
@rule cot(~x) => cos(~x) * sin(~x)^-1

@rule tanh(~x) => sinh(~x) * cosh(~x)^-1
@rule csch(~x) => sinh(~x)^-1
@rule sech(~x) => cosh(~x)^-1
@rule coth(~x) => cosh(~x) * sinh(~x)^-1


@rule 𝛷(asinh(~x)) => one(~x) + asinh(~x) + ~x*asinh(~x) + sqrt(~x*~x + 1)
@rule 𝛷(acosh(~x)) => one(~x) + acosh(~x) + ~x*acosh(~x) + sqrt(~x*~x - 1)
@rule 𝛷(atanh(~x)) => one(~x) + atanh(~x) + ~x*atanh(~x) + log(~x + 1)
@rule 𝛷(acsch(~x)) => one(~x) + acsch(~x)
@rule 𝛷(asech(~x)) => one(~x) + asech(~x)
@rule 𝛷(acoth(~x)) => one(~x) + acoth(~x) + ~x*acot(~x) + log(~x + 1)
