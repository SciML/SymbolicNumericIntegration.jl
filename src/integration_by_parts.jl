function try_integration_by_parts(eq, x; kwargs...)
    f = factors(eq, x)
    if length(f) <= 2 return zero(x), eq, Inf end

    D = Differential(x)
    ϵ₀ = Inf

    for u in f
        v′ = eq / u
        if !is_number(u) && !is_number(v′)
            v, r, ϵ = integrate_term(v′, x; kwargs...)
            if isequal(r, 0)
                uv = expand_derivatives(v*D(u))
                s, r, ϵ = integrate_sum(uv, x; kwargs...)
                if isequal(r, 0)
                    return expand(u*v - s), 0, ϵ
                else
                    zero(x), eq, ϵ
                    # ϵ₀ = min(ϵ, ϵ₀)
                end
            end
        end
    end

    return zero(x), eq, ϵ₀
end

factors(eq, x) = isdependent(eq, x) ? [one(x), eq] : [one(x)]

function factors(eq::SymbolicUtils.Pow, x)
    p, k = arguments(eq)
    [p^(i*sign(k)) for i=0:abs(k)]
end

function factors(eq::SymbolicUtils.Mul, x)
    terms = [factors(q,x) for q in arguments(eq)]
    n = length(terms)

    l = Any[one(x)]

    for j = 1:n
        m = length(l)
        for t in terms[j]
            for k = 1:m
                push!(l, l[k]*t)
            end
        end
    end

    unique(l)
end
