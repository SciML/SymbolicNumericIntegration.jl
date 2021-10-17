function try_integration_by_parts(eq, x, l; kwargs...)
    args = Dict(kwargs)
    verbose = args[:verbose]

    f = factors(eq, x)
    if length(f) <= 2 return zero(x), eq, Inf end

    D = Differential(x)
    ϵ₀ = Inf

    kw = copy(kwargs)
    kw[:bypart] = false

    k = 1

    for u in f
        v′ = eq * inverse(u)
        if !is_number(u) && !is_number(v′)
            attempt(l, "Integrating by parts (trial $k). u is ", u)
            inform(l, "v' is ", v′)

            v, r, ϵ = integrate_term(v′, x, l; kwargs...)
            if isequal(r, 0)
                inform(l, "v is ", v)
                u′v = expand(v*expand_derivatives(D(u)))
                inform(l, "u′*v is ", u′v)
                s, r, ϵ = integrate_sum_fast(u′v, x, l; kw...)
                if isequal(r, 0)
                    y = expand(u*v - s)
                    result(l, "Success", y)
                    return y, 0, ϵ
                end
            end
            result(l, "Trial $k failed")
            k += 1
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
