using DataStructures

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x; homotopy = true)
    # if homotopy return generate_homotopy(eq, x) end
    eq = expand(expr(eq))
    S = 0 #Set{Any}()
    for t in terms(eq)
        q = equivalent(t, x)
        f = kernel(q)

        C₁ = closure(f, x)
        p = q * inverse(f)

        if homotopy
            if isdependent(p, x)
                C₂ = generate_homotopy(p, x)
            else
                C₂ = 1
            end
        else
            C₂ = find_candidates(p, x)
        end

        S += sum(c₁ * c₂ for c₁ in C₁ for c₂ in C₂)
    end
    return cache.(unique([one(x); [equivalent(t, x) for t in terms(S)]]))
end

function expand_basis(basis, x)
    b = sum(expr.(basis))
    δb = sum(deriv!.(basis, x))
    eq = (1 + x) * (b + δb)
    eq = expand(eq)
    S = Set{Any}()
    enqueue_expr!(S, eq, x)
    return cache.([one(x); [s for s in S]])
end

function closure(eq, x; max_terms = 50)
    if !isdependent(eq, x)
        return [one(x)]
    end
    D = Differential(x)
    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr!(S, q, eq, x)

    while !isempty(q) && length(S) < max_terms
        y = dequeue!(q)
        enqueue_expr!(S, q, expand_derivatives(D(y)), x)
    end
    unique([one(x); [s for s in S]; [s * x for s in S]])
end

function find_candidates(eq, x)
    eq = apply_d_rules(eq)
    D = Differential(x)

    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr!(S, q, eq, x)

    for y in q
        ∂y = expand_derivatives(D(y))
        enqueue_expr!(S, ∂y, x)
    end

    return unique([one(x); [s for s in S]])
end

function candidate_pow_minus(p, k; abstol = 1e-8)
    if isnan(poly_deg(p))
        return [p^k, p^(k + 1), log(p)]
    end

    x = var(p)
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

    # ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
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

    # return [[p^k, p^(k+1)]; candidates(q₁, x)]
    if k ≈ -1
        return [[p^k]; q]
    else
        return [[p^k, p^(k + 1)]; q]
    end
end

function candidate_sqrt(p, k)
    x = var(p)

    h = Any[p^k, p^(k + 1)]

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

    Δ = expand_derivatives(Differential(x)(p))
    push!(h, log(0.5 * Δ + sqrt(p)))
    return h
end

###############################################################################

function enqueue_expr!(S, q, eq::SymbolicUtils.Add, x)
    for t in arguments(eq)
        enqueue_expr!(S, q, t, x)
    end
end

function enqueue_expr!(S, q, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x)
        enqueue!(q, y)
        push!(S, y)
    end
end

function enqueue_expr!(S, eq::SymbolicUtils.Add, x)
    for t in arguments(eq)
        enqueue_expr!(S, t, x)
    end
end

function enqueue_expr!(S, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x)
        push!(S, y)
    end
end
