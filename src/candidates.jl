using DataStructures

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x)
    eq = expand(eq)
    S = Set{Any}()
    for t in terms(eq)
        q = t / coef(t, x)
        f = kernel(q)
        C₁ = closure(f, x) # find_candidates(f, x)
        C₂ = find_candidates_nonsolvable(q * inverse(f), x)
        # C₂ = generate_by_parts(q * inverse(f), x)

        for c₁ in C₁
            enqueue_expr_ex!(S, c₁, x)
        end

        for c₂ in C₂
            enqueue_expr_ex!(S, c₂, x)
        end

        for c₁ in C₁
            for c₂ in C₂
                enqueue_expr_ex!(S, expand(c₁*c₂), x)
            end
        end
    end
    return unique([one(x); [s for s in S]])
end

function closure(eq, x; max_terms=50)
    if !isdependent(eq, x) return [one(x)] end
    D = Differential(x)
    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr_ex!(S, q, eq, x)

    while !isempty(q) && length(S) < max_terms
        y = dequeue!(q)
        enqueue_expr_ex!(S, q, expand_derivatives(D(y)), x)
    end
    unique([one(x); [s for s in S]; [s*x for s in S]])
end

function find_candidates_nonsolvable(eq, x)
    eq = apply_d_rules(eq)
    D = Differential(x)

    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr_ex!(S, q, eq, x)

    for y in q
        ∂y = expand_derivatives(D(y))
        enqueue_expr_ex!(S, ∂y, x)
    end

    return unique([one(x); [s for s in S]])
end

function candidate_pow_minus(p, k; abstol=1e-8)
    if isnan(poly_deg(p))
        return [p^k, p^(k+1), log(p)]
    end

    x = var(p)
    r, s = find_roots(p, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    # ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
    q = Any[log(x - u) for u in r]
    for i in eachindex(s)
        β = s[i]
        if abs(imag(β)) > abstol
            push!(q, atan((x - real(β))/imag(β)))
            push!(q, (1+x)*log(x^2 - 2*real(β)*x + abs2(β)))
        else
            push!(q, log(x - real(β)))
        end
    end
    q = unique(q)

    # return [[p^k, p^(k+1)]; candidates(q₁, x)]
    if k ≈ -1
        return [[p^k]; q]
    else
        return [[p^k, p^(k+1)]; q]
    end
end

function candidate_sqrt(p, k)
    x = var(p)

    h = Any[p^k, p^(k+1)]

    if poly_deg(p) == 2
        r, s = find_roots(p, x)
        l = leading(p, x)
        if length(r) == 2
            if sum(r) ≈ 0
                r₁ = abs(r[1])
                if l > 0
                    push!(h, acosh(x/r₁))
                else
                    push!(h, asin(x/r₁))
                end
            end
        elseif real(s[1]) ≈ 0
            push!(h, asinh(x/imag.(s[1])))
        end
    end

    Δ = expand_derivatives(Differential(x)(p))
    push!(h, log(0.5*Δ + sqrt(p)))
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
    if y ∉ S && isdependent(y, x) # && all(u->u>=0, extract_power(y))
        enqueue!(q, y)
        push!(S, y)
    end
end

function enqueue_expr_ex!(S, q, eq::SymbolicUtils.Add, x)
    for t in arguments(eq)
        enqueue_expr_ex!(S, q, t, x)
    end
end

function enqueue_expr_ex!(S, q, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x)
        enqueue!(q, y)
        push!(S, y)
    end
end

function enqueue_expr_ex!(S, eq::SymbolicUtils.Add, x)
    for t in arguments(eq)
        enqueue_expr_ex!(S, t, x)
    end
end

function enqueue_expr_ex!(S, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x)
        push!(S, y)
    end
end

###############################################################################

extract_power(eq::SymbolicUtils.Pow) = [arguments(eq)[2]]
extract_power(eq::SymbolicUtils.Term) = [1]
extract_power(eq::SymbolicUtils.Mul) = union([extract_power(t) for t in arguments(eq)]...)
extract_power(eq) = []
