using DataStructures

# this is the main heurisctic used to find the test fragments
function generate_basis(eq, x, h=[]; use_closure=true, use_rules=false)
    if use_closure
        f = kernel(eq, x)
        C₂ = find_candidates(f, x)
    else
        f = 1
        C₂ = [one(x)]
    end

    g = eq / f

    if use_rules
        C₁ = find_candidates_nonsolvable(g, x)
    else
        h , _ = collect_hints(g, x)
        H = prod(h; init=one(x))
        Δg = expand_derivatives(Differential(x)(g))
        kers = expand(g + Δg + H)
        C₁ = [one(x); candidates(kers, x)]
    end

    # println("C₁ = ", C₁)
    # println("C₂ = ", C₂)

    return [c₁*c₂ for c₁ in C₁ for c₂ in C₂]
end

function find_candidates_nonsolvable(eq, x)
    eq = apply_d_rules(eq)
    # println(">>> ", eq)
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

"""
    candidates returns a list of candidate expressions to form the integration
    basis
"""
candidates(eq, x) = isdependent(eq,x) ? [eq] : []
candidates(eq::Num, x) = candidates(value(eq), x)

# the candidates of an Add is the union of the candidates of the terms
# ∫ Σᵢ fᵢ(x) dx = Σᵢ ∫ fᵢ(x) dx
candidates(eq::SymbolicUtils.Add, x) = unique(∪([candidates(t,x) for t in arguments(eq)]...))

# the candidates of a Mul is the outer product of the candidates of the terms
# d(uv)/dx = u dv/dx + v + du/dx
function candidates(eq::SymbolicUtils.Mul, x)
    terms = [candidates(q,x) for q in arguments(eq)]
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

    unique(l[2:end])    # removing the initial 1
end

function candidates_ex(eq::SymbolicUtils.Mul, x)
    terms = [sum(candidates(q,x); init=1) for q in arguments(eq)]
    # println(">> ", terms)
    eq = expand(prod(terms; init=1))
    # println(">>> ", eq)
    D = Differential(x)
    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr_ex!(S, q, eq, x)

    for y in q
        ∂y = expand_derivatives(D(y))
        enqueue_expr_ex!(S, ∂y, x)
    end

    return [one(x); [s for s in S]]
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

function enqueue_expr_ex!(S, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x)
        push!(S, y)
    end
end


# the candidates of a Pow encode different integration rules
function candidates(eq::SymbolicUtils.Pow, x)
    if !isdependent(eq,x) return [one(x)] end

    p = arguments(eq)[1]    # eq = p ^ k
    k = arguments(eq)[2]

    # if k < 0 && k ≈ round(k)
    if k ≈ -1
        return candidate_pow_minus(p, k, x)
    elseif k ≈ 0.5 || k ≈ -0.5
        if check_poly(p,x) == :real_poly && leading(p,x) < 0 p = -p end
        # ∫ √f(x) dx = ... + c * log(df/dx + √f) if deg(f) == 2
        Δ = expand_derivatives(Differential(x)(p))
        return [p^k, p^(k+1), log(0.5*Δ + sqrt(p))]
    elseif k < 0
        return [p^i for i=k:0 if i<0]
    end

    # ∫ p^k dp = c * p^(k+1)
    # return [p^k, p^(k+1)]
    return [p^i for i=k+1:-1:0 if i>0]
end

nice_abs2(u) = abs2(u)     # nice_parameter(abs2(u))

function candidate_pow_minus(p, k, x)
    check = check_poly(p, x)
    if check == :not_poly || check == :complex_poly
        return [p^k, p^(k+1), log(p)]
        # Δp = expand_derivatives(Differential(x)(p))
        # return [p^k, p^(k+1), log(p); Δp]
    end

    r, s = find_roots(p, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    # ∫ 1 / ((x-z₁)(x-z₂)) dx = ... + c₁ * log(x-z₁) + c₂ * log(x-z₂)
    q = [[log(x - u) for u in r];
          [atan((x - real(u))/imag(u)) for u in s];
          [log(x^2 - 2*real(u)*x + nice_abs2(u)) for u in s]
         ]

    # return [[p^k, p^(k+1)]; candidates(q₁, x)]
    if k ≈ -1
        return [[p^k]; q]
    else
        return [[p^k, p^(k+1)]; q]
    end
end

###############################################################################


function enqueue_expr!(S, q, eq::SymbolicUtils.Add, x)
    for t in arguments(eq)
        enqueue_expr!(S, q, t, x)
    end
end

function enqueue_expr!(S, q, eq, x)
    y = eq / coef(eq, x)
    if y ∉ S && isdependent(y, x) && all(u->u>=0, extract_power(y))
        enqueue!(q, y)
        push!(S, y)
    end
end

function closure(eq, x; max_terms=50)
    if !isdependent(eq, x) return [one(x)] end
    D = Differential(x)
    S = Set{Any}()
    q = Queue{Any}()
    enqueue_expr!(S, q, eq, x)

    while !isempty(q) && length(S) < max_terms
        y = dequeue!(q)
        enqueue_expr!(S, q, expand_derivatives(D(y)), x)
    end
    [one(x); [s for s in S]]
end

function find_candidates(eq, x)
    P, Q = split_frac(eq)
    PP = closure(P, x)
    if isequal(Q, 1) return PP end
    QQ = generate_basis(1/Q, x)

    S = Set()

    for p in PP
        for q in QQ
            push!(S, p*q)
        end
    end
    [s for s in S]
end

###############################################################################

extract_power(eq::SymbolicUtils.Pow) = [arguments(eq)[2]]
extract_power(eq::SymbolicUtils.Term) = [1]
extract_power(eq::SymbolicUtils.Mul) = union([extract_power(t) for t in arguments(eq)]...)
extract_power(eq) = []

function split_frac(eq::SymbolicUtils.Pow)
    if arguments(eq)[2] >=0
        return [eq,1]
    else
        return [1,1/eq]
    end
end

split_frac(eq) = [eq,1]

function split_frac(eq::SymbolicUtils.Mul)
    P = 1
    Q = 1
    for t in arguments(eq)
        p, q = split_frac(t)
        P *= p
        Q *= q
    end
    [P, Q]
end

extract_factors(eq::SymbolicUtils.Pow) = [arguments(eq)[1]]
extract_factors(eq::SymbolicUtils.Mul) = unique(union([extract_factors(t) for t in arguments(eq)]...))
extract_factors(eq) = [eq]

###############################################################################

kernel(eq::SymbolicUtils.Add, x) = sum(kernel(t,x) for t in arguments(eq); init=0)
kernel(eq::SymbolicUtils.Mul, x) = prod(kernel(t,x) for t in arguments(eq); init=1)

function kernel(eq::SymbolicUtils.Pow, x)
    p = arguments(eq)[1]    # eq = p ^ k
    k = arguments(eq)[2]
    if isinteger(k) && k >= 0
        kernel(p, x)^k
    else
        return 1
    end
end

function kernel(eq::SymbolicUtils.Term, x)
    op = operation(eq)
    p = arguments(eq)[1]

    if is_linear_poly(p,x) && (op == sin || op == cos || op == exp || op == sinh || op == cosh)
        return eq
    else
        return 1
    end
end

kernel(eq, x) = isequal(eq, x) ? x : 1
