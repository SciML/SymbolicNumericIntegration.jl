"""
    generate_basis(eq, x, try_kernel = true)

Generate candidate antiderivative terms for a symbolic-numeric integration problem.

## Arguments

  - `eq`: Scalar symbolic integrand.
  - `x`: Independent symbolic variable.
  - `try_kernel`: Whether to factor candidate kernels before generating terms.

## Returns

A vector of internal cached symbolic expressions. Pass the result to `best_hints`
or use it to inspect the candidate basis; it is not a stable cache representation.

## Examples

```julia
julia> using Symbolics, SymbolicNumericIntegration

julia> @variables x

julia> !isempty(generate_basis(x^2, x))
true
```
"""
function generate_basis(eq, x, try_kernel = true)
    if !try_kernel
        S = sum(generate_homotopy(expr(eq), x))
        return cache.(unique([equivalent(t, x) for t in terms(S)]))
    end

    S = 0
    eq = expand(expr(eq))

    for t in terms(eq)
        q = equivalent(t, x)
        f = kernel(q)
        p = q / f

        if isdependent(p, x)
            C₂ = generate_homotopy(p, x)
        else
            C₂ = 1
        end

        C₁ = closure(f, x)

        S += sum(c₁ * c₂ for c₁ in C₁ for c₂ in C₂)
    end
    return cache.(unique([equivalent(t, x) for t in terms(S)]))
end

function expand_basis(basis, x; Kmax = 1000)
    if isempty(basis)
        return basis, false
    end

    b = sum(expr.(basis))

    Kb = complexity(b) # Kolmogorov complexity
    if is_proper(Kb) && Kb > Kmax
        return basis, false
    end

    δb = sum(deriv!.(basis, x))
    eq = (1 + x) * (b + δb)
    eq = expand(eq)
    S = Set{Any}()
    enqueue_expr!(S, eq, x)
    return cache.([s for s in S]), true
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
        y = popfirst!(q)
        enqueue_expr!(S, q, expand_derivatives(D(y)), x)
    end
    return unique([[s for s in S]; [s * x for s in S]])
end

###############################################################################

enqueue_expr!(S, q, eq, x) = enqueue_expr!!(S, q, ops(eq)..., x)

function enqueue_expr!!(S, q, ::Add, eq, x)
    for t in arguments(eq)
        enqueue_expr!(S, q, t, x)
    end
    return
end

function enqueue_expr!!(S, q, ::Any, eq, x)
    y = eq / term_coefficient(eq, x)
    return if y ∉ S && isdependent(y, x)
        push!(q, y)
        push!(S, y)
    end
end

enqueue_expr!(S, eq, x) = enqueue_expr!!(S, ops(eq)..., x)

function enqueue_expr!!(S, ::Add, eq, x)
    for t in arguments(eq)
        enqueue_expr!(S, t, x)
    end
    return
end

function enqueue_expr!!(S, ::Any, eq, x)
    y = eq / term_coefficient(eq, x)
    return if y ∉ S && isdependent(y, x)
        push!(S, y)
    end
end
