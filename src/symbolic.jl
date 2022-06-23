function try_symbolic(T, eq, x, basis, Δbasis; kwargs...)
    eq isa Num || return 0, 0
    n = length(basis)
    # @syms θ[1:n]
    @variables θ[1:n]

    Δeq = sum(θ[j] * Δbasis[j] for j in 1:n) - eq
    Δeq = expand(Δeq)

    terms = collect_terms(Δeq, x)
    eqs = collect(values(terms))
    # eqs = filter(p->length(get_variables(p)) > 0, eqs)

    sol = solve_symbolic(eqs)

    for i in 1:n
        if !haskey(sol, θ[i])
            sol[θ[i]] = 0
        end
    end

    p = substitute(expand(sum(θ[j] * basis[j] for j in 1:n)), sol)
    return p, 0
end

mutable struct Fragment
    eq::Any
    lhs::Any
end

function solve_symbolic(eqs)
    n = length(eqs)
    solved = Set()
    unfinished = true
    frags = Fragment[]
    k = 1

    while unfinished
        unfinished = false
        for eq in eqs
            δf = [v for v in get_variables(eq) if v ∉ solved]
            if length(δf) == 1
                push!(solved, δf[1])
                push!(frags, Fragment(eq, δf[1]))
                unfinished = true
            end
        end
    end

    sol = Dict()
    for f in frags
        eq = substitute(f.eq, sol)
        u = Symbolics.solve_for(eq ~ 0, f.lhs)
        sol[f.lhs] = nice_parameter(u)
    end

    sys = []
    vars = Set()

    for eq in eqs
        δf = [v for v in get_variables(eq) if v ∉ solved]
        if length(δf) > 1
            for v in δf
                push!(vars, v)
            end
            q = substitute(eq, sol)
            push!(sys, q)
        end
    end

    sys = unique(sys)
    vars = [v for v in vars]

    if !isempty(vars) && length(vars) == length(sys)
        try
            vals = nice_parameters.(Symbolics.solve_for(sys .~ 0, vars))
            # vals = Symbolics.solve_for(sys .~ 0, vars)
            for (v, u) in zip(vars, vals)
                sol[v] = u
            end
        catch e
            # println("from symbolic: ", e)
        end
    end

    return sol
end

function collect_terms(eq::SymbolicUtils.Add, x)
    println(eq)
    d = Dict{Any,Any}(1 => 0)
    for t in arguments(eq)
        if isdependent(t, x)
            for s in collect_terms(t, x)
                v, c = first(s), last(s)
                if haskey(d, v)
                    d[v] += c
                else
                    d[v] = c
                end
            end
        else
            d[1] += t
        end
    end
    return d
end

function collect_terms(eq, x)
    c = coef(eq, x)
    return Dict(eq / c => c)
end
