mutable struct ExprCache
    eq::Any# the primary expression
    f::Any# compiled eq
    δeq::Any# the symbolic derivative of eq
    δf::Any# compiled δeq
end

cache(eq::ExprCache) = eq
cache(eq) = ExprCache(eq, nothing, nothing, nothing)

expr(c::ExprCache) = c.eq
expr(c) = c

function deriv!(c::ExprCache, xs...)
    if c.δeq == nothing
        c.δeq = expand_derivatives(Differential(xs[1])(expr(c)))
    end
    return c.δeq
end

function deriv!(c, xs...)
    return expand_derivatives(Differential(xs[1])(c))
end

function fun!(c::ExprCache, xs...)
    if c.f == nothing
        c.f = build_function(expr(c), xs...; expression = false)
    end
    return c.f
end

function fun!(c, xs...)
    return build_function(c, xs...; expression = false)
end

function deriv_fun!(c::ExprCache, xs...)
    if c.δf == nothing
        c.δf = build_function(deriv!(c, xs...), xs...; expression = false)
    end
    return c.δf
end

function deriv_fun!(c, xs...)
    return build_function(deriv!(c, xs...), xs...; expression = false)
end

Base.show(io::IO, c::ExprCache) = show(expr(c))
