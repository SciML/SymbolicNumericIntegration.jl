const DEBUG_CACHE = true

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

function deriv!(c::ExprCache, x)
    if c.δeq == nothing
        c.δeq = expand_derivatives(Differential(x)(expr(c)))
    end
    return c.δeq
end

function deriv!(c, x)
    if DEBUG_CACHE
        error("ExprCache object expected")
    end
    return expand_derivatives(Differential(x)(c))
end

function fun!(c::ExprCache, x)
    if c.f == nothing
        c.f = build_function(expr(c), x; expression = false)
    end
    return c.f
end

function fun!(c, x)
    if DEBUG_CACHE
        error("ExprCache object expected")
    end
    return build_function(c, x; expression = false)
end

function deriv_fun!(c::ExprCache, x)
    if c.δf == nothing
        c.δf = build_function(deriv!(c, x), x; expression = false)
    end
    return c.δf
end

function deriv_fun!(c, x)
    if DEBUG_CACHE
        error("ExprCache object expected")
    end
    return build_function(deriv!(c, x), x; expression = false)
end

Base.show(io::IO, c::ExprCache) = show(expr(c))
