"""
    solve_newton is a symbolic Newton-Ralphson solver
    f is a symbolic equation to be solved (f ~ 0)
    x is the variable to solve
    x₀ is the initial guess
"""
function solve_newton(T, p, ∂p, x, x₀, zs; abstol = 1e-10, maxiter = 50, s = 1)
    d = Dict(x => x₀)
    xₙ = x₀

    for i in 1:maxiter
        d[x] = xₙ
        f = T(substitute(p, d))
        f′ = T(substitute(∂p, d))
        ρ = sum(1 / (xₙ - z) for z in zs; init = 0)
        xₙ₊₁ = xₙ - s * f / (f′ - s * ρ * f)

        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return nothing
end

function find_roots(T, p, x; abstol = 1e-8, num_roots = 0)
    n = (num_roots == 0 ? poly_deg(p) : num_roots)
    abstol = T(abstol)

    zs = Complex{T}[]
    r = T[]
    s = Complex{T}[]

    while !isequal(p, 0)
        q = expand(p * x^-1)
        if !is_poly(q)
            break
        end
        push!(r, 0)
        p = q
        n -= 1
    end

    if isequal(p, 0)
        return r, s
    end

    ∂p = expand_derivatives(Differential(x)(p))

    while length(zs) < n
        z = solve_newton(Complex{T}, p, ∂p, x, cis(2π * rand()), zs; abstol)
        if z != nothing
            if abs(imag(z)) < abstol
                push!(zs, Complex(real(z)))
                push!(r, real(z))
            else
                if abs(real(z)) < abstol
                    z = Complex(0, imag(z))
                end
                push!(zs, z)
                push!(s, z)

                z = conj(z)
                if abs(Complex{T}(substitute(p, Dict(x => z)))) < abstol
                    push!(zs, z)
                    push!(s, z)
                end
            end
        else
            break
        end
    end

    return sort(r), s
end

find_roots(p, x; abstol = 1e-8, num_roots = 0) = find_roots(Float64, p, x; num_roots, abstol)

function find_poles(T, p, x; abstol = 1e-10)
    abstol = Float64(abstol)
    ∂p = expand_derivatives(Differential(x)(p))

    zs = Complex{T}[]
    r = Float64[]
    s = Complex{T}[]

    while true
        z = solve_newton(Complex{T}, p, ∂p, x, cis(2π * rand()), zs; abstol, s = -1)
        if z != nothing
            if abs(imag(z)) < abstol
                push!(zs, Complex(real(z)))
                push!(r, real(z))
            else
                if abs(real(z)) < abstol
                    z = Complex(0, imag(z))
                end
                push!(zs, z)
                push!(s, z)

                z = conj(z)
                if abs(Complex{Float64}(substitute(p, Dict(x => z)))) < abstol
                    push!(zs, z)
                    push!(s, z)
                end
            end
        else
            break
        end
    end

    return sort(r), s
end

find_poles(p, x; abstol = 1e-10) = find_poles(Float64, p, x; abstol)
