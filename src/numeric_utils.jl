
rms(x) = sqrt(sum(x .^ 2) / length(x))

"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol = 1e-3)
    Q, R = qr(A)
    abs.(diag(R)) .> abstol
end

function test_point(complex_plane, radius)
    if complex_plane
        return radius * sqrt(rand()) * cis(2π * rand())
    else
        return Complex(radius * (2 * rand() - 1))
    end
end

function accept_solution(eq, x, sol, radius)
    try
        x₀ = test_point(true, radius)
        Δ = substitute(expand_derivatives(Differential(x)(sol) - eq), Dict(x => x₀))
        return abs(Δ)
    catch e
        #
    end
    return Inf
end

"""
    converts float to int or small rational numbers
"""

function nice_parameter(u::Complex{T}; abstol = 1e-3, M = 10) where {T <: Real}
    α = nice_parameter(real(u))
    β = nice_parameter(imag(u))
    return β ≈ 0 ? α : Complex(α, β)
end

nice_parameter(x; abstol=1e-7) = small_rational(x; abstol)

nice_parameters(p; abstol = 1e-3) = nice_parameter.(p)

function small_rational(x::T; abstol=1e-6, M=20) where {T <: Real}
    # c = lcm(collect(1:M)...)
   	a = floor(Int, x)
	r = x - a
    for den in 2:M
        try
            if abs(round(r * den) - r * den) < abstol
                s = a + round(Int, r * den) // den
                return (denominator(s) == 1 ? numerator(s) : s)
            end
        catch e
        end
    end
    return x
end




