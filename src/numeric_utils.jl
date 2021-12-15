
rms(x) = sqrt(sum(x.^2) / length(x))

"""
    returns a list of the indices of a linearly independent subset of the columns of A
"""
function find_independent_subset(A; abstol=1e-3)
    Q, R = qr(A)
    abs.(diag(R)) .> abstol
end

function test_point(complex_plane, radius)
    if complex_plane
        return radius * sqrt(rand()) * cis(2π*rand())
    else
        return Complex(radius * (2*rand() - 1))
    end
end

function accept_solution(eq, x, sol, radius)
    try
        x₀ = test_point(true, radius)
        Δ = substitute(expand_derivatives(Differential(x)(sol)-eq), Dict(x => x₀))
        return abs(Δ)
    catch e
        #
    end
    return Inf
end

"""
    converts float to int or small rational numbers
"""
function nice_parameters(p; abstol=1e-3)
    c = lcm(collect(1:10)...)
    n = length(p)
    q = Array{Any}(undef, n)
    for i = 1:n
        den = 1
        while den < 10
            if abs(round(p[i]*den) - p[i]*den) < abstol
                a = round(Int, p[i]*den) // den
                q[i] = (denominator(a) == 1 ? numerator(a) : a)
                den = 10
            else
                q[i] = Float64(p[i])
            end
            den += 1
        end
    end
    q
end

function nice_parameter(u::T; abstol=1e-3, M=10) where T<:Real
    c = lcm(collect(1:M)...)
    for den = 1:M
        try
            if abs(round(u*den) - u*den) < abstol
                a = round(Int, u*den) // den
                return (denominator(a) == 1 ? numerator(a) : a)
            end
        catch e
        end
    end
    return u
end

function nice_parameter(u::Complex{T}; abstol=1e-3, M=10) where T<:Real
    α = nice_parameter(real(u))
    β = nice_parameter(imag(u))
    return β ≈ 0 ? α : Complex(α, β)
end
