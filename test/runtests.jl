using SymbolicNumericIntegration
using SymbolicNumericIntegration: value
using Symbolics

using SymbolicUtils
using SymbolicUtils.Rewriters

using Test

include("axiom.jl")

##############################################################################

@variables x a b c β

"""
    a list of basic standard integral tests
    based on http://integral-table.com/ with modifications
"""
basic_integrals = [
    # Basic Forms
    1,
    x^2,
    4x^3,
    # Integrals of Rational Functions
    1 / x,
    1 / (2x + 5),
    1 / (x + 1)^2,
    (x + 3)^3,
    x * (x - 2)^4,
    1 / (1 + x^2),
    1 / (9 + x^2),
    x / (4 + x^2),
    x^2 / (16 + x^2),
    x^3 / (1 + x^2),
    1 / (x^2 - 5x + 6),
    1 / (x^2 + x + 1),
    x / (x + 4)^2,
    x / (x^2 + x + 1),
    # Integrals with Roots
    sqrt(x - 2),
    1 / sqrt(x - 1),
    1 / sqrt(x + 1),
    1 / sqrt(4 - x),
    x * sqrt(x - 3),
    sqrt(2x + 5),
    (3x - 1)^1.5,
    x / sqrt(x - 1),
    x / sqrt(x + 1),
    sqrt(x / (4 - x)),
    sqrt(x / (4 + x)),
    x * sqrt(2x + 3),
    sqrt(x * (x + 2)),
    sqrt(x^3 * (x + 3)),
    sqrt(x^2 + 4),
    sqrt(x^2 - 4),
    sqrt(4 - x^2),
    x * sqrt(x^2 + 9),
    x * sqrt(x^2 - 9),
    1 / sqrt(x^2 + 4),
    1 / sqrt(x^2 - 4),
    1 / sqrt(4 - x^2),
    x / sqrt(x^2 + 4),
    x / sqrt(x^2 - 4),
    x / sqrt(4 - x^2),
    x^2 / sqrt(x^2 + 4),
    x^2 / sqrt(x^2 - 4),
    sqrt(x^2 - 5x + 6),
    x * sqrt(x^2 - 5x + 6),
    1 / sqrt(x^2 - 5x + 6),
    1 / (4 + x^2)^1.5,
    # Integrals with Logarithms
    log(x),
    x * log(x),
    x^2 * log(x),
    log(2x) / x,
    log(x) / x^2,
    log(2x + 1),
    log(x^2 + 4),
    log(x^2 - 4),
    log(x^2 - 5x + 6),
    x * log(x + 2),
    x * log(9 - 4x^2),
    log(x)^2,
    log(x)^3,
    x * log(x)^2,
    x^2 * log(x)^2,
    # Integrals with Exponentials
    exp(x),
    sqrt(x) * exp(x),
    x * exp(x),
    x * exp(3x),
    x^2 * exp(x),
    x^2 * exp(5x),
    x^3 * exp(x),
    x^3 * exp(2x),
    exp(x^2),
    x * exp(x^2),
    # Integrals with Trigonometric Functions
    sin(4x),
    sin(x)^2,
    sin(x)^3,
    cos(3x),
    cos(x)^2,
    cos(2x)^3,
    sin(x) * cos(x),
    sin(3x) * cos(5x),
    sin(x)^2 * cos(x),
    sin(3x)^2 * cos(x),
    sin(x) * cos(x)^2,
    sin(x) * cos(5x)^2,
    sin(x)^2 * cos(x),
    sin(x)^2 * cos(x)^2,
    sin(4x)^2 * cos(4x)^2,
    tan(x),
    tan(7x),
    tan(x)^2,
    tan(x)^3,
    sec(x),
    sec(x) * tan(x),
    sec(x)^2 * tan(x),
    csc(x),
    sec(x) * csc(x),
    # Products of Trigonometric Functions and Monomials
    x * cos(x),
    x * cos(3x),
    x^2 * cos(x),
    x^2 * cos(5x),
    x * sin(x),
    x * sin(3x),
    x^2 * sin(x),
    x^2 * sin(5x),
    x * cos(x)^2,
    x * sin(x)^2,
    x * tan(x)^2,
    x * sec(x)^2,
    x^3 * sin(x),
    x^4 * cos(2x),
    sin(x)^2 * cos(x)^3,
    # Products of Trigonometric Functions and Exponentials
    exp(x) * sin(x),
    exp(3x) * sin(2x),
    exp(x) * cos(x),
    exp(2x) * cos(7x),
    x * exp(x) * sin(x),
    x * exp(x) * cos(x),
    # Integrals of Hyperbolic Functions
    cosh(x),
    exp(x) * cosh(x),
    sinh(3x),
    exp(2x) * sinh(3x),
    tanh(x),
    exp(x) * tanh(x),
    cos(x) * cosh(x),
    cos(x) * sinh(x),
    sin(x) * cosh(x),
    sin(x) * sinh(x),
    sinh(x) * cosh(x),
    sinh(3x) * cosh(5x),
    # Misc
    exp(x) / (1 + exp(x)),
    cos(exp(x)) * sin(exp(x)) * exp(x),
    cos(exp(x))^2 * sin(exp(x)) * exp(x),
    1 / (x * log(x)),
    1 / (exp(x) - 1),
    1 / (exp(x) + 5),
    sqrt(x) * log(x),
    log(log(x)) / x,
    x^3 * exp(x^2),
    sin(log(x)),
    x * cos(x) * exp(x),
    log(x - 1)^2,
    1 / (exp(2x) - 1),
    exp(x) / (exp(2x) - 1),
    x / (exp(2x) - 1),
    # derivative-divide examples (Lamangna 7.10.2)
    exp(x) * exp(exp(x)),
    exp(sqrt(x)) / sqrt(x),
    log(log(x)) / (x * log(x)),
    log(cos(x)) * tan(x),
    # rothstein-Trager examples (Lamangna 7.10.9)
    1 / (x^3 - x),
    1 / (x^3 + 1),
    1 / (x^2 - 8),
    (x + 1) / (x^2 + 1),
    x / (x^4 - 4),
    x^3 / (x^4 + 1),
    1 / (x^4 + 1),
    # exponential/trigonometric/logarithmic integral functions
    exp(2x) / x,
    exp(x + 1) / (x + 1),
    x * exp(2x^2 + 1) / (2x^2 + 1),
    sin(3x) / x,
    sin(x + 1) / (x + 1),
    cos(5x) / x,
    x * cos(x^2 - 1) / (x^2 - 1),
    1 / log(3x - 1),
    1 / (x * log(log(x))),
    x / log(x^2),
    # from Geddes & Stefanus
    exp(1 / (1 + log(x))) * ((2 * log(x) + 1) / x),
    1 / (1 + exp(x)),
    sin(2x) * exp(x),
    -10 * exp(x^-10) / x^11,
    # bypass = true
    β,      # turn of bypass = true
    (log(x - 1) + (x - 1)^-1) * log(x),
    exp(x) / x - exp(x) / x^2,
    cos(x) / x - sin(x) / x^2,
    1 / log(x) - 1 / log(x)^2,
]

sym_integrals = [
    # Basic Forms
    a * x^2,
    a * x + b * x^2 - c * x^3,
    (3(x^2) + 3(a^2) - 4) / (4 * a * b),
    a / x,
    1 / (a * x + 5),
    1 / (x + a),
    x / (x + a),
    1 / (x + a)^2,
    x / (x + a)^2,
    (x + a)^3,
    x * (x - a)^4,
    1 / (a + x^2),
    sqrt(x - a),
    1 / sqrt(a * x - 1),
    x * sqrt(a * x + b),
    log(a * x),
    x * log(a * x),
    x^2 * log(a * x),
    log(a * x) / x,
    log(x^2 - a * x + b),
    log(a * x)^2,
    log(a + x),
    log(a + x^2) * x,
    # x^2 * log(a * x + b)^2,
    exp(a * x),
    x * exp(a * x),
    x^2 * exp(a * x),
    x * exp(a * x^2),
    sin(a * x),
    sin(a * x)^2,
    cos(a * x + b)^2,
    sin(a * x) * cos(a * x),
    sin(a * x) * cos(b * x),
    tan(a * x),
    sec(a * x),
    x * cos(a * x),
    x^2 * cos(a * x),
    exp(a * x) * sin(b * x),
    x * exp(a * x) * sin(a * x),
    x * exp(a * x) * cos(b * x),
    cosh(a * x),
    exp(a * x) * cosh(b * x),
    cos(a * x) * cosh(b * x),
    sin(a * x) * cos(b * x) * exp(c * x),
    sin(a * x) * sinh(b * x) * exp(c * x),
    sec(a * x)^2 * tan(a * x),
    exp(a * x) / (1 + exp(a * x)),
    exp(a * x) / exp(b * x),
    cos(exp(a * x)) * sin(exp(a * x)) * exp(a * x),
    1 / (x * log(a * x)),
    log(log(a * x)) / x,
    log(a + log(x)) / x,
    sin(log(a * x)),
    # x / (exp(a * x) - b),
    exp(a * x) / (b * exp(a * x) + c),
    exp(a * x) * exp(exp(a * x)),
    log(cos(a * x)) * tan(a * x),
    1 / (x^3 + a),
    exp(a * x + b) / x,
    sin(x + a) / (x + a),
    cos(a * x) / x,
    x / log(a * x^2),
    # bypass = true
    β,      # turn of bypass = true
    exp(a * x) / x - exp(a * x) / x^2,
    cos(a * x) / x - sin(a * x) / x^2,
]

function test_integrals(basic = true, subs = nothing; kw...)
    args = isempty(kw) ? Dict() : Dict(kw)
    args[:detailed] = false
    misses = []
    k = 1

    integrals = basic ? basic_integrals : sym_integrals
    args[:symbolic] = !basic

    for (i, eq) in enumerate(integrals)
        if isequal(eq, β)
            printstyled("**** bypass on ****\n"; color = :red)
            bypass = true
            args[:bypass] = true
        else
            if subs != nothing
                eq = substitute(eq, subs)
            end

            printstyled(k, ": "; color = :blue)
            k += 1
            printstyled(eq, " =>\n"; color = :green)
            sol = SymbolicNumericIntegration.integrate(eq, x; args...)
            if sol == nothing
                printstyled("\t<no solution>\n"; color = :red)
                push!(misses, eq)
            else
                printstyled('\t', sol, '\n'; color = :cyan)
            end
        end
    end

    n = length(misses)
    if n > 0
        println("**** missess (n=$n) *****")
    end
    for eq in misses
        printstyled(eq, '\n'; color = :red)
    end
    return n
end

@testset "integral" begin
    n = test_integrals(;
        symbolic = false, verbose = false, homotopy = true, num_steps = 2,
        num_trials = 10
    )
    @test n > 0
end

@testset "vector expression error handling" begin
    @variables x α

    # Test that vector expressions throw an appropriate error
    @test_throws ErrorException("Vector expressions are not supported. Please use element-wise integration with `integrate.([expr1, expr2, ...], x)` instead.") integrate([x])
    @test_throws ErrorException("Vector expressions are not supported. Please use element-wise integration with `integrate.([expr1, expr2, ...], x)` instead.") integrate(
        [1, 2 * α], α)

    # Test that scalar integration still works
    @test integrate(x) == ((1 // 2) * (x^2), 0, 0)
    @test integrate(2 * α, α) == (α^2, 0, 0)

    # Test that element-wise integration works
    results = integrate.([1, 2 * α], α)
    @test length(results) == 2
    @test results[1] == (α, 0, 0)
    @test results[2] == (α^2, 0, 0)
end

@testset "definite integral with infinite bounds" begin
    @variables x

    # Test exp(-x) from 0 to Inf - should give 1
    result1 = integrate(exp(-x), (x, 0, Inf); symbolic = true, detailed = false)
    @test result1 ≈ 1.0

    # Test x*exp(-x) from 0 to Inf - should give 1 (Gamma(2) = 1!)
    result2 = integrate(x * exp(-x), (x, 0, Inf); symbolic = true, detailed = false)
    @test result2 ≈ 1.0

    # Test x^2*exp(-x) from 0 to Inf - should give 2 (Gamma(3) = 2!)
    result3 = integrate(x^2 * exp(-x), (x, 0, Inf); symbolic = true, detailed = false)
    @test result3 ≈ 2.0

    # Test exp(x) from -Inf to 0 - should give 1
    result4 = integrate(exp(x), (x, -Inf, 0); symbolic = true, detailed = false)
    @test result4 ≈ 1.0

    # Note: x*exp(x) from -Inf to 0 is skipped due to a bug in SymbolicLimits.jl
    # with limits at -Inf involving products

    # Finite bounds should still work
    result6 = integrate(x, (x, 0, 1); symbolic = false, detailed = false)
    @test result6 == 1 // 2
end
