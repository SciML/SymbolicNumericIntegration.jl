using SymbolicNumericIntegration
using Symbolics

using SymbolicUtils
using SymbolicUtils.Rewriters

using Test
using PyCall, SymPy


include("axiom.jl")

##############################################################################

@variables x β

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
    sin(4x)^2 * cos(4x^2),
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
    (log(x - 1) + (x - 1)^-1) * log(x),
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
    # bypass = true
    β,      # turn of bypass = true
    exp(x) / x - exp(x) / x^2,
    cos(x) / x - sin(x) / x^2,
    1 / log(x) - 1 / log(x)^2,
]

function test_integrals(; kw...)
    args = Dict(kw)
    misses = []
    k = 1

    for eq in basic_integrals
        if isequal(eq, β)
            printstyled("**** bypass on ****\n"; color=:red)
            bypass = true
            args[:bypass] = true
        else
            printstyled(k, ": "; color=:blue)
            k += 1
            printstyled(eq, " =>\n"; color=:green)
            solved, unsolved = SymbolicNumericIntegration.integrate(eq; args...)
            printstyled('\t', solved; color=:cyan)
            if isequal(unsolved, 0)
                println()
            else
                printstyled(" + ∫ ", unsolved, '\n'; color=:red)
                push!(misses, eq)
            end
        end
    end

    n = length(misses)
    if n > 0
        println("**** missess (n=$n) *****")
    end
    for eq in misses
        printstyled(eq, '\n'; color=:red)
    end
    return n
end

@testset "integral" begin
    n = test_integrals(; symbolic=false, verbose=false, homotopy=true, num_steps=2, num_trials=10)
    @test n > 0
end
