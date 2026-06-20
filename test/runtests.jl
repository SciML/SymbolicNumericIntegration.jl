using Test

using SymbolicNumericIntegration
using SymbolicNumericIntegration: value
using Symbolics

using SymbolicUtils
using SymbolicUtils.Rewriters

using SafeTestsets
using SciMLTesting

run_tests(;
    core = () -> begin
        @safetestset "integral" begin
            include("integral_tests.jl")
        end
        @safetestset "vector expression error handling" begin
            include("vector_expression_tests.jl")
        end
        @safetestset "definite integral with infinite bounds" begin
            include("definite_integral_tests.jl")
        end
    end,
    qa = (; env = joinpath(@__DIR__, "qa"), body = joinpath(@__DIR__, "qa", "qa.jl")),
)
