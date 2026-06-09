using SymbolicNumericIntegration
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(SymbolicNumericIntegration)
end

@testset "JET" begin
    JET.test_package(SymbolicNumericIntegration; target_defined_modules = true)
end
