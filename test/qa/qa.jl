using SafeTestsets

@safetestset "Aqua" begin
    using SymbolicNumericIntegration
    using Aqua
    using Test

    Aqua.test_all(SymbolicNumericIntegration; piracies = false)
    @test_broken false  # Aqua piracy: Base.signbit(::Complex)/signbit(::SymbolicUtils.Sym) in src/integral.jl + DataDrivenSparse.active_set! in src/sparse.jl — see https://github.com/SciML/SymbolicNumericIntegration.jl/issues/156
end

@safetestset "JET" begin
    using SymbolicNumericIntegration
    using JET
    using Test

    JET.test_package(SymbolicNumericIntegration; target_defined_modules = true)
end
