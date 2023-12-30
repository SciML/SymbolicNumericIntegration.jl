using SymbolicNumericIntegration, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(SymbolicNumericIntegration)
    Aqua.test_ambiguities(SymbolicNumericIntegration, recursive = false)
    Aqua.test_deps_compat(SymbolicNumericIntegration)
    Aqua.test_piracies(SymbolicNumericIntegration, broken = true)
    Aqua.test_project_extras(SymbolicNumericIntegration)
    Aqua.test_stale_deps(SymbolicNumericIntegration)
    Aqua.test_unbound_args(SymbolicNumericIntegration)
    Aqua.test_undefined_exports(SymbolicNumericIntegration)
end
