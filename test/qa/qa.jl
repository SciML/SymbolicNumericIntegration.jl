using SymbolicNumericIntegration, Test
using JET
using SciMLTesting

run_qa(
    SymbolicNumericIntegration;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    # Aqua piracy: Base.signbit(::Complex)/signbit(::SymbolicUtils.Sym) in
    # src/integral.jl + DataDrivenSparse.active_set! in src/sparse.jl
    # https://github.com/SciML/SymbolicNumericIntegration.jl/issues/156
    aqua_broken = (:piracies,),
    ei_kwargs = (;
        # Re-exports: accessed from the re-exporting module rather than the owner.
        # coef -> StatsAPI (via DataDrivenSparse); scalarize/toexpr/unwrap ->
        # SymbolicUtils[.Code] (via Symbolics).
        all_qualified_accesses_via_owners = (;
            ignore = (:coef, :scalarize, :toexpr, :unwrap),
        ),
        # Non-public names of upstream packages accessed via qualification.
        # AbstractDataDrivenAlgorithm: DataDrivenDiffEq; active_set!, coef:
        # DataDrivenSparse; Sym: SymbolicUtils; derivative, scalarize, toexpr,
        # unwrap, value: Symbolics.
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractDataDrivenAlgorithm, :Sym, :active_set!, :coef,
                :derivative, :scalarize, :toexpr, :unwrap, :value,
            ),
        ),
        # Non-public names explicitly imported from upstream packages.
        # BasicSymbolic, issym: SymbolicUtils; get_variables, value: Symbolics.
        all_explicit_imports_are_public = (;
            ignore = (:BasicSymbolic, :issym, :get_variables, :value),
        ),
    ),
    # Heavy `using` of the Symbolics/SymbolicUtils/DataDriven stacks (macros +
    # module bindings used for qualified access) — tracked for a follow-up
    # explicit-import pass rather than a risky mass refactor.
    # https://github.com/SciML/SymbolicNumericIntegration.jl/issues/164
    ei_broken = (:no_implicit_imports,),
)
