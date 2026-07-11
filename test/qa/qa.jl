using SymbolicNumericIntegration, Test
using JET
using SciMLTesting

run_qa(
    SymbolicNumericIntegration;
    explicit_imports = true,
    api_docs_kwargs = (; rendered = true),
    jet_kwargs = (; target_defined_modules = true),
    # Aqua piracy: Base.signbit(::Complex)/signbit(::SymbolicUtils.Sym) in
    # src/integral.jl + DataDrivenSparse.active_set! in src/sparse.jl
    # https://github.com/SciML/SymbolicNumericIntegration.jl/issues/156
    aqua_broken = (:piracies,),
    ei_kwargs = (;
        # Re-exports: accessed from the re-exporting module rather than the owner.
        # coef -> StatsAPI (via DataDrivenSparse); toexpr -> SymbolicUtils.Code
        # (via Symbolics). Neither owner is a SciML make-public target.
        all_qualified_accesses_via_owners = (;
            ignore = (:coef, :toexpr),
        ),
        # Non-public names of upstream packages accessed via qualification.
        # AbstractDataDrivenAlgorithm: DataDrivenDiffEq; active_set!, coef:
        # DataDrivenSparse; toexpr: Symbolics.
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractDataDrivenAlgorithm, :active_set!, :coef, :toexpr,
            ),
        ),
    ),
    # Heavy `using` of the Symbolics/SymbolicUtils/DataDriven stacks (macros +
    # module bindings used for qualified access) — tracked for a follow-up
    # explicit-import pass rather than a risky mass refactor.
    # https://github.com/SciML/SymbolicNumericIntegration.jl/issues/164
    ei_broken = (:no_implicit_imports,),
)
