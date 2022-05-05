using Documenter
using SymbolicNumericIntegration

makedocs(
    sitename="SymbolicNumericIntegration",
    format=Documenter.HTML(),
    modules=[SymbolicNumericIntegration]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="https://github.com/SciML/SymbolicNumericIntegration.jl";
    push_preview=true
)
