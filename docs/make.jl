using Documenter, SymbolicNumericIntegration

include("pages.jl")

makedocs(
    sitename = "SymbolicNumericIntegration.jl",
    authors = "Shahriar Iravanian",
    modules = [SymbolicNumericIntegration],
    clean = true,
    doctest = true,
    checkdocs = :all,
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/SymbolicNumericIntegration/stable/"
    ),
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/SymbolicNumericIntegration.jl.git";
    push_preview = true
)
