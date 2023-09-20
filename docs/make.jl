using Documenter, SymbolicNumericIntegration

# cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
# cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "SymbolicNumericIntegration.jl",
    authors = "Shahriar Iravanian",
    modules = [SymbolicNumericIntegration],
    clean = true, doctest = false, linkcheck = true,
    warnonly = true,
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/SymbolicNumericIntegration/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/SymbolicNumericIntegration.jl.git";
    push_preview = true)
