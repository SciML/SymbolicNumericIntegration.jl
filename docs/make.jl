using Documenter, SymbolicNumericIntegration

include("pages.jl")

makedocs(sitename = "SymbolicNumericIntegration.jl",
         authors = "Shahriar Iravanian",
         modules = [SymbolicNumericIntegration],
         clean = true, doctest = false,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://symbolicnumericintegration.sciml.ai/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/SymbolicNumericIntegration.jl.git";
           push_preview = true)
