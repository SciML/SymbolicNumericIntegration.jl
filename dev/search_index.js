var documenterSearchIndex = {"docs":
[{"location":"symbolicnumericintegration/#API","page":"API","title":"API","text":"","category":"section"},{"location":"symbolicnumericintegration/","page":"API","title":"API","text":"integrate","category":"page"},{"location":"symbolicnumericintegration/#SymbolicNumericIntegration.integrate","page":"API","title":"SymbolicNumericIntegration.integrate","text":"integrate(eq, x; kwargs...)\n\nis the main entry point to integrate a univariate expression eq with respect to `x' (optional). \n\njulia> using Symbolics, SymbolNumericIntegration\n\njulia> @variables x a\n\njulia> integrate(x * sin(2x))\n((1//4)*sin(2x) - (1//2)*x*cos(2x), 0, 0)\n\njulia> integrate(x * sin(a*x), x; symbolic=true, detailed=false)\n(sin(a*x) - a*x*cos(a*x)) / (a^2)\n\njulia> integrate(x * sin(a*x), (x, 0, 1); symbolic=true, detailed=false)\n(sin(a) - a*cos(a)) / (a^2)\n\nArguments:\n\neq: a univariate expression\nx: independent variable (optional if eq is univariate) or a tuple        of (independent variable, lower bound, upper bound) for definite integration.\n\nKeyword Arguments:\n\nabstol (default: 1e-6): the desired tolerance\nnum_steps (default: 2): the number of different steps with expanding basis to be tried\nnum_trials (default: 10): the number of trials in each step (no changes to the basis)\nshow_basis (default: false): if true, the basis (list of candidate terms) is printed\nbypass (default: false): if true do not integrate terms separately but consider all at once\nsymbolic (default: false): try symbolic integration first (will be forced if eq has symbolic constants)\nmax_basis (default: 100): the maximum number of candidate terms to consider\nverbose (default: false): print a detailed report\ncomplex_plane (default: true): generate random test points on the complex plane (if false, the points will be on real axis)\nradius (default: 1.0): the radius of the disk in the complex plane to generate random test points\nopt (default: STLSQ(exp.(-10:1:0))): the sparse regression optimizer (from DataDrivenSparse)\nhomotopy (default: true): use the homotopy algorithm to generate the basis (deprecated, will be removed in a future version)\nuse_optim (default: false): use Optim.jl minimize function instead of the STLSQ algorithm (experimental)\ndetailed (default: true): (solved, unsolved, err) output format. If detailed=false, only the final integral is returned. \n\nOutput:\n\nsolved: the solved integral \nunsolved: the residual unsolved portion of the input\nerr: the numerical error in reaching the solution\n\n\n\n\n\nDefinite integral\n\n\n\n\n\n","category":"function"},{"location":"#SymbolicNumericIntegration.jl","page":"Home","title":"SymbolicNumericIntegration.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SymbolicNumericIntegration.jl is a hybrid symbolic/numerical integration package that works on the Julia Symbolics expressions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"SymbolicNumericIntegration.jl uses a randomized algorithm based on a hybrid of the method of undetermined coefficients and sparse regression and can solve a large subset of basic standard integrals (polynomials, exponential/logarithmic, trigonometric and hyperbolic, inverse trigonometric and hyperbolic, rational and square root). The symbolic part of the algorithm is similar (but not identical) to Risch-Bronstein's poor man's integrator and generates a list of ansatzes (candidate terms). The numerical part uses sparse regression adopted from the Sparse identification of nonlinear dynamics (SINDy) algorithm to prune down the ansatzes and find the corresponding coefficients. The basis of how it works and the theory of integration using the Symbolic-Numeric methods refer to Basis of Symbolic-Numeric Integration.","category":"page"},{"location":"","page":"Home","title":"Home","text":"hyint is the python counterpart of SymbolicNumericIntegration.jl that works with sympy expressions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Originally, SymbolicNumericIntegration.jl expected only univariate expression with constant real or complex coefficients as input. As of version 1.2, integrate function accepts symbolic constants for a subset of solvable integrals.","category":"page"},{"location":"","page":"Home","title":"Home","text":"integrate returns a tuple with three values. The first one is the solved integral, the second one is the sum of the unsolved terms, and the third value is the residual error. If integrate is successful, the unsolved portion is reported as 0. If we pass detailed=false to integrate', the output is simplified to only the resulting integrals. In this case,nothing` is returned if the integral cannot be solved. Note that the simplified output will become the default in a future version.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install SymbolicNumericIntegration.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"SymbolicNumericIntegration\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Examples:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using SymbolicNumericIntegration\njulia> using Symbolics\n\njulia> @variables x a b\n\n# if `detailed = true` (default), the output is a tuple of (solution, unsolved portion, err)\n\njulia> integrate(3x^3 + 2x - 5)\n(x^2 + (3//4)*(x^4) - (5x), 0, 0)\n\njulia> integrate((5 + 2x)^-1)\n((1//2)*log((5//2) + x), 0, 0.0)\n\n# `detailed = false` simplifies the output to just the resulting integral\n\njulia> integrate(x^2 / (16 + x^2); detailed = false)\nx + 4atan((-1//4)*x)\n\njulia> integrate(x^2 * log(x); detailed = false)\n(1//3)*(x^3)*log(x) - (1//9)*(x^3)\n\njulia> integrate(sec(x) * tan(x); detailed = false)\nsec(x)\n\n# Symbolic integration. Here, a is a symbolic constant; therefore, we need \n# to explicitly define the independent variable (say, x). Also, we set\n# `symbolic = true` to force using the symbolic solver\n\njulia> integrate(sin(a * x), x; detailed = false, symbolic = true)\n(-cos(a*x)) / a\n\njulia> integrate(x^2 * cos(a * x), x; detailed = false, symbolic = true)\n((a^2)*(x^2)*sin(a*x) + 2.0a*x*cos(a*x) - 2.0sin(a*x)) / (a^3)\n\njulia> integrate(log(log(a * x)) / x, x; detailed = false, symbolic = true)\nlog(a*x)*log(log(a*x)) - log(a*x)\n\n# multiple symbolic constants\n\njulia> integrate(cosh(a * x) * exp(b * x), x; detailed = false, symbolic = true)\n(a*sinh(a*x)*exp(b*x) - b*cosh(a*x)*exp(b*x)) / (a^2 - (b^2))\n\n# definite integration, passing a tuple of (x, lower bound, higher bound) in the \n# second argument\n\njulia> integrate(x * sin(a * x), (x, 0, 1); symbolic = true, detailed = false)\n(sin(a) - a*cos(a)) / (a^2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"SymbolicNumericIntegration.jl exports some special integral functions (defined over Complex numbers) and uses them in solving integrals:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ei: exponential integral (define as ∫ exp(x) / x dx)\nSi: sine integral (define as ∫ sin(x) / x dx)\nCi: cosine integral (define as ∫ cos(x) / x dx)\nLi: logarithmic integral (define as ∫ 1 / log(x) dx)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For examples:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> integrate(exp(x + 1) / (x + 1))\n(SymbolicNumericIntegration.Ei(1 + x), 0, 0.0)\n\njulia> integrate(x * cos(a*x^2 - 1) / (a*x^2 - 1), x; detailed=false, symbolic=true)\n((1//2)*SymbolicNumericIntegration.Ci(a*(x^2) - 1)) / a\n\njulia> integrate(1 / (x*log(log(x))), x; detailed=false, symbolic=true)\nSymbolicNumericIntegration.Li(log(x))","category":"page"},{"location":"","page":"Home","title":"Home","text":"integrate(eq, x; kwargs...)","category":"page"},{"location":"#SymbolicNumericIntegration.integrate-Tuple{Any, Any}","page":"Home","title":"SymbolicNumericIntegration.integrate","text":"integrate(eq, x; kwargs...)\n\nis the main entry point to integrate a univariate expression eq with respect to `x' (optional). \n\njulia> using Symbolics, SymbolNumericIntegration\n\njulia> @variables x a\n\njulia> integrate(x * sin(2x))\n((1//4)*sin(2x) - (1//2)*x*cos(2x), 0, 0)\n\njulia> integrate(x * sin(a*x), x; symbolic=true, detailed=false)\n(sin(a*x) - a*x*cos(a*x)) / (a^2)\n\njulia> integrate(x * sin(a*x), (x, 0, 1); symbolic=true, detailed=false)\n(sin(a) - a*cos(a)) / (a^2)\n\nArguments:\n\neq: a univariate expression\nx: independent variable (optional if eq is univariate) or a tuple        of (independent variable, lower bound, upper bound) for definite integration.\n\nKeyword Arguments:\n\nabstol (default: 1e-6): the desired tolerance\nnum_steps (default: 2): the number of different steps with expanding basis to be tried\nnum_trials (default: 10): the number of trials in each step (no changes to the basis)\nshow_basis (default: false): if true, the basis (list of candidate terms) is printed\nbypass (default: false): if true do not integrate terms separately but consider all at once\nsymbolic (default: false): try symbolic integration first (will be forced if eq has symbolic constants)\nmax_basis (default: 100): the maximum number of candidate terms to consider\nverbose (default: false): print a detailed report\ncomplex_plane (default: true): generate random test points on the complex plane (if false, the points will be on real axis)\nradius (default: 1.0): the radius of the disk in the complex plane to generate random test points\nopt (default: STLSQ(exp.(-10:1:0))): the sparse regression optimizer (from DataDrivenSparse)\nhomotopy (default: true): use the homotopy algorithm to generate the basis (deprecated, will be removed in a future version)\nuse_optim (default: false): use Optim.jl minimize function instead of the STLSQ algorithm (experimental)\ndetailed (default: true): (solved, unsolved, err) output format. If detailed=false, only the final integral is returned. \n\nOutput:\n\nsolved: the solved integral \nunsolved: the residual unsolved portion of the input\nerr: the numerical error in reaching the solution\n\n\n\n\n\n","category":"method"},{"location":"#Testing","page":"Home","title":"Testing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"test/runtests.jl contains a test suite of 170 easy to moderate test integrals (can be run by calling test_integrals). Currently, SymbolicNumericIntegration.jl solves more than 95% of its test suite.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Additionally, 12 test suites from the Rule-based Integrator (Rubi) are included in the /test directory. For example, we can test the first one as below. Axiom refers to the format of the test files)","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SymbolicNumericIntegration\ninclude(\"test/axiom.jl\")  # note, you may need to use the correct path\n\nL = convert_axiom(:Apostle)   # can also use L = convert_axiom(1)  \ntest_axiom(L, false; bypass = false, verbose = false, homotopy = true)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The test suites description based on the header of the files in the Rubi directory are","category":"page"},{"location":"","page":"Home","title":"Home","text":"name id comment\n:Apostle 1 Tom M. Apostol - Calculus, Volume I, Second Edition (1967)\n:Bondarenko 2 Vladimir Bondarenko Integration Problems\n:Bronstein 3 Manuel Bronstein - Symbolic Integration Tutorial (1998)\n:Charlwood 4 Kevin Charlwood - Integration on Computer Algebra Systems (2008)\n:Hearn 5 Anthony Hearn - Reduce Integration Test Package\n:Hebisch 6 Waldek Hebisch - email May 2013\n:Jeffrey 7 David Jeffrey - Rectifying Transformations for Trig Integration (1997)\n:Moses 8 Joel Moses - Symbolic Integration Ph.D. Thesis (1967)\n:Stewart 9 James Stewart - Calculus (1987)\n:Timofeev 10 A. F. Timofeev - Integration of Functions (1948)\n:Welz 11 Martin Welz - posts on Sci.Math.Symbolic\n:Webster 12 Michael Wester","category":"page"},{"location":"#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use SymbolicNumericIntegration.jl, please cite Symbolic-Numeric Integration of Univariate Expressions based on Sparse Regression:","category":"page"},{"location":"","page":"Home","title":"Home","text":"@article{Iravanian2022,   \n   author = {Shahriar Iravanian and Carl Julius Martensen and Alessandro Cheli and Shashi Gowda and Anand Jain and Julia Computing and Yingbo Ma and Chris Rackauckas},\n   doi = {10.48550/arxiv.2201.12468},\n   month = {1},\n   title = {Symbolic-Numeric Integration of Univariate Expressions based on Sparse Regression},\n   url = {https://arxiv.org/abs/2201.12468v2},\n   year = {2022},\n}","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nSee the SciML Style Guide for common coding practices and other style decisions.\nThere are a few community forums:\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Slack\nThe #diffeq-bridged and #sciml-bridged channels in the Julia Zulip\nOn the Julia Discourse forums\nSee also SciML Community page","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(; mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nusing Markdown\nversion = TOML.parse(read(\"../../Project.toml\", String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\", String))[\"name\"]\nlink_manifest = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n                \"/assets/Manifest.toml\"\nlink_project = \"https://github.com/SciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n               \"/assets/Project.toml\"\nMarkdown.parse(\"\"\"You can also download the\n[manifest]($link_manifest)\nfile and the\n[project]($link_project)\nfile.\n\"\"\")","category":"page"}]
}
