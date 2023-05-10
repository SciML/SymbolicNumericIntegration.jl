# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://docs.sciml.ai/Symbolics/stable/) expressions.

**SymbolicNumericIntegration.jl** uses a randomized algorithm based on a hybrid of the *method of undetermined coefficients* and *sparse regression* and can solve a large subset of basic standard integrals (polynomials, exponential/logarithmic, trigonometric and hyperbolic, inverse trigonometric and hyperbolic, rational and square root).
The basis of how it works and the theory of integration using the Symbolic-Numeric methods refer to [Basis of Symbolic-Numeric Integration](https://github.com/SciML/SymbolicNumericIntegration.jl/blob/main/docs/theory.ipynb).

Function `integrate` returns the integral of a univariate expression with *constant* real or complex coefficients. `integrate` returns a tuple with three values. The first one is the solved integral, the second one is the sum of the unsolved terms, and the third value is the residual error. If `integrate` is successful, the unsolved portion is reported as 0.

## Installation

To install SymbolicNumericIntegration.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SymbolicNumericIntegration")
```

Examples:

```julia
using Symbolics
using SymbolicNumericIntegration

@variables x
```

```julia
julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)

julia> integrate((5 + 2x)^-1)
((1//2)*log((5//2) + x), 0, 0.0)

julia> integrate(1 / (6 + x^2 - (5x)))
(log(x - 3) - log(x - 2), 0, 3.339372764128952e-16)

julia> integrate(1 / (x^2 - 16))
((1//8)*log(x - 4) - ((1//8)*log(4 + x)), 0, 1.546926788028958e-16)

julia> integrate(x^2 / (16 + x^2))
(x + 4atan((-1//4)*x), 0, 1.3318788420751984e-16)

julia> integrate(x^2 / sqrt(4 + x^2))
((1//2)*x*((4 + x^2)^0.5) - ((2//1)*log(x + sqrt(4 + x^2))), 0, 8.702422633074313e-17)

julia> integrate(x^2 * log(x))
((1//3)*log(x)*(x^3) - ((1//9)*(x^3)), 0, 0)

julia> integrate(x^2 * exp(x))
(2exp(x) + exp(x)*(x^2) - (2x*exp(x)), 0, 0)

julia> integrate(tan(2x))
((-1//2)*log(cos(2x)), 0, 0)

julia> integrate(sec(x) * tan(x))
(cos(x)^-1, 0, 0)

julia> integrate(cosh(2x) * exp(x))
((2//3)*exp(x)*sinh(2x) - ((1//3)*exp(x)*cosh(2x)), 0, 7.073930088880992e-8)

julia> integrate(cosh(x) * sin(x))
((1//2)*sin(x)*sinh(x) - ((1//2)*cos(x)*cosh(x)), 0, 4.8956233716268386e-17)

julia> integrate(cosh(2x) * sin(3x))
(0.153845sinh(2x)*sin(3x) - (0.23077cosh(2x)*cos(3x)), 0, 4.9807620877373405e-6)

julia> integrate(log(log(x)) * (x^-1))
(log(x)*log(log(x)) - log(x), 0, 0)

julia> integrate(exp(x^2))
(0, exp(x^2), Inf)    # as expected!
```

SymbolicNumericIntegration.jl exports some special integral functions (defined over Complex numbers) and uses them in solving integrals:

  - `Ei`: exponential integral (define as ∫ exp(x) / x dx)
  - `Si`: sine integral (define as ∫ sin(x) / x dx)
  - `Ci`: cosine integral (define as ∫ cos(x) / x dx)
  - `Li`: logarithmic integral (define as ∫ 1 / log(x) dx)

For examples:

```
julia> integrate(exp(x + 1) / (x + 1))
(SymbolicNumericIntegration.Ei(1 + x), 0, 1.1796119636642288e-16)

julia> integrate(x * cos(x^2 - 1) / (x^2 - 1))
((1//2)*SymbolicNumericIntegration.Ci(x^2 - 1), 0, 2.7755575615628914e-17)

julia> integrate(1 / (x*log(log(x))))
(SymbolicNumericIntegration.Li(log(x)), 0, 1.1102230246251565e-16)
```

```@docs
integrate(eq, x; kwargs...)
```

## Testing

`test/runtests.jl` contains a test suite of 160 easy to moderate test integrals (can be run by calling `test_integrals`). Currently, **SymbolicNumericIntegration.jl** solves more than 90% of its test suite.

Additionally, 12 test suites from the *Rule-based Integrator* ([Rubi](https://rulebasedintegration.org/)) are included in the `/test` directory. For example, we can test the first one as below ([Axiom](http://www.axiom-developer.org/) refers to the format of the test files)

```julia
using SymbolicNumericIntegration
include("test/axiom.jl")  # note, you may need to use the correct path

L = convert_axiom(:Apostle)   # can also use L = convert_axiom(1)  
test_axiom(L, false; bypass = false, verbose = false, homotopy = true)
```

The test suites description based on the header of the files in the Rubi directory are

| name        | id | comment                                                                |
|:----------- |:-- |:---------------------------------------------------------------------- |
| :Apostle    | 1  | Tom M. Apostol - Calculus, Volume I, Second Edition (1967)             |
| :Bondarenko | 2  | Vladimir Bondarenko Integration Problems                               |
| :Bronstein  | 3  | Manuel Bronstein - Symbolic Integration Tutorial (1998)                |
| :Charlwood  | 4  | Kevin Charlwood - Integration on Computer Algebra Systems (2008)       |
| :Hearn      | 5  | Anthony Hearn - Reduce Integration Test Package                        |
| :Hebisch    | 6  | Waldek Hebisch - email May 2013                                        |
| :Jeffrey    | 7  | David Jeffrey - Rectifying Transformations for Trig Integration (1997) |
| :Moses      | 8  | Joel Moses - Symbolic Integration Ph.D. Thesis (1967)                  |
| :Stewart    | 9  | James Stewart - Calculus (1987)                                        |
| :Timofeev   | 10 | A. F. Timofeev - Integration of Functions (1948)                       |
| :Welz       | 11 | Martin Welz - posts on Sci.Math.Symbolic                               |
| :Webster    | 12 | Michael Wester                                                         |

## Citation

If you use **SymbolicNumericIntegration.jl**, please cite [Symbolic-Numeric Integration of Univariate Expressions based on Sparse Regression](https://arxiv.org/abs/2201.12468):

```
@article{Iravanian2022,   
   author = {Shahriar Iravanian and Carl Julius Martensen and Alessandro Cheli and Shashi Gowda and Anand Jain and Julia Computing and Yingbo Ma and Chris Rackauckas},
   doi = {10.48550/arxiv.2201.12468},
   month = {1},
   title = {Symbolic-Numeric Integration of Univariate Expressions based on Sparse Regression},
   url = {https://arxiv.org/abs/2201.12468v2},
   year = {2022},
}
```

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the 
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Manifest.toml"
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Project.toml"
```

```@raw html
">project</a> file.
```
