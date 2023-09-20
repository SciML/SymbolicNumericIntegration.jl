# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://docs.sciml.ai/Symbolics/stable/) expressions.

**SymbolicNumericIntegration.jl** uses a randomized algorithm based on a hybrid of the *method of undetermined coefficients* and *sparse regression* and can solve a large subset of basic standard integrals (polynomials, exponential/logarithmic, trigonometric and hyperbolic, inverse trigonometric and hyperbolic, rational and square root). The symbolic part of the algorithm is similar (but not identical) to Risch-Bronstein's poor man's integrator and generates a list of ansatzes (candidate terms). The numerical part uses sparse regression adopted from the Sparse identification of nonlinear dynamics (SINDy) algorithm to prune down the ansatzes and find the corresponding coefficients. The basis of how it works and the theory of integration using the Symbolic-Numeric methods refer to [Basis of Symbolic-Numeric Integration](https://github.com/SciML/SymbolicNumericIntegration.jl/blob/main/docs/theory.ipynb).

[hyint](https://github.com/siravan/hyint) is the python counterpart of **SymbolicNumericIntegration.jl** that works with **sympy** expressions.

Originally, **SymbolicNumericIntegration.jl** expected only univariate expression with *constant* real or complex coefficients as input. As of version 1.2, `integrate` function accepts symbolic constants for a subset of solvable integrals.

`integrate` returns a tuple with three values. The first one is the solved integral, the second one is the sum of the unsolved terms, and the third value is the residual error. If `integrate` is successful, the unsolved portion is reported as 0. If we pass `detailed=false` to `integrate', the output is simplified to only the resulting integrals. In this case, `nothing` is returned if the integral cannot be solved. Note that the simplified output will become the default in a future version.

## Installation

To install SymbolicNumericIntegration.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("SymbolicNumericIntegration")
```

Examples:

```julia
julia> using SymbolicNumericIntegration
using Symbolics

julia> @variables x a b

julia> integrate(3x^3 + 2x - 5)

julia> integrate((5 + 2x)^-1)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)

julia> integrate(x^2 / (16 + x^2); detailed = false)
((1//2)*log((5//2) + x), 0, 0.0)

# `detailed = false` simplifies the output to just the resulting integral

julia> integrate(x^2 * log(x); detailed = false)
x + 4atan((-1//4)*x)

julia> integrate(sec(x) * tan(x); detailed = false)
(1//3)*(x^3)*log(x) - (1//9)*(x^3)

julia> integrate(sin(a * x), x; detailed = false, symbolic = true)
sec(x)

# Here, a is a symbolic constant; therefore, we need to explicitly
# define the independent variable (say, x). Also, we set
# `symbolic = true` to force using the symbolic solver

julia> integrate(x^2 * cos(a * x), x; detailed = false, symbolic = true)
(-cos(a*x)) / a

julia> integrate(cosh(a * x) * exp(b * x), x; detailed = false, symbolic = true)
((a^2)*(x^2)*sin(a*x) + 2.0a*x*cos(a*x) - 2.0sin(a*x)) / (a^3)

# multiple symbolic constants

julia> integrate(log(log(a * x)) / x, x; detailed = false, symbolic = true)
(a*sinh(a*x)*exp(b*x) - b*cosh(a*x)*exp(b*x)) / (a^2 - (b^2))

julia> integrate(x * sin(a * x), (x, 0, 1); symbolic = true, detailed = false)
log(a*x)*log(log(a*x)) - log(a*x)

# definite integration, passing a tuple of (x, lower bound, higher bound) in the 
# second argument
```

SymbolicNumericIntegration.jl exports some special integral functions (defined over Complex numbers) and uses them in solving integrals:

  - `Ei`: exponential integral (define as ∫ exp(x) / x dx)
  - `Si`: sine integral (define as ∫ sin(x) / x dx)
  - `Ci`: cosine integral (define as ∫ cos(x) / x dx)
  - `Li`: logarithmic integral (define as ∫ 1 / log(x) dx)

For examples:

```
julia> integrate(exp(x + 1) / (x + 1))
(SymbolicNumericIntegration.Ei(1 + x), 0, 0.0)

julia> integrate(x * cos(a*x^2 - 1) / (a*x^2 - 1), x; detailed=false, symbolic=true)
((1//2)*SymbolicNumericIntegration.Ci(a*(x^2) - 1)) / a

julia> integrate(1 / (x*log(log(x))), x; detailed=false, symbolic=true)
SymbolicNumericIntegration.Li(log(x))
```

```@docs
integrate(eq, x; kwargs...)
```

## Testing

`test/runtests.jl` contains a test suite of 170 easy to moderate test integrals (can be run by calling `test_integrals`). Currently, **SymbolicNumericIntegration.jl** solves more than 95% of its test suite.

Additionally, 12 test suites from the *Rule-based Integrator* ([Rubi](https://rulebasedintegration.org/)) are included in the `/test` directory. For example, we can test the first one as below. *Axiom* refers to the format of the test files)

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
