# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

**SymbolicNumericIntegration.jl** uses a randomized algorithm based on a hybrid of the *method of undetermined coefficients* and *sparse regression* and is able to solve a large subset of basic standard integrals (polynomials, exponential/logarithmic, trigonometric and hyperbolic, inverse trigonometric and hyperbolic, rational and square root).
The basis of how it works and the theory of integration using the Symbolic-Numeric methods refer to [Basis of Symbolic-Numeric Integration](docs/theory.ipynb).

Function `integrate` returns the integral of a univariate expression with *constant* real or complex coefficients. `integrate` returns a tuple with three values. The first one is the solved integral, the second one is the sum of the unsolved terms, and the third value is the residual error. If `integrate` is successful, the unsolved portion is reported as 0.

```julia
using SymbolicUtils
using SymbolicNumericIntegration

@syms x

julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)

julia> integrate((5 + 2x)^-1)
((1//2)*log((5//2) + x), 0, 0.0)

julia> integrate(1 / (6 + x^2 - (5x)))
(log(x - 3) - log(x - 2), 0, 3.339372764128952e-16)

y = integrate(1 / (x^2 - 16))
((1//8)*log(x - 4) - ((1//8)*log(4 + x)), 0, 1.546926788028958e-16)

julia> integrate(x^2/(16 + x^2))
(x + 4atan((-1//4)*x), 0, 1.3318788420751984e-16)

julia> integrate(x^2/sqrt(4 + x^2))
((1//2)*x*((4 + x^2)^0.5) - ((2//1)*log(x + sqrt(4 + x^2))), 0, 8.702422633074313e-17)

julia> integrate(x^2*log(x))
((1//3)*log(x)*(x^3) - ((1//9)*(x^3)), 0, 0)

julia> integrate(x^2*exp(x))
(2exp(x) + exp(x)*(x^2) - (2x*exp(x)), 0, 0)

julia> integrate(tan(2x))
((-1//2)*log(cos(2x)), 0, 0)

julia> integrate(sec(x)*tan(x))
(cos(x)^-1, 0, 0)

julia> integrate(cosh(2x)*exp(x))
((2//3)*exp(x)*sinh(2x) - ((1//3)*exp(x)*cosh(2x)), 0, 7.073930088880992e-8)

julia> integrate(cosh(x)*sin(x))
((1//2)*sin(x)*sinh(x) - ((1//2)*cos(x)*cosh(x)), 0, 4.8956233716268386e-17)

julia> integrate(cosh(2x)*sin(3x))
(0.153845sinh(2x)*sin(3x) - (0.23077cosh(2x)*cos(3x)), 0, 4.9807620877373405e-6)

julia> integrate(log(log(x))*(x^-1))
(log(x)*log(log(x)) - log(x), 0, 0)

julia> integrate(exp(x^2))
(0, exp(x^2), Inf)    # as expected!
```

`integrate` has the form `integrate(y; kw...)` or `integrate(y, x; kw...)`, where `y` is the integrand and the optional `x` is the variable of integration. The keyword parameters are:

* `abstol` (default `1e-6`): the error tolerance to accept a solution.
* `symbolic` (default `true`): if true, pure symbolic integration is attempted first.
* `bypass` (default `false`): if true, the whole expression is considered at once and not per term.
* `num_steps` (default `2`): one plus the number of expanded basis to check (if `num_steps` is 1, only the main basis is checked).
* `num_trials` (default `5`): the number of attempts to solve the integration numerically for each basis set.
* `show_basis` (default `false`): print the basis set, useful for debugging. Only works if `verbose` is also set.
* `homotopy` (default: `true` as of version 0.7.0): uses the continuous Homotopy operators to generate the integration candidates.
* `verbose` (default `false`): if true, prints extra (and voluminous!) debugging information.
* `radius` (default `1.0`): the starting radius to generate random test points.
* `opt` (default `STLSQ(exp.(-10:1:0))`): the optimizer passed to `sparse_regression!`.
* `max_basis` (default `110`): the maximum number of expression in the basis.
* `complex_plane` (default `true`): random test points are generated on the complex plane (only over the real axis if `complex_plane` is `false`).

## Testing

`test/runtests.jl` contains a test suite of 160 easy to moderate test integrals (can be run by calling `test_integrals`). Currently, **SymbolicNumericIntegration.jl** solves more than 90% of its test suite.

Additionally, 12 test suites from the *Rule-based Integrator* ([Rubi](https://rulebasedintegration.org/)) are included in the `/test` directory. For example, we can test the first one as below ([Axiom](http://www.axiom-developer.org/) refers to the format of the test files)

```julia
  using SymbolicNumericIntegration
  include("test/axiom.jl")  # note, you may need to use the correct path

  L = convert_axiom(:Apostle)   # can also use L = convert_axiom(1)  
  test_axiom(L, false; bypass=false, verbose=false, homotopy=true)
```

The test suites description based on the header of the files in the Rubi directory are

| name        | id | comment                                  |
|-------------|----|------------------------------------------|
|:Apostle     | 1  | Tom M. Apostol - Calculus, Volume I, Second Edition (1967) |
|:Bondarenko  | 2  | Vladimir Bondarenko Integration Problems |
|:Bronstein   | 3  | Manuel Bronstein - Symbolic Integration Tutorial (1998) |
|:Charlwood   | 4  | Kevin Charlwood - Integration on Computer Algebra Systems (2008) |
|:Hearn       | 5  | Anthony Hearn - Reduce Integration Test Package |
|:Hebisch     | 6  | Waldek Hebisch - email May 2013 |
|:Jeffrey     | 7  | David Jeffrey - Rectifying Transformations for Trig Integration (1997) |
|:Moses       | 8  | Joel Moses - Symbolic Integration Ph.D. Thesis (1967) |
|:Stewart     | 9  | James Stewart - Calculus (1987) |
|:Timofeev    | 10 | A. F. Timofeev - Integration of Functions (1948) |
|:Welz        | 11 | Martin Welz - posts on Sci.Math.Symbolic |
|:Webster     | 12 | Michael Wester |
