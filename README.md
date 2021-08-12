# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a symbolic and semi-numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

Function `integrate` returns the integral of an expression with *constant* coefficients. It uses a randomized algorithm based on a hybrid of the *method of undetermined coefficients* and *sparse regression* and is able to solve a large subset of basic standard integrals (polynomials, exponential/logarithmic, trigonometric and hyperbolic, inverse trigonometric and hyperbolic, rational and squarer root).

`integrate` returns a tuple with three values. The first one is the solved integral, the second one is the sum of the unsolved terms, and the third value is the residual error. If `integrate` is successful, the unsolved portion is reported as 0.

```julia
using SymbolicUtils
using SymbolicIntegration

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

`test/runtests.jl` contains a test suite of hundreds of expressions (function `test_integrals`).

`integrate` has the form `integrate(y; kw...)` or `integrate(y, x; kw...)`, where `y` is the integrand and the optional `x` is the variable of integration. The keyword parameters are:

* `abstol` (default `1e-6`): the error tolerance to accept a solution.
* `symbolic` (default `true`): if true, symbolic integration is attempted first.
* `bypass` (default `false`): if true, the whole expression is considered at once and not per term.
* `bypart` (default `true`): if true, integration by parts is tried.
* `num_steps` (default `2`): the number of different basis sets checked.
* `num_trials` (default `10`): the number of attempts for each basis set.
* `show_basis` (default `false`): print the basis set, useful for debugging.
* `verbose` (default `false`): if true, prints extra debugging information.
* `radius` (default `1.0`): the starting radius to generate random test points.
* `opt` (default `STLSQ(exp.(-10:1:0))`): the optimizer passed to `sparse_regression!`.
* `attempt_ratio` (default `5`): the maximum number of random points generated.
* `max_basis` (default `110`): the maximum number of expression in the basis.
