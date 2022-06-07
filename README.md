# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://symbolicnumericintegration.sciml.ai/stable/). Use the
[in-development documentation](https://symbolicnumericintegration.sciml.ai/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

```julia
using Symbolics
using SymbolicNumericIntegration

@variables x

julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)
```
