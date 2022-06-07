# SymbolicNumericIntegration.jl

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

```julia
using Symbolics
using SymbolicNumericIntegration

@variables x

julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)
```
