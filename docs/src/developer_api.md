# Developer API

## Extension Boundary

SymbolicNumericIntegration.jl does not expose an abstract extension interface.
Its `Op` and `IntegrationAlgorithm` hierarchies are implementation details that
depend on private expression and solver representations. External packages must
use the documented public functions instead of subtyping or dispatching on these
types.

```@docs
SymbolicNumericIntegration.Op
SymbolicNumericIntegration.IntegrationAlgorithm
SymbolicNumericIntegration.NumericalPlan
```
