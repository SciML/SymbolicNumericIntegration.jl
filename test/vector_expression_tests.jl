using SymbolicNumericIntegration
using Symbolics
using Test

@variables x α

# Test that vector expressions throw an appropriate error
@test_throws ErrorException("Vector expressions are not supported. Please use element-wise integration with `integrate.([expr1, expr2, ...], x)` instead.") integrate([x])
@test_throws ErrorException("Vector expressions are not supported. Please use element-wise integration with `integrate.([expr1, expr2, ...], x)` instead.") integrate(
    [1, 2 * α], α
)

# Test that scalar integration still works
@test integrate(x) == ((1 // 2) * (x^2), 0, 0)
@test integrate(2 * α, α) == (α^2, 0, 0)

# Test that element-wise integration works
results = integrate.([1, 2 * α], α)
@test length(results) == 2
@test results[1] == (α, 0, 0)
@test results[2] == (α^2, 0, 0)
