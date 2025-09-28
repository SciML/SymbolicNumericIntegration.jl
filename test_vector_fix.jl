using Pkg
Pkg.activate(".")
using Symbolics, SymbolicNumericIntegration

println("Testing vector expression error handling:")
println()

# Test case from issue #58
@variables α
exp2 = [1, 2 * α]

println("Test 1: integrate([1, 2α], α)")
try
    result = integrate(exp2, α)
    println("ERROR: Should have thrown an error but got: $result")
catch e
    println("✓ Correctly caught error: $(e)")
end

println()
println("Test 2: Verify element-wise integration still works")
try
    result = integrate.(exp2, α)
    println("✓ Element-wise integration works: $result")
catch e
    println("ERROR: Element-wise should work but got: $(e)")
end

println()
# Test case from issue #106
@variables x
println("Test 3: integrate([x])")
try
    result = integrate([x])
    println("ERROR: Should have thrown an error but got: $result")
catch e
    println("✓ Correctly caught error: $(e)")
end

println()
println("Test 4: Verify scalar integration still works")
try
    result = integrate(2 * α, α)
    println("✓ Scalar integration works: $result")
catch e
    println("ERROR: Scalar should work but got: $(e)")
end