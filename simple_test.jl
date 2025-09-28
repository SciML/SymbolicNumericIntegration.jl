using Pkg
Pkg.activate(".")

# Force load the package
using SymbolicNumericIntegration
using Symbolics

println("Testing vector expression error handling...")

# Test case 1: Vector with symbolic variables
@variables x
try
    result = integrate([x])
    println("❌ Test 1 FAILED: Should have thrown an error but got: $result")
catch e
    if occursin("Vector expressions are not supported", string(e))
        println("✅ Test 1 PASSED: Vector input correctly rejected")
    else
        println("❌ Test 1 FAILED: Wrong error message: $(e)")
    end
end

# Test case 2: Scalar integration should still work
try
    result = integrate(x)
    println("✅ Test 2 PASSED: Scalar integration works: $result")
catch e
    println("❌ Test 2 FAILED: Scalar integration should work but got: $(e)")
end

println("\nAll tests completed!")