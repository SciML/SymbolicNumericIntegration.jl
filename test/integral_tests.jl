using SymbolicNumericIntegration
using Symbolics
using Test

include("integral_setup.jl")

n = test_integrals(;
    symbolic = false, verbose = false, homotopy = true, num_steps = 2,
    num_trials = 10
)
@test n > 0
