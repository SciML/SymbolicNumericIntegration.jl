module SymbolicIntegration

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters
using Symbolics

using DataDrivenDiffEq.Optimize

include("utils.jl")

export deg, check_poly

include("roots.jl")

export solve_newton, find_roots, find_poles

include("integral.jl")

export integrate, generate_basis, test_integrals

include("symbolic.jl")
include("integration_by_parts.jl")

end # module
