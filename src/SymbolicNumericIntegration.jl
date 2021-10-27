# module SymbolicNumericIntegration

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters
using Symbolics

using DataDrivenDiffEq.Optimize

include("utils.jl")
include("roots.jl")

export find_roots

include("rules.jl")
include("candidates.jl")
include("integral.jl")

export integrate, generate_basis

include("symbolic.jl")
include("integration_by_parts.jl")

include("logger.jl")
include("prune.jl")

# end # module
