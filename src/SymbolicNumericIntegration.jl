module SymbolicNumericIntegration

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters
using Symbolics

using DataDrivenDiffEq

include("utils.jl")
include("roots.jl")

export find_roots

include("rules.jl")
include("candidates.jl")
include("homotopy.jl")

include("numeric_utils.jl")
include("integral.jl")

export integrate, generate_basis

include("symbolic.jl")
include("logger.jl")

end # module
