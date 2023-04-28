module SymbolicNumericIntegration

using SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics
using Symbolics: value, get_variables, expand_derivatives
using SymbolicUtils.Rewriters

using DataDrivenDiffEq

include("utils.jl")
include("special.jl")

include("cache.jl")
include("roots.jl")

include("rules.jl")
include("candidates.jl")
include("homotopy.jl")

export Ei, Si, Ci, Li

include("numeric_utils.jl")
include("sparse.jl")
include("optim.jl")
include("integral.jl")

export integrate, generate_basis

include("symbolic.jl")
include("logger.jl")

end # module
