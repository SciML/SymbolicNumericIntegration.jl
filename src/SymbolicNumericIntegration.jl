module SymbolicNumericIntegration

using ForwardDiff
using TermInterface: iscall
using SymbolicUtils
using SymbolicUtils: operation, arguments
using Symbolics
using Symbolics: value, get_variables, expand_derivatives, coeff, Equation
using SymbolicUtils.Rewriters
using SymbolicUtils: issym, BasicSymbolic

using SymbolicLimits: limit

using DataDrivenDiffEq, DataDrivenSparse

struct NumericalPlan
    abstol::Float64
    radius::Float64
    complex_plane::Bool
    opt::DataDrivenDiffEq.AbstractDataDrivenAlgorithm
end

default_plan() = NumericalPlan(1.0e-6, 5.0, true, STLSQ(exp.(-10:1:0)))

include("utils.jl")
include("tree.jl")
include("special.jl")

include("cache.jl")
include("roots.jl")

include("rules.jl")
include("candidates.jl")
include("homotopy.jl")

export Ei, Si, Ci, Li

include("numeric_utils.jl")
include("sparse.jl")
include("integral.jl")

export integrate, generate_basis, best_hints

include("symbolic.jl")

end # module
