module SymbolicNumericIntegration

using DataDrivenDiffEq: DataDrivenCommonOptions
using DataDrivenSparse: STLSQ, SparseLinearSolver
using DataStructures: Queue
using LinearAlgebra: diag, qr
using SpecialFunctions: cosint, erfi, expint, sinint
using StatsAPI: coef
using Statistics: mean, std
using SymbolicLimits: limit
using SymbolicUtils: @acrule, @rule, @syms, BasicSymbolic, arguments, expand, issym, operation,
                   simplify, substitute
using SymbolicUtils.Code: toexpr
using SymbolicUtils.Rewriters: Chain, Fixpoint, PassThrough, Prewalk
import Symbolics
using Symbolics: @register_symbolic, @variables, Differential, Equation, Num, build_function,
                 expand_derivatives, get_variables, scalarize, unwrap, value
import Symbolics: derivative
using TermInterface: iscall

"""
    NumericalPlan

Internal configuration for the numerical verification stage of symbolic-numeric
integration. This type and its fields are implementation details, not an extension
interface or part of the public API.
"""
struct NumericalPlan
    abstol::Float64
    radius::Float64
    complex_plane::Bool
    opt
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
