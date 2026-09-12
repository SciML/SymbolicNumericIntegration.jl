using SymbolicNumericIntegration, BenchmarkTools
using Symbolics

const SUITE = BenchmarkGroup()

@variables x

integrate = SymbolicNumericIntegration.integrate

# =============================================================================
# Indefinite integration — representative integrands across the table:
# polynomial, rational, root, and transcendental forms
# =============================================================================

SUITE["indefinite"] = BenchmarkGroup()

integrands = Dict(
    "poly" => 4x^3,
    "rational" => 1 / (x^2 + x + 1),
    "rational_root" => x / sqrt(x^2 + 1),
    "exp" => x * exp(-x^2),
    "trig" => sin(x)^2 * cos(x),
)

for (name, eq) in integrands
    SUITE["indefinite"][name] = @benchmarkable $integrate(
        $eq, $x; symbolic = false, homotopy = true, num_steps = 2, num_trials = 4,
        detailed = false, verbose = false
    ) seconds = 120
end

# =============================================================================
# Symbolic-plan path and definite integration
# =============================================================================

SUITE["symbolic_definite"] = BenchmarkGroup()

SUITE["symbolic_definite"]["exp_neg_x"] = @benchmarkable $integrate(
    exp(-$x), ($x, 0, Inf); symbolic = true, detailed = false
) seconds = 120
SUITE["symbolic_definite"]["poly_01"] = @benchmarkable $integrate(
    $x, ($x, 0, 1); symbolic = false, detailed = false
) seconds = 120
