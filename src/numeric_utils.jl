# A list of miscellaneous numerical utility functions

rms(x) = sqrt(sum(x .^ 2) / length(x))

# returns a list of the indices of a linearly independent
# subset of the columns of A
function find_independent_subset(A; abstol = 1.0e-3)
    Q, R = qr(A)
    return abs.(diag(R)) .> abstol
end

function test_point(complex_plane, radius)
    if complex_plane
        return radius * sqrt(rand()) * cis(2π * rand())
    else
        return Complex(radius * (2 * rand() - 1))
    end
end

accept_solution(eq::ExprCache, x, sol, radius) = accept_solution(expr(eq), x, sol, radius)

function accept_solution(eq, x, sol; plan = default_plan())
    try
        # x₀ = test_point(plan.complex_plane, plan.radius)
        # Δ = substitute(diff(sol, x) - expr(eq), Dict(x => x₀))
        S = subs_symbols(eq, x; include_x = true, plan.radius)
        # Also substitute any symbols in sol that may not be in eq
        # Use value() to ensure consistent comparison with x (which may be BasicSymbolic)
        x_val = value(x)
        for v in get_variables(value(sol))
            if !haskey(S, v) && !isequal(v, x_val)
                S[v] = Complex(randn())
            end
        end
        Δ = substitute(diff(sol, x) - expr(eq), S)

        # Check if Δ is zero - handle both Num and BasicSymbolic cases
        # For Num, isequal works; for BasicSymbolic, we need to extract the value
        if isequal(Δ, 0)
            return 0.0
        end
        # Try to extract numeric value for BasicSymbolic zero
        try
            Δ_val = Symbolics.value(Num(Δ))
            if Δ_val isa Real || Δ_val isa Complex
                return abs(Δ_val)
            end
        catch
        end

        result = abs(Δ)

        # Note: Num <: Number, so we must check for Num BEFORE Number
        # First try to extract numeric value from Num type
        if result isa Num
            inner = Symbolics.unwrap(result)
            # Check if unwrapped value is a concrete number (not symbolic)
            if inner isa Real || inner isa Complex
                return abs(inner)
            end
            # Check if inner is structurally zero (e.g., abs(0))
            if isequal(inner, 0) || (
                    SymbolicUtils.iscall(inner) && SymbolicUtils.operation(inner) === abs &&
                        let arg = SymbolicUtils.arguments(inner)[1]
                        isequal(arg, 0) || (arg isa Real && arg == 0)
                    end
                )
                return 0.0
            end
            # Try value extraction for symbolic wrapper (SymbolicUtils v4+)
            try
                val = Symbolics.value(result)
                if val isa Real || val isa Complex
                    return abs(val)
                end
            catch
            end
            # If still symbolic, return Inf
            return Inf
        end
        # For non-Num concrete numbers
        if result isa Real || result isa Complex
            return abs(result)
        end
        # Try to extract numeric value from symbolic wrapper (SymbolicUtils v4+)
        try
            val = Symbolics.value(Num(result))
            if val isa Real || val isa Complex
                return abs(val)
            end
        catch
        end
        # Fallback: try unwrap
        try
            val = Symbolics.unwrap(result)
            if val isa Real || val isa Complex
                return abs(val)
            end
        catch
        end
        return Inf
    catch e
        #
    end
    return Inf
end

# converts float to int or small rational numbers
function nice_parameter(u::Complex{T}; abstol = 1.0e-6) where {T <: Real}
    α = nice_parameter(real(u); abstol)
    β = nice_parameter(imag(u); abstol)
    return β ≈ 0 ? α : Complex(α, β)
end

nice_parameter(x; abstol = 1.0e-6) = small_rational(x; abstol)

nice_parameters(p; abstol = 1.0e-6) = nice_parameter.(p; abstol)

function small_rational(x::T; abstol = 1.0e-6, M = 20) where {T <: Real}
    a = floor(Int, x)
    r = x - a
    for den in 2:M
        try
            if abs(round(r * den) - r * den) < abstol
                s = a + round(Int, r * den) // den
                return (denominator(s) == 1 ? numerator(s) : s)
            end
        catch e
        end
    end
    return x
end
