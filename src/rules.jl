########################## Predicates #########################################

is_pos_int(x) = is_proper(x) && isinteger(x) && x > 0
is_neg_int(x) = is_proper(x) && isinteger(x) && x < 0
is_int_gt_one(x) = is_proper(x) && isinteger(x) && x > 1
is_pos(x) = is_proper(x) && x > 0
is_neg(x) = is_proper(x) && x < 0
is_neg_one(x) = is_proper(x) && (x ≈ -1)
is_pos_half(x) = is_proper(x) && (x ≈ 0.5)
is_neg_half(x) = is_proper(x) && (x ≈ -0.5)
is_abs_half(x) = is_proper(x) && (x ≈ 0.5 || x ≈ -0.5)

###############################################################################

terms(eq::SymbolicUtils.Add) = [t for t in arguments(eq)]
terms(eq) = [eq]

factors(eq, x) = isdependent(eq, x) ? [one(x), eq] : [one(x)]

function factors(eq::SymbolicUtils.Pow, x)
    p, k = arguments(eq)
    return [p^(i * sign(k)) for i in 0:abs(k)]
end

function factors(eq::SymbolicUtils.Mul, x)
    terms = [factors(q, x) for q in arguments(eq)]
    n = length(terms)

    l = Any[one(x)]

    for j in 1:n
        m = length(l)
        for t in terms[j]
            for k in 1:m
                push!(l, l[k] * t)
            end
        end
    end

    return unique(l)
end

extract_power(eq::SymbolicUtils.Pow) = [arguments(eq)[2]]
extract_power(eq::SymbolicUtils.Term) = [1]
extract_power(eq::SymbolicUtils.Mul) = union([extract_power(t) for t in arguments(eq)]...)
extract_power(eq) = []

###############################################################################

is_var(eq) = isequal(eq, var(eq))

@syms Ω(x)

# rules to calculate the degree of a polynomial expression or NaN if not a poly
p_rules = [@rule Ω(+(~~xs)) => max(map(Ω, ~~xs)...)
           @rule Ω(~x * ~y) => Ω(~x) + Ω(~y)
           @rule Ω(^(~x, ~k::is_proper)) => is_pos_int(~k) ? ~k * Ω(~x) : NaN
           @rule Ω(~x::is_number) => 0
           @rule Ω(~x::is_var) => 1
           @rule Ω(~x) => NaN]

function poly_deg(eq)
    ω = Prewalk(Chain(p_rules))(Ω(value(eq)))
    return substitute(ω, Dict())
end

is_poly(eq) = !isnan(poly_deg(eq))
is_linear_poly(eq) = poly_deg(eq) == 1

# rules to extract the kernel (the solvable portion) of an expression
s_rules = [@rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
           @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
           @rule Ω(~x / ~y) => Ω(~x) * Ω(^(~y, -1))
           @rule Ω(^(~x::is_linear_poly, ~k::is_pos_int)) => ^(~x, ~k)
           @rule Ω(~x::is_var) => ~x
           @rule Ω(~x::is_linear_poly) => ~x
           @rule Ω(sin(~x::is_linear_poly)) => sin(~x)
           @rule Ω(cos(~x::is_linear_poly)) => cos(~x)
           @rule Ω(sinh(~x::is_linear_poly)) => sinh(~x)
           @rule Ω(cosh(~x::is_linear_poly)) => cosh(~x)
           @rule Ω(exp(~x::is_linear_poly)) => exp(~x)
           @rule Ω(^(sin(~x::is_linear_poly), ~k::is_pos_int)) => ^(sin(~x), ~k)
           @rule Ω(^(cos(~x::is_linear_poly), ~k::is_pos_int)) => ^(cos(~x), ~k)
           @rule Ω(^(sinh(~x::is_linear_poly), ~k::is_pos_int)) => ^(sinh(~x), ~k)
           @rule Ω(^(cosh(~x::is_linear_poly), ~k::is_pos_int)) => ^(cosh(~x), ~k)
           @rule Ω(^(exp(~x::is_linear_poly), ~k::is_pos_int)) => ^(exp(~x), ~k)
           @rule Ω(^(~x, ~k)) => 1
           @rule Ω((~f)(~x)) => 1
           @rule Ω(~x) => 1]

kernel(eq) = Prewalk(Chain(s_rules))(Ω(value(eq)))

coef(eq::SymbolicUtils.Mul, x) = prod(t for t in arguments(eq) if !isdependent(t, x); init = 1)
coef(eq::SymbolicUtils.Add, x) = minimum(abs(coef(t, x)) for t in arguments(eq))
coef(eq, x) = is_number(eq) ? eq : 1

equivalent(eq, x) = eq / coef(eq, x)

########################## Transformation Rules ###############################

trigs_rules = [@rule tan(~x) => sin(~x) * cos(~x)^-1
               @rule sec(~x) => one(~x) * cos(~x)^-1
               @rule csc(~x) => one(~x) * sin(~x)^-1
               @rule cot(~x) => cos(~x) * sin(~x)^-1
               @rule sin(~n::is_int_gt_one * ~x) => sin((~n - 1) * ~x) * cos(~x) +
                                                    cos((~n - 1) * ~x) * sin(~x)
               @rule cos(~n::is_int_gt_one * ~x) => cos((~n - 1) * ~x) * cos(~x) -
                                                    sin((~n - 1) * ~x) * sin(~x)
               @rule tan(~n::is_int_gt_one * ~x) => (tan((~n - 1) * ~x) + tan(~x)) *
                                                    (1 - tan((~n - 1) * ~x) * tan(~x))^-1
               @rule csc(~n::is_int_gt_one * ~x) => sec((~n - 1) * ~x) * sec(~x) *
                                                    csc((~n - 1) * ~x) * csc(~x) *
                                                    (sec((~n - 1) * ~x) * csc(~x) +
                                                     csc((~n - 1) * ~x) * sec(~x))^-1
               @rule sec(~n::is_int_gt_one * ~x) => sec((~n - 1) * ~x) * sec(~x) *
                                                    csc((~n - 1) * ~x) * csc(~x) *
                                                    (csc((~n - 1) * ~x) * csc(~x) -
                                                     sec((~n - 1) * ~x) * sec(~x))^-1
               @rule cot(~n::is_int_gt_one * ~x) => (cot((~n - 1) * ~x) * cot(~x) - 1) *
                                                    (cot((~n - 1) * ~x) + cot(~x))^-1
               @rule ^(sin(~x), ~k::is_neg) => ^(csc(~x), -~k)
               @rule ^(cos(~x), ~k::is_neg) => ^(sec(~x), -~k)
               @rule ^(tan(~x), ~k::is_neg) => ^(cot(~x), -~k)
               @rule ^(csc(~x), ~k::is_neg) => ^(sin(~x), -~k)
               @rule ^(sec(~x), ~k::is_neg) => ^(cos(~x), -~k)
               @rule ^(cot(~x), ~k::is_neg) => ^(tan(~x), -~k)
               @rule sin(~x + ~y) => sin(~x) * cos(~y) + cos(~x) * sin(~y)
               @rule cos(~x + ~y) => cos(~x) * cos(~y) - sin(~x) * sin(~y)
               @rule tan(~x + ~y) => (tan(~x) + tan(~y)) * (1 - tan(~x) * tan(~y))^-1
               @rule csc(~x + ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) *
                                     (sec(~x) * csc(~y) + csc(~x) * sec(~y))^-1
               @rule sec(~x + ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) *
                                     (csc(~x) * csc(~y) - sec(~x) * sec(~y))^-1
               @rule cot(~x + ~y) => (cot(~x) * cot(~y) - 1) * (cot(~x) + cot(~y))^-1
               @rule sin(~x - ~y) => sin(~x) * cos(~y) - cos(~x) * sin(~y)
               @rule cos(~x - ~y) => cos(~x) * cos(~y) + sin(~x) * sin(~y)
               @rule tan(~x - ~y) => (tan(~x) - tan(~y)) * (1 + tan(~x) * tan(~y))^-1
               @rule csc(~x - ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) *
                                     (sec(~x) * csc(~y) - csc(~x) * sec(~y))^-1
               @rule sec(~x - ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) *
                                     (csc(~x) * csc(~y) + sec(~x) * sec(~y))^-1
               @rule cot(~x - ~y) => (cot(~x) * cot(~y) + 1) * (cot(~x) - cot(~y))^-1

               # @rule sin(2*~x) => 2*sin(~x)*cos(~x)
               # @rule cos(2*~x) => 2*cos(~x)^2 - 1
               # @rule tan(2*~x) => 2*tan(~x) / (1 - tan(~x)^2)
               # @rule cot(2*~x) => (cot(~x)^2 - 1) / (2*cot(~x))
               # @rule sec(2*~x) => sec(~x)^2 / (2 - sec(~x)^2)
               # @rule csc(2*~x) => sec(~x)*csc(~x) / 2

               # @rule sin(3*~x) => 3*sin(~x) - 4*sin(~x)^3
               # @rule cos(3*~x) => 4*cos(~x)^3 - 3*cos(~x)
               # @rule tan(3*~x) => (3*tan(~x) - tan(~x)^3) / (1 - 3*tan(~x)^2)
               # @rule cot(3*~x) => (3*cot(~x) - cot(~x)^3) / (1 - 3*cot(~x)^2)
               # @rule sec(3*~x) => sec(~x)^3 / (4 - 3*sec(~x)^2)
               # @rule csc(3*~x) => csc(~x)^3 / (3*csc(~x)^2 - 4)

               @rule tanh(~x) => sinh(~x) * cosh(~x)^-1
               @rule sech(~x) => one(~x) * cosh(~x)^-1
               @rule csch(~x) => one(~x) * sinh(~x)^-1
               @rule coth(~x) => cosh(~x) * sinh(~x)^-1
               @rule sqrt(~x) => ^(~x, 0.5)
               @acrule exp(~x) * exp(~y) => exp(~x + ~y)]

simplify_trigs(eq) = expand(Fixpoint(Prewalk(PassThrough(Chain(trigs_rules))))(value(eq)))

##############################################################################

function factor_poly(p)
    x = var(p)
    r, s = find_roots(p, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)
    return [[(x - u) for u in r];
            [(x^2 - 2 * real(u) * x + abs2(u)) for u in s]]
end

function decompose_rational(eq)
    if poly_deg(eq) == 1
        return inverse(eq)
    end
    x = var(eq)

    d = poly_deg(eq)

    for j in 1:10  # will try 10 times to find the roots
        r, s = find_roots(eq, x)
        if length(r) + length(s) >= d
            break
        end
    end

    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    F = Any[(x^2 - 2 * real(u) * x + abs2(u))^-1 for u in s] ∪
        [x * (x^2 - 2 * real(u) * x + abs2(u))^-1 for u in s]
    for i in eachindex(r)
        μ = sum(r[1:i] .== r[i])
        push!(F, (x - r[i])^-μ)
    end
    F = unique(F)

    n = length(F)
    A = zeros(Complex, (n, n))
    b = zeros(Complex, n)

    for i in 1:n
        x₀ = test_point(false, 1.0)
        d = Dict(x => x₀)
        b[i] = 1 / substitute(eq, d)
        for j in 1:n
            A[i, j] = substitute(F[j], d)
        end
    end

    q₀ = A \ b
    q = nice_parameter.(q₀)
    p = sum(F[i] * q[i] for i = 1:n if q[i] != 0; init = zero(x))
    return p
end

rational_rules = [@rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
                  @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
                  @rule Ω(^(~x::is_poly, ~k::is_neg_int)) => expand(^(decompose_rational(~x), -~k))
                  @rule Ω(~x / ^(~y::is_poly, ~k::is_pos_int)) => expand(~x *
                                                                         ^(decompose_rational(~y),
                                                                           ~k))
                  @rule Ω(~x / ~y::is_poly) => expand(~x * decompose_rational(~y))
                  # @rule Ω(^(~x,~k)) => ^(~x, ~k)
                  @rule Ω(sqrt(~x::is_poly)) => prod(sqrt(f) for f in factor_poly(~x))
                  @rule Ω(log(~x::is_poly)) => sum(log(f) for f in factor_poly(~x))
                  @rule Ω(~x) => ~x]

factor_rational(eq) = Prewalk(PassThrough(Chain(rational_rules)))(Ω(value(eq)))

###############################################################################

div_rule = @rule ~x / ~y => ~x * ^(~y, -1)

apply_div_rule(eq) = Prewalk(PassThrough(div_rule))(value(eq))

###############################################################################

inv_rules = [@rule Ω(1 / ~x) => ~x
             @rule Ω(~x / ~y) => Ω(~x) * ~y
             @rule Ω(^(~x, -1)) => ~x
             @rule Ω(^(~x, ~k)) => ^(Ω(~x), ~k)
             @rule Ω(~x) => ^(~x, -1)]

inverse(eq) = Prewalk(PassThrough(Chain(inv_rules)))(Ω(value(eq)))

###############################################################################

# the integration rules for non-solvable portion of an expression

@syms 𝛷(x)

d_rules = [@rule 𝛷(sin(~x)) => one(~x) + sin(~x) + cos(~x)
           @rule 𝛷(cos(~x)) => one(~x) + cos(~x) + sin(~x)
           @rule 𝛷(tan(~x)) => one(~x) + tan(~x) + log(cos(~x))
           @rule 𝛷(csc(~x)) => one(~x) + csc(~x) + log(sin(~x)^-1 - cos(~x) * sin(~x)^-1)
           @rule 𝛷(sec(~x)) => one(~x) + sec(~x) + log(cos(~x)^-1 + sin(~x) * cos(~x)^-1)
           @rule 𝛷(cot(~x)) => one(~x) + cot(~x) + log(sin(~x))
           @rule 𝛷(sinh(~x)) => one(~x) + sinh(~x) + cosh(~x)
           @rule 𝛷(cosh(~x)) => one(~x) + cosh(~x) + sinh(~x)
           @rule 𝛷(tanh(~x)) => one(~x) + tanh(~x) + log(cosh(~x))
           @rule 𝛷(csch(~x)) => one(~x) + csch(~x) + log(sinh(~x)^-1 - cosh(~x) * sinh(~x)^-1)
           @rule 𝛷(sech(~x)) => one(~x) + sech(~x) + log(cosh(~x)^-1 + sinh(~x) * cosh(~x)^-1)
           @rule 𝛷(coth(~x)) => one(~x) + coth(~x) + log(sinh(~x))
           @rule 𝛷(asin(~x)) => one(~x) + asin(~x) + ~x * asin(~x) + sqrt(1 - ~x * ~x)
           @rule 𝛷(acos(~x)) => one(~x) + acos(~x) + ~x * acos(~x) + sqrt(1 - ~x * ~x)
           @rule 𝛷(atan(~x)) => one(~x) + atan(~x) + ~x * atan(~x) + log(~x * ~x + 1)
           @rule 𝛷(acsc(~x)) => one(~x) + acsc(~x)
           @rule 𝛷(asec(~x)) => one(~x) + asec(~x)
           @rule 𝛷(acot(~x)) => one(~x) + acot(~x) + ~x * acot(~x) + log(~x * ~x + 1)
           @rule 𝛷(asinh(~x)) => one(~x) + asinh(~x) + ~x * asinh(~x) + sqrt(~x * ~x + 1)
           @rule 𝛷(acosh(~x)) => one(~x) + acosh(~x) + ~x * acosh(~x) + sqrt(~x * ~x - 1)
           @rule 𝛷(atanh(~x)) => one(~x) + atanh(~x) + ~x * atanh(~x) + log(~x + 1)
           @rule 𝛷(acsch(~x)) => one(~x) + acsch(~x)
           @rule 𝛷(asech(~x)) => one(~x) + asech(~x)
           @rule 𝛷(acoth(~x)) => one(~x) + acoth(~x) + ~x * acot(~x) + log(~x + 1)
           @rule 𝛷(log(~x)) => one(~x) + log(~x) + ~x + ~x * log(~x) + 𝛷(inverse(~x))
           @rule 𝛷(sqrt(~x)) => sum(candidate_sqrt(~x, 0.5); init = one(~x));
           @rule 𝛷(^(sqrt(~x), -1)) => 𝛷(^(~x, -0.5))
           @rule 𝛷(cbrt(~x)) => one(~x) + cbrt(~x) + ^(~x, 4 / 3)
           @rule 𝛷(^((~f)(~x), ~k::is_int_gt_one)) => var(~x) * ^((~f)(~x), ~k) +
                                                      𝛷(var(~x) * ^((~f)(~x), ~k - 1))
           @rule 𝛷(^(~x, ~k::is_abs_half)) => sum(candidate_sqrt(~x, ~k); init = 𝛷(~x));
           @rule 𝛷(^(~x, ~k::is_neg)) => sum(candidate_pow_minus(~x, ~k); init = 𝛷(~x));
           @rule 𝛷(^(~x::is_linear_poly, ~k::is_pos)) => one(~x) + ^(~x, ~k) + ^(~x, ~k + 1)
           @rule 𝛷(^(~x, ~k::is_pos)) => one(~x) + ^(~x, ~k) + ^(~x, ~k + 1)
           @rule 𝛷(exp(~x)) => one(~x) + exp(~x)
           @rule 𝛷(+(~~xs)) => sum(map(𝛷, ~~xs))
           @rule 𝛷(*(~~xs)) => prod(map(𝛷, ~~xs))
           @rule 𝛷(~x / ~y) => 𝛷(~x * inverse(~y))
           @rule 𝛷(~x) => one(~x) + ~x
           @rule 𝛷(-~x) => -𝛷(~x)]

apply_d_rules(eq) = expand(Fixpoint(Prewalk(PassThrough(Chain(d_rules))))(𝛷(value(eq))))
