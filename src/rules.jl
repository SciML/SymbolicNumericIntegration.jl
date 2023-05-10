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

terms(eq) = terms(ops(eq)...)
terms(::Add, eq) = [t for t in arguments(eq)]
terms(::Any, eq) = [eq]

factors(eq, x) = factors(ops(eq)..., x)
factors(::Any, eq, x) = isdependent(eq, x) ? [one(x), eq] : [one(x)]

function factors(::Pow, eq, x)
    p, k = arguments(eq)
    [p^(i * sign(k)) for i in 0:abs(k)]
end

function factors(::Mul, eq, x)
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

extract_power(::Pow, eq) = [arguments(eq)[2]]
extract_power(::Term, eq) = [1]
extract_power(::Mul, eq) = union([extract_power(t) for t in arguments(eq)]...)
extract_power(::Any, eq) = []
extract_power(eq) = extract_power(ops(eq)...)


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
    substitute(ω, Dict())
end

is_poly(eq) = !isnan(poly_deg(eq))
is_linear_poly(eq) = poly_deg(eq) == 1

# rules to extract the kernel (the solvable portion) of an expression
s_rules = [@rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
           @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
           # @rule Ω(~x / ~y) => Ω(~x) * Ω(^(~y, -1))
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

function coef(::Mul, eq, x)
    prod(t for t in arguments(eq) if !isdependent(t, x); init = 1)
end
coef(::Add, eq, x) = minimum(abs(coef(t, x)) for t in arguments(eq))
coef(::Any, eq, x) = is_number(eq) ? eq : 1
coef(eq, x) = coef(ops(eq)..., x)

equivalent(eq, x) = eq / coef(eq, x)

########################## Transformation Rules ###############################

trigs_rules = [@rule tan(~x) => sin(~x) / cos(~x)
               @rule sec(~x) => one(~x) / cos(~x)
               @rule csc(~x) => one(~x) / sin(~x)
               @rule cot(~x) => cos(~x) / sin(~x)
               @rule sin(~n::is_int_gt_one * ~x) => sin((~n - 1) * ~x) * cos(~x) +
                                                    cos((~n - 1) * ~x) * sin(~x)
               @rule cos(~n::is_int_gt_one * ~x) => cos((~n - 1) * ~x) * cos(~x) -
                                                    sin((~n - 1) * ~x) * sin(~x)
               @rule tan(~n::is_int_gt_one * ~x) => (tan((~n - 1) * ~x) + tan(~x)) / (1 - tan((~n - 1) * ~x) * tan(~x))
               @rule csc(~n::is_int_gt_one * ~x) => sec((~n - 1) * ~x) * sec(~x) *
                                                    csc((~n - 1) * ~x) * csc(~x) / (sec((~n - 1) * ~x) * csc(~x) + csc((~n - 1) * ~x) * sec(~x))
               @rule sec(~n::is_int_gt_one * ~x) => sec((~n - 1) * ~x) * sec(~x) *
                                                    csc((~n - 1) * ~x) * csc(~x) / (csc((~n - 1) * ~x) * csc(~x) - sec((~n - 1) * ~x) * sec(~x))
               @rule cot(~n::is_int_gt_one * ~x) => (cot((~n - 1) * ~x) * cot(~x) - 1) / (cot((~n - 1) * ~x) + cot(~x))

               @rule 1 / sin(~x) => csc(~x)
               @rule 1 / cos(~x) => sec(~x)
               @rule 1 / tan(~x) => cot(~x)
               @rule 1 / csc(~x) => sin(~x)
               @rule 1 / sec(~x) => cos(~x)
               @rule 1 / cot(~x) => tan(~x)

               @rule 1 / ^(sin(~x), ~k) => ^(csc(~x), ~k)
               @rule 1 / ^(cos(~x), ~k) => ^(sec(~x), ~k)
               @rule 1 / ^(tan(~x), ~k) => ^(cot(~x), ~k)
               @rule 1 / ^(csc(~x), ~k) => ^(sin(~x), ~k)
               @rule 1 / ^(sec(~x), ~k) => ^(cos(~x), ~k)
               @rule 1 / ^(cot(~x), ~k) => ^(tan(~x), ~k)
               
               @rule sin(~x + ~y) => sin(~x) * cos(~y) + cos(~x) * sin(~y)
               @rule cos(~x + ~y) => cos(~x) * cos(~y) - sin(~x) * sin(~y)
               @rule tan(~x + ~y) => (tan(~x) + tan(~y)) / (1 - tan(~x) * tan(~y))
               @rule csc(~x + ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) / (sec(~x) * csc(~y) + csc(~x) * sec(~y))
               @rule sec(~x + ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) / (csc(~x) * csc(~y) - sec(~x) * sec(~y))
               @rule cot(~x + ~y) => (cot(~x) * cot(~y) - 1) / (cot(~x) + cot(~y))
               @rule sin(~x - ~y) => sin(~x) * cos(~y) - cos(~x) * sin(~y)
               @rule cos(~x - ~y) => cos(~x) * cos(~y) + sin(~x) * sin(~y)
               @rule tan(~x - ~y) => (tan(~x) - tan(~y)) / (1 + tan(~x) * tan(~y))
               @rule csc(~x - ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) / (sec(~x) * csc(~y) - csc(~x) * sec(~y))
               @rule sec(~x - ~y) => sec(~x) * sec(~y) * csc(~x) * csc(~y) / (csc(~x) * csc(~y) + sec(~x) * sec(~y))
               @rule cot(~x - ~y) => (cot(~x) * cot(~y) + 1) / (cot(~x) - cot(~y))

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

               @rule tanh(~x) => sinh(~x) / cosh(~x)
               @rule sech(~x) => one(~x) / cosh(~x)
               @rule csch(~x) => one(~x) / sinh(~x)
               @rule coth(~x) => cosh(~x) / sinh(~x)
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
    p = sum(F[i] * q[i] for i in 1:n if q[i] != 0; init = zero(x))
    return p
end

rational_rules = [@rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
                  @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
                  @rule Ω(^(~x::is_poly, ~k::is_neg_int)) => expand(^(decompose_rational(~x),
                                                                      -~k))
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

inv_rules = [@rule Ω(1 / ~x) => ~x
             @rule Ω(~x / ~y) => Ω(~x) * ~y
             @rule Ω(^(~x, -1)) => ~x
             @rule Ω(^(~x, ~k)) => ^(Ω(~x), ~k)
             @rule Ω(~x) => 1 / ~x]

inverse(eq) = Prewalk(PassThrough(Chain(inv_rules)))(Ω(value(eq)))

###############################################################################

@syms ω

is_sym(x) = first(ops(x)) isa Sym

h_rules = [
			@rule +(~~xs) => ω + sum(~~xs)
			@rule *(~~xs) => ω + sum(~~xs)
			@rule ~x / ~y => ω + ~x + ~y
			@rule ^(~x, ~y) => ω + ~x + ~y
			@rule (~f)(~x) => ω + ~x
			@rule ~x::is_sym => ω
		  ]			

# complexity returns a measure of the complexity of an equation
# it is roughly similar ro kolmogorov complexity
function complexity(eq)
	_, eq = ops(eq)	
	h = Prewalk(PassThrough(Chain(h_rules)))(eq)
	return substitute(h, Dict(ω => 1))
end

