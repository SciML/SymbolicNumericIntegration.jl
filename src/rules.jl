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

hints(eq::SymbolicUtils.Add, x, k) = map(t->hints(t,x,k), arguments(eq))
hints(eq::SymbolicUtils.Mul, x, k) = map(t->hints(t,x,k), arguments(eq))
hints(eq::SymbolicUtils.Pow, x, k) = hints(arguments(eq)[1],x,k)

function hints(eq::SymbolicUtils.Term, x, k)
    s = Symbol(operation(xq))
    u = arguments(eq)[1]

    if !isequal(u, x)
        _, _, islinear = Symbolics.linear_expansion(u, x)
        if islinear
            push!(k, u)
        end
    end
end

function collect_hints(eq, x)
    k = []  # kernels
    hints(eq, x, k)
    k
end

##############################################################################

# the integration rules for non-solvable portion of an expression

@syms 𝛷(x)

d_rule_g1 = @rule 𝛷(sin(~x)) => one(~x) + sin(~x) + cos(~x)
d_rule_g2 = @rule 𝛷(cos(~x)) => one(~x) + cos(~x) + sin(~x)
d_rule_g3 = @rule 𝛷(tan(~x)) => one(~x) + tan(~x) + log(cos(~x))
d_rule_g4 = @rule 𝛷(csc(~x)) => one(~x) + csc(~x) + log(1/sin(~x) - cos(~x)/sin(~x))
d_rule_g5 = @rule 𝛷(sec(~x)) => one(~x) + sec(~x) + log(1/cos(~x) + sin(~x)/cos(~x))
d_rule_g6 = @rule 𝛷(cot(~x)) => one(~x) + cot(~x) + log(sin(~x))

d_rule_h1 = @rule 𝛷(sinh(~x)) => one(~x) + sinh(~x) + cosh(~x)
d_rule_h2 = @rule 𝛷(cosh(~x)) => one(~x) + cosh(~x) + sinh(~x)
d_rule_h3 = @rule 𝛷(tanh(~x)) => one(~x) + tanh(~x) + log(cosh(~x))
d_rule_h4 = @rule 𝛷(csch(~x)) => one(~x) + csch(~x) + log(1/sinh(~x) - cosh(~x)/sinh(~x))
d_rule_h5 = @rule 𝛷(sech(~x)) => one(~x) + sech(~x) + log(1/cosh(~x) + sinh(~x)/cosh(~x))
d_rule_h6 = @rule 𝛷(coth(~x)) => one(~x) + coth(~x) + log(sinh(~x))

d_rule_i1 = @rule 𝛷(asin(~x)) => one(~x) + asin(~x) + ~x*asin(~x) + sqrt(1 - ~x*~x)
d_rule_i2 = @rule 𝛷(acos(~x)) => one(~x) + acos(~x) + ~x*acos(~x) + sqrt(1 - ~x*~x)
d_rule_i3 = @rule 𝛷(atan(~x)) => one(~x) + atan(~x) + ~x*atan(~x) + log(~x*~x + 1)
d_rule_i4 = @rule 𝛷(acsc(~x)) => one(~x) + acsc(~x)
d_rule_i5 = @rule 𝛷(asec(~x)) => one(~x) + asec(~x)
d_rule_i6 = @rule 𝛷(acot(~x)) => one(~x) + acot(~x) + ~x*acot(~x) + log(~x*~x + 1)

d_rule_j1 = @rule 𝛷(asinh(~x)) => one(~x) + asinh(~x) + ~x*asinh(~x) + sqrt(~x*~x + 1)
d_rule_j2 = @rule 𝛷(acosh(~x)) => one(~x) + acosh(~x) + ~x*acosh(~x) + sqrt(~x*~x - 1)
d_rule_j3 = @rule 𝛷(atanh(~x)) => one(~x) + atanh(~x) + ~x*atanh(~x) + log(~x + 1)
d_rule_j4 = @rule 𝛷(acsch(~x)) => one(~x) + acsch(~x)
d_rule_j5 = @rule 𝛷(asech(~x)) => one(~x) + asech(~x)
d_rule_j6 = @rule 𝛷(acoth(~x)) => one(~x) + acoth(~x) + ~x*acot(~x) + log(~x + 1)

d_rule_l1 = @rule 𝛷(log(~x)) => one(~x) + log(~x) + ~x + ~x * log(~x) + 𝛷(1/~x)
d_rule_l2 = @rule 𝛷(sqrt(~x)) => sum(candidate_sqrt(~x,0.5); init=one(~x))
d_rule_l3 = @rule 𝛷(^(sqrt(~x),-1)) => 𝛷(^(~x,-0.5))
d_rule_l4 = @rule 𝛷(cbrt(~x)) => one(~x) + cbrt(~x) + ^(~x, 4/3)

d_rule_p1 = @rule 𝛷(^(~x, ~k::is_abs_half)) => sum(candidate_sqrt(~x,~k); init=𝛷(~x))
d_rule_p2 = @rule 𝛷(^(~x, ~k::is_neg)) => sum(candidate_pow_minus(~x,~k); init=𝛷(~x))
d_rule_p3 = @rule 𝛷(^(~x, ~k::is_pos)) => one(~x) + ^(~x,~k+1) + 𝛷(~x)

d_rule_e1 = @rule 𝛷(exp(~x)) => one(~x) + exp(~x)
d_rule_e2 = @rule 𝛷(+(~~xs)) => sum(map(𝛷, ~~xs))
d_rule_e3 = @rule 𝛷(*(~~xs)) => prod(map(𝛷, ~~xs))
d_rule_e4 = @rule 𝛷(~x) => one(~x) + ~x
d_rule_e5 = @rule 𝛷(-~x) => -𝛷(~x)


d_rules = [
    d_rule_g1,
    d_rule_g2,
    d_rule_g3,
    d_rule_g4,
    d_rule_g5,
    d_rule_g6,

    d_rule_h1,
    d_rule_h2,
    d_rule_h3,
    d_rule_h4,
    d_rule_h5,
    d_rule_h6,

    d_rule_i1,
    d_rule_i2,
    d_rule_i3,
    d_rule_i4,
    d_rule_i5,
    d_rule_i6,

    d_rule_l1,
    d_rule_l2,
    d_rule_l3,
    # d_rule_l4,

    d_rule_p1,
    d_rule_p2,
    d_rule_p3,

    d_rule_e1,
    d_rule_e2,
    d_rule_e3,
    d_rule_e4,
    d_rule_e5,
]

apply_d_rules(eq) = expand(Fixpoint(Prewalk(PassThrough(Chain(d_rules))))(𝛷(value(eq))))


###############################################################################

is_var(eq) = isequal(eq, var(eq))

@syms Ω(x)

# rules to calculate the degree of a polynomial expression or NaN if not a poly
p_rules = [
    @rule Ω(+(~~xs)) => max(map(Ω, ~~xs)...)
    @rule Ω(~x * ~y) => Ω(~x) + Ω(~y)
    @rule Ω(^(~x, ~k::is_proper)) => is_pos_int(~k) ? ~k * Ω(~x) : NaN
    @rule Ω(~x::is_number) => 0
    @rule Ω(~x::is_var) => 1
    @rule Ω(~x) => NaN
]

function poly_deg(eq)
    ω = Prewalk(Chain(p_rules))(Ω(value(eq)))
    substitute(ω, Dict())
end

is_poly(eq) = !isnan(poly_deg(eq))
is_linear_poly(eq) = poly_deg(eq) == 1

# rules to extract the kernel (the solvable portion) of an expression
s_rules = [
    @rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
    @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
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
    @rule Ω(~x) => 1
]

kernel(eq) = Prewalk(Chain(s_rules))(Ω(value(eq)))

coef(eq::SymbolicUtils.Mul, x) = prod(t for t in arguments(eq) if !isdependent(t,x); init=1)
coef(eq::SymbolicUtils.Add, x) = minimum(abs(coef(t,x)) for t in arguments(eq))
coef(eq, x) = 1

########################## Transformation Rules ###############################

int_rules = [
    @rule tan(~x) => sin(~x) / cos(~x)
    @rule sec(~x) => one(~x) / cos(~x)
    @rule csc(~x) => one(~x) / sin(~x)
    @rule cot(~x) => cos(~x) / sin(~x)
    @rule sin(~n::is_int_gt_one*~x) => sin((~n - 1)*~x)*cos(~x) + cos((~n - 1)*~x)*sin(~x)
    @rule cos(~n::is_int_gt_one*~x) => cos((~n - 1)*~x)*cos(~x) - sin((~n - 1)*~x)*sin(~x)
    @rule tan(~n::is_int_gt_one*~x) => (tan((~n - 1)*~x) + tan(~x)) / (1 - tan((~n - 1)*~x)*tan(~x))
    @rule csc(~n::is_int_gt_one*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (sec((~n - 1)*~x)*csc(~x) + csc((~n - 1)*~x)*sec(~x))
    @rule sec(~n::is_int_gt_one*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (csc((~n - 1)*~x)*csc(~x) - sec((~n - 1)*~x)*sec(~x))
    @rule cot(~n::is_int_gt_one*~x) => (cot((~n - 1)*~x)*cot(~x) - 1) / (cot((~n - 1)*~x) + cot(~x))

    @rule ^(sin(~x), ~k::is_neg) => ^(csc(~x), -~k)
    @rule ^(cos(~x), ~k::is_neg) => ^(sec(~x), -~k)
    @rule ^(tan(~x), ~k::is_neg) => ^(cot(~x), -~k)
    @rule ^(csc(~x), ~k::is_neg) => ^(sin(~x), -~k)
    @rule ^(sec(~x), ~k::is_neg) => ^(cos(~x), -~k)
    @rule ^(cot(~x), ~k::is_neg) => ^(tan(~x), -~k)

    @rule sin(~x + ~y) => sin(~x)*cos(~y) + cos(~x)*sin(~y)
    @rule cos(~x + ~y) => cos(~x)*cos(~y) - sin(~x)*sin(~y)
    @rule tan(~x + ~y) => (tan(~x) + tan(~y)) / (1 - tan(~x)*tan(~y))
    @rule csc(~x + ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (sec(~x)*csc(~y) + csc(~x)*sec(~y))
    @rule sec(~x + ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (csc(~x)*csc(~y) - sec(~x)*sec(~y))
    @rule cot(~x + ~y) => (cot(~x)*cot(~y) - 1) / (cot(~x) + cot(~y))

    @rule sin(~x - ~y) => sin(~x)*cos(~y) - cos(~x)*sin(~y)
    @rule cos(~x - ~y) => cos(~x)*cos(~y) + sin(~x)*sin(~y)
    @rule tan(~x - ~y) => (tan(~x) - tan(~y)) / (1 + tan(~x)*tan(~y))
    @rule csc(~x - ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (sec(~x)*csc(~y) - csc(~x)*sec(~y))
    @rule sec(~x - ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (csc(~x)*csc(~y) + sec(~x)*sec(~y))
    @rule cot(~x - ~y) => (cot(~x)*cot(~y) + 1) / (cot(~x) - cot(~y))

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

    @rule sqrt(~x) => ^(~x, 0.5)
    @acrule exp(~x) * exp(~y) => exp(~x + ~y)
]

apply_integration_rules(eq) = expand(Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq)))

##############################################################################

function factor_poly(p)
    x = var(p)
    r, s = find_roots(p, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)
    return [[(x - u) for u in r];
            [(x^2 - 2*real(u)*x + abs2(u)) for u in s]
           ]
end

function decompose_rational(eq, k=1)
    if poly_deg(1/eq) == 1 return eq^k end
    x = var(eq)
    r, s = find_roots(1/eq, x)
    s = s[1:2:end]
    r = nice_parameter.(r)
    s = nice_parameter.(s)

    F = [1/(x^2 - 2*real(u)*x + abs2(u)) for u in s] ∪
        [x/(x^2 - 2*real(u)*x + abs2(u)) for u in s]
    for i in eachindex(r)
        μ = sum(r[1:i] .== r[i])
        push!(F, 1/(x-r[i])^μ)
    end
    F = unique(F)
    
    n = length(F)
    A = zeros(Complex, (n,n))
    b = zeros(Complex, n)

    for i = 1:n
        x₀ = test_point(false, 1.0)
        d = Dict(x => x₀)
        b[i] = substitute(eq, d)
        for j = 1:n
            A[i,j] = substitute(F[j], d)
        end
    end

    q₀ = A \ b
    q = nice_parameter.(q₀)
    p = sum(F[i]*q[i] for i=1:n if q[i] != 0; init=zero(x))
    return expand(p^k)
end

q_rules = [
    @rule Ω(+(~~xs)) => sum(map(Ω, ~~xs))
    @rule Ω(*(~~xs)) => prod(map(Ω, ~~xs))
    @rule Ω(^(~x::is_poly, ~k::is_neg_int)) => decompose_rational(1/~x, -~k)
    # @rule Ω(^(~x,~k)) => ^(~x, ~k)
    @rule Ω(sqrt(~x::is_poly)) => prod(sqrt(f) for f in factor_poly(~x))
    @rule Ω(log(~x::is_poly)) => sum(log(f) for f in factor_poly(~x))
    @rule Ω(~x) => ~x
]

apply_q_rules(eq) = Prewalk(PassThrough(Chain(q_rules)))(Ω(value(eq)))
