########################## Predicates #########################################

is_pos_int(x) = isinteger(x) && x > 0
is_pos(x) = x > 0
is_neg(x) = x < 0
is_neg_one(x) = (x ≈ -1)
is_pos_half(x) = (x ≈ 0.5)
is_neg_half(x) = (x ≈ -0.5)
is_abs_half(x) = (x ≈ -0.5 || x ≈ -0.5)

###############################################################################

function U(u...)
    u = map(x -> x isa AbstractArray ? x : [], u)
    return union(u...)
end

hints(eq::SymbolicUtils.Add, x, h, k) = map(t->hints(t,x,h,k), arguments(eq))
hints(eq::SymbolicUtils.Mul, x, h, k) = map(t->hints(t,x,h,k), arguments(eq))
hints(eq::SymbolicUtils.Pow, x, h, k) = hints(arguments(eq)[1],x,h,k)

function hints(eq::SymbolicUtils.Term, x, h, k)
    s = Symbol(operation(eq))
    u = arguments(eq)[1]

    if !isequal(u, x)
        _, _, islinear = Symbolics.linear_expansion(u, x)
        if islinear
            push!(k, u)
        end
    end

    if s == :sec
        push!(h, log(1/cos(u) + sin(u)/cos(u)))
    elseif s == :csc
        push!(h, log(1/sin(u) - cos(u)/sin(u)))
    elseif s == :tan
        push!(h, log(cos(u)))
    elseif s == :cot
        push!(h, log(sin(u)))
    elseif s == :tanh
        push!(h, log(cosh(u)))
    # elseif s == :log
    #     if check_poly(u, x) == :real_poly && deg(u, x) == 2
    #         r, s = find_roots(u, x)
    #         if !isempty(r)
    #             push!(h, log(x - r[1]))
    #             push!(h, log(x - r[2]))
    #         else
    #             push!(h, log(u))
    #         end
    #     end
    end
end

function hints(eq, x, h, k)
end

function collect_hints(eq, x)
    h = []  # hints
    k = []  # kernels
    hints(eq, x, h, k)
    h, k
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
    @rule Ω(^(~x, ~k)) => ~k >= 0 ? ~k * Ω(~x) : NaN
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

n_rules = [
    @rule Ω(+(~~xs)) => sum(map(Ω,~~xs))
    @rule Ω(*(~~xs)) => union(map(Ω,~~xs))
    @rule Ω(^(~x, ~k)) => Ω(~x)
    @rule Ω((~f)(~x::is_linear_poly)) => [~x]
    @rule Ω(~x) => []
]

inners(eq) = Fixpoint(Prewalk(PassThrough(Chain(n_rules))))(Ω(value(eq)))

coef(eq::SymbolicUtils.Mul, x) = prod(t for t in arguments(eq) if !isdependent(t,x); init=1)
coef(eq::SymbolicUtils.Add, x) = minimum(abs(coef(t,x)) for t in arguments(eq))
coef(eq, x) = 1

########################## Transformation Rules ###############################

is_int_gt_one(x) = isinteger(x) && x > 1

trig_rule1 = @rule tan(~x) => sin(~x) / cos(~x)
trig_rule2 = @rule sec(~x) => one(~x) / cos(~x)
# trig_rule2 = @rule sec(~x) => (tan(~x)/cos(~x) + 1/cos(~x)^2) / (tan(~x) + 1/cos(~x))
trig_rule3 = @rule csc(~x) => one(~x) / sin(~x)
# trig_rule3 = @rule csc(~x) => (cot(~x)/sin(~x) + 1/sin(~x)^2) / (cot(~x) + 1/sin(~x))
trig_rule4 = @rule cot(~x) => cos(~x) / sin(~x)
trig_rule5 = @rule sin(~n::is_int_gt_one*~x) => sin((~n - 1)*~x)*cos(~x) + cos((~n - 1)*~x)*sin(~x)
trig_rule6 = @rule cos(~n::is_int_gt_one*~x) => cos((~n - 1)*~x)*cos(~x) - sin((~n - 1)*~x)*sin(~x)
trig_rule7 = @rule tan(~n::is_int_gt_one*~x) => (tan((~n - 1)*~x) + tan(~x)) / (1 - tan((~n - 1)*~x)*tan(~x))
trig_rule8 = @rule csc(~n::is_int_gt_one*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (sec((~n - 1)*~x)*csc(~x) + csc((~n - 1)*~x)*sec(~x))
trig_rule9 = @rule sec(~n::is_int_gt_one*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (csc((~n - 1)*~x)*csc(~x) - sec((~n - 1)*~x)*sec(~x))
trig_rule10 = @rule cot(~n::is_int_gt_one*~x) => (cot((~n - 1)*~x)*cot(~x) - 1) / (cot((~n - 1)*~x) + cot(~x))

trig_rule11 = @rule sin(~x + ~y) => sin(~x)*cos(~y) + cos(~x)*sin(~y)
trig_rule12 = @rule cos(~x + ~y) => cos(~x)*cos(~y) - sin(~x)*sin(~y)
trig_rule13 = @rule tan(~x + ~y) => (tan(~x) + tan(~y)) / (1 - tan(~x)*tan(~y))
trig_rule14 = @rule csc(~x + ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (sec(~x)*csc(~y) + csc(~x)*sec(~y))
trig_rule15 = @rule sec(~x + ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (csc(~x)*csc(~y) - sec(~x)*sec(~y))
trig_rule16 = @rule cot(~x + ~y) => (cot(~x)*cot(~y) - 1) / (cot(~x) + cot(~y))

trig_rule17 = @rule sin(~x - ~y) => sin(~x)*cos(~y) - cos(~x)*sin(~y)
trig_rule18 = @rule cos(~x - ~y) => cos(~x)*cos(~y) + sin(~x)*sin(~y)
trig_rule19 = @rule tan(~x - ~y) => (tan(~x) - tan(~y)) / (1 + tan(~x)*tan(~y))
trig_rule20 = @rule csc(~x - ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (sec(~x)*csc(~y) - csc(~x)*sec(~y))
trig_rule21 = @rule sec(~x - ~y) => sec(~x)*sec(~y)*csc(~x)*csc(~y) / (csc(~x)*csc(~y) + sec(~x)*sec(~y))
trig_rule22 = @rule cot(~x - ~y) => (cot(~x)*cot(~y) + 1) / (cot(~x) - cot(~y))

# trig_rule5 = @rule sin(2*~x) => 2*sin(~x)*cos(~x)
# trig_rule6 = @rule cos(2*~x) => 2*cos(~x)^2 - 1
# trig_rule7 = @rule tan(2*~x) => 2*tan(~x) / (1 - tan(~x)^2)
# trig_rule8 = @rule cot(2*~x) => (cot(~x)^2 - 1) / (2*cot(~x))
# trig_rule9 = @rule sec(2*~x) => sec(~x)^2 / (2 - sec(~x)^2)
# trig_rule10 = @rule csc(2*~x) => sec(~x)*csc(~x) / 2
#
# trig_rule11 = @rule sin(3*~x) => 3*sin(~x) - 4*sin(~x)^3
# trig_rule12 = @rule cos(3*~x) => 4*cos(~x)^3 - 3*cos(~x)
# trig_rule13 = @rule tan(3*~x) => (3*tan(~x) - tan(~x)^3) / (1 - 3*tan(~x)^2)
# trig_rule14 = @rule cot(3*~x) => (3*cot(~x) - cot(~x)^3) / (1 - 3*cot(~x)^2)
# trig_rule15 = @rule sec(3*~x) => sec(~x)^3 / (4 - 3*sec(~x)^2)
# trig_rule16 = @rule csc(3*~x) => csc(~x)^3 / (3*csc(~x)^2 - 4)

trig_rules = [trig_rule1, trig_rule2, trig_rule3, trig_rule4,
              trig_rule5, trig_rule6, trig_rule7, trig_rule8, trig_rule9, trig_rule10,
              trig_rule11, trig_rule12, trig_rule13, trig_rule14, trig_rule15, trig_rule16,
              trig_rule17, trig_rule18, trig_rule19, trig_rule20, trig_rule21, trig_rule22]

hyper_rule1 = @rule tanh(~x) => sinh(~x) / cosh(~x)
hyper_rule2 = @rule sech(~x) => one(~x) / cosh(~x)
hyper_rule3 = @rule csch(~x) => one(~x) / sinh(~x)
hyper_rule4 = @rule coth(~x) => cosh(~x) / sinh(~x)

hyper_rules = [hyper_rule1, hyper_rule2, hyper_rule3, hyper_rule4]

misc_rule1 = @rule sqrt(~x) => ^(~x, 0.5)
misc_rule2 = @acrule exp(~x) * exp(~y) => exp(~x + ~y)

misc_rules = [misc_rule1]

int_rules = [trig_rules; hyper_rules; misc_rules]
# int_rules = misc_rules

apply_integration_rules(eq) = expand(Fixpoint(Prewalk(PassThrough(Chain(int_rules))))(value(eq)))
