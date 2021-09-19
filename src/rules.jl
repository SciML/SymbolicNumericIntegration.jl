########################## Transformation Rules ###############################

trig_rule1 = @rule tan(~x) => sin(~x) / cos(~x)
trig_rule2 = @rule sec(~x) => one(~x) / cos(~x)
# trig_rule2 = @rule sec(~x) => (tan(~x)/cos(~x) + 1/cos(~x)^2) / (tan(~x) + 1/cos(~x))
trig_rule3 = @rule csc(~x) => one(~x) / sin(~x)
# trig_rule3 = @rule csc(~x) => (cot(~x)/sin(~x) + 1/sin(~x)^2) / (cot(~x) + 1/sin(~x))
trig_rule4 = @rule cot(~x) => cos(~x) / sin(~x)

is_pos_int(x::T) where T<:Integer = (x >= 2)
is_pos_int(x) = false

trig_rule5 = @rule sin(~n::is_pos_int*~x) => sin((~n - 1)*~x)*cos(~x) + cos((~n - 1)*~x)*sin(~x)
trig_rule6 = @rule cos(~n::is_pos_int*~x) => cos((~n - 1)*~x)*cos(~x) - sin((~n - 1)*~x)*sin(~x)
trig_rule7 = @rule tan(~n::is_pos_int*~x) => (tan((~n - 1)*~x) + tan(~x)) / (1 - tan((~n - 1)*~x)*tan(~x))
trig_rule8 = @rule csc(~n::is_pos_int*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (sec((~n - 1)*~x)*csc(~x) + csc((~n - 1)*~x)*sec(~x))
trig_rule9 = @rule sec(~n::is_pos_int*~x) => sec((~n - 1)*~x)*sec(~x)*csc((~n - 1)*~x)*csc(~x) / (csc((~n - 1)*~x)*csc(~x) - sec((~n - 1)*~x)*sec(~x))
trig_rule10 = @rule cot(~n::is_pos_int*~x) => (cot((~n - 1)*~x)*cot(~x) - 1) / (cot((~n - 1)*~x) + cot(~x))

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

########################## Expansion Rules ####################################

x_rule_g1 = @rule sin(~x) => (exp(im * ~x) - exp(-im * ~x)) / 2im
x_rule_g2 = @rule cos(~x) => (exp(im * ~x) + exp(-im * ~x)) / 2
x_rule_g3 = @rule tan(~x) => -im * (exp(2im * ~x) - 1) / (exp(2im * ~x) + 1)
x_rule_g4 = @rule csc(~x) => -2im / (exp(im * ~x) - exp(-im * ~x))
x_rule_g5 = @rule sec(~x) => 2 / (exp(im * ~x) + exp(im * ~x))
x_rule_g6 = @rule cot(~x) => im * (exp(2im * ~x) + 1) / (exp(2im * ~x) - 1)

x_rule_h1 = @rule sinh(~x) => (exp(~x) - exp(-~x)) / 2
x_rule_h2 = @rule cosh(~x) => (exp(~x) + exp(-~x)) / 2
x_rule_h3 = @rule tanh(~x) => (exp(2 * ~x) - 1) / (exp(2 * ~x) + 1)
x_rule_h4 = @rule csch(~x) => 2 / (exp(~x) - exp(-~x))
x_rule_h5 = @rule sech(~x) => 2 / (exp(~x) + exp(-~x))
x_rule_h6 = @rule coth(~x) => (exp(2 * ~x) + 1) / (exp(2 * ~x) - 1)

x_rule_i1 = @rule asin(~x) => -im * log(sqrt(1-x^2) + im*x)
x_rule_i2 = @rule acos(~x) => -im * log(x + im*sqrt(1-x^2))
x_rule_i3 = @rule atan(~x) => im/2 * log((x + im) / (-x + im))
x_rule_i4 = @rule acsc(~x) => -im * log(sqrt(1-1.0/x^2) + im/x)
x_rule_i5 = @rule asec(~x) => -im * log(1.0/x + im*sqrt(1-1.0/x^2))
x_rule_i6 = @rule acot(~x) => im/2 * log((1.0/x + im) / (-1.0/x + im))

x_rule_j1 = @rule asin(~x) => log(sqrt(1+x^2) + x)
x_rule_j2 = @rule acos(~x) => log(x + sqrt(1-x^2))
x_rule_j3 = @rule atan(~x) => 1/2 * log((x + 1) / (-x + 1))
x_rule_j4 = @rule acsc(~x) => log(sqrt(1+1.0/x^2) + 1.0/x)
x_rule_j5 = @rule asec(~x) => log(1.0/x + sqrt(1-1.0/x^2))
x_rule_j6 = @rule acot(~x) => 1/2 * log((1.0/x + 1) / (-1.0/x + 1))

x_rule_e1 = @acrule exp(~x) * exp(~y) => exp(~x + ~y)

expansion_rules = [
    x_rule_g1,
    x_rule_g2,
    x_rule_g3,
    x_rule_g4,
    x_rule_g5,
    x_rule_g6,

    x_rule_h1,
    x_rule_h2,
    x_rule_h3,
    x_rule_h4,
    x_rule_h5,
    x_rule_h6,

    x_rule_i1,
    x_rule_i2,
    x_rule_i3,
    x_rule_i4,
    x_rule_i5,
    x_rule_i6,

    x_rule_e1,
]

apply_expansion(eq) = Fixpoint(Prewalk(PassThrough(Chain(expansion_rules))))(value(eq))

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
    else s == :tanh
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
