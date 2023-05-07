@syms 𝑥
@syms u[20]

transformer(eq) = transformer(ops(eq)...)

function transformer(::Mul, eq)
    return vcat([transformer(t) for t in arguments(eq)]...)
end

function transformer(::Div, eq)
	a = transformer(arguments(eq)[1])
	b = transformer(arguments(eq)[2])
	b = [(1/q, k) for (q, k) in b]
    return [a; b]
end

function transformer(::Pow, eq)
    y, k = arguments(eq)
    if is_number(k)
    	r = nice_parameter(k)
    	if denominator(r) == 1
	    	return [(y, k)]
	    else
	    	return [(y^(1/denominator(r)), numerator(r))]
	    end
    else
    	return [(eq, 1)]
	end
end

function transformer(::Any, eq)
    return [(eq, 1)]    
end

function transform(eq, x)
    eq = substitute(eq, Dict(x => 𝑥))
    p = transformer(eq)
    p = p[isdependent.(first.(p),  𝑥)]
    
    return p
end

function rename_factors(p)
	n = length(p)

	q = 1
	sub = Dict()
	ks = Int[]
	
	for (i,(y,k)) in enumerate(p)
		μ = u[i]
		q *= μ ^ k
		sub[μ] = y
		push!(ks, k)
	end
	
	return q, sub, ks
end

##############################################################################

Symbolics.@register_symbolic Ei(z)
Symbolics.@register_symbolic Si(z)
Symbolics.@register_symbolic Ci(z)
Symbolics.@register_symbolic Li(z)

Symbolics.derivative(::typeof(Ei), args::NTuple{1, Any}, ::Val{1}) = exp(args[1]) / args[1]
Symbolics.derivative(::typeof(Si), args::NTuple{1, Any}, ::Val{1}) = sin(args[1]) / args[1]
Symbolics.derivative(::typeof(Ci), args::NTuple{1, Any}, ::Val{1}) = cos(args[1]) / args[1]
Symbolics.derivative(::typeof(Li), args::NTuple{1, Any}, ::Val{1}) = 1 / log(args[1])

@syms si(𝑥) ci(𝑥) ei(𝑥) li(𝑥)

##############################################################################

function substitute_x(eq, x, sub)
    eq = substitute(eq, sub)
    return substitute(eq, Dict(𝑥 => x))
end

guard_zero(x) = isequal(x, 0) ? one(x) : x

function generate_homotopy(eq, x)
    eq = eq isa Num ? eq.val : eq
    x = x isa Num ? x.val : x

	p = transform(eq, x)
    q, sub, ks = rename_factors(p)
    S = 0

    for i in 1:length(sub)
		μ = u[i]
		h₁, ∂h₁ = apply_partial_int_rules(sub[μ])
		h₁ = substitute(h₁, Dict(si => Si, ci => Ci, ei => Ei, li => Li))		    
	    h₁ = substitute_x(h₁, x, sub)
		
    	for j = 1:ks[i]
		    h₂ = substitute_x((q / μ^j) / ∂h₁, x, sub)
		    S += expand((1 + h₁) * guard_zero(1 + h₂))
		end
    end    
    
    ζ = [x^k for k=1:maximum(ks)+1]
    
    unique([one(x); ζ; [equivalent(t, x) for t in terms(S)]])
end

##############################################################################

function ∂(x)
    d = expand_derivatives(Differential(𝑥)(x))
    return isequal(d, 0) ? 1 : d
end

@syms 𝛷(x)

partial_int_rules = [
                     # trigonometric functions
                     @rule 𝛷(sin(~x)) => (cos(~x) + si(~x), ∂(~x))
                     @rule 𝛷(cos(~x)) => (sin(~x) + ci(~x), ∂(~x))
                     @rule 𝛷(tan(~x)) => (log(cos(~x)), ∂(~x))
                     @rule 𝛷(csc(~x)) => (log(csc(~x) + cot(~x)), ∂(~x))
                     @rule 𝛷(sec(~x)) => (log(sec(~x) + tan(~x)), ∂(~x))
                     @rule 𝛷(cot(~x)) => (log(sin(~x)), ∂(~x))
                     # hyperbolic functions
                     @rule 𝛷(sinh(~x)) => (cosh(~x), ∂(~x))
                     @rule 𝛷(cosh(~x)) => (sinh(~x), ∂(~x))
                     @rule 𝛷(tanh(~x)) => (log(cosh(~x)), ∂(~x))
                     @rule 𝛷(csch(~x)) => (log(tanh(~x / 2)), ∂(~x))
                     @rule 𝛷(sech(~x)) => (atan(sinh(~x)), ∂(~x))
                     @rule 𝛷(coth(~x)) => (log(sinh(~x)), ∂(~x))
                     # 1/trigonometric functions
                     @rule 𝛷(1 / sin(~x)) => (log(csc(~x) + cot(~x)) + log(sin(~x)), ∂(~x))
                     @rule 𝛷(1 / cos(~x)) => (log(sec(~x) + tan(~x)) + log(cos(~x)), ∂(~x))
                     @rule 𝛷(1 / tan(~x)) => (log(sin(~x)) + log(tan(~x)), ∂(~x))
                     @rule 𝛷(1 / csc(~x)) => (cos(~x) + log(csc(~x)), ∂(~x))
                     @rule 𝛷(1 / sec(~x)) => (sin(~x) + log(sec(~x)), ∂(~x))
                     @rule 𝛷(1 / cot(~x)) => (log(cos(~x)) + log(cot(~x)), ∂(~x))
                     # 1/hyperbolic functions
                     @rule 𝛷(1 / sinh(~x)) => (log(tanh(~x / 2)) + log(sinh(~x)), ∂(~x))
                     @rule 𝛷(1 / cosh(~x)) => (atan(sinh(~x)) + log(cosh(~x)), ∂(~x))
                     @rule 𝛷(1 / tanh(~x)) => (log(sinh(~x)) + log(tanh(~x)), ∂(~x))
                     @rule 𝛷(1 / csch(~x)) => (cosh(~x) + log(csch(~x)), ∂(~x))
                     @rule 𝛷(1 / sech(~x)) => (sinh(~x) + log(sech(~x)), ∂(~x))
                     @rule 𝛷(1 / coth(~x)) => (log(cosh(~x)) + log(coth(~x)), ∂(~x))
                     # inverse trigonometric functions
                     @rule 𝛷(asin(~x)) => (~x * asin(~x) + sqrt(1 - ~x * ~x), ∂(~x))
                     @rule 𝛷(acos(~x)) => (~x * acos(~x) + sqrt(1 - ~x * ~x), ∂(~x))
                     @rule 𝛷(atan(~x)) => (~x * atan(~x) + log(~x * ~x + 1), ∂(~x))
                     @rule 𝛷(acsc(~x)) => (~x * acsc(~x) + atanh(1 - ^(~x, -2)), ∂(~x))
                     @rule 𝛷(asec(~x)) => (~x * asec(~x) + acosh(~x), ∂(~x))
                     @rule 𝛷(acot(~x)) => (~x * acot(~x) + log(~x * ~x + 1), ∂(~x))
                     # inverse hyperbolic functions
                     @rule 𝛷(asinh(~x)) => (~x * asinh(~x) + sqrt(~x * ~x + 1), ∂(~x))
                     @rule 𝛷(acosh(~x)) => (~x * acosh(~x) + sqrt(~x * ~x - 1), ∂(~x))
                     @rule 𝛷(atanh(~x)) => (~x * atanh(~x) + log(~x + 1), ∂(~x))
                     @rule 𝛷(acsch(~x)) => (acsch(~x), ∂(~x))
                     @rule 𝛷(asech(~x)) => (asech(~x), ∂(~x))
                     @rule 𝛷(acoth(~x)) => (~x * acot(~x) + log(~x + 1), ∂(~x))
                     # logarithmic and exponential functions
                     @rule 𝛷(log(~x)) => (~x + ~x * log(~x) +
                                          sum(candidate_pow_minus(~x, -1); init = one(~x)),
                                          ∂(~x))
                     @rule 𝛷(1 / log(~x)) => (log(log(~x)) + li(~x), ∂(~x))
                     @rule 𝛷(exp(~x)) => (exp(~x) + ei(~x), ∂(~x))
                     @rule 𝛷(^(exp(~x), ~k::is_neg)) => (^(exp(-~x), -~k), ∂(~x))
                     # square-root functions
                     @rule 𝛷(^(~x, ~k::is_abs_half)) => (sum(candidate_sqrt(~x, ~k);
                                                             init = one(~x)), 1);
                     @rule 𝛷(sqrt(~x)) => (sum(candidate_sqrt(~x, 0.5); init = one(~x)), ∂(~x));
                     @rule 𝛷(1 / sqrt(~x)) => (sum(candidate_sqrt(~x, -0.5); init = one(~x)), ∂(~x));
                     # rational functions                                                              
                     @rule 𝛷(1 / ^(~x::is_poly, ~k::is_pos_int)) => (sum(candidate_pow_minus(~x, -~k);
                                                                 init = one(~x)), 1)
                     @rule 𝛷(1 / ~x::is_poly) => (sum(candidate_pow_minus(~x, -1);
                                                                 init = one(~x)), 1)
                     @rule 𝛷(^(~x, -1)) => (log(~x), ∂(~x))
                     @rule 𝛷(^(~x, ~k::is_neg_int)) => (sum(^(~x, i) for i in (~k + 1):-1),
                                                        ∂(~x))
                     @rule 𝛷(1 / ~x) => (log(~x), ∂(~x))
                     @rule 𝛷(^(~x, ~k::is_pos_int)) => (sum(^(~x, i+1) for i=1:~k+1), ∂(~x))
                     @rule 𝛷(1) => (𝑥, 1)
                     @rule 𝛷(~x) => ((~x + ^(~x, 2)), ∂(~x))]

function apply_partial_int_rules(eq)
    expand(Fixpoint(Prewalk(Chain(partial_int_rules))))(𝛷(value(eq)))
end
