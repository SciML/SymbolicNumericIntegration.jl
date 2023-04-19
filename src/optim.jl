using Optim

function solve_optim(T, eq, x, basis, radius, rounds=2; kwargs...)
    args = Dict(kwargs)
    abstol, complex_plane, verbose = args[:abstol], args[:complex_plane], args[:verbose]
    n = length(basis)
    λ = 1e-9
    
    A = zeros(Complex{T}, (2n, n))
    X = zeros(Complex{T}, 2n)
    init_basis_matrix!(T, A, X, x, eq, basis, radius, complex_plane; abstol)    
    # modify_basis_matrix!(T, A, X, x, eq, basis, radius; abstol)
    
    l = find_independent_subset(A; abstol) 
    # l .&= rand(length(l)) .> 0.5   
    A, basis = A[:, l], basis[l]
   	q, ϵ = find_minimizer(A, λ)
   	
   	if ϵ > abstol
   		return 0, ϵ
   	end
   	
   	qa = q
   	μ = maximum(abs.(qa))
    
    for ρ in exp10.(-1:-1:-8)    
		l = abs.(qa) .> ρ * μ
		q, ϵ = find_minimizer(A[:, l], λ)
		if ϵ < abstol
			q = nice_parameter.(q)    
			basis = basis[l]
			y = sum(q[i] * expr(basis[i]) for i=1:length(basis))
			return y, ϵ
		end
    end
   	
   	return 0, ϵ
end


function find_minimizer(A, λ)
	n, m = size(A)
    B = real.(A' * A)
    b = real.(A' * ones(n))
	
	f = function(q)
		l = λ * sum(abs.(q))		# L1 norm
		# l = λ * sqrt(sum(q .^ 2))	# L2 norm
		l += sum(abs2.(A * q .- 1))
		return l
	end
	
	g! = function(G, q)		
		G .= 2 * (B * q .- b) .+ λ * sign.(q)
	end
	
	q0 = randn(m)	
	res = Optim.optimize(f, g!, q0, LBFGS())	
	q = Optim.minimizer(res)
	ϵ = Optim.minimum(res)
	
	return q, ϵ
end

