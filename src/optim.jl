using Optim

function solve_optim(eq, x, basis, radius; kwargs...)
    args = Dict(kwargs)
    abstol, complex_plane, verbose = args[:abstol], args[:complex_plane], args[:verbose]

    n = length(basis)
    λ = 1e-9

    A, X = init_basis_matrix(eq, x, basis, radius, complex_plane; abstol)

    l = find_independent_subset(A; abstol)
    A, basis, n = A[:, l], basis[l], sum(l)
    p = rank_basis(A, basis)

    qm = zeros(n)
    ϵm = 1e6
    lm = qm .> 0

    for i in 1:n
        l = (1:n .<= i)
        q, ϵ = find_minimizer(A[:, l], λ)
        nz = sum(abs.(q) .> abstol)

        # println(i, '\t', ϵ, '\t', nz)

        if ϵ * nz < ϵm
            ϵm = ϵ * nz
            qm = q
            lm = l
        end
    end

    if ϵm < abstol
        return reconstruct(qm, basis[lm]), ϵm
    else
        return 0, ϵm
    end
end

function reconstruct(q, basis)
    q = nice_parameter.(q)
    y = sum(q[i] * expr(basis[i]) for i in 1:length(basis))
    return y
end

# returns a vector of indices of basis elems from the most important to the least
function rank_basis(A, basis)
    n, m = size(A)
    q = A \ ones(n)
    w = [abs(q[i]) * norm(A[:, i]) for i in 1:m]
    p = sortperm(w; rev = true)
    return p
end

clamp_down(x, η) = abs(x) < η ? 0 : x

function find_minimizer(A, λ)
    n, m = size(A)

    f = function (q)
        q .= clamp_down.(q, maximum(abs.(q)) * 1e-6)
        t = A * q .- 1
        l = sum(t' * t)
        l += λ * norm(q, 1)
        return l
    end

    g! = function (G, q)
        t = A * q .- 1
        G .= 2 * real(A' * t)
    end

    q0 = A \ ones(n) # randn(m)	
    res = Optim.optimize(f, g!, q0, LBFGS())
    q = Optim.minimizer(res)
    ϵ = Optim.minimum(res)

    return q, ϵ
end
