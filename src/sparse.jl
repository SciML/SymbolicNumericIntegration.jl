

function solve_sparse(eq, x, basis; plan = default_plan())    
    abstol = plan.abstol

    A, X = init_basis_matrix(eq, x, basis; plan)

    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, basis = A[l, l], basis[l]

    y₁, ε₁ = sparse_fit(A, basis; plan)
    if ε₁ < abstol
        return y₁, ε₁, basis
    end

    rank = sum(l)

    if rank == 1
        y₂, ε₂ = find_singlet(A, basis; abstol)
        if ε₂ < abstol
            return y₂, ε₂
        end
    elseif rank < 8
        y₃, ε₃ = find_dense(A, basis; abstol)
        if ε₃ < abstol
            return y₃, ε₃
        end
    end

    # moving toward the poles
    modify_basis_matrix!(A, X, eq, x, basis)
    y₄, ε₄ = sparse_fit(A, basis; plan)

    if ε₄ < abstol || ε₄ < ε₁
        return y₄, ε₄
    else
        return y₁, ε₁
    end
end

function prune_basis(eq, x, basis; plan = default_plan())
    A, X = init_basis_matrix(eq, x, basis; plan)
    l = find_independent_subset(A; plan.abstol)
    return basis[l]
end

function init_basis_matrix(eq, x, basis; plan = default_plan())
    n = length(basis)
    eq = value(eq)
    basis = value.(basis)

    # A is an nxn matrix holding the values of the fragments at n random points
    A = zeros(Complex{Float64}, (n, n))
    X = zeros(Complex{Float64}, n)

    S = subs_symbols(eq, x)
    if !isempty(S)
        eq = substitute(eq, S)
        basis = [substitute(y, S) for y in basis]
    end

    eq_fun = fun!(eq, x)
    Δbasis_fun = deriv_fun!.(basis, x)

    k = 1
    l = 10 * n# max attempt

    while k <= n && l > 0
        try
            x₀ = test_point(plan.complex_plane, plan.radius)
            X[k] = x₀
            b₀ = eq_fun(X[k])

            if is_proper(b₀)
                for j in 1:n
                    A[k, j] = Δbasis_fun[j](X[k]) / b₀
                end
                if all(is_proper, A[k, :])
                    k += 1
                end
            end
        catch e
              println("Error from init_basis_matrix!: ", e)
        end
        l -= 1
    end

    return A, X
end

function modify_basis_matrix!(A, X, eq, x, basis)
    n, m = size(A)
    eq_fun = fun!(eq, x)
    Δeq_fun = deriv_fun!(eq, x)
    Δbasis_fun = deriv_fun!.(basis, x)

    for k in 1:n
        # One Newton iteration toward the poles
        # note the + sign instead of the usual - in Newton-Raphson's method. This is
        # because we are moving toward the poles and not zeros.

        X[k] += eq_fun(X[k]) / Δeq_fun(X[k])
        b₀ = eq_fun(X[k])
        for j in 1:m
            A[k, j] = Δbasis_fun[j](X[k]) / b₀
        end
    end
end

function DataDrivenSparse.active_set!(idx::BitMatrix, p::SoftThreshold,
    x::Matrix{ComplexF64}, λ::Float64)
    DataDrivenSparse.active_set!(idx, p, abs.(x), λ)
end

function sparse_fit(A, basis; plan = default_plan())
    n, m = size(A)

    try
        b = ones((1, n))
        solver = SparseLinearSolver(plan.opt,
            options = DataDrivenCommonOptions(verbose = false,
                maxiters = 1000))
        res, _... = solver(A', b)
        q₀ = DataDrivenSparse.coef(first(res))

        ε = rms(A * q₀' .- 1)
        q = nice_parameter.(q₀)
        if sum(iscomplex.(q)) > 2
            return nothing, Inf
        end   # eliminating complex coefficients
        return sum(q[i] * expr(basis[i]) for i in 1:length(basis) if q[i] != 0;
            init = 0),
        abs(ε)
    catch e
        println("Error from sparse_fit", e)
        return nothing, Inf
    end
end


function find_singlet(A, basis; abstol = 1e-6)
    σ = vec(std(A; dims = 1))
    μ = vec(mean(A; dims = 1))
    l = (σ .< abstol) .* (abs.(μ) .> abstol)
    if sum(l) == 1
        k = findfirst(l)
        return nice_parameter(1 / μ[k]) * expr(basis[k]), σ[k]
    else
        return nothing, Inf
    end
end


function find_dense(A, basis; abstol = 1e-6)
    n = size(A, 1)
    b = ones(n)

    try
        q = A \ b
        if minimum(abs.(q)) > abstol
            ε = maximum(abs.(A * q .- b))
            if ε < abstol
                y = sum(nice_parameter.(q) .* expr.(basis))
                return y, ε
            end
        end
    catch e
        #
    end
    return nothing, Inf
end

###############################################################################

function hints(eq, x, basis; plan = default_plan())    
    abstol = plan.abstol
    A, X = init_basis_matrix(eq, x, basis; plan)
    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, basis = A[l, l], basis[l]

    n, m = size(A)

    try
        b = ones((1, n))
        solver = SparseLinearSolver(plan.opt,
            options = DataDrivenCommonOptions(verbose = false,
                maxiters = 1000))
        res, _... = solver(A', b)
        q = DataDrivenSparse.coef(first(res))
        err = abs(rms(A * q' .- 1))
        if err < abstol
            sel = abs.(q) .> abstol
            h = [basis[i] for i in 1:length(basis) if sel[i]]
        else
            h = []
        end
        return h, err
    catch e
        # println("Error from hints: ", e)        
    end

    return 0, Inf
end

function best_hints(eq, x, basis; plan = default_plan(), num_trials = 10)
    H = []
    L = Int[]

    for _ in 1:num_trials
        try
            h, err = hints(eq, x, basis; plan)
            push!(H, h)
            push!(L, err < plan.abstol ? length(h) : length(basis))
        catch e
            #
        end
    end

    if isempty(L)
        return nothing
    end

    _, idx = findmin(L)
    return H[idx]
end
