

function solve_sparse(T, eq, x, basis, radius; kwargs...)
    args = Dict(kwargs)
    abstol, opt, complex_plane, verbose = args[:abstol], args[:opt], args[:complex_plane],
                                          args[:verbose]
    n = length(basis)

    # A is an nxn matrix holding the values of the fragments at n random points
    A = zeros(Complex{T}, (n, n))
    X = zeros(Complex{T}, n)

    init_basis_matrix!(T, A, X, x, eq, basis, radius, complex_plane; abstol)

    y₁, ϵ₁ = sparse_fit(T, A, x, basis, opt; abstol)
    if ϵ₁ < abstol
        return y₁, ϵ₁
    end

    y₂, ϵ₂ = find_singlet(T, A, basis; abstol)
    if ϵ₂ < abstol
        return y₂, ϵ₂
    end

    if n < 8    # here, 8 is arbitrary and signifies a small basis
        y₃, ϵ₃ = find_dense(T, A, basis; abstol)
        if ϵ₃ < abstol
            return y₃, ϵ₃
        end
    end

    # moving toward the poles
    modify_basis_matrix!(T, A, X, x, eq, basis, radius; abstol)
    y₄, ϵ₄ = sparse_fit(T, A, x, basis, opt; abstol)

    if ϵ₄ < abstol || ϵ₄ < ϵ₁
        return y₄, ϵ₄
    else
        return y₁, ϵ₁
    end
end

function init_basis_matrix!(T, A, X, x, eq, basis, radius, complex_plane; abstol = 1e-6)
    n, m = size(A)
    k = 1
    i = 1

    eq_fun = fun!(eq, x)
    Δbasis_fun = deriv_fun!.(basis, x)

    while k <= n
        try
            x₀ = test_point(complex_plane, radius)
            X[k] = x₀ # move_toward_roots_poles(x₀, x, eq)
            b₀ = eq_fun(X[k])

            if is_proper(b₀)
                for j in 1:m
                    A[k, j] = Δbasis_fun[j](X[k]) / b₀
                end
                if all(is_proper, A[k, :])
                    k += 1
                end
            end
        catch e
            println("Error from init_basis_matrix!: ", e)
        end
    end
end

function move_toward_roots_poles(z, x, eq; n = 1, max_r = 100.0)
    eq_fun = fun!(eq, x)
    Δeq_fun = deriv_fun!(eq, x)
    is_root = rand() < 0.5
    z₀ = z
    for i in 1:n
        dz = eq_fun(z) / Δeq_fun(z)
        if is_root
            z -= dz
        else
            z += dz
        end
        if abs(z) > max_r
            return z₀
        end
    end
    return z
end

function modify_basis_matrix!(T, A, X, x, eq, basis, radius; abstol = 1e-6)
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


function sparse_fit(T, A, x, basis, opt; abstol = 1e-6)
    n = length(basis)
    # find a linearly independent subset of the basis
    l = find_independent_subset(A; abstol)
    A, basis, n = A[l, l], basis[l], sum(l)

    try
        b = ones(n)
        # q₀ = A \ b
        q₀ = DataDrivenDiffEq.init(opt, A, b)
        @views sparse_regression!(q₀, A, permutedims(b)', opt, maxiter = 1000)
        ϵ = rms(A * q₀ - b)
        q = nice_parameter.(q₀)
        if sum(iscomplex.(q)) > 2
            return nothing, Inf
        end   # eliminating complex coefficients
        return sum(q[i] * expr(basis[i]) for i in 1:length(basis) if q[i] != 0;
                   init = zero(x)),
               abs(ϵ)
    catch e
        println("Error from sparse_fit", e)
        return nothing, Inf
    end
end

function find_singlet(T, A, basis; abstol)
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

function find_dense(T, A, basis; abstol = 1e-6)
    n = size(A, 1)
    b = ones(T, n)

    try
        q = A \ b
        if minimum(abs.(q)) > abstol
            ϵ = maximum(abs.(A * q .- b))
            if ϵ < abstol
                y = sum(nice_parameter.(q) .* expr.(basis))
                return y, ϵ
            end
        end
    catch e
        #
    end
    return nothing, Inf
end
