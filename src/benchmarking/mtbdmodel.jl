using SparseArrays
using LinearAlgebra

function construct_mtbd_band_matrices(N)
    I = spdiagm(0 => ones(N+1))  # identity matrix
    A = spdiagm(1 => 1:N)
    B = spdiagm(0 => 0:N)
    C = spdiagm(-1 => 0:N-1)
    D = spdiagm(-1 => ones(N))
    E = spdiagm(-1 => 1:N)
    return I, A, B, C, D, E
end


function Q_mul!(out::Vector{Float64}, v::Vector{Float64}, N::Int, α::Vector{Float64}, λ::Matrix{Float64}, Q_bands::Tuple{SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Int64, Int64}})
    # Reshape vec(v) into matrix V
    V = reshape(v, N+1, N+1)  # columns = infected, rows = susceptible

    # Preallocate output matrix
    QV = zeros(N+1, N+1)

    # Load banded matrices
    I, A, B, C, D, E = Q_bands

    QV .+= α[1] * (I * V * A')  # Death inflow
    QV .-= α[1] * (I * V * B')  # Death outflow

    QV .+= α[2] * (A * V * I')  # Death inflow
    QV .-= α[2] * (B * V * I')  # Death outflow

    QV .+= λ[1, 1] * (I * V * C')  # Infection inflow
    QV .-= λ[1, 1] * (I * V * B')  # Infection outflow

    QV .+= λ[1, 2] * (B * V * D') # Infection inflow
    QV .-= λ[1, 2] * (B * V * I') # Infection outflow

    QV .+= λ[2, 1] * (D * V * B') # Infection inflow
    QV .-= λ[2, 1] * (I * V * B') # Infection outflow

    QV .+= λ[2, 2] * (C * V * I') # Infection inflow
    QV .-= λ[2, 2] * (B * V * I') # Infection outflow

    # Flatten result into output vector
    out .= vec(QV)
    return out
end


function uniformize_leaky(p0::Vector{Float64}, t::Float64, γ::Float64,
                          Q_mul!::Function, Q_args...; ε=1e-12, n_max=1000)

    p = zeros(length(p0))
    q = copy(p0)
    n = 0
    w = 1.0  # Poisson weight = (γt)^n / n!
    expfactor = exp(-γ * t)

    while w * expfactor > ε && n < n_max
        p .+= w * expfactor * q

        temp = similar(q)
        Q_mul!(temp, q, Q_args...)
        q .+= (1 / γ) .* temp

        n += 1
        w *= γ * t / n
    end

    return p
end

N = 3
α = [0.1, 0.2]
λ = [0.1 0.2; 0.3 0.4]
Q_bands = construct_mtbd_band_matrices(N)
p0 = zeros((N+1)^2)
γ = N * (sum(λ) + sum(α))
I10 = 1
I20 = 0
p0[I10 + I20 * (N+1) + 1] = 1.0

p1 = uniformize(p0, 1.0, γ, Q_mul!, N, α, λ, Q_bands)
p2 = uniformize(p1, 1.0, γ, Q_mul!, N, α, λ, Q_bands)
pcheck = uniformize(p0, 2.0, γ, Q_mul!, N, α, λ, Q_bands)
[p0 p1 p2 pcheck]