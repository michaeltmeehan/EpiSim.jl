using SparseArrays
using LinearAlgebra

function construct_sir_band_matrices(N)
    I = spdiagm(0 => ones(N+1))  # identity matrix

    # S+inf: superdiag(1:N)
    S_plus_inf = spdiagm(1 => 1:N)

    # S−inf: diag(0:N)
    S_minus_inf = spdiagm(0 => 0:N)

    # I+inf: subdiag(0:N-1)
    I_plus_inf = spdiagm(-1 => 0:N-1)

    # I−inf: diag(0:N-1, 0)
    I_minus_inf = spdiagm(0 => [0:N-1; 0])  # pad with 0 for size N+1

    # S+rec and S−rec are just identity
    S_plus_rec = I
    S_minus_rec = I

    # I+rec: superdiag(1:N)
    I_plus_rec = spdiagm(1 => 1:N)

    # I−rec: diag(0:N)
    I_minus_rec = spdiagm(0 => 0:N)

    return (
        S_plus_inf, S_minus_inf,
        I_plus_inf, I_minus_inf,
        S_plus_rec, S_minus_rec,
        I_plus_rec, I_minus_rec
    )
end


function Q_mul!(out::Vector{Float64}, v::Vector{Float64}, N::Int, α::Float64, β::Float64, Q_bands::Tuple{SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Int64, Int64}, SparseMatrixCSC{Int64, Int64}})
    # Reshape vec(v) into matrix V
    V = reshape(v, N+1, N+1)  # columns = infected, rows = susceptible

    # Preallocate output matrix
    QV = zeros(N+1, N+1)

    # Load banded matrices
    S⁺ᵢ, S⁻ᵢ, I⁺ᵢ, I⁻ᵢ, S⁺ʳ, S⁻ʳ, I⁺ʳ, I⁻ʳ = Q_bands

    # Infection inflow: (β/N) * (S⁺ ⊗ I⁺)
    QV .+= (β/N) * (I⁺ᵢ * V * S⁺ᵢ')

    # Recovery inflow: α * (S⁺ ⊗ I⁺)
    QV .+= α * (I⁺ʳ * V * S⁺ʳ')

    # Infection outflow (negative diag terms)
    QV .-= (β/N) * (I⁻ᵢ * V * S⁻ᵢ')

    # Recovery outflow (negative diag terms)
    QV .-= α * (I⁻ʳ * V * S⁻ʳ')

    # Flatten result into output vector
    out .= vec(QV)
    return out
end


function uniformize(p0::Vector{Float64}, t::Float64, γ::Float64,
                    Q_mul!::Function, Q_args...; ε=1e-10, n_max=1000)

    p = zeros(length(p0))
    q = copy(p0)
    w = 1.0
    n = 0
    norm_deficit = 1.0

    while norm_deficit > ε && n < n_max
        p .+= exp(-γ * t) * w * q

        # Apply Q·q → temp, then P·q = q + (1/γ) * Q·q
        temp = similar(q)
        Q_mul!(temp, q, Q_args...)
        q .+= (1 / γ) .* temp

        # Update Poisson weight
        n += 1
        w *= γ * t / n

        norm_deficit = 1 - sum(p)  # L1 mass deficit
    end

    return p
end


β = 2.
α = 1.
N = 4
S0 = 3
I0 = 1
Q_bands = construct_sir_band_matrices(N)
p0 = zeros((N+1)^2)
p0[I0 + S0 * (N+1) + 1] = 1.0
y = similar(p0)
Q_mul!(y, p0, N, α, β, Q_bands)


γ = (N-1)*α + max((N-1)*β, α)
p1 = uniformize(p0, 1.0, γ, Q_mul!, N, α, β, Q_bands)
p2 = uniformize(p1, 1.0, γ, Q_mul!, N, α, β, Q_bands)
pcheck = uniformize(p0, 2.0, γ, Q_mul!, N, α, β, Q_bands)
[p0 p1 p2 pcheck]