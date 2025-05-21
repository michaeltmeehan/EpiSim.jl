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


β = 2.
α = 1.
N = 4
S0 = 3
I0 = 1
Q_bands = construct_sir_band_matrices(N)
p0 = zeros((N+1)^2)
p0[I0 + S0 * (N+1) + 1] = 1.0
y = similar(p0)


γ = (N-1)*α + max((N-1)*β, α)

model = SIRModel(N=N, I=I0, transmission_rate=β, recovery_rate=α, sampling_rate=0.0)
rng = Random.MersenneTwister(1234)

ens = simulate(rng, model, 100_000, stop_condition = (state) -> state.t > 5.)

tvec = collect(0.:1.:5.)
empirical_prevalence = get_prevalence(ens, tvec)
n_max = maximum(empirical_prevalence)
empirical_pop_distribution = get_distribution(empirical_prevalence, n_max)
theoretical_pop_distribution = reduce(vcat, [aggregate_p(uniformize(p0, t, γ, Q_mul!, N, α, β, Q_bands), N)' for t in tvec])