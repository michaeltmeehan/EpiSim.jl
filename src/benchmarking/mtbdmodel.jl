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

N = 5
λ = [2. 2.; 2. 2.]
μ = [0.5, 0.5]
ψ = [0.5, 0.5]
α = μ .+ ψ
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

model = MTBDModel(birth_rate=birth_rate, death_rate=death_rate, sampling_rate=sampling_rate)
rng = Random.MersenneTwister(1234)

ens = simulate(rng, model, 100_000, stop_condition = (state) -> state.t > 0.5, max_iter=500)

tvec = collect(0.:0.1:0.5)
I1_prevalence = fill(0, length(tvec), 100_000)
I2_prevalence = fill(0, length(tvec), 100_000)
for (j, out) in enumerate(ens.replicates)
    for (i, t) in enumerate(tvec)
        idx = searchsortedlast(out.state_log.t, t)
        I1_prevalence[i, j] = out.state_log.I_1[idx]
        I2_prevalence[i, j] = out.state_log.I_2[idx]
    end
end

I1_max = maximum(I1_prevalence)
I1_empirical_distribution = get_distribution(I1_prevalence, I1_max)
reduce(vcat, [aggregate_p(uniformize(p0, t, γ, Q_mul!, N, α, λ, Q_bands), N)' for t in tvec])

I2_max = maximum(I2_prevalence)
I2_empirical_distribution = get_distribution(I2_prevalence, I2_max)
reduce(hcat, [aggregate_p(uniformize(p0, t, γ, Q_mul!, N, α, λ, Q_bands), N, dims=1)' for t in tvec])'