# ------------------------------------------------------------
# SIR generator (matrix-free)
# ------------------------------------------------------------

struct SIRGenerator <: AbstractGenerator
    N::Int
    α::Float64
    β::Float64
    γ::Float64
    bands::NTuple{8, SparseMatrixCSC{Float64, Int}}
end

# ------------------------------------------------------------
# Construct band matrices
# ------------------------------------------------------------

function construct_sir_band_matrices(N::Int)

    I = spdiagm(0 => ones(Float64, N+1))

    S_plus_inf  = spdiagm(1 => Float64.(1:N))
    S_minus_inf = spdiagm(0 => Float64.(0:N))

    I_plus_inf  = spdiagm(-1 => Float64.(0:N-1))
    I_minus_inf = spdiagm(0 => Float64.([0:N-1; 0]))

    S_plus_rec  = I
    S_minus_rec = I

    I_plus_rec  = spdiagm(1 => Float64.(1:N))
    I_minus_rec = spdiagm(0 => Float64.(0:N))

    return (
        S_plus_inf, S_minus_inf,
        I_plus_inf, I_minus_inf,
        S_plus_rec, S_minus_rec,
        I_plus_rec, I_minus_rec
    )
end

# ------------------------------------------------------------
# Compute exact uniformization rate γ
# ------------------------------------------------------------

function compute_gamma(N::Int,
                       α::Float64,
                       β::Float64)::Float64

    γ = 0.0
    for S in 0:N
        for I in 0:N
            rate = (β * S * I / N) + (α * I)
            γ = max(γ, rate)
        end
    end
    return γ
end

# ------------------------------------------------------------
# Constructor
# ------------------------------------------------------------

function SIRGenerator(N::Int,
                      α::Float64,
                      β::Float64)

    bands = construct_sir_band_matrices(N)
    γ = compute_gamma(N, α, β)

    return SIRGenerator(N, α, β, γ, bands)
end

# ------------------------------------------------------------
# Uniformization rate
# ------------------------------------------------------------

uniformization_rate(gen::SIRGenerator) = gen.γ

# ------------------------------------------------------------
# Q multiplication
# ------------------------------------------------------------

function Q_mul!(out::Vector{Float64},
                v::Vector{Float64},
                gen::SIRGenerator)

    N = gen.N
    α = gen.α
    β = gen.β

    S⁺ᵢ, S⁻ᵢ,
    I⁺ᵢ, I⁻ᵢ,
    S⁺ʳ, S⁻ʳ,
    I⁺ʳ, I⁻ʳ = gen.bands

    V = reshape(v, N+1, N+1)
    QV = zeros(Float64, N+1, N+1)

    # Infection inflow
    QV .+= (β / N) .* (I⁺ᵢ * V * S⁺ᵢ')

    # Recovery inflow
    QV .+= α .* (I⁺ʳ * V * S⁺ʳ')

    # Infection outflow
    QV .-= (β / N) .* (I⁻ᵢ * V * S⁻ᵢ')

    # Recovery outflow
    QV .-= α .* (I⁻ʳ * V * S⁻ʳ')

    out .= vec(QV)
    return out
end