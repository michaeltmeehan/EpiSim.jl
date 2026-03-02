# =========================================
# SEIR Generator (DU benchmark)
# =========================================

struct SEIRGenerator <: AbstractGenerator
    N::Int
    α::Float64
    β::Float64
    σ::Float64
    index_state::Vector{NTuple{3,Int}}
    state_index::Dict{NTuple{3,Int},Int}
    γ::Float64
end


# -----------------------------------------
# Constructor
# -----------------------------------------

function SEIRGenerator(N::Int, α::Float64, β::Float64, σ::Float64)

    index_state = NTuple{3,Int}[]
    state_index = Dict{NTuple{3,Int},Int}()

    # Enumerate all (S,E,I) such that S+E+I ≤ N
    for S in 0:N
        for E in 0:(N-S)
            for I in 0:(N-S-E)
                push!(index_state, (S,E,I))
                state_index[(S,E,I)] = length(index_state)
            end
        end
    end

    # Compute γ (max exit rate)
    γ = 0.0
    for (S,E,I) in index_state
        rate = (β * S * I / N) + (σ * E) + (α * I)
        γ = max(γ, rate)
    end

    return SEIRGenerator(N, α, β, σ,
                         index_state,
                         state_index,
                         γ)
end


# -----------------------------------------
# Uniformization rate
# -----------------------------------------

uniformization_rate(gen::SEIRGenerator) = gen.γ


# -----------------------------------------
# Generator multiplication
# -----------------------------------------

function Q_mul!(out::Vector{Float64},
                v::Vector{Float64},
                gen::SEIRGenerator)

    fill!(out, 0.0)

    N = gen.N
    α = gen.α
    β = gen.β
    σ = gen.σ

    index_state = gen.index_state
    state_index = gen.state_index

    @inbounds for k in eachindex(index_state)

        S,E,I = index_state[k]
        vk = v[k]
        λ = 0.0

        # Infection: S -> S-1, E -> E+1
        if S > 0 && I > 0
            rate = β * S * I / N
            k2 = state_index[(S-1, E+1, I)]
            out[k2] += rate * vk
            λ += rate
        end

        # Activation: E -> E-1, I -> I+1
        if E > 0
            rate = σ * E
            k2 = state_index[(S, E-1, I+1)]
            out[k2] += rate * vk
            λ += rate
        end

        # Recovery: I -> I-1
        if I > 0
            rate = α * I
            k2 = state_index[(S, E, I-1)]
            out[k2] += rate * vk
            λ += rate
        end

        # Diagonal
        out[k] -= λ * vk
    end

    return out
end