struct MTBD <: AbstractModel
    λ::Matrix{Float64}      # Birth rate
    μ::Vector{Float64}      # Death rate
    ψ::Vector{Float64}      # Sampling rate
    ε::Vector{Float64}      # Seeding rate
end


mutable struct MTBDCount <: AbstractCount
    I::Vector{Int}
end


mutable struct MTBDAgent <: AbstractAgent
    I::Vector{Int}
    active::Vector{Vector{Int}}
    cum_infected::Int
    cum_sampled::Int
end


const MTBDEVENTS = [Transmission, Recovery, Sampling, Seeding]


get_events(model::MTBD) = MTBDEVENTS


capture(state::Union{MTBDCount, MTBDAgent}) = MTBDCount(Vector{Int}(state.I))


function update_rates!(λ::Vector{Float64}, model::MTBD, state::AbstractState)
    I = state.I

    λ[1] = zero(Float64)
    @inbounds for col in eachcol(model.λ)
        λ[1] += sum(col .* I)
    end     # Birth
    λ[2] = model.μ ⋅ I      # Death
    λ[3] = model.ψ ⋅ I      # Sampling
    λ[4] = zero(Float64)
    @inbounds for rate in model.ε
        λ[4] += rate
    end      # Seeding
    return nothing
end


# For each event in MTBDEVENTS, define how it updates the state
### Count states ###
@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDCount, 
                               ::Type{Transmission})::Transmission
    type = wsampleindex(rng, model.λ * state.I)
    state.I[type] += 1
    return Transmission(0, 0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDCount, 
                               ::Type{Recovery})::Recovery
    type = wsampleindex(rng, model.μ .* state.I)
    state.I[type] -= 1
    return Recovery(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDCount, 
                               ::Type{Sampling})::Sampling
    type = wsampleindex(rng, model.ψ .* state.I)
    state.I[type] -= 1
    return Sampling(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDCount, 
                               ::Type{Seeding})::Seeding
    type = wsampleindex(rng, model.ε)
    state.I[type] += 1
    return Seeding(0)
end


### Agent states ###
@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDAgent, 
                               ::Type{Transmission})::Transmission
    infectee_type, infector_type = wsampleindex_cols(rng, model.λ, state.I)
    state.I[infectee_type] += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    infector = rand(rng, state.active[infector_type])
    push!(state.active[infectee_type], infectee)
    return Transmission(infector, infectee)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDAgent, 
                               ::Type{Recovery})::Recovery
    type = wsampleindex(rng, model.μ .* state.I)
    state.I[type] -= 1
    recovered = popr!(rng, state.active[type])
    return Recovery(recovered)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDAgent, 
                               ::Type{Sampling})::Sampling
    type = wsampleindex(rng, model.ψ .* state.I)
    state.I[type] -= 1
    state.cum_sampled += 1
    sampled = popr!(rng, state.active[type])
    return Sampling(sampled)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MTBD,
                               state::MTBDAgent, 
                               ::Type{Seeding})::Seeding
    type = wsampleindex(rng, model.ε)
    state.I[type] += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    push!(state.active[type], infectee)
    return Seeding(infectee)
end


function SSBD(R::Float64, ρ::Float64, c::Float64, μ::Vector{Float64}, ψ::Vector{Float64}, ε::Vector{Float64})
    δ = μ .+ ψ
    λ = R / (1. - (1. - ρ) * (1. - c)) * [c, 1-c] * [δ[1] ρ * δ[2]]
    return MTBD(λ, μ, ψ, ε)
end


SSBDCount(I::Vector{Int}) = MTBDCount(I)
SSBDAgent(I::Vector{Int}) = MTBDAgent(I, [(i == 1) ? [1] : Int[] for _ in 1:length(I)], 1, 0)