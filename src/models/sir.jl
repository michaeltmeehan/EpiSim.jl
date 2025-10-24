struct SIR <: AbstractModel
    β::Float64      # Transmission rate
    α::Float64      # Recovery rate
    ψ::Float64      # Sampling rate
    ε::Float64      # Seeding rate
end


mutable struct SIRCount <: AbstractCount
    S::Int
    I::Int
    R::Int
    N::Int
end


mutable struct SIRAgent <: AbstractAgent
    S::Int
    I::Int
    R::Int
    N::Int
    active::Vector{Int}
    cum_infected::Int
    cum_sampled::Int
end


const SIREVENTS = [Transmission, Recovery, Sampling, Seeding]


get_events(model::SIR) = SIREVENTS


capture(state::Union{SIRCount, SIRAgent}) = SIRCount(state.S, state.I, state.R, state.N)


function update_rates!(λ::Vector{Float64}, model::SIR, state::AbstractState)
    S = state.S
    I = state.I
    N = state.N

    λ[1] = model.β * S * I / N      # Transmission
    λ[2] = model.α * I              # Recovery
    λ[3] = model.ψ * I              # Sampling
    λ[4] = model.ε                  # Seeding
    return nothing
end


# For each event in SIREvents, define how it updates the state
### Count states ###
@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRCount, 
                               ::Type{Transmission})::Transmission
    state.S -= 1
    state.I += 1
    return Transmission(0, 0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRCount, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    state.R += 1
    return Recovery(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRCount, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    state.R += 1
    return Sampling(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRCount, 
                               ::Type{Seeding})::Seeding
    state.I += 1
    state.N += 1
    return Seeding(0)
end


### Agent states ###
@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRAgent, 
                               ::Type{Transmission})::Transmission
    state.S -= 1
    state.I += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    infector = rand(rng, state.active)
    push!(state.active, infectee)
    return Transmission(infector, infectee)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRAgent, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    state.R += 1
    recovered = popr!(rng, state.active)
    return Recovery(recovered)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRAgent, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    state.R += 1
    sampled = popr!(rng, state.active)
    state.cum_sampled += 1
    return Sampling(sampled)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SIR,
                               state::SIRAgent, 
                               ::Type{Seeding})::Seeding
    state.I += 1
    state.N += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    push!(state.active, infectee)
    return Seeding(infectee)
end