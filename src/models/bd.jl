struct BD <: AbstractModel
    λ::Float64      # Birth rate
    μ::Float64      # Death rate
    ψ::Float64      # Sampling rate
    ε::Float64      # Seeding rate
end


mutable struct BDCount <: AbstractCount
    I::Int
end


mutable struct BDAgent <: AbstractAgent
    I::Int
    active::Vector{Int}
    cum_infected::Int
    cum_sampled::Int
end


const BDEVENTS = [Transmission, Recovery, Sampling, Seeding]


get_events(model::BD) = BDEVENTS


capture(state::Union{BDCount, BDAgent}) = BDCount(state.I)


function update_rates!(λ::Vector{Float64}, model::BD, state::AbstractState)
    I = state.I

    λ[1] = model.λ * I      # Birth
    λ[2] = model.μ * I      # Death
    λ[3] = model.ψ * I      # Sampling
    λ[4] = model.ε          # Seeding
    return nothing
end


function initialize_rates(model::BD, state::AbstractState)
    I = state.I

    λ = Vector{Float64}(undef, 4)

    λ[1] = model.λ * I      # Birth
    λ[2] = model.μ * I      # Death
    λ[3] = model.ψ * I      # Sampling
    λ[4] = model.ε          # Seeding
    return λ
end


function update!(rng::AbstractRNG,
                 model::BD,
                 state::AbstractState,
                 λ::Vector{Float64},
                 event_type::Type{Seeding})

    # Update state and return event
    state.I += 1

    # Update rates
    λ[1] += model.λ      # Birth
    λ[2] += model.μ      # Death
    λ[3] += model.ψ      # Sampling

    return Seeding(0)
end


function update!(rng::AbstractRNG,
                 model::BD,
                 state::AbstractState,
                 λ::Vector{Float64},
                 event_type::Type{Transmission})

    # Update state and return event
    state.I += 1

    # Update rates
    λ[1] += model.λ      # Birth
    λ[2] += model.μ      # Death
    λ[3] += model.ψ      # Sampling

    return Transmission(0, 0)
end


function update!(rng::AbstractRNG,
                 model::BD,
                 state::AbstractState,
                 λ::Vector{Float64},
                 event_type::Type{Recovery})

    # Update state and return event
    state.I -= 1

    # Update rates
    λ[1] -= model.λ      # Birth
    λ[2] -= model.μ      # Death
    λ[3] -= model.ψ      # Sampling

    return Recovery(0)
end


function update!(rng::AbstractRNG,
                 model::BD,
                 state::AbstractState,
                 λ::Vector{Float64},
                 event_type::Type{Sampling})

    # Update state and return event
    state.I -= 1

    # Update rates
    λ[1] -= model.λ      # Birth
    λ[2] -= model.μ      # Death
    λ[3] -= model.ψ      # Sampling

    return Sampling(0)
end

# For each event in BDEVENTS, define how it updates the state
### Count states ###
@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDCount, 
                               ::Type{Transmission})::Transmission
    state.I += 1
    return Transmission(0, 0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDCount, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    return Recovery(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDCount, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    return Sampling(0)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDCount, 
                               ::Type{Seeding})::Seeding
    state.I += 1
    return Seeding(0)
end


### Agent states ###
@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDAgent, 
                               ::Type{Transmission})::Transmission
    state.I += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    infector = rand(rng, state.active)
    push!(state.active, infectee)
    return Transmission(infector, infectee)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDAgent, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    recovered = popr!(rng, state.active)
    return Recovery(recovered)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDAgent, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    state.cum_sampled += 1
    sampled = popr!(rng, state.active)
    return Sampling(sampled)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BD,
                               state::BDAgent, 
                               ::Type{Seeding})::Seeding
    state.I += 1
    state.cum_infected += 1
    infectee = state.cum_infected
    push!(state.active, infectee)
    return Seeding(infectee)
end