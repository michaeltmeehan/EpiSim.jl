struct MTBDParameters <: AbstractParameters
    birth_rate::Matrix{Float64}
    death_rate::Vector{Float64}
    sampling_rate::Vector{Float64}
end


const MTBD_EVENT_TYPES = [Transmission, Recovery, Sampling]


mutable struct AgenticMTBDState <: AgenticState
    t::Float64
    I::Vector{Int}
    currently_infected::Vector{Vector{Int}}
    n_sampled::Int
    n_cumulative::Int
end


function AgenticMTBDState(; t::Float64=0., 
                            I::Vector{Int}=[1, 0])
    currently_infected = [collect(1:I[1]), collect(1:I[2])]
    n_sampled = 0
    n_cumulative = sum(I) 
    return AgenticMTBDState(t, I, currently_infected, n_sampled, n_cumulative)
end


mutable struct AggregateMTBDState <: AggregateState
    t::Float64
    I::Vector{Int}
end


AggregateMTBDState(; t::Float64=0., I::Vector{Int}=[1, 0]) = AggregateMTBDState(t, I)


MTBDState = Union{AgenticMTBDState, AggregateMTBDState}


function capture(state::MTBDState)
    labels = [:t, Symbol.("I_", 1:length(state.I))...]
    values = (state.t, state.I...)
    return (; zip(labels, values)...)
end


function MTBDModel(; birth_rate::Matrix{Float64}=[1. 0.1; 20. 2.],
                    death_rate::Vector{Float64}=[0.9, 0.9], 
                    sampling_rate::Vector{Float64}=[0.1, 0.1], 
                    I::Vector{Int}=[1, 0], 
                    agentic::Bool=true)
    par = MTBDParameters(birth_rate, death_rate, sampling_rate)
    state = agentic ?
        AgenticMTBDState(; I=I) :
        AggregateMTBDState(; I=I)
    return Model(par, MTBD_EVENT_TYPES, state)
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     par::MTBDParameters, 
                                     state::MTBDState)
    λ = par.birth_rate
    μ = par.death_rate
    ψ = par.sampling_rate
    Λ = sum(λ, dims=2)[:]  # Total birth rate for each type

    I = state.I

    event_rates[1] = Λ ⋅ I # Infection rate
    event_rates[2] = μ ⋅ I # Recovery rate
    event_rates[3] = ψ ⋅ I # Sampling rate
end


@inline function sample_types(rng::AbstractRNG,
                              w::Matrix{Float64})::Tuple{Int, Int}
    w ./= sum(w)
    n = size(w, 1)
    r = rand(rng)
    cumulative_probability = w[1]
    idx = 1
    while r > cumulative_probability && idx < n^2
        idx += 1
        cumulative_probability += w[idx]
    end
    return mod(idx - 1, n) + 1, div(idx - 1, n) + 1
end


@inline function sample_type(rng::AbstractRNG,
                             w::Vector{Float64})::Int
    w ./= sum(w)
    n = length(w)
    r = rand(rng)
    cumulative_probability = w[1]
    idx = 1
    while r > cumulative_probability && idx < n
        idx += 1
        cumulative_probability += w[idx]
    end
    return idx
end


# TODO: Perhaps still return an event log aggregate models
@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AggregateMTBDState, 
                               ::Type{Transmission})::Nothing
    weights = par.birth_rate * diagm(state.I)
    child_type, _ = sample_types(rng, weights)
    state.I[child_type] += 1
    return Transmission(nothing, nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AgenticMTBDState, 
                               ::Type{Transmission})::Transmission
    weights = par.birth_rate * diagm(state.I)
    child_type, parent_type = sample_types(rng, weights)
    state.I[child_type] += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected[parent_type])
    push!(state.currently_infected[child_type], infectee)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AggregateMTBDState, 
                               ::Type{Recovery})::Nothing
    # Recovery event
    weights = par.death_rate .* state.I
    recovery_type = sample_type(rng, weights)
    state.I[recovery_type] -= 1
    return Recovery(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AgenticMTBDState, 
                               ::Type{Recovery})::Recovery
    # Recovery event
    weights = par.death_rate .* state.I
    recovery_type = sample_type(rng, weights)
    state.I[recovery_type] -= 1
    recovered = pop_random!(rng, state.currently_infected[recovery_type])
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AggregateMTBDState, 
                               ::Type{Sampling})::Nothing
    # Recovery event
    weights = par.sampling_rate .* state.I
    sampled_type = sample_type(rng, weights)
    state.I[sampled_type] -= 1
    return Sampling(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::MTBDParameters,
                               state::AgenticMTBDState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    weights = par.sampling_rate .* state.I
    sampled_type = sample_type(rng, weights)
    state.I[sampled_type] -= 1
    sampled = pop_random!(rng, state.currently_infected[sampled_type])
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end